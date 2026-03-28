"""VariantAgent CLI — command-line interface for variant interpretation."""

from __future__ import annotations

import json
import re
import sys

import typer
from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.tree import Tree

app = typer.Typer(
    name="variantagent",
    help="Multi-agent clinical variant interpretation system",
    no_args_is_help=True,
)
console = Console()


def _parse_variant_string(variant_str: str) -> dict:
    """Parse a variant string like 'chr17:7674220 G>A' into components."""
    # Pattern: chr17:7674220 G>A  or  17:7674220:G:A  or  chr17:7674220 G A
    patterns = [
        r"(?P<chr>(?:chr)?\d+):(?P<pos>\d+)\s+(?P<ref>[ACGT]+)>(?P<alt>[ACGT]+)",
        r"(?P<chr>(?:chr)?\d+):(?P<pos>\d+):(?P<ref>[ACGT]+):(?P<alt>[ACGT]+)",
        r"(?P<chr>(?:chr)?\d+):(?P<pos>\d+)\s+(?P<ref>[ACGT]+)\s+(?P<alt>[ACGT]+)",
    ]
    for pattern in patterns:
        match = re.match(pattern, variant_str.strip(), re.IGNORECASE)
        if match:
            return match.groupdict()
    raise typer.BadParameter(
        f"Could not parse variant: '{variant_str}'. "
        f"Expected format: 'chr17:7674220 G>A' or '17:7674220:G:A'"
    )


@app.command()
def analyze(
    variant_str: str = typer.Argument(
        help="Variant to analyze (e.g., 'chr17:7674220 G>A')"
    ),
    gene: str | None = typer.Option(None, "--gene", "-g", help="Gene symbol (e.g., TP53)"),
    sample_id: str | None = typer.Option(None, "--sample", "-s", help="Sample identifier"),
    batch_id: str | None = typer.Option(None, "--batch", "-b", help="Batch identifier"),
    output: str | None = typer.Option(None, "--output", "-o", help="Output file path (JSON)"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Show detailed provenance"),
) -> None:
    """Analyze a genetic variant using the multi-agent system."""
    from variantagent.agents.orchestrator import analyze_variant
    from variantagent.models.variant import Variant

    # Parse variant
    parsed = _parse_variant_string(variant_str)

    variant = Variant(
        chromosome=parsed["chr"],
        position=int(parsed["pos"]),
        reference=parsed["ref"],
        alternate=parsed["alt"],
        gene=gene,
    )

    console.print(
        Panel(
            f"[bold]Variant:[/bold] {variant.variant_id}\n"
            f"[bold]Gene:[/bold] {gene or 'not specified'}\n"
            f"[bold]Sample:[/bold] {sample_id or 'not specified'}",
            title="[blue]VariantAgent[/blue]",
            border_style="blue",
        )
    )

    with console.status("[bold blue]Running multi-agent analysis...", spinner="dots"):
        report = analyze_variant(variant, sample_id=sample_id, batch_id=batch_id)

    # --- Display results ---

    # Classification
    if report.classification:
        cls = report.classification
        color = {
            "Pathogenic": "red",
            "Likely Pathogenic": "red",
            "Uncertain Significance": "yellow",
            "Likely Benign": "green",
            "Benign": "green",
        }.get(cls.classification.value, "white")

        console.print(Panel(
            f"[bold {color}]{cls.classification.value}[/bold {color}]\n"
            f"Confidence: {report.overall_confidence:.0%}\n"
            f"Rule: {cls.classification_rule}\n"
            f"Evidence codes: {', '.join(cls.applied_codes_summary) or 'none'}",
            title="[bold]Classification[/bold]",
            border_style=color,
        ))

    # QC Assessment
    if report.qc_assessment:
        qc = report.qc_assessment
        qc_color = {"pass": "green", "warn": "yellow", "fail": "red"}.get(qc.overall_status.value, "white")
        console.print(f"\n[bold]QC Status:[/bold] [{qc_color}]{qc.overall_status.value.upper()}[/{qc_color}]")
        if qc.issues:
            for issue in qc.issues:
                console.print(f"  [{qc_color}]•[/{qc_color}] {issue.metric}: {issue.description}")

    # Annotation summary
    if report.annotation:
        ann = report.annotation
        table = Table(title="Annotation Summary", show_header=True)
        table.add_column("Source", style="cyan")
        table.add_column("Result")

        if ann.clinvar.found:
            table.add_row("ClinVar", f"{ann.clinvar.clinical_significance} ({ann.clinvar.review_stars}★)")
        else:
            table.add_row("ClinVar", "[dim]Not found[/dim]")

        if ann.gnomad.found and ann.gnomad.overall_af is not None:
            table.add_row("gnomAD AF", f"{ann.gnomad.overall_af:.6f}")
        else:
            table.add_row("gnomAD", "[dim]Not found[/dim]")

        if ann.ensembl_vep.found:
            vep_parts = [ann.ensembl_vep.consequence_type or ""]
            if ann.ensembl_vep.sift_prediction:
                vep_parts.append(f"SIFT: {ann.ensembl_vep.sift_prediction}")
            if ann.ensembl_vep.polyphen_prediction:
                vep_parts.append(f"PolyPhen: {ann.ensembl_vep.polyphen_prediction}")
            table.add_row("VEP", " | ".join(vep_parts))
        else:
            table.add_row("VEP", "[dim]Not found[/dim]")

        console.print(table)

    # Reviewer concerns
    concerns = [f for f in report.reviewer_findings if f.concern]
    if concerns:
        console.print(f"\n[bold yellow]Reviewer Concerns ({len(concerns)}):[/bold yellow]")
        for finding in concerns:
            risk_color = {"high": "red", "medium": "yellow", "low": "green"}.get(
                finding.hallucination_risk, "white"
            )
            console.print(f"  [{risk_color}]⚠[/{risk_color}] {finding.concern}")

    # Summary
    console.print(f"\n[dim]{report.natural_language_summary}[/dim]")

    # Provenance (verbose only)
    if verbose:
        console.print("\n[bold]Provenance Trail:[/bold]")
        tree = Tree("[bold]Analysis[/bold]")
        for entry in report.provenance:
            node = tree.add(f"[cyan]{entry.agent}[/cyan]: {entry.action}")
            node.add(f"Input: {entry.input_summary}")
            node.add(f"Output: {entry.output_summary}")
            if entry.duration_ms is not None:
                node.add(f"Duration: {entry.duration_ms}ms")
            if entry.error:
                node.add(f"[red]Error: {entry.error}[/red]")
        console.print(tree)

    # JSON output
    if output:
        report_json = report.model_dump_json(indent=2)
        with open(output, "w") as f:
            f.write(report_json)
        console.print(f"\n[green]Report saved to {output}[/green]")

    # Trace ID for reference
    console.print(f"\n[dim]Trace ID: {report.trace_id}[/dim]")


@app.command()
def serve(
    host: str = typer.Option("0.0.0.0", help="Host to bind to"),
    port: int = typer.Option(8000, help="Port to bind to"),
) -> None:
    """Start the FastAPI server."""
    import uvicorn

    uvicorn.run("variantagent.api.app:app", host=host, port=port, reload=True)


@app.command()
def version() -> None:
    """Show version information."""
    from variantagent import __version__

    console.print(f"VariantAgent v{__version__}")


if __name__ == "__main__":
    app()
