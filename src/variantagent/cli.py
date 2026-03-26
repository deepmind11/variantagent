"""VariantAgent CLI — command-line interface for variant interpretation."""

import typer
from rich.console import Console
from rich.panel import Panel

app = typer.Typer(
    name="variantagent",
    help="Multi-agent clinical variant interpretation system",
    no_args_is_help=True,
)
console = Console()


@app.command()
def analyze(
    variant: str = typer.Argument(
        help="Variant to analyze (e.g., 'chr17:7674220 G>A' or path to VCF file)"
    ),
    sample_id: str | None = typer.Option(None, "--sample", "-s", help="Sample identifier"),
    batch_id: str | None = typer.Option(None, "--batch", "-b", help="Batch identifier"),
    output: str | None = typer.Option(None, "--output", "-o", help="Output file path (JSON)"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Show detailed provenance"),
) -> None:
    """Analyze a genetic variant using the multi-agent system."""
    console.print(
        Panel(
            f"[bold]Analyzing variant:[/bold] {variant}",
            title="VariantAgent",
            border_style="blue",
        )
    )
    # TODO: Wire up to orchestrator
    console.print("[yellow]Not yet implemented — scaffold only[/yellow]")


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
