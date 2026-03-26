# ADR-002: MCP (Model Context Protocol) for Database Tool Interfaces

## Status
Accepted

## Context
The Annotation Agent needs to query multiple biological databases (ClinVar, gnomAD, Ensembl VEP). We need to decide how to expose these database interfaces to the agent.

## Decision
Implement database interfaces as **standalone MCP servers** that can be used independently of VariantAgent.

## Rationale
1. **Reusability.** MCP servers are standalone tools that any MCP-compatible client can use. A ClinVar MCP server has value beyond this project.
2. **Community contribution.** Anthropic has a life sciences MCP ecosystem. Contributing ClinVar, gnomAD, and Ensembl VEP servers adds to this ecosystem.
3. **Separation of concerns.** Database interaction logic (rate limiting, caching, error handling) is isolated from agent reasoning logic.
4. **Industry adoption.** MCP has 97M+ monthly SDK downloads and appears in 6/10 agentic AI job descriptions in biotech, trending upward.
5. **Testability.** MCP servers can be tested independently with mock data, without running the full agent system.

## Consequences
- Each MCP server must conform to the MCP specification (tools, resources, prompts).
- Additional setup step for users (MCP server must be running or configured).
- Docker Compose handles this automatically for the full-stack deployment.
