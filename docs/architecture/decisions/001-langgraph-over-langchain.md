# ADR-001: LangGraph over basic LangChain

## Status
Accepted

## Context
We need an agent orchestration framework for a multi-agent variant interpretation system. Options considered:
- LangChain (AgentExecutor)
- LangGraph
- OpenAI Agents SDK
- CrewAI

## Decision
Use **LangGraph** as the primary orchestration framework.

## Rationale
1. **Stateful workflows.** Variant interpretation is a multi-step process where later agents depend on earlier results. LangGraph's state machine model makes this explicit and inspectable.
2. **Dynamic routing.** We need conditional edges (e.g., if QC fails, skip annotation). LangGraph's graph structure supports this natively.
3. **Production adoption.** LangGraph has 38M+ monthly PyPI downloads and is the most common framework in enterprise agent deployments (as of early 2026).
4. **Observability.** LangSmith integration provides trace-level visibility into every agent decision, which is critical for a system making clinical-adjacent assessments.
5. **Industry demand.** LangChain/LangGraph appears in 8/10 agentic AI job descriptions in biotech (based on analysis of 10+ JDs).

## Alternatives Considered
- **OpenAI Agents SDK:** Simpler API but lacks stateful orchestration and is tied to OpenAI's ecosystem. May add as a secondary implementation for comparison.
- **CrewAI:** Good for quick prototyping but teams frequently outgrow it for production use.
- **Basic LangChain AgentExecutor:** Insufficient for multi-agent routing with conditional logic.

## Consequences
- Developers must learn LangGraph's state machine concepts (nodes, edges, conditional edges).
- The system is not locked to a specific LLM provider — LangGraph works with any LangChain-compatible model.
