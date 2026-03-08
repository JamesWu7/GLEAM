# GLEAM API Migration: pathway -> signature

GLEAM now uses **signature** as the canonical umbrella term in the public API.
Legacy `*_pathway*` names remain available as compatibility wrappers.

## Canonical names

- `score_signature()` (preferred) | compatibility: `score_pathway()`
- `aggregate_signature()` (preferred) | compatibility: `aggregate_pathway()`
- `test_signature()` (preferred) | compatibility: `test_pathway()`
- `differential_signature()` (preferred) | compatibility: `differential_pathway()`
- `test_signature_spatial()` (preferred) | compatibility: `test_pathway_spatial()`
- `differential_signature_spatial()` (preferred) | compatibility: `differential_pathway_spatial()`
- `test_signature_trajectory()` (preferred) | compatibility: `test_pathway_trajectory()`
- `differential_signature_trajectory()` (preferred) | compatibility: `differential_pathway_trajectory()`

## Workflow entry points

- Canonical: `run_gleam()`
- Compatibility wrapper: `run_scpathway()`

## Notes

- Plotting functions remain separate by design.
- Backward-compatible wrappers are kept to reduce migration friction.
