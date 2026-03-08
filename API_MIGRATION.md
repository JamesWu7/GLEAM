# GLEAM API Migration: pathway -> signature (breaking cleanup)

GLEAM now uses **signature** as the canonical umbrella term in the public API.
Legacy `*_pathway*` names are no longer exported in v0.2+.

## Canonical names

- `score_signature()`
- `aggregate_signature()`
- `test_signature()`
- `test_signature_spatial()`
- `test_signature_trajectory()`

## Workflow entry points

- `run_gleam()`

## Notes

- Plotting functions remain separate by design.
- Differential entrypoints were merged into the `test_signature*` family.
- Use `compare_scoring_methods()` (the `compare_methods()` alias is removed from exports).
