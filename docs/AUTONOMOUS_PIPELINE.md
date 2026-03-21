# Autonomous Development Pipeline (GitHub Review Driven)

This repository provides a semi-automated maintenance pipeline that follows:

`analyze -> fix -> validate -> PR -> GitHub human review -> merge`

## Safety Model

- No direct push to `main` in the autonomous pipeline.
- Only one stage (`fix_engineer_agent`) is allowed to generate code changes.
- Validation runs before PR creation.
- PR is created as **draft** and requires human approval on GitHub.

## Agent Topology

The workflow `.github/workflows/autonomous-dev-pipeline.yaml` maps requested agent stages to dedicated jobs:

1. `explorer_agent`
2. `reviewer_agent`
3. `ci_explorer_agent`
4. `pkgdown_explorer_agent`
5. `bug_analyzer_agent`
6. `fix_engineer_agent`
7. `test_runner_agent`
8. `final_reviewer_agent`
9. `pr_writer_agent`
10. `supervisor_agent` (creates branch + draft PR)

## Triggering

- Manual: GitHub Actions -> `autonomous-dev-pipeline` -> `Run workflow`
- Scheduled: daily (UTC) by cron in workflow file.

## Artifacts

Each stage uploads artifacts so validation and diagnostics are visible on GitHub:

- root-cause report
- fix patch
- check/pkgdown logs
- final review summary
- generated PR metadata

## Human Review Checklist (GitHub)

1. Open PR checks and inspect uploaded artifacts.
2. Verify docs links and pkgdown pages.
3. Verify tutorials and plots (UMAP/spatial/dot-bar/violin/ridge).
4. Confirm API consistency and no regressions.
5. Merge only after explicit approval.

## Notes

- Optional `monocle3` checks are kept in a manual workflow (`optional-monocle3.yaml`).
- `pkgdown` now builds/deploys via GitHub Pages artifact deployment (no commit back to `main`).
