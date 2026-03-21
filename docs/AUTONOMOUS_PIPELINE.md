# Autonomous Development Pipeline (Resident Mode)

This repository uses a streamlined autonomous maintenance pipeline:

`analyze -> root-cause -> fix -> validate -> review -> auto-push PR`

## Safety Model

- No direct push to `main`/`master`.
- Only `fix_engineer_agent` is allowed to generate code changes.
- Validation runs before auto-push.
- Changes are pushed to a `codex/*` branch and PR is created automatically.

## Agent Topology (Simplified)

The workflow `.github/workflows/autonomous-dev-pipeline.yaml` keeps only non-duplicated roles:

1. `ci_explorer_agent`
2. `pkgdown_explorer_agent`
3. `bug_analyzer_agent`
4. `fix_engineer_agent`
5. `test_runner_agent`
6. `final_reviewer_agent`
7. `pr_writer_agent`
8. `supervisor_agent`

## Triggering

- Manual: GitHub Actions -> `autonomous-dev-pipeline` -> `Run workflow`
- Scheduled: daily (UTC) by cron in workflow file.

## Auto-Push and Auto-Repair

- After review + validation pass, supervisor auto-pushes fix branch and opens PR.
- If PR checks fail (`R-CMD-check` / `pkgdown`), `.github/workflows/auto-repair-check-failure.yaml`:
  - pulls failed logs via GitHub API,
  - reruns root-cause + fix scripts,
  - validates locally,
  - pushes follow-up repair commit to the same branch.

## Artifacts

Each stage uploads artifacts so diagnostics stay visible on GitHub:

- root-cause report
- fix patch
- check/pkgdown logs
- final review summary
- generated PR metadata

## Notes

- Optional `monocle3` checks remain in manual workflow (`optional-monocle3.yaml`).
- `pkgdown` deploy remains via GitHub Pages artifact deployment (no docs commit to main).
