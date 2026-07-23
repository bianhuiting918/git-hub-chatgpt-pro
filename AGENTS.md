# Repository working agreements

## Default coding surface

- Perform code, documentation, configuration, branch, commit, and pull-request work through this GitHub repository and Codex Web by default.
- Unless the user explicitly authorizes a specific local write, local workspaces are read-only: do not create or modify project files, worktrees, branches, commits, intermediate results, or caches.
- When writes are required, use a GitHub branch and let Codex Web implement and verify the change. If GitHub or Codex Web is unavailable, stop and report the blocker; do not silently fall back to local edits.
- Scientific computations explicitly assigned to a remote compute server still run in that server's project directory and must retain reproducible scripts, a `RUNBOOK.md`, and append-only run records.
- A user authorization for one specific local write is limited to that target and does not grant continuing local-write permission.

## Project context

- Keep verified PASS artifacts immutable unless the user explicitly requests recomputation.
- Separate topology, geometry, resource, convergence, and scientific-result gates. Technical completion alone is not a scientific PASS.
- Do not store credentials, tokens, private keys, or passwords in the repository, scripts, logs, or chat output.
