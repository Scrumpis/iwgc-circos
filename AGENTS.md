# iwgc-circos

This is a small shell-based genomics visualization tool. Human workflows are in
[README.md](README.md) and [docs/getting-started.md](docs/getting-started.md).
The portable end-user skill is [skills/iwgc-circos/SKILL.md](skills/iwgc-circos/SKILL.md).
Read it when installing, running, or interpreting the tool. It is not required
for human use.

For development, read [docs/validation/audit.md](docs/validation/audit.md) and
[docs/validation/verification.md](docs/validation/verification.md). Keep changes
small. Preserve original data and existing results. Derive scientific test
expectations independently of candidate output; distinguish execution,
regression, numerical agreement, and independent review. Record commands,
versions, failures and limitations. Do not silently change coordinate contracts,
transformations, filters or defaults: present before/after evidence for
maintainer review. Update the human docs and skill when behavior changes.

Fast checks: `python3 -m unittest discover -s tests -v` and
`python3 tools/validate_skill.py`. Rendering: `python3 tests/smoke.py` in the
supported dependency environment. See the verification record for exact
container commands. Work on development branches. Do not push, merge, tag,
publish, alter upstream repositories, or delete existing results without
explicit authorization. No new orchestration framework is needed.
