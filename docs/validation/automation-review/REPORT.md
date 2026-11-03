# Actions and Docker Hub convention/safety audit

Branch: `audit/release-candidate-0.4.0`, still based on
`fe224c81c41e4aee5905e587ed01e0c1829ca0a3`, with uncommitted work preserved.
This is a bounded automation follow-up to the
[scientific release review](../final-review/GO-NO-GO.md), not a new scientific
audit. **GO for the first authorized candidate CI run; NO-GO for publication.**
The candidate remains `0.4.0-rc.1`, which the stable publisher intentionally
rejects. No commits, staging, push, merge, branch creation, tag, publication or
live settings changes were performed by this follow-up.

## Preservation and reference

Before editing, all incoming candidate files, their hashes and the full tracked
binary diff were backed up outside the repository. The scientific runtime,
Dockerfile, font assets, existing scientific tests, skill and all earlier
validation records remain unchanged, including the 15 baseline observations,
original figures, numerical evidence and original source manifest. That
manifest and the final-review inventory describe the preceding review state;
they intentionally retain the old workflow filename and hashes. The
[supplemental manifest](source-manifest.json) identifies this automation update.
No historical evidence was rewritten to appear to validate these later edits.

Inspected every local workflow, the Docker context/recipe, shared container
checks, skill validator, rendering smoke and release instructions. Read dndsR's
`publish.yaml`, `container.yaml`, `docs/dev/GITHUB-SETUP.md` and
`docs/dev/RELEASE.md` at local reference commit
`5983405287ed8ecc62feeaedb15f176ef09b9790`. Its unrelated pre-existing working
tree edit was untouched. Public GitHub ruleset detail GETs supplement the local
reference; [responses are preserved](live-rulesets.json). No upstream repo was
modified.

## Findings and bounded corrections

| Finding | Evidence / severity | Resolution or remaining action |
|---|---|---|
| Publisher filename differed from convention | `.github/workflows/publish.yml`; low | Renamed to **`.github/workflows/publish.yaml`**, without a duplicate publisher. |
| Prerelease tags were accepted | Previous regex permitted `-rc.N`; medium policy mismatch | Stable `vX.Y.Z` only, no leading zeroes, prerelease or build suffix. This does not change scientific behavior. |
| Approval was only optional in instructions | Environment name alone cannot enforce approval; high if enabled without reviewers | Setup now mandates Required reviewers and selected release tags. Actual environment protection remains unverified; keep publication disabled until verified. |
| Immutable publication was overstated by a tag lookup | HTTP check is not atomic; high if another writer races | Repeat lookup before login; require Docker Hub immutable version tags and remove competing publishers. Registry setting remains unverified. |
| Main history was checked but review enforcement incomplete | Live `protect_main` has PR protection, no required checks; high release setup gap | Add `fast` and `build-and-test` after first CI. Explicit fresh main-ref fetch and checkout/event-SHA equality added. Ancestry is not independent proof of review. |
| Required checks need stable, unconditional PR execution | Previous names were implicit; push path filters could omit builds | Explicit names retained; both workflows run for every PR to main/dev and every canonical-branch push, without path filters. |
| Shell syntax check covered only its first script | `bash -n script1 script2` treats the rest as arguments; low validation gap | Explicit loop checks all five shell files in fast CI and shared container checks. No calculation or test expectation changed. |

Checkout credentials are no longer persisted. Added a small offline check of
the actual workflow shell steps using temporary real Git history and fake
registry/Docker commands. No new runtime dependency or orchestration framework
was introduced; PyYAML was already installed by fast CI.

## Exact configuration and safety assessment

| File | Workflow name | Job ID / check name |
|---|---|---|
| `.github/workflows/validate.yml` | `Validate` | `fast` |
| `.github/workflows/container.yml` | `Clean container` | `build-and-test` |
| `.github/workflows/publish.yaml` | `Publish versioned image` | `publish` |

Required PR checks are **`fast`** and **`build-and-test`**, not the workflow
display names. `publish` is tag-only and must not be a required PR check.
Deployment environment: **`dockerhub-release`**. Docker Hub image:
**`scrumpis/iwgc-circos`**, with only **`vX.Y.Z`** published.

| Requested property | Verified scope / limit |
|---|---|
| Clean PR build and real end-to-end tests | Workflow checks out GitHub's PR test merge in a fresh hosted runner; builds with `--pull --no-cache`; executes full suite, independent arithmetic, actual Circos/SVG/PNG, gap fidelity and human example in that image. Defined and statically checked; runner execution pending. |
| Approved main history | HEAD must equal event SHA and be an ancestor of freshly fetched main. Offline approved/unapproved-history cases pass. Main required-check settings still need adding; deployment approval is separately required. |
| Exact tagged source rebuilt/tested | Publisher checks tag/VERSION/skill identity, builds once, tests and pushes that exact local image; no rebuild after testing. Actual Docker execution pending. |
| Approved deployment | Publishing job references `dockerhub-release`. Required reviewers must be enabled in settings; this cannot be guaranteed by naming the environment. Not live-verified. |
| Disabled by default | Only `vars.DOCKERHUB_PUBLISH_ENABLED == 'true'` permits the job. Live repository variable value is unverified. |
| No overwrite/latest | Existing tags and lookup errors fail before build and again before login. Only the event version tag is pushed, with per-tag workflow concurrency. Registry immutability is mandatory to prevent external races; live setting unverified. |
| Secrets confined to publisher | Only publish job's push step references `DOCKERHUB_USERNAME` and `DOCKERHUB_TOKEN`. No checking workflow references secrets or environment. Confirm secrets exist only in the environment, without repository/org duplicates. |
| Versioned skill | Actual metadata validator runs before build and inside the full image; release artifacts include `iwgc-circos-skill-vX.Y.Z.tar.gz` and manifest SHA256. Artifact creation/upload/published digest untested until real release. |
| Stable check names | Explicit `fast` and `build-and-test`; no matrices, conditional jobs or path filters that omit PR checks. Select the actual GitHub Actions check contexts after first run. |

The shared full-image checker preserves the existing numerical tests and the
known SVG/native gap comparison unchanged. A failure stops publication. It
does not convert an incomplete pipeline into apparent success. Logs/artifacts
use fresh output directories; source is mounted read-only. The image is still
a dependency environment with licensed font, used with a matching checkout.

## Live settings observed and still required

Read-only public inspection found only `main` at the original base SHA. There
is no remote `dev` yet. Visible target rulesets:

- [`protect_main`, 26406147](https://github.com/Scrumpis/iwgc-circos/rules/26406147)
  is active for main, prevents deletion/force push and requires PRs with zero
  approvals. **No required status checks are configured.** Add `fast` and
  `build-and-test` from GitHub Actions with strict up-to-date enforcement after
  the first run makes the check names selectable.
- [`protect_tags`, 26406246](https://github.com/Scrumpis/iwgc-circos/rules/26406246)
  is active for `v*`, restricting creation, updates and deletion. Its bypass is
  a repository administrative role, broader than dndsR's single-user bypass;
  review membership and narrow to intended release maintainers if needed.
- `protect_dev` is absent; create alongside the permanent branch after the
  initial release merge, preventing deletion and force push as in dndsR.

The connected reader supports public rulesets but not administration-only
environment/secret/variable settings. No secret values were requested or read.
Required environment reviewers, tag-only deployment policy, environment secret
placement, repository variable value and Docker Hub immutability/other build
rules therefore remain **unverified**, not implicitly satisfied.

Follow the exact setup table in the [release process](../../release-process.md):
required environment reviewers (manual self-approval is possible for one
maintainer), selected `v*` tags, no approval bypass, environment-only
`DOCKERHUB_USERNAME`/`DOCKERHUB_TOKEN`, immutable Docker version tags and no
competing publisher. Enable the repository variable **last**. GitHub's
[environment protection](https://docs.github.com/en/actions/reference/workflows-and-actions/deployments-and-environments)
and Docker's [immutable tag documentation](https://docs.docker.com/docker-hub/repos/manage/hub-images/immutable-tags/)
explain the distinct service-side guarantees. No workflow edit can install
those settings by itself.

## Executed checks and readiness

Exact commands, all exit **0**:

```bash
/usr/bin/python3 tests/check_workflows.py
actionlint -shellcheck='' -pyflakes=''
for script in iwgc_circos_tracks.sh create_configs.sh add_legend.sh preprocess/minimap2_self_alignment.sh tools/check_container.sh; do
  bash -n "$script" || exit
done
git diff --check
```

The [offline workflow record](workflow-checks.json) has **17 passing cases**:
stable approved source accepted; malformed/prerelease tags, VERSION/skill/SHA
mismatches, unapproved Git history, existing tags and registry/network errors
rejected; a fake Docker confirms only the exact version is passed to push.
These are guard/control tests, not a registry publication. [Check metadata](checks.json)
records actionlint 1.7.12; shellcheck/pyflakes integrations were disabled locally.

The earlier 62-test regression and 24-comparison scientific records are
preserved, not rerun or newly claimed to validate these workflow edits. Docker
is still unavailable locally, so **no new image ID, clean-container outcome or
GitHub Actions success is claimed**. The old SVG/native mismatch is still a
release gate, and independent scientific review remains distinct from the
same-author arithmetic checks.

The candidate is ready to be frozen and, **only after authorization**, pushed
for its first non-publishing CI run. Required deployment settings are not
prerequisites for that checking-only run. They and all scientific/rendering CI
gates must be satisfied before a stable main merge/tag/publication sequence.
The first CI may expose build or rendering failures; readiness to run is not
prediction of a passing result.

## Permanent dev strategy

After the initial release is merged into main, create `dev` directly from that
merged main commit and push it only when authorized. Install `protect_dev` and
keep main as the stable default branch. Future topic branches start at dev;
review them into dev, then release via dev-to-main PRs with both required checks.
Preserve permanent-branch ancestry with merge commits; synchronize main back
to dev by fast-forward where possible or an ordinary reviewed merge. No reset,
rebase, force push, initial extra merge or rewritten history is needed. Exact
recommended commands are in the human release process; none were executed.
