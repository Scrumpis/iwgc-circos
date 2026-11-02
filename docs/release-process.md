# Automated validation and Docker Hub releases

Stable releases live on `main`; ongoing development lives on `dev` after the
initial release. Temporary work uses `feature/*`, `fix/*` or `audit/*`, including
the current `audit/release-candidate-0.4.0` branch. The default branch stays `main`.
This follows dndsR's branch, environment and versioned publication conventions.

First CI update: the initial push of `a5940bc` ran both checking workflows;
both failed on missing environment dependencies. See the
[first-run diagnosis and bounded corrections](validation/ci-first-run.md).
The corrected clean Docker build still needs its next authorized push/run.

The workflows are prepared locally. No push, merge, tag, publication or live
settings change was made by the automation audit. The paragraphs below record
that pre-CI setup state. The candidate remains `0.4.0-rc.1`:
**ready for its first CI run after an authorized freeze/push, not for release**.
The clean image has not been built here; the recorded SVG/native rendering
mismatch remains a required failing gate until resolved in the new image.
See the [scientific release review](validation/final-review/GO-NO-GO.md) and
[automation audit](validation/automation-review/REPORT.md).

## Workflows and required checks

| File | Workflow name | Job ID / check name | Trigger and purpose |
|---|---|---|---|
| `.github/workflows/validate.yml` | `Validate` | `fast` | Native tests, skill, workflow guards and synthetic rendering |
| `.github/workflows/container.yml` | `Clean container` | `build-and-test` | Uncached full Docker build and real end-to-end validation |
| `.github/workflows/publish.yaml` | `Publish versioned image` | `publish` | Approved stable version-tag build/test/publication |

Both checking workflows run on pushes to `main`, `dev`, `feature/**`, `fix/**`
and `audit/**`; PRs targeting `main` or `dev`; and manual dispatch. There are
no path filters: even documentation-only PRs produce both required check names.
PR checkout uses GitHub's test merge commit, in a fresh hosted runner. Checkout
does not persist credentials. Neither checking job references Docker Hub secrets
or the deployment environment. There is no `pull_request_target` workflow.

`container.yml` builds with `--pull --no-cache --platform linux/amd64`. It tests
that same image with the checkout mounted read-only and a fresh output directory.
`tools/check_container.sh` requires the full dependency set, verifies the licensed
font, checks each shell script's syntax, validates the skill, runs the complete
regression suite and independent numerical checks, renders actual Circos/SVG and
labeled PNG output, checks the SVG/native histogram gaps, and runs the human
genes/GC example. The rendering gate is retained without changed expectations.
It records dependency versions if preceding checks pass. A failed check stops
publication; it is not evidence of a complete successful run.

Artifacts `smoke` and `clean-container-evidence` retain diagnostics (the latter
for 30 days). Retain release evidence separately before artifact expiry. The
job/check names `fast` and `build-and-test` are explicit and must remain stable;
select these actual GitHub Actions checks after the first run, rather than the
workflow display names. Do not require the tag-only `publish` check on PRs.

## Required GitHub and Docker Hub settings

These settings are not installed by workflow YAML. Live public rulesets were
read without modification; environment approval, credential placement, the
repository variable and Docker Hub immutability remain unverified here.
**Keep publication disabled until the settings below are verified.**

| Setting | Required configuration / current observation |
|---|---|
| `protect_main` | Active on `main`; require PRs, prevent deletion/force push, require `fast` and `build-and-test` from GitHub Actions with strict up-to-date checks. Currently active with PR protection but **no required status checks**; add them after the first CI run. Match dndsR's single-maintainer policy of 0 PR approvals if appropriate; deployment approval below is still mandatory. |
| `protect_dev` | Active on `dev`, prevent deletion/force push. Create when introducing `dev`; currently absent. Like dndsR, PRs need not be mandatory for this development branch. |
| `protect_tags` | Active on `v*`, restrict creation, updates and deletion to authorized release maintainers. Already active with an administrative-role bypass; review who holds that role. Never reuse a release tag. |
| Environment | Name exactly `dockerhub-release`. **Enable Required reviewers** and select the releasing maintainer/team. Restrict deployment to selected **tags** `v*`, not branches. Disable protection-rule bypass where available. |
| Solo-maintainer approval | A maintainer may trigger and manually approve the deployment, as in dndsR. Leave Prevent self-review off for this case; enable it only when a second reviewer can approve. Approval is still explicit. |
| Environment secrets | Store `DOCKERHUB_USERNAME` and `DOCKERHUB_TOKEN` only in `dockerhub-release`; do not duplicate these names as repository/organization secrets accessible to checking jobs. Use a suitably scoped Docker Hub Read & Write token for `scrumpis/iwgc-circos`. |
| Repository variable | `DOCKERHUB_PUBLISH_ENABLED`; only the exact value `true` enables the publishing job. Leave absent/false through initial CI and configuration review; enable last. |
| Actions permissions | Allow the checkout/artifact actions; `GITHUB_TOKEN` needs only `contents: read`. The workflows do not request repository write access. |
| Docker Hub image | `scrumpis/iwgc-circos`; enable immutable version tags, for example the RE2 rule `^v[0-9]+\.[0-9]+\.[0-9]+$`. This is required for immutable publications: preflight HTTP lookups alone cannot prevent external publisher races. Verify availability/enforcement before enabling publication. |
| Other publishers | Disable competing Docker Hub autobuild rules or other automation that could push these version tags or `latest` outside the reviewed release process. |

An environment name in YAML alone does not require approval. Required reviewers
must exist in GitHub settings. GitHub documents this gate in
[deployment environments](https://docs.github.com/en/actions/reference/workflows-and-actions/deployments-and-environments).
Registry immutability is separate from Git tag protection; see
[Docker Hub immutable tags](https://docs.docker.com/docker-hub/repos/manage/hub-images/immutable-tags/).
If these controls cannot be enabled, leave publishing disabled and resolve that
limitation explicitly rather than claiming immutable/approved publication.

## Stable versioned publication

Only stable semantic tags `vX.Y.Z` are accepted, without leading zeroes,
prerelease suffixes or build metadata. The Actions tag glob only selects a run;
a strict shell check rejects other forms before building or reading credentials.
The current `0.4.0-rc.1` source is deliberately ineligible for publication.
Do not retag it as stable without updating/reviewing VERSION, skill metadata and
release documentation consistently.

`.github/workflows/publish.yaml` performs these steps in its single `publish`
job, after the configured environment approval:

1. Require `DOCKERHUB_PUBLISH_ENABLED == 'true'`, an exact VERSION/tag match and
   validated matching skill metadata. Verify the checked-out HEAD equals the
   event SHA; fetch current `main` explicitly and require that SHA to be its
   ancestor. This is a history check, not independent proof of review: enforce
   `protect_main` and authorize the merge as well.
2. Require Docker Hub's version-tag lookup to return 404; existing tags and
   lookup/network errors stop the job.
3. Build the exact tagged source with `--pull --no-cache`, linux/amd64, and
   version/revision/source labels. Save image inspection.
4. Run the full shared checks **inside that newly built image**. No skipped
   scientific checks, rendering exceptions or replacement build is permitted.
5. Recheck tag absence immediately before login; read the two environment
   secrets only in that push step, and push only `scrumpis/iwgc-circos:vX.Y.Z`.
   Use temporary Docker authentication and log out afterward. Registry
   immutability closes the race between lookup and push. Never move `latest`.
6. Upload the versioned skill tarball, logs, image inspection and
   `release-manifest.json` with commit, pullable image digest and skill SHA256.
   Artifact `release-evidence` retains these for 90 days; preserve them with
   the authorized release before expiry. No GitHub release is created by CI.

The image contains dependencies and the licensed font. Runtime scripts are
mounted from the matching versioned checkout. Consumers should retain/use the
manifest's immutable image digest. The first real publication can only be
verified during a separately authorized tag/publication run.

## First CI and initial release sequence

1. Freeze the candidate plus these bounded workflow/documentation corrections
   on `audit/release-candidate-0.4.0`. Preserve the prior scientific validation
   records; the automation audit supplements their historical source manifest.
2. Obtain authorization to push that branch. Its push runs both checking jobs
   without Docker Hub credentials or publication enabled. Open a PR to `main`
   when authorized; both checks run for the PR too.
3. Preserve actual CI evidence and resolve genuine failures narrowly. In
   particular, do not remove the SVG/native gap check or substitute an old image.
4. Add passing `fast` and `build-and-test` as required checks in `protect_main`.
   Review the qualified telomere/font limitations and independent scientific
   review status. A first CI run is not a release approval.
5. Before stable release, update VERSION, skill metadata and release notes to
   `0.4.0`, rerun required CI and authorize the PR merge. Configure and verify
   the mandatory release environment and Docker Hub settings above, enable the
   variable last, then separately authorize creation/push of `v0.4.0`.
6. Approve the deployment and retain its actual manifest, image digest and skill
   artifact. Publishing a GitHub release remains a separate authorized action.

## Introduce permanent dev after the initial merge

Create `dev` from the newly merged `origin/main`, after the first release merge
and from a clean checkout. Recommended commands, **not executed here**:

```bash
git fetch origin
git switch -c dev origin/main
git push -u origin dev
```

Install `protect_dev`; keep `main` as the stable default branch. Subsequent
`feature/*`, `fix/*` and `audit/*` branches start from `dev`, return through PRs
to `dev`, and releases use a reviewed `dev` → `main` PR. Prefer merge commits
between the permanent branches to preserve ancestry; after release, bring
`main` back into `dev` by fast-forward when possible, otherwise an ordinary
reviewed merge. Do not reset, rebase or force-push either permanent branch.
The initial `dev` branch needs no merge or history rewriting because it starts
at the already merged main commit. Creating/pushing it still needs authorization.

Git history remains intact. Historical Git objects contain the old Arial binary;
freeze its existing candidate deletion before creating distributable archives.
`.gitattributes` omits original render outputs from archives while preserving
those outputs and the original defect evidence in the repository.
