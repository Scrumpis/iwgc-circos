# First CI failures and dependency corrections

The maintainer pushed candidate commit
`a5940bcddcc31150371e035fcb173b5dd2f5c3c3` on
`audit/release-candidate-0.4.0`. Both first workflows executed and failed; the
previous records describing unexecuted CI are retained as historical evidence.
This follow-up changes dependencies and failure diagnostics, not scientific
calculations, defaults, rendering expectations or release protections.

## Observed failures

| Run / job | Actual outcome | Evidence and cause |
|---|---|---|
| [Validate, 41546692988](https://github.com/Scrumpis/iwgc-circos/actions/runs/41546692988), `fast` / 125623877650 | Workflow guards and skill validator passed; 62 unit tests discovered, 60 passed and 2 optional detector tests skipped; rendered smoke failed | Circos printed `missing List::MoreUtils` and exited **0** without PNG/SVG. The harness then raised FileNotFoundError for the SVG. |
| [Clean container, 41546693013](https://github.com/Scrumpis/iwgc-circos/actions/runs/41546693013), `build-and-test` / 125623878118 | Clean Docker build failed with **1** during RIdeogram installation; container checks were skipped | Source dependencies could not compile: `make: not found`, missing `pkg-config`, XML/rsvg development requirements; downstream R namespace errors cascaded from these failures. The final `requireNamespace` assertion correctly stopped the build. |

The smoke artifact `13188169169` contains command statuses and captured Circos
stdout, confirming that this was a missing Perl module rather than a successful
render. Its ZIP SHA256 matches GitHub's upload log:
`324168e455000c6b0635d407ff232ea4d9c6cfb15556c65862c67a6ebcf4b59d`.
It was downloaded through the GitHub connector and inspected without modifying
the original artifact. Relevant command records are preserved in
[ci-first-run.json](ci-first-run.json), with runner-specific checkout paths
normalized to `<repo>`. Generated figures, ZIPs, access links and full runner
logs are not added to the software distribution.

Both environments installed Ubuntu Circos `0.69.9+dfsg-2`; installing `circos`
alone did not supply the required List::MoreUtils module in these actual runs.
Neither workflow reached a successful qualified end-to-end result. No candidate
image ID or registry digest exists for this failed full build. Publication was
not triggered or attempted.

## Corrections and scope

- Explicitly install `liblist-moreutils-perl` in fast CI and the Docker core
  stage. Its absence was reproduced by the independently provisioned CI run
  and confirmed in the saved artifact's subprocess stdout.
- Add `r-base-dev`, `libxml2-dev` and `librsvg2-dev` to the full stage, retaining
  the existing CRAN install and mandatory RIdeogram load check. Ubuntu's
  [r-base-dev dependencies](https://packages.ubuntu.com/jammy/r-base-dev)
  include make/compiler tooling, pkg-config and PNG/JPEG/ICU headers. The extra
  XML/rsvg headers address the corresponding source dependency requirements.
  This retains quarTeT plotting compatibility; no detector or R version/default
  is deliberately changed. Transitive package versions still float as previously
  disclosed.
- Require nonempty PNG/SVG output immediately after Circos in the smoke helper.
  A misleading zero status now raises an error including the original stdout
  and stderr, before XML parsing. All existing SVG geometry, native PNG size,
  legend checks and the full-image SVG/native gap check remain in place.
- Add one focused regression using a real subprocess that exits zero and emits
  the missing-module diagnostic without producing output. It verifies that the
  helper rejects the false success while retaining the diagnostic and status.

The missing build dependencies were my recipe omission. The first real clean
CI build exposed them; the earlier SIF results did not establish Docker build
correctness. No failed scientific expectation was relaxed to obtain a pass.

## Local verification after correction

All listed commands exited **0**:

```bash
python3 -m unittest discover -s tests -p test_smoke.py -v
/usr/bin/python3 tests/check_workflows.py
actionlint -shellcheck='' -pyflakes=''
singularity exec --cleanenv "$IWGC_SIF" python3 tests/smoke.py --output /tmp/iwgc-ci-dependency-fix-smoke
git diff --check
```

The new regression passed; all 17 offline workflow guard cases passed;
actionlint 1.7.12 passed with shellcheck/pyflakes integrations disabled locally.
`IWGC_SIF` is the previously identified existing image, SHA256
`5067065f346b5233dc708bf032053a3287a38b7150ca12ca0900cb8419db7cb7`.
The real Circos auxiliary smoke passed with six tracks, two chromosomes,
rendered link and legend. This image lacks the SVG rasterizer: **no new
SVG-to-PNG fidelity result or clean Docker validation is claimed**. The full
unit suite was not rerun locally for dependency/test-harness-only changes; the
next CI suite will discover 63 tests, including the new regression.

Docker remains unavailable locally. The dependency additions are not yet
execution-validated in a newly built image. A new push must run both workflows
on the correction; rerunning the old SHA would repeat the old recipe. Preserve
its actual build logs, image identity and test outcomes. The known rendering
fidelity gate and all remaining scientific/maintainer release gates still
apply: **NO-GO for release pending successful corrected CI and review**.
No push, tag, registry publication, upstream edit or history rewrite was made
by this diagnostic follow-up.
