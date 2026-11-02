# Initial verification record

Date: 2026-11-02. Source baseline `fe224c81c41e4aee5905e587ed01e0c1829ca0a3`,
local branch `audit/release-candidate-0.4.0`, candidate `0.4.0-rc.1`.
This records the initial pass. The [pre-review manifest](final-review/pre-review-manifest.json)
identifies its runtime/test/skill bytes; the [follow-up review](final-review/GO-NO-GO.md)
records font removal, additional numerical checks and the new rendering concern.
The current candidate remains uncommitted and its bytes are in
[source-manifest.json](source-manifest.json). No push, merge, tag or release.

The same Codex session authored changes, derived expected values and ran tests.
The numerical oracle uses explicit base counts and interval arithmetic outside
the implementation, but the evaluator is not an independent scientific reviewer.
A separate fresh agent performed the explicitly loaded skill trial described
below. This is bounded evidence, not certification of biological correctness.

## Environment and exact commands

Linux amd64 HPC, Singularity CE 4.1.2-jammy. Existing dependency SIF SHA256:
`5067065f346b5233dc708bf032053a3287a38b7150ca12ca0900cb8419db7cb7`.
Tools: Python 3.13.2, GNU awk 5.1.0, samtools 1.13, bedtools 2.30.0,
Circos 0.69-8, quarTeT 1.2.5, tidk 0.2.0. This is not a registry digest or a
candidate Dockerfile build. See [container assessment](container.md).

From the checkout, set `IWGC_SIF` to that existing image's absolute path:

```bash
singularity inspect --json "$IWGC_SIF"
sha256sum "$IWGC_SIF"
bash -n iwgc_circos_tracks.sh create_configs.sh add_legend.sh preprocess/minimap2_self_alignment.sh
python3 tools/validate_skill.py
singularity exec --cleanenv "$IWGC_SIF" python3 -m unittest discover -s tests -v
singularity exec --cleanenv "$IWGC_SIF" python3 tests/smoke.py --output /tmp/iwgc-smoke-approved-candidate
singularity exec --cleanenv "$IWGC_SIF" python3 tests/reproduce_baseline.py --output /tmp/iwgc-baseline-final-capture
git diff --check
git diff --exit-code HEAD -- iwgc_circos/tmp
test ! -e fonts/ArialBold.ttf
```

These temporary output directories were fresh during execution; choose new
paths when rerunning. All container calls used the existing local image, no
pull. Container execution needed scoped permission after restricted startup
failed with missing FUSE/setuid extraction. The failing startup is recorded
in the fresh-agent trial. No runtime security controls were disabled.

Also executed the available skill-creator validator:

```bash
/usr/bin/python3 "$SKILL_CREATOR/scripts/quick_validate.py" skills/iwgc-circos
```

Here `SKILL_CREATOR` refers to the installed system skill-creator folder, not a
repository dependency. `/usr/bin/python3` had PyYAML 5.4.1; the repository's
validator and tests need only the Python standard library. Workflow YAML was
parsed with that PyYAML and its three-job inventory checked. This is not a
GitHub Actions run or validation of every GitHub expression.

## Final outcomes

| Check | Observed result | Claim supported |
|---|---|---|
| Unit/integration suite | **62 passed**, no skips, 11.895 s in final captured run | Known-answer transformations, input errors, optional detector controls and config invariants in this environment. [Full output](tests.txt). |
| Real Circos smoke | Passed: six numeric tracks, two chromosomes (including `scaffold_T1`), rendered link, SVG and 3000×3000 PNG, seven ring letters | Actual four-stage execution and meaningful SVG geometry/structure; [commands and logs](smoke-commands.json). |
| Legend angle override | Actual 0° SVG has common y and right-side radial x | Corrects previous x/y disagreement. |
| Labeled PNG in existing image | **Not exercised**: converter unavailable; SVG still succeeds | No labeled-PNG claim; explicit requested conversion failure is now nonzero. |
| Frozen baseline reproducer | 15 observations captured, including incorrect scientific values, timeout and harmless shell marker | Independently reproduced original defects; [observations](baseline-reproduced.json). It reports failures rather than treating them as candidate passes. |
| Skill metadata/portable links/version/symlink | Passed repository validator and skill-creator validator | File format, local reference portability and link target; not automatic discovery. |
| Bash syntax / whitespace | Passed | Narrow engineering checks, not semantic correctness. |
| Existing figures/font integrity at initial pass | Original figures and old font matched baseline then | Figures still match; the later removal of the unused font is an authorized follow-up change. |
| Clean Docker build / Docker execution / GitHub Actions | **Not run** | Docker unavailable; CI definitions are prospective. |

Final source hashes and synthetic fixture hashes are in the manifest. Smoke
artifacts remain under `/tmp/iwgc-smoke-approved-candidate` in this workspace;
command/exit-status records are stored in this repository. Large generated
figures are not duplicated into version control. Original figure files were
parsed and visually inspected, not regenerated from their absent biological
inputs.

## Independent numerical expectations

Tolerance for normalized tracks is absolute **3×10⁻⁶**, preselected to cover
historical six-digit output and final six-decimal LTR serialization, not to
mask a change in window/base inclusion. Interval coordinates, counts, class
membership, selected links and error/nonerror status are compared exactly.
Near-constant 0/1 regressions have exact endpoint expectations. No output from
iwgc-circos was copied as a numerical expected fixture.

| Synthetic case | Oracle / expectation | Actual result |
|---|---|---|
| Chr1 genes 2..4, 5..8, 10..11; two 12 bp chromosomes; window 4 | Covered bases 3,4,2,0,0,0; square-root fractions and min/max 0/1 give `0.8660254038,1,0.7071067812,0,0,0`. mRNA attributes must not add coverage. | Agrees within tolerance. Baseline coordinate-only gave `0.816497,1,0.577351,0,0,0`; attribute trap also altered the third value. |
| Overlapping repeats 2..4 and 3..4 | Union is 3 bases, not 5; remaining intervals give same fractions above | Agrees; independently enumerated BED fixture also gives raw `.75,1,.5`. |
| Nonzero minimum coverage | Transform `.25,.5,1` then subtract `.5` and divide `.5` | Middle display `(sqrt(.5)-.5)/.5 = 0.4142135624`; agrees. |
| Boundary and one-base genes | GFF3 4..4 and 5..5 occupy adjacent 4 bp windows once each | Exactly two positive windows; agrees. |
| Fixed window tail | Length 10/window 4 => [0,4), [4,8), [8,10) | Exact coordinates and length-based GC denominator agree. |
| Sliding window tail | Length 9/window 4/step .5 => [0,4), [2,6), [4,8), [6,9); drop [8,9) | Exact coordinates agree. Zero/rounded-zero step and no surviving windows fail. |
| GC ambiguity and lowercase | `gcNN`, `ATNN`, `GCGC`, `AAAA` => .5,0,1,0 | Raw counts/fractions and normalized values agree. Ns remain in denominator. |
| Constant full coverage / constant GC | min=max => display 0 | All zeros, no NaN; agrees. This does not imply biological absence. |
| Near-full 10 Mb windows | Fractions 9,999,999/10,000,000 and 1 => normalized 0,1 | Exactly 0,1 after precision fix; baseline 0,0. GC small-difference count case also retains 0,1. |
| LTR means 0,10,20,missing | Invert present min-max to 2,1.5,1; missing 0 | Exactly 2,1.5,1,0,0,0. Reverse interval and zero age handled. |
| LTR overlap weighting | Ages 0 and 20 both overlap window, with different overlap lengths => arithmetic mean 10 | Raw mean 10, not length-weighted 4; scaled result agrees. |
| Identical/empty LTR ages | Present all-zero age =>1.5; missing=>0. Empty file=>all missing | Agrees; bedtools empty-map version quirk handled explicitly. |
| MUMmer reverse query | 1..4 and 10..7 => [0,4), [6,10) | Exact six columns agree for 11 and 13-column inputs. |
| Link merging and selection | Both-side overlap/touch, extend minima/maxima; disjoint second side stays separate; totals 2,4,6 with top2 =>4,6 | Exact spans, threshold inclusion, asc/dsc and pre-selection chromosome filtering agree. |
| PAF helper | PAF starts 0 and 6 stay 0 and 6; exact diagonal excluded | Exact columns agree using a fixed synthetic aligner output. Real minimap2 sensitivity is not tested by this stub. |
| Telomere markers | Stubbed left/right status => bounded end markers; real two-chromosome 30 kb fixture with >50 plant motifs at each end => four markers of width 300 | Expected marker coordinates agree. Single-chromosome detector case fails and is reported; this is a dependency limitation, not a validated absence. |
| Layout | 0–13 track counts: ordered positive nonoverlapping radii; compact 1–8 inner radius≈.5; links inner-.005 | All invariants agree, including 9–12 track fallback. No whole-config snapshots used. |
| Unusual chromosome name | Valid `scaffold_T1` is a chromosome, `band` rows are not | Both IDs in config/tick counts; real smoke draws link to that chromosome. |

## Failures, recovery and limits

The first 49-test run failed five checks: empty LTR output and link-radius
invariants for 9–12 tracks. Both were reproduced independently on baseline and
fixed. The next 49-test run passed. The expanded 54-test run exposed a failed
single-chromosome telomere positive expectation; a later 57-test run still
reproduced it. These failures were retained in
[initial-test-failures.txt](initial-test-failures.txt) and
[detector-test-failure.txt](detector-test-failure.txt), not overwritten as passes.

The detector investigation tested one-record motifs at 100 and 1000 copies and
on longer inputs; neither nominated a motif in tidk 0.2.0. A two-record positive
control succeeded. Upstream source inspection suggests candidate aggregation
and terminal-run handling limitations, but the entire detector was not audited.
The failed case is now an explicit **error-propagation regression**, with a
separate two-record sensitivity positive control; changing the expectation does
not establish that the original single-record biological claim was correct.
No upstream dependency was modified or silently upgraded.

An initial skill-creator validation rejected the optional `compatibility`
frontmatter field, despite the current [Agent Skills specification](https://agentskills.io/specification)
allowing it. Requirements were moved into the skill body to work with the
available older validator; both validators then passed. This did not remove
runtime requirements or establish automatic discovery.

The available environment lacked Docker and an SVG converter. Shellcheck was
not available for this local execution; Bash syntax plus runtime/negative tests
were performed. No tests were run on macOS, Windows, ARM or Apptainer. Annotation
producer accuracy, age estimation/model assumptions, orthology and full-scale
performance remain outside these small known-answer cases.

## Fresh-agent usability

A separate agent with no conversation history was given only the software
checkout, an explicit skill path, synthetic FASTA/GFF3 locations and existing
image path. It read the skill/references, human getting-started guide, fixtures,
help and metadata; it did not read source implementation, audit conclusions or
numerical expected tests. It wrote to a new private workspace with the source
mounted read-only.

It completed genes+GC tracks with 4 bp windows, configuration, Circos SVG/PNG
and three ring letters. The parent independently re-read the resulting tracks
and XML against enumerated expected base counts. See
[trial report](skill-trial.md), [script/image hashes](skill-trial-identities.sha256)
and [exact trial commands](skill-trial-commands.txt).

The trial initially had a restricted container startup failure, then succeeded
with scoped execution permission. Its prerequisite run was **53/54**, with the
then-unresolved optional detector failure; it disclosed that limitation and
continued only the unrelated genes/GC request after parent confirmation. The
trial was on an earlier candidate, identified by its recorded hashes, not an
unchanged final code snapshot. Final production/tests were separately verified
in the 62-test and real-render runs above.

The agent identified ambiguous labeled-output paths, unused-optional-dependency
failure policy and center-fallback guidance. Docs/skill now provide the full
path, distinguish core failures from documented unused optional failures, and
explain root-center fallback on symmetric ideograms. Automatic skill discovery,
Claude runtime loading of the symlink and cross-client discovery were **not**
verified. Explicit loading is the tested path; copying the self-contained skill
folder is structurally supported, not a runtime discovery claim.
