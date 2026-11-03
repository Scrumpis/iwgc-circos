# 0.4.0-rc.1 final validation review — NO-GO for release

Review date: 2026-11-02. Branch `audit/release-candidate-0.4.0`, base
`fe224c81c41e4aee5905e587ed01e0c1829ca0a3`; changes remain uncommitted at
review. **GO to freeze the candidate for CI review; NO-GO to release/publish.**
A clean candidate Docker build has not executed, and the old runtime has an
observed SVG/native rendering mismatch that the new image must resolve or
explicitly qualify. Nothing was pushed, merged, tagged or published.

## Candidate preservation and changes

Before follow-up edits, the full tracked binary diff, all 44 untracked files,
and per-file hashes were backed up outside this repository in a private scratch
workspace. The pre-review content manifest is preserved verbatim as
[pre-review-manifest.json](pre-review-manifest.json). The original
[15 observations](../baseline-reproduced.json), [initial capture](../baseline-probes.json),
and `tests/reproduce_baseline.py` are byte-identical to their incoming versions.
Original figure files under `iwgc_circos/tmp/` remain unchanged.

Reviewed the tracked diff and all new source/configuration/test/skill files,
fixtures, notices and retained evidence. Incoming work was not reset or
replaced. No new feature or scientific default was introduced in this follow-up.
The small follow-up changes are:

- Remove unused `fonts/ArialBold.ttf`; retain the unmodified licensed Liberation
  font and full OFL notice; add attribution. Update the SVG fallback family to
  Liberation Sans, and install/cache that exact font in the candidate image.
- Correct the Docker help example to the candidate's documented working
  directory/tag. Keep the scientific processing/CLI interfaces unchanged.
- Add bounded independent arithmetic and rendering-gap checks as release
  validation. The new gap check records a failure; it does not hide one.
- Automate clean full-image build/test and version-tag publication with the
  exact same tested image. No credentials or live repository settings changed.
- Restrict Docker build context to the Dockerfile and licensed font files;
  ignore local example run directories and omit historical render outputs from
  source archives, while preserving them in Git.

The complete reviewable candidate is `git diff --binary HEAD` plus the new
files listed in [candidate-files.json](candidate-files.json); hashes and scope
are in the current [source manifest](../source-manifest.json). No current-HEAD
archive should be published before freezing the deletion: historical HEAD
still contains Arial, as required by the no-history-rewrite policy.

Final hygiene completion inventories 87 existing candidate paths; the current
manifest hashes 86 entries and excludes itself. The three original renders are
preserved byte-for-byte in the working tree and excluded from source archives
by `.gitattributes`. The proposed archive was constructed in a temporary
repository, without staging or committing in this repository. Baseline evidence
and independently checked scientific source hashes still match their saved
records. No completed numerical/regression tests were restarted during this
final inventory/report completion.

## Executed evidence

| Check | Result / exit status | Scope |
|---|---|---|
| Existing automated suite in identified SIF | **62 passed**, no skips; 0 | Auxiliary regression, 12.469 s; [log](regression-tests.txt). Not a clean build. |
| Additional independent numerical check | **24 passed**; 0 | Literal windows, enumerated one-based base sets, 45-digit Decimal arithmetic, hand-constructed link unions; [record](independent-check.json). |
| Single-chromosome telomere controls | Expected missing-report failure 1; direct search and two-record control 0 | [Detailed evidence](telomere-check.json), below. |
| Native old/new-font Circos render + lettering | Both 0 | Same plot paths and seven letter positions; different glyph metrics; [record](font-rendering.json). |
| Host SVG conversion | Commands returned 0; **visual fidelity not accepted** | ImageMagick also misplaces/clips text; librsvg and ImageMagick expose colored gap mismatch. CLI success is not figure validation. |
| PNG gap probe | Native/native control 0; native/SVG comparison **1** | Native expected `[255,255,255]`; SVG `[239,59,44]` and `[253,141,60]`. New publishing check must pass. |
| Clean Docker build attempt | **127**, `docker: command not found` | No image ID/digest exists for this attempted build. No old image substituted. |
| Docker availability/rootless prerequisites | No client/daemon/socket/module; no subordinate UID/GID ranges | Read-only host checks; no host configuration changed. No connected remote execution session was available. |
| Workflow static validation | actionlint **1.7.12**, 0; YAML parsed | Official release checksum verified. Shellcheck/pyflakes integration disabled for this local lint. Not GitHub execution. |
| Bash syntax, portable skill checks | 0 | Existing CLI/help paths, metadata, version, portable links, Claude symlink and Arial absence. |
| Candidate distribution hygiene | **Passed**; [record](hygiene.json), [inventory](candidate-files.json) | 87 candidate paths; source archive excludes historical renders, Arial, generated run outputs and private paths. |

Implementer and scientific rechecker are the **same Codex session**. Expected
numerical values were constructed independently of program outputs, but that
does not constitute independent human/model scientific review. The earlier
fresh agent tested end-user usability, not scientific correctness. Its report
is preserved; no new automatic-discovery claim is made. A second fresh trial
was not required for the limited font/CI/documentation changes.

### Numerical scope

`tests/independent_check.py` does not import implementation code or use bedtools
as the oracle. It enumerates inclusive genomic base positions, compares literal
fixed/sliding BED windows, computes union coverage and GC from input bases,
then uses independent Decimal square-root/min-max arithmetic. It checks genes,
repeat overlap union, intact classes, missing/zero LTR ages and unweighted means,
constant tracks and a near-constant cancellation case. Link cases use explicit
expected spans for containment/overlap, reverse closed MUMmer endpoints, 11/13
columns, longest-N/order, minimum-length inclusion, filtering and PAF single-base
conversion/self-diagonal removal.

Ordinary coverage/GC tolerance is 1e-12; near-constant cancellation uses 1e-9
because a ~1e-6 range amplifies double roundoff; LTR output uses 1e-6 for its
six-decimal serialization and bedtools mean precision. Coordinates, window
membership and links compare exactly. Source/input hashes, actual/expected
values, commands and statuses are retained. Passing this finite designed set
supports these transformations, not annotation accuracy, orthology, age-model
validity, whole-genome performance or universal numerical correctness.

### Telomere mechanism and expected behavior

The 860 bp minimal discovery fixture has 60 `TTTAGGG` repeats, 20 intervening
bases and 60 `CCCTAAA` repeats. tidk 0.2.0 exits 0 but its candidate summary is
header-only. A 30 kb one-record control has 100 motifs at each end and separate
10 kb bins. Known-motif `tidk search` finds exactly 100 forward copies in the
first bin and 100 reverse copies in the last; discovery still nominates none.
The same sequence under two record IDs produces a candidate and four bands.

This is consistent with [tidk 0.2.0's source](https://github.com/tolkit/telomeric-identifier/blob/0.2.0/src/explore.rs):
`get_telomeric_repeat_estimates` accumulates matching pairs of candidate entries;
a single entry has no pairs. The location table records a leading tract while
the terminal-run handling also limits the entries. This source interpretation
explains the control outcome; it is not proof of the exact binary build source.
quarTeT treats the empty candidate table as exit 0 without a report. The
candidate wrapper correctly turns that missing report into exit 1.

Biologically the one-record control contains terminal motifs and would require
two marker bands if detection succeeded. Its current supported failure behavior
is **unavailable detection with nonzero status**, never biological absence.
This does not imply every single-chromosome genome fails. Do not duplicate
biological input to force a successful result. No dependency upgrade, fixed-motif
fallback, threshold/clade change, or new dating method was applied. Such a
scientific change requires separate review; the qualified limitation remains.

### Font licensing and rendering

The candidate tree contains only unmodified Liberation Sans Bold 2.1.5 as its
bundled font; embedded metadata names the Google/Red Hat copyrights and OFL 1.1.
Full notice and attribution ship in `fonts/`. Keep those with font distributions;
the font is not relicensed as MIT. See [font notice](../../../fonts/NOTICE.md).
Git history is preserved; a legacy SVG can reference Arial by name/path, but
no embedded Arial font bytes were found in the candidate figures.

The synthetic old/new-font comparison has identical plot paths and seven ring
letter coordinates. Three native SVG text anchors change: two tick labels by
about 0.5 px, and the long `scaffoldT1` label from (55.7,1523.5) to
(63.7,1523.3). Native PNG differs at 15,875/9,000,000 pixels (0.1764%); this is
text/glyph evidence, not proof of pixel equivalence. Original figures were
not regenerated from absent biological inputs.

A host ImageMagick conversion clipped/misplaced chromosome text. A separate
librsvg API render with the exact font registered in that process preserves
readable text, but both SVG backends paint two intended white gaps red/orange
where the native PNG is white. The gap points are (707,1500) and (808,1500) in
the fixed 3000 px six-track fixture, midway in the 0.02r gaps on the zero-TE
chromosome, outside the link radius. `tests/check_render_gaps.py` detects the
mismatch independently of file existence/dimensions. It is included in the
new-image release check. No upstream renderer or plot defaults were patched.

## Clean build and automation status

Exact attempted local command (exit 127):

```bash
docker build --pull --no-cache --platform linux/amd64 -t iwgc-circos:clean-review .
```

The host has no Docker endpoint and cannot meet documented rootless
prerequisites without host changes. None were made. The user's preferred
GitHub route is prepared following the read-only dndsR workflow reference
`5983405287ed8ecc62feeaedb15f176ef09b9790`:

- automatic uncached full builds/tests on relevant updates and PRs;
- reusable full-image checks, including human example, PNG/SVG gap fidelity and
  dependency inventory;
- enabled version-tag build/test/push, tested-image identity, manifest and skill
  artifact; no `latest` movement or existing version overwrite;
- environment credentials and one-time activation documented in the
  [release process](../../release-process.md).

A local image ID is a configuration digest, not automatically a pullable
registry digest. Neither exists for the failed local build. CI will save
`docker image inspect`; publication records the registry digest after push.
No value is invented from the existing SIF hash. Workflow lint exit 0 does not
resolve the build blocker. GitHub Actions execution remains **pending an
explicitly authorized branch push**; publishing also needs setup and a separate
authorized merge/tag/publication sequence.

## Outstanding gates and freeze recommendation

1. Run the prepared clean-container workflow on the frozen source. Preserve the
   build log/status, image ID, full suite/independent/example outcomes and
   dependency inventory. A clean full build is required, not an old Hub/SIF run.
2. Require the new image's SVG/native gap check to pass. If it fails, request
   review of a narrowly reproduced rendering correction; do not waive the
   failure by changing its expected colors or claiming structural smoke suffices.
3. Maintainer must accept the qualified telomere limitation and review the
   font/layout differences. Independent scientific review remains outstanding.
4. Configure/enable the publishing environment only when ready. No production
   publication is authorized by a static workflow definition or this report.

One atomic freeze commit is appropriate for this incoming uncommitted candidate;
it keeps interfaces, metadata, skill, tests and evidence aligned. These are
**recommended commands, not executed actions**:

```bash
git switch audit/release-candidate-0.4.0
git add -A -- Dockerfile README.md VERSION AGENTS.md CLAUDE.md \
  .gitignore .gitattributes .dockerignore .github/workflows \
  .claude/skills/iwgc-circos iwgc_circos_tracks.sh create_configs.sh add_legend.sh \
  iwgc_circos_template.config iwgc_circos/iwgc_circos.conf \
  preprocess tests tools skills fonts docs
git diff --cached --check
git diff --cached --name-status
git commit -m "Prepare 0.4.0-rc.1 for clean-container validation"
```

After an authorized push and completed CI, a second evidence-only commit should
be `Record clean container validation for 0.4.0-rc.1`, containing the actual CI
run URL, frozen commit, image ID/dependency inventory, gap-check outcome and
updated GO/NO-GO decision under `docs/validation/`. Do not create that commit
with hypothetical results. Stable 0.4.0 version/doc updates and any subsequent
merge/tag/publication require their own approved release step.
