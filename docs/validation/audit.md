# Release-candidate audit

Date: 2026-11-02. Candidate: **0.4.0-rc.1**, unpublished. Work branch:
`audit/release-candidate-0.4.0`. Starting tree was clean at
`fe224c81c41e4aee5905e587ed01e0c1829ca0a3` (`main` and local `origin/main`).
Latest local release tag inspected: v0.3.0,
`7e327f66903d396a07cbde65c0a71eff86451db0`. HEAD already differs from that tag
in square-root transforms, docs, legend offset and file modes. No remote release
or registry state is certified by those local refs.

All 23 original tracked files were inventoried: three main shell scripts,
template, all nine config fragments/master, two preprocessing examples,
Dockerfile, README, MIT license, Arial font, placeholder and three representative
figure files. Both SVGs were parsed, PNG/SVG structure inspected, and the
original and synthetic PNG visually reviewed. Original figure inputs are not
bundled, so the existing biological figure was **not scientifically reproduced**.
No original results were changed or deleted. Original files under
`iwgc_circos/tmp/` remain byte-identical to baseline. The initial audit retained
the historical font; the follow-up request authorized removing that unused
binary. It is now deleted in the candidate diff. Current review/gates are in
[GO-NO-GO.md](final-review/GO-NO-GO.md); initial captures remain preserved.

## Foundation and reference use

Read agentG 0.1.0 at `ca1632977b2399a51bffb56f9f1b1ca243cdcc6c`:
constitution, all four roles, initialization/adoption skill, independent
validation skill, architecture, validation evidence and evidence-level lesson.
Used its distinction between execution, regression, independent numerical
expectation and independent review. Did not run the initializer, copy its
instruction bundle, alter agentG, or introduce orchestration.

Read dndsR's local development checkout at
`4e2161d2e22e78213e385955be4a85916395060b` (local `origin/dev` was
`5983405287ed8ecc62feeaedb15f176ef09b9790`). That checkout had unrelated local
changes, which were preserved. Inspected AGENTS/CLAUDE guidance, portable skill
and verification reference, numerical/GFF tests, test-agent workflow, validation
criteria, installation and getting-started docs. These are architectural
examples, not evidence that iwgc-circos works. No dN/dS thresholds or domain
requirements were imported. Liberation Sans Bold and its license were copied
as separately licensed font assets, with hashes in the source manifest.

Main implementation, probes and scientific expectations were authored/executed
by the same Codex session. A separate fresh agent tested end-user usability;
it did not perform independent scientific review. No human scientific sign-off
or blinded comparative review is claimed.

## Findings and resolutions

Severity reflects impact: **high** can silently alter scientific inputs/results
or execute unintended commands; **medium** materially affects figure completeness,
reliability or reproducibility; **low** affects documentation/usability.
Confirmed original defects below have an executed discriminator, except the
Dockerfile grammar issue, whose duplicated `FROM`/pasted prompt is visible in
the original Git blob. That build was not executed and is explicitly unverified.

| ID | Severity / classification | Location, evidence and impact | Candidate status |
|---|---|---|---|
| S1 | High, confirmed error | Track GFF3-to-BED conversion passed starts through unchanged. Coordinate-only probe returns `0.816497,1,0.577351` rather than `0.8660254,1,0.7071068` for Chr1. Single-base/boundary coverage is also affected. | Approved correction: subtract one for gene/repeat/intact/LTR starts; independent fixture tests. |
| S2 | High, confirmed error | Gene `grep gene` selected mRNA/exon/attribute text. The fixture's mRNA changes the third display value from intended `0.7071068` to `1`. | Exact column-3 `gene` selection, approved and tested. |
| S3 | High, confirmed error | MUMmer reverse swap maps 10..7 to 10..10, shortening to zero; one-based MUMmer spans and PAF helper starts were not consistently converted. | Approved corrected swap/conversions. Six-column contract now explicitly zero-based half-open; strand still discarded. |
| S4 | High, confirmed error | Merge probe with disjoint second intervals merges them and drops their lower bound. | Approved both-side overlap/touch requirement and min/max bounds; greedy sorted merging retained. |
| S5 | Medium, confirmed error | Default ascending `-top-n 2` chooses totals 2,4 instead of 4,6; filtering happened after truncation. | Approved longest-N eligible merged links, followed by requested draw order; filter against retained IDs first. |
| S6 | High, confirmed numerical error | Two nearly full 10 Mb windows collapse to `0,0` because intermediate six-digit output loses differences. | Separately approved precision correction; 17-digit intermediates, integer coverage/GC counts. Regression expects `0,1`. |
| E1 | High, confirmed security error | FASTA filename `genome$(touch${IFS}IWGC_MARKER).fa` reaches quarTeT shell interpolation; frozen baseline creates a harmless marker in its private scratch dir and returns 0 after detector failure. | Fixed relative `input.fa` / `iwgc` dependency arguments in a unique private subdirectory; no user path is interpolated by the detector. Regression checks the boundary. |
| E2 | High, confirmed error | `-window 0` produces empty GC with exit 0 after bedtools error. Empty LTR mapping also returns 0 without rows in bedtools 2.30.0. Missing `-gene` value loops until timeout. | Strict shell status/pipe handling, argument validation, explicit empty-LTR mapping, complete numeric row-count checks. Errors remain nonzero. Partial outputs can remain and must not be rendered. |
| E3 | Medium, confirmed error | Original intact splitting can reuse a pre-existing class intermediate in the output directory; whitespace splits file lists. | Private per-run class files, quoted arrays, streaming validation plus external sort. Final old output files are deliberately not deleted; use fresh run directories. |
| E4 | Medium, confirmed presentation error | With nine leading numeric tracks, baseline keeps link radius 0.305r instead of 0.515r (innermost ring 0.52r). Similar gap for 10–12 tracks. | Correct innermost-radius fallback; layout invariants tested at all 0–13 counts. |
| E5 | Medium, confirmed ambiguity risk | Original glob resolver always selects lexicographically first match, possibly silently choosing a different genome. | Ambiguous nonempty matches fail. Empty tracks/links omitted. Custom outdir gets missing include files. |
| E6 | Medium, build defect by inspection; build unverified | Dockerfile contains two recipes and an un-commented shell prompt; uses disabled TLS validation, unpinned latest installer/quarTeT/CPAN, HTTP CRAN and lowercase `rideogram`. | Replaced duplicate recipe with core/full dependency stages; signed Ubuntu Circos, TLS, pinned/checksummed installer, pinned quarTeT/tidk, correct `RIdeogram` name. Clean Docker build remains required. |
| E7 | Medium, confirmed presentation error | `--theta-deg` overrides y but still uses fixed gap midpoint x. | Override now uses radial x too; actual SVG 0° geometry tested. Malformed/repeated SVG and identical input/output refused; explicit PNG failure is nonzero. |
| E8 | Medium, confirmed silent omission | A valid `scaffold_T1` chromosome is removed from display and tick statistics by a redundant name heuristic, despite its `chr` record type. Frozen baseline writes only `Chr1` for a two-chromosome input. | Separately approved removal of name exclusion; all valid `chr` records retained, bands excluded by type. Regression and real render use the unusual ID. |
| L1 | High release gate, unresolved redistribution evidence | Arial metadata identifies Monotype copyright and product-license-dependent use; no redistribution grant is present. Repo MIT license does not establish font rights. | Liberation Sans Bold/OFL is the default. Follow-up validation removes the unused Arial binary from the candidate, retains font notice/attribution and preserves original figures/history. Current font/render evidence is in the final review. |
| M1 | Medium, confirmed dependency limitation | quarTeT 1.2.5 / tidk 0.2.0 fails to nominate a motif for a designed single-chromosome positive sequence; also reproduced on a longer one-record control. Two-chromosome positive control succeeds. | No upstream detector change. Wrapper fails when report is absent. Separate tests cover positive detection and error propagation; do not equate the latter with sensitivity validation. |
| M2 | Methodological concern, retained | Relative min-max scales, zero-variance = 0, overlap-dependent LTR means and orientation-free greedy links can be misread as absolute density, ages or orthology. | Definitions/limits documented; original formulas and orientation-free default retained. No biological estimator/orthology validation claimed. |
| D1 | Low, confirmed docs inconsistencies | README/help disagree on repeat cube vs actual square-root; gap spacing and ideogram default path wrong; some container commands refer to old names. | Corrected prose and added tested human guide, exact input reference and portable skill. Historical image tag identity remains explicitly unverified. |

Review and approval evidence for scientific/default changes is in
[proposed-behavior-changes.md](proposed-behavior-changes.md). Before/after values
and the exact test commands are in [verification.md](verification.md).
[Baseline observations](baseline-reproduced.json) come from
`tests/reproduce_baseline.py` running original Git objects. The initial
[probe capture](baseline-probes.json) is also preserved.

## Remaining limits and release gates

1. Arial removal is implemented in the candidate; freeze that deletion before
   making release archives. Retain the Liberation notice. Font metrics and
   SVG/native rendering differences require the review recorded below.
2. Execute clean candidate **core and full Docker builds**, tests and rendering;
   then record dependency inventory, immutable registry digest and pull route.
   Docker/BuildKit were unavailable. SIF execution proves only that existing
   image, not candidate build correctness. CI jobs have not run on GitHub. The follow-up adds automatic full-image
   validation and tested-version publishing; execution awaits an authorized push.
3. Keep the single-chromosome telomere limitation explicit; either qualify
   supported use or separately validate an approved dependency upgrade before
   claiming that route is generally validated. Plant-only assumptions and
   genome-scale telomere sensitivity remain unvalidated.
4. The follow-up found an SVG/native rendering mismatch in two histogram gaps
   under the old 0.69-8 runtime. The new image must pass the added gap check
   before publication. No plot default was changed to conceal it.
5. Obtain independent maintainer/scientific review before calling this a
   professionally validated release. Hand-derived oracle agreement is bounded
   evidence, not a certification of all annotations, EDTA ages or genome reports.

Other untested conditions: macOS/Windows/ARM, Apptainer, gzip FASTA, user-supplied
custom template geometries, arbitrary edited/translated SVG, labels/highlights
combined with automatic ring lettering, very large annotations under measured
resource budgets, and labeled-PNG rendering in this existing image. The input
converter streams features and delegates large sorts to system `sort`; no
whole-annotation Python list is retained. This is engineering design, not a
benchmark. Existing indexes are trusted apart from mtime/schema checks; a
plausible but wrong `.fai` needs separate user verification.

LTR means still use bedtools map's arithmetic and numeric output precision;
near-equal ages at extreme magnitudes were not independently validated.
Coverage scaling and final LTR output do not carry statistical uncertainty.

## Reusable lessons proposed for agentG later

- Test the status **and** expected data inventory when external scientific
  tools can return success after printing an error. Reproduced with bedtools
  empty mapping and quarTeT missing reports; do not universalize it to unrelated tools.
- Include a one-base boundary and a near-constant numerical probe alongside a
  normal example; normalization can conceal coordinate/precision errors.
- A portable skill should retain its required references when copied alone;
  repository-only relative links are not a standalone installation strategy.
- Freeze source bytes for agent trials and distinguish explicit skill loading
  from automatic discovery. A fresh agent found documentation issues; it did
  not independently certify scientific methods.

These are proposals only. No agentG/dndsR repository, catalog, policy, hooks,
remote refs, releases or registries were modified.
