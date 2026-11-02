# iwgc-circos 0.4.0-rc.1 — candidate, unpublished

Recommend **0.4.0** after the release gates below. This pre-1.0 minor release
makes coordinate contracts and stricter failure behavior explicit and adds a
portable skill. No release/tag has been created.

## Changes

- Correct GFF3/LTR one-based closed to BED zero-based half-open conversion and
  select exact `gene` features. Regenerate annotation-derived tracks.
- Correct MUMmer reverse intervals and coordinate conversion; preserve PAF
  starts in the self-alignment helper. Existing six-column inputs must be
  zero-based half-open; do not silently reuse old helper outputs.
- Correct merged link bounds, chromosome filtering and longest-N selection;
  retain orientation-free links and existing draw-order options.
- Retain numeric precision through coverage/GC transforms. Two near-full
  windows can now retain their intended 0/1 contrast instead of becoming all zero.
- Keep square-root gene/total repeat and min-max transforms, constant-track zero,
  LTR missing zero / constant ages 1.5, unweighted mean ages and sliding tail cutoff.
- Fail on malformed/unknown/out-of-range intervals, bad options, stale/duplicate
  indexes, dependency errors and incomplete numeric output. Empty links and
  empty annotation coverage remain valid. Runs can leave partial files on error;
  preserve logs and retry in a fresh output directory.
- Protect the quarTeT filename boundary with fixed private names. Clamp marker
  widths to chromosome lengths. No change to plant detection threshold or clade.
- Retain all valid chromosome records, including names containing `_T1`; distinguish
  telomere bands by record type rather than excluding chromosome names.
- Correct dynamic link radius for partial layouts; reject ambiguous globs and
  malformed karyotypes; omit empty tracks. Preserve the three main shell interfaces.
- Use licensed Liberation Sans Bold by default, embedding the chosen font in
  labeled SVG. Text metrics differ from historical figures; inspect captions/labels.
- Add small synthetic known-answer tests, real Circos smoke, CI checks, human
  getting-started guide, project instructions, and portable versioned Agent Skill.
- Prepare automatic clean full-image validation and version-tag Docker Hub
  publishing. The exact tested image is pushed only after all checks pass; CI
  execution and registry publication remain pending authorization/setup.
- Repair container recipe and dependency provenance; core/full stages preserve
  an optional telomere workflow without imposing its environment on core use.

Scientific/default corrections were reviewed and authorized in the audit
conversation; see [the proposal/evidence](validation/proposed-behavior-changes.md).
See [verification](validation/verification.md) for actual numerical results and
executed commands. A rendered plot is not scientific certification.

## Gates before release

- The unused historical Arial binary is removed from the candidate; keep the
  licensed Liberation font, notice and attribution in distributions.
- Perform a clean full Docker build (including its core stage) and CI runs; publish no image until its
  identity and tests are recorded. No reproducible-registry digest is claimed here.
- Review the known single-chromosome tidk 0.2.0 limitation and the qualified
  telomere workflow. Current wrapper fails on missing reports.
- Resolve or qualify the old 0.69-8 SVG/native histogram-gap mismatch in the
  clean candidate image; the publishing check now detects it.
- Review figures/font metrics and commission independent scientific/maintainer
  review. Automatic skill discovery and other container engines/platforms are
  not established by the explicit fresh-agent trial.

All existing representative results were preserved. No push, merge, tag,
release publication or upstream changes were made.

Current decision and exact proposed freeze commits: [GO/NO-GO report](validation/final-review/GO-NO-GO.md).
