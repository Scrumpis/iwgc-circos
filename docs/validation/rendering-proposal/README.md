# Review proposal: zero-width Circos SVG fill correction

Status update: the maintainer approved applying this one-line correction after
review of the evidence below. Its canonical copy is now under `patches/` and
the Docker/fast-CI definitions apply it. The following records the pre-approval
proposal and temporary experiment; it does not claim the corrected CI has run.

After independently reproducing the Ubuntu include failure, the corrected
probe exposed a second Ubuntu behavior: its global ideogram/tick files silently
shadow unqualified project includes. Using those defaults shrank chromosomes
to narrow arcs and made the gap pixels pass vacuously. The candidate now uses
explicit project includes and the smoke checks the independently calculated
ideogram radius, `1500 * 0.825 = 1237.5` pixels, so that substitution is rejected.

With the intended configuration restored, the original documented SVG/native
gap mismatch persists. The exact Ubuntu package's `Circos/SVG.pm` explains it:
`draw_slice` has a special branch for `radius_from == radius_to`. Its comment
says a zero-width slice cannot have fill, but the emitted arc retains the fill
style. SVG implicitly closes the arc with a chord, painting a circular segment.
The native PNG treats this as a polyline and does not fill it. These are false
colored areas for zero-height histogram bins, not real TE coverage.

Proposed resolution: the supplied [one-line patch](circos-zero-width-svg.patch)
sets `fill:none` **only in that equal-radius SVG branch**. It preserves arc
geometry, stroke styling, positive-width slice fill, input tracks, all numerical
calculations and the native PNG. If approved, apply it to the supported image
and fast-CI Circos package, record it as a local dependency correction, and
document the requirement for equivalent manual installations. Do not change
expected gap colors or remove the fidelity gate.

The patch was exercised only in a temporary copy of the exact Ubuntu
`circos_0.69.9+dfsg-2_all.deb` source/configuration, mounted into the previously
identified SIF along with its Ubuntu font dependencies. It is not a clean
Docker build or a claim about all SVG viewers. The temporary installation and
original images remain outside this repository; original project results and
all initial defect records remain untouched.

[Comparison evidence](comparison.json):

- All SVG path geometry is identical; 17 zero-width path fill styles change.
- Positive-width slices retain their styles. Both native PNGs are byte-identical,
  SHA256 `c8f39e3c43bbf49d7ac8ec82320ead685648a2c2e95618a04875c868f4eef7f6`.
- At (707,1500) and (808,1500), native PNG is white. Before correction, SVG
  rasterization yields red `[239,59,44]` and orange `[253,141,60]`; after correction
  both are `[255,255,255]`. ImageMagick independently sampled these values.
- Both unmodified and proposed-renderer figures pass the same structural
  six-track/two-chromosome/link/legend smoke, including the new project-radius
  invariant. Numerical and nonzero-height geometry were not changed.
- Rasterization used the saved librsvg API reproducer with the exact Liberation
  font registered. An actual corrected-image `rsvg-convert` run remains required.

The correction is approved and configured locally. A successful clean image
build, unchanged fidelity gate and full CI execution are still release gates.
