# Skill and software verification

Candidate: software 0.4.0-rc.1, skill 1.0.0-rc.1, based on iwgc-circos
`fe224c81c41e4aee5905e587ed01e0c1829ca0a3`, after tag v0.3.0.
This candidate is not published. The release evidence lives in the checkout's
`docs/validation/verification.md`, `audit.md`, `source-manifest.json` and test
sources. Do not assume a skill copy describes an older installed executable.

The executing and implementing Codex session is the same evaluator. Known
answers use independent base/interval arithmetic; that does not make the
reviewer independent or certify biological interpretations.

Tested: Linux amd64, Singularity CE 4.1.2; existing dependency SIF SHA256
`5067065f346b5233dc708bf032053a3287a38b7150ca12ca0900cb8419db7cb7`.
It identifies its origin as `scrumpis/iwgc-circos-tracks:latest` and supplies
samtools 1.13, bedtools 2.30.0, Circos 0.69-8, GNU awk, Python 3.13.2,
quarTeT 1.2.5 and tidk 0.2.0. This is a local SIF hash, not a pullable registry
digest or evidence that it was built from the candidate Dockerfile.

The synthetic suite checks GFF3/BED/alignment coordinates, overlapping coverage,
window tails, scaling, missing/zero LTR ages, GC ambiguity, link parsing,
merging/filters/order, input failures and config geometry. The end-to-end smoke
actually generates tracks/config, renders PNG/SVG with Circos, and adds labels.
It checks six nonempty numeric SVG rings, ideogram paths, a rendered link,
letter count/order, finite coordinates, angle override, and PNG dimensions.
Numerical accuracy is checked separately from rendering. The final suite
passed 62 tests with no skips in this image, including a real two-chromosome
telomere positive control and a check that the known single-chromosome detector
limitation is surfaced as failure. This last check validates error propagation,
not the detector's biological sensitivity. See the checkout's
verification record for final counts and exact commands.

The local image lacks an SVG-to-PNG converter, so labeled PNG conversion is
not established by its smoke run. Circos' own PNG output is checked. Candidate
Dockerfile includes librsvg; it still needs an actual clean Docker build.

Unverified unless a later release record explicitly says otherwise: clean
Docker build, Docker execution, GitHub Actions, Apptainer, macOS, Windows,
ARM, general browser font portability, genome-scale performance and biological
annotation/dating accuracy. The Dockerfile pins some upstream artifacts but
is not a dependency lock or bit-for-bit reproducible build.

Skill metadata/links are checked by `python3 tools/validate_skill.py` and the
available skill-creator validator. A fresh agent explicitly loaded this skill, completed the genes/GC workflow,
and recorded the optional detector failure and missing labeled-PNG converter.
It used an earlier candidate revision; script hashes and documentation changes
are recorded in the checkout. Automatic discovery was not tested. Runtime
discovery and fresh-agent trial outcomes are recorded in the checkout's verification record; a valid directory
alone does not demonstrate automatic discovery. Claude has a repository-local
symlink; other clients can explicitly load this SKILL.md or install the whole
folder using their supported convention. The managed `.agents` directory in
the audit environment was read-only and was not modified.

The historical Arial binary is removed from the candidate tree in the follow-up
release review; Git history and original figures are retained. The licensed
Liberation Sans Bold asset/notice are shipped, and the Dockerfile now installs
that exact font for its SVG rasterizer. Native renders confirm unchanged plot
paths and legend positions, with documented text-metric differences. The clean
Docker build remains unexecuted here. Current review evidence and CI activation
steps are in the checkout's `docs/validation/final-review/GO-NO-GO.md` and
`docs/release-process.md`. The original 15 defect observations are preserved.


The follow-up native/SVG raster comparison found color in two intended white
histogram gaps in the old 0.69-8 SVG under both tested backends. Native PNG
probes are white. `tests/check_render_gaps.py` detects this difference; the
new-image CI/publishing checks must pass it. Structural smoke success alone did
not establish SVG/native fidelity. No plot/scientific defaults were changed to
hide this observation; clean-image verification remains pending.
