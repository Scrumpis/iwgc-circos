---
name: iwgc-circos
description: Install, run, troubleshoot, and interpret iwgc-circos genome visualization from FASTA, GFF3, LTR age tables, and alignment coordinates using native tools or containers. Use for iwgc-circos track generation, Circos configuration, rendering, and SVG track labels.
license: MIT
metadata:
  version: "1.0.0-rc.1"
  software-version: "0.4.0-rc.1"
  source: "https://github.com/Scrumpis/iwgc-circos"
---

# Run iwgc-circos

Requires Bash, GNU awk, Python 3.9+, samtools, bedtools and Circos; optional quarTeT/tidk for plant telomeres and an SVG converter for labeled PNG. Linux amd64 with an existing Singularity image was tested. Downloads need network access. Docker builds and other platforms are not yet verified.

Create genome overview figures, not biological significance tests. Read
[references/verification.md](references/verification.md) for executed checks and
limits, and [references/inputs-and-meaning.md](references/inputs-and-meaning.md)
before preparing annotations or interpreting values. The skill folder is
self-contained; runtime scripts come from the software checkout.

## Locate and verify the software

Find the user's checkout, inspect `VERSION` and `git rev-parse HEAD`, and run
`bash iwgc_circos_tracks.sh --help` there. Candidate `0.4.0-rc.1` is unpublished;
do not invent a release tag or assume the public default branch has its fixes.
If no checkout exists, clone `https://github.com/Scrumpis/iwgc-circos` into a
user-approved working location, inspect its revision and version, and use that
revision's documentation. Versions through the audited `fe224c8` contain known
coordinate and link defects; report this and obtain the candidate checkout
before relying on corrected scientific results. Never overwrite an existing
checkout or results to install a newer copy.

Check `command -v singularity apptainer docker samtools bedtools circos gawk python3`.
An existing SIF is the tested route. The repository documents pulling
`docker://scrumpis/iwgc-circos:latest`, but that mutable tag has not been verified
against this candidate. Prefer the user's already verified image or build the
candidate Dockerfile, record its identity, and run the checks below. The image
holds dependencies; execute the scripts from the mounted checkout. Do not run
an older copy of the scripts installed inside an image.

## Check the environment

From the checkout, with an absolute `IWGC_SIF` path supplied by the user:

```bash
singularity exec --cleanenv "$IWGC_SIF" python3 -m unittest discover -s tests -v
singularity exec --cleanenv "$IWGC_SIF" python3 tests/smoke.py
```

The second command runs all four stages on synthetic data and checks numeric
and SVG invariants. It uses a private temporary workspace. A failure means the
installation is not fully verified; examine logs before analyzing user data.
If a documented failure is confined to an optional dependency that the requested
workflow does not use, record it and proceed only with the unaffected workflow.
A core track/config/render failure blocks plotting. Passing
these checks supports the tested cases, not scientific correctness for every
genome. `apptainer` is a documented substitute for `singularity`, untested here.

For native operation install the listed dependencies using the platform's
package manager, then omit `singularity exec --cleanenv "$IWGC_SIF"` from each
command. Docker users can build `docker build --pull --no-cache --platform linux/amd64 -t iwgc-circos:candidate .` and
replace that prefix with `docker run --rm -v "$PWD":/data -w /data iwgc-circos:candidate`.
Maintainer CI is configured to clean-build/test source updates and publish tested
version tags after setup and authorization; it has not executed yet. See the
checkout's `docs/release-process.md`.
The full image includes telomere tools; `--target core` builds only core tracks
and rendering. Builds download dependencies and are not fully locked. Do not
publish images as part of a plotting request.

## Run the four stages

Use a fresh run directory so wildcard matching cannot mix genomes. Copy the
checkout's `iwgc_circos_template.config`, `iwgc_circos/*.conf`, and the default `fonts/LiberationSans-Bold.ttf` plus its `LICENSE.LiberationSans`
into that directory, preserving the same paths; create `iwgc_circos_data` and
`iwgc_circos/tmp`. Invoke the scripts by absolute path in the checkout. Human
step-by-step commands are in its `docs/getting-started.md`.

Paths outside the current directory must be visible in the container. For
Singularity add `-B /absolute/input:/input:ro` before the image and use `/input/...`
paths. FASTA indexing writes `<FASTA>.fai` beside the input: provide a current
index before mounting read-only, or copy the FASTA into the run directory.
Docker additionally needs the input bind mount. For the simplest run, place
inputs in the fresh run directory and use relative names without spaces.

```bash
singularity exec --cleanenv -B "$IWGC_REPO:$IWGC_REPO:ro" -B "$PWD:$PWD" "$IWGC_SIF" bash "$IWGC_REPO/iwgc_circos_tracks.sh" genome.fa \
  -gene genes.gff3 -repeat repeats.gff3 -intact intact.gff3 \
  -ltr-dating genome.pass.list -links genome.coords -gc \
  -window 1000000 -out iwgc_circos_data
singularity exec --cleanenv -B "$IWGC_REPO:$IWGC_REPO:ro" -B "$PWD:$PWD" "$IWGC_SIF" bash "$IWGC_REPO/create_configs.sh" -gap
singularity exec --cleanenv "$IWGC_SIF" circos \
  -conf iwgc_circos/iwgc_circos.conf -outputdir iwgc_circos/tmp
singularity exec --cleanenv -B "$IWGC_REPO:$IWGC_REPO:ro" -B "$PWD:$PWD" "$IWGC_SIF" bash "$IWGC_REPO/add_legend.sh"
```

Set `IWGC_REPO` to the absolute checkout path, not the installed skill folder.
Omit options for unavailable tracks. FASTA alone generates a karyotype; add
`-gc` for a minimal numeric plot. `-intact` needs an EDTA intact annotation,
separate from the total repeat annotation. `-telomere` is optional and runs
quarTeT plant-mode detection at minimum 50 repeats. Use `-keep-temp` to retain
raw intermediate coverage/ages; use `-sliding -step 0.5` for half-window steps.

Outputs are `<FASTA basename>_*.circos` tracks,
`iwgc_circos/iwgc_circos.conf`, `iwgc_circos/ticks.conf`,
`iwgc_circos/tmp/iwgc_circos.svg` and `.png`. The optional label step writes
`iwgc_circos/tmp/iwgc_circos.labeled.svg`, and PNG if a converter is available. Labels `a.`,
`b.`, etc. identify rings, not track names or a color scale. Build the figure
caption from the actual generated plot order, not a hard-coded example.

Stop on nonzero status. A failed multi-track run may leave partial files; do
not render them or interpret them as a completed run. Preserve logs and retry
in a fresh output directory. An empty links file is valid; configuration omits
empty tracks/links. Ambiguous wildcard matches fail rather than selecting a
random genome. Check retained chromosome IDs and counts before rendering.

## Interpret and report

Gene and total-repeat tracks are square-root transformed covered-base fractions,
then min-max scaled over all retained windows. Intact class coverage and GC
are min-max scaled without that power transform. Constant tracks become zero,
even with full coverage. These values are relative within each track/run:
a value of 1 does not imply 100% coverage and values cannot be compared directly
between runs. LTR values invert the mean supplied age into [1,2], with missing
windows 0 and constant present ages 1.5. Zero age is present, not missing.

Links lose orientation and are greedily merged spans, not validated orthology
or proof of rearrangements. Template rules hide intrachromosomal and <20 kb
links. Read the input reference for exact thresholds and coordinate contracts.

Report the code version/commit, image digest or SIF SHA256, input hashes,
command options, retained chromosomes, window/step, normalization, generated
files, failed/skipped steps, and the scope of validation. Keep source data and
existing figures intact. Do not imply independent peer review or a validated
new biological conclusion from a successful render.

## Troubleshoot

- Missing dependency: use the container or install the named tool; GNU awk is required.
- Stale `.fai`: rebuild explicitly with `samtools faidx genome.fa` on a writable copy.
- Malformed input: use tab-separated nine-column GFF3; links are 6/11/13 columns,
  not raw PAF or a decorated MUMmer report. See the input reference.
- No chromosomes/windows: inspect `-filter-chrs` regex and the half-window sliding tail cutoff.
- Missing intact rings: check `Classification=type/class`; only template classes render automatically.
- Single-chromosome telomere detection: tidk 0.2.0 fails on the designed one-record
  positive input; no report now causes failure. Do not interpret that as absence.
- Center fallback warning: symmetric ideogram chords may be parallel, so root-box
  center is used. Inspect letter placement; edited/translated SVG requires review.
- PNG warning: SVG may still be valid. Prefer `rsvg-convert` in the qualified image. The old
  0.69-8 SVG differs from its native PNG at two expected white gaps in both
  tested backends, and ImageMagick also misplaces text. Do not accept such
  conversion as validated output; see the current review record. Explicit `-png FILE` requires converter success.
- Legend geometry error: use original Circos SVG with ideogram/plot groups and at least
  one numeric ring. `-conf`/`-ideo` are retained compatibility options; geometry comes from SVG.
- FUSE/Singularity startup failure: report the runtime error; try a locally supported
  container engine or native environment. Do not disable security controls or claim a render ran.
