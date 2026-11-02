# Run a genome overview

You can use iwgc-circos directly from a terminal. No AI agent is needed.
Candidate version is `0.4.0-rc.1`; it is not a published release. See
[the audit](validation/audit.md) before reusing tracks from earlier releases.

## Dependencies

Use an existing tested Singularity image, or build the candidate Dockerfile.
The historical pull command is:

```bash
singularity pull iwgc-circos.sif docker://scrumpis/iwgc-circos:latest
```

That tag is mutable and has not been verified against this candidate. Do not
assume it contains candidate scripts. Run scripts from your checkout, mounted
into the image. Record `sha256sum iwgc-circos.sif` and `git rev-parse HEAD`.
The verification record identifies the existing image actually tested.

Docker build (not executed in this audit; Docker is unavailable):

```bash
docker build --pull --no-cache --platform linux/amd64 -t iwgc-circos:candidate .
```

For automated clean builds and tested-image Docker Hub releases, see
[the release process](release-process.md). Workflows are prepared but have not
executed on GitHub yet.

The default image includes optional quarTeT/tidk. For all other tracks and
rendering, `docker build --target core -t iwgc-circos:core .` avoids the
additional R/conda tools. The image supplies dependencies and a system-installed licensed Liberation font; it
does not copy or install the runtime scripts from this checkout. `CMD` opens Bash.

Native Linux requirements: Bash, GNU awk, Python 3.9+, samtools, bedtools,
Circos and its modules; optional minimap2 for preprocessing, quarTeT/tidk for
plant telomeres, `rsvg-convert` for labeled PNG. On Ubuntu, the core packages
are `gawk python3 samtools bedtools circos minimap2 librsvg2-bin`. Other
platforms are not tested. Apptainer may substitute for Singularity but was not
exercised here.

## Verify the installation

In the checkout, set an absolute path to your image:

```bash
IWGC_SIF=/absolute/path/to/iwgc-circos.sif
singularity exec --cleanenv "$IWGC_SIF" python3 -m unittest discover -s tests -v
singularity exec --cleanenv "$IWGC_SIF" python3 tests/smoke.py
```

The smoke test makes a private temporary genome, generates tracks/config,
renders SVG/PNG, adds ring letters, checks structure and geometry, and cleans
its temporary outputs. To preserve them, append `--output /absolute/new/directory`.
Passing is evidence for these cases, not proof of correctness on every genome.
Inspect failures before continuing. A documented optional-dependency failure
may leave an unrelated genes/GC workflow usable; report that limitation. A core
track, config, or rendering failure blocks use of the figure.
Native users omit the container prefix. Docker users replace it with:

```bash
docker run --rm -v "$PWD":/data -w /data iwgc-circos:candidate python3 tests/smoke.py
```

## Make a fresh run directory

From the checkout:

```bash
IWGC_REPO="$PWD"
RUN_DIR=$(mktemp -d "$PWD/run.XXXXXX")
mkdir -p "$RUN_DIR/iwgc_circos/tmp" "$RUN_DIR/iwgc_circos_data"
cp iwgc_circos_template.config "$RUN_DIR/"
cp iwgc_circos/*.conf "$RUN_DIR/iwgc_circos/"
mkdir "$RUN_DIR/fonts"
cp fonts/LiberationSans-Bold.ttf fonts/LICENSE.LiberationSans "$RUN_DIR/fonts/"
cd "$RUN_DIR"
```

Copy your FASTA and desired annotation/alignment files into this directory.
For a minimal synthetic example, these commands provide FASTA and GFF3:

```bash
cp "$IWGC_REPO/tests/fixtures/genome.fa" genome.fa
cp "$IWGC_REPO/tests/fixtures/genes.gff3" genes.gff3
```

This very small example demonstrates commands; the smoke test also checks
larger intervals and visible links. Genomic IDs must agree exactly across
inputs. FASTA indexing creates `genome.fa.fai` beside your input. A stale index
fails: rebuild it explicitly with `samtools faidx` on a writable copy.

## Create tracks, configure and render

The minimal example uses 4 bp windows; for a real genome choose a suitable
window size (default 1 Mb) and supply your annotation files.

```bash
singularity exec --cleanenv -B "$IWGC_REPO:$IWGC_REPO:ro" -B "$PWD:$PWD" "$IWGC_SIF" bash "$IWGC_REPO/iwgc_circos_tracks.sh" genome.fa \
  -gene genes.gff3 -gc -window 4 -keep-temp -out iwgc_circos_data
singularity exec --cleanenv -B "$IWGC_REPO:$IWGC_REPO:ro" -B "$PWD:$PWD" "$IWGC_SIF" bash "$IWGC_REPO/create_configs.sh" -gap
singularity exec --cleanenv "$IWGC_SIF" circos \
  -conf iwgc_circos/iwgc_circos.conf -outputdir iwgc_circos/tmp
singularity exec --cleanenv -B "$IWGC_REPO:$IWGC_REPO:ro" -B "$PWD:$PWD" "$IWGC_SIF" bash "$IWGC_REPO/add_legend.sh"
```

For data outside the run directory, add `-B /absolute/inputs:/input:ro` to
Singularity before the image name, and use `/input/...` inside it. Ensure the
FASTA already has its current `.fai`, or copy it into the writable run folder.
Docker needs both the run directory mounted at `/data` and any input bind.
When the run directory differs from the checkout, mount the source separately;
for example, the corresponding track command is:

```bash
docker run --rm -v "$IWGC_REPO":/opt/iwgc-circos:ro -v "$PWD":/data -w /data \
  iwgc-circos:candidate bash /opt/iwgc-circos/iwgc_circos_tracks.sh genome.fa \
  -gene genes.gff3 -gc -window 4 -out iwgc_circos_data
```

Use a container-visible absolute path to the checkout when invoking its scripts.
Avoid spaces in figure run/config paths: track generation accepts them but
Circos configuration parsing has not been validated for them.

Additional track options:

```text
-repeat repeats.gff3             total repeat coverage
-intact intact.gff3               per-class EDTA intact coverage
-ltr-dating genome.pass.list      supplied LTR insertion ages
-links genome.coords             syntenic spans (6/11/13 columns)
-telomere                        plant-mode quarTeT telomere markers
-sliding -step 0.5               half-window steps
-filter-chrs '^Chr[0-9]+$'       case-insensitive full-ID regex
```

Read [input contracts and interpretation](../skills/iwgc-circos/references/inputs-and-meaning.md)
before converting coordinates or interpreting colors. GFF3 is one-based
inclusive; generated BED and six-column links are zero-based half-open.
Raw PAF is not link input: use the preprocessing helper.

## Inspect outputs and failures

Final tracks: `iwgc_circos_data/genome.fa_*.circos`. Configuration:
`iwgc_circos/iwgc_circos.conf`, `ticks.conf`, and gap-patched `ideogram.conf`.
Figures: `iwgc_circos/tmp/iwgc_circos.svg` and `.png`. Optional labels:
`iwgc_circos/tmp/iwgc_circos.labeled.svg` and, when a converter is available, the corresponding `.png`.
Legend letters identify rings; supply track names and transformation meanings
in your caption. They do not generate a color key.

Check every command's exit code. Multi-track failures may leave partial files:
preserve logs, start a new run directory and regenerate, rather than rendering
an incomplete run. Empty links mean no eligible alignments and are omitted.
Ambiguous wildcard matches fail: keep only the intended inputs in the data
directory or combine compatible tracks explicitly. Unknown intact classes need
a corresponding custom template plot to render.

Gene and total repeat values are square-root transformed fractions, followed
by per-track min-max scaling. Intact and GC tracks use min-max alone. Constant
tracks are zero even at complete coverage. LTR values invert mean supplied age
into [1,2]; missing is 0 and constant present ages are 1.5. These display scales
are not physical percentages or biological hypothesis tests. See the input
reference for full formulas and limitations.

Record code/image identities, input hashes, commands, window/step, filters and
validation evidence with your figure. Keep original data and figures intact.


The tested tidk 0.2.0 fails to nominate a telomere motif on the designed
single-chromosome positive example; the two-chromosome positive control passes.
The wrapper now fails when no detector report exists, rather than treating an
upstream success status as evidence of telomere absence. This limitation is not
fixed upstream here. The legend may warn about root-center fallback on
symmetric ideograms; inspect the letter positions, especially on edited SVGs.


The single-record detector limitation is explained with a minimal fixture,
known-motif positive search and expected error behavior in the
[current review](validation/final-review/GO-NO-GO.md). It does not mean the
chromosome lacks telomeres. Do not duplicate sequences merely to force motif
nomination; record unavailable detection or omit the optional track.


For labeled PNG, prefer `rsvg-convert` (included in the candidate Dockerfile).
In the old Circos 0.69-8 runtime, both tested SVG rasterization backends paint
two intended white ring gaps differently from the native Circos PNG; ImageMagick
also misplaced chromosome text. Those conversions are not accepted as figure
validation. The clean Docker image must pass the added gap-fidelity check; its
CLI conversion still awaits CI execution. The font comparison records a roughly 8 px shift
for the long chromosome label and sub-pixel tick-label shifts; check figures
when migrating fonts. Plot paths and automatic ring-letter positions were
unchanged in the comparison.
