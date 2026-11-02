# Inputs, coordinates and figure meaning

Applies to candidate 0.4.0-rc.1. Earlier releases have different and sometimes
incorrect coordinate handling; regenerate derived tracks for this candidate.

## Data contracts

| Input | Required representation | Conversion / use |
|---|---|---|
| FASTA | Uncompressed, nonempty sequence records; unique IDs; current samtools `.fai` | IDs are first header tokens. Supported ID characters: letters, numbers, `_ . : -`. Index beside input is created if absent; stale index is rejected. Copy inputs if the source directory is read-only. |
| Gene GFF3 | Nine tab-separated fields; one-based inclusive starts/ends, start ≤ end | Column 3 exactly `gene`; start minus one, end unchanged. Other feature types are validated but excluded. Comments begin with `#`; `##FASTA` ends feature parsing. |
| Repeat GFF3 | Same coordinates and columns | All feature records contribute their union of bases; nested/duplicate features do not multiply coverage. |
| EDTA intact GFF3 | Same, plus `Classification=type/class` attribute | One coverage file per second slash component (`Gypsy`, `Copia`, etc.). Class token must be safe as a filename. Missing classification fails. |
| LTR pass.list | Tab-separated, at least 12 columns; column 1 `chr:start..end`; column 12 finite nonnegative age | Reverse endpoints are ordered, start minus one. No age estimation occurs here; units and substitution-rate assumptions are inherited from the producer. |
| Six-column links | `chrA startA endA chrB startB endB`, whitespace-separated | Zero-based half-open intervals. Reverse endpoints are ordered; strand is discarded. |
| MUMmer links | Headerless `show-coords` tables with 11 or 13 columns, IDs last | First four fields are one-based inclusive endpoints, other non-ID fields numeric. Convert each ordered start minus one. Decorated pipe-separated reports and raw PAF are rejected. |
| PAF preprocessing | `preprocess/minimap2_self_alignment.sh -f genome.fa -p prefix -t 8` | Runs minimap2 `-x asm5` self-alignment; converts PAF starts/ends without shifting; removes exact same-ID/same-interval diagonal hits. Drops strand. |

All intervals must name a FASTA ID and fall within its length. Known FASTA IDs
excluded by `-filter-chrs` are omitted, rather than treated as unknown. IDs in
GFF3 percent-encoded differently from FASTA must be reconciled first; the tool
does not decode sequence IDs or invent mappings. Gzip FASTA is not supported by
the current input/index naming contract. Annotation BED is not accepted in the
GFF3 options: convert it explicitly before use. Empty annotation tables mean
zero coverage; an annotation without exact `gene` records likewise gives zero.

## Windows and transformations

Default fixed windows tile each retained chromosome in 1,000,000 bp intervals,
including a shorter final window. `-sliding -step 0.5` advances half a window;
step is rounded to the nearest integer bp by AWK (ties depend on the runtime).
Only sliding windows at least half the requested length survive. A step of zero
or one rounding to zero is an error. Empty retained chromosome/window sets fail.

Coverage is the number of distinct bases covered by one or more features divided
by actual window length. It is not the number of features or the sum of their
lengths. For gene and total repeat coverage let `t = sqrt(coverage)`; for intact
classes let `t = coverage`. Output is `(t - min(t))/(max(t)-min(t))`, calculated
across all retained windows for that track, with constant tracks set to zero.

GC uses `(G + C) / actual window length`, case-insensitive. Ambiguous/N bases
stay in the denominator. It is scaled by the same min-max rule; constant GC,
including all-N windows, becomes zero. Displayed GC values are not raw GC
percentages. Use `-keep-temp` and `gc_content.bed` for raw nucleotide counts.

LTR age per window is the unweighted arithmetic mean of every overlapping
insertion's supplied age. A single base of overlap counts that insertion once;
its age contributes fully, not weighted by overlap length. Normalize only
windows having insertions, then display `2 - (mean - min)/(max - min)`.
Identical present means (including zero) produce 1.5. Windows without insertions
produce 0 and are painted white. No mutation rate, sequence divergence, age
confidence interval, or unit conversion is calculated by iwgc-circos.

Windows overlap under sliding mode and are not independent observations.
Scaling removes absolute magnitudes, so different tracks, classes, or genomes
cannot be compared as physical percentages from their plotted colors alone.

## Links and filtering

Empty links are a valid no-alignments result. Malformed noncomment records
fail, including a bad record mixed into an otherwise valid table. Both endpoint
intervals are ordered; the output carries no strand or twist information.

Rows are sorted by first chromosome/start, second chromosome/start. Consecutive
rows with the same chromosome pair are greedily merged when both intervals
overlap or touch; merged bounds are minima/maxima. This is not general connected
component clustering, and intervening chromosome pairs can prevent merging.

`-min-tlen` filters on the sum of the two merged interval lengths (inclusive
threshold), not on either length separately. `-top-n N` keeps the longest N
eligible merged spans, after chromosome filtering; 0 means all. `-link-order asc`
emits selected links from short to long; `dsc` from long to short. Ties follow
C-locale sort order. The input does not establish orthology; self-alignments can
include repeats and paralogs. Reciprocal alignments are not deduplicated across
swapped chromosome pairs.

The default Circos template then hides intrachromosomal links and links with
`max(size1,size2) < 20 kb`. Colors are blue by default, green if *both* sizes
exceed 50 kb, yellow above 100 kb, orange above 200 kb and red above 1 Mb.
Rules use strict greater-than thresholds; they operate on rendered merged spans.

`-filter-chrs` is case-insensitive and applies a GNU AWK extended regex to the
ID as a substring. Omitted pattern is `Chr0?[1-9][0-9]?`, deliberately retained
for compatibility; it is not a complete nuclear chromosome classifier. Supply
anchors (e.g. `^Chr[0-9]+$`) if a full-ID match is required.

## Configuration, labels and optional telomeres

Default template recognizes 13 dynamic tracks: gene, total repeats, LTR age,
Gypsy, Copia, unknown, DTA, DTC, DTH, DTM, DTT, Helitron, GC. Other intact classes
can produce data files but need an explicit template plot to appear. Every
wildcard must resolve to at most one nonempty file; combine species inputs
explicitly or run them in separate directories. Configuration does not infer
scientific compatibility from filenames.

For 1–8 numeric tracks, rings are repacked from outer radius 0.98r toward
inner radius about 0.50r with 0.02r spacing; widths expand equally. For larger
sets, present leading template radii are retained and tracks after a missing
predecessor are repacked with 0.04r width and 0.03r spacing. Link radius is 0.005r
inside the innermost numeric ring. No numeric rings means template link radius.
Tick spacing is a display heuristic based on total length and chromosome count,
with a minimum major tick of 1 Mb; tiny synthetic genomes can have no internal
ticks. `-outdir` writes the master/ticks and copies missing ideogram includes;
it does not change where the template looks for track data.

`-gap` patches the gap between last and first karyotype chromosomes. The copied
example ideogram has species-specific pairwise spacing until `-gap` is run.
Legend letters follow ideogram then SVG `plotN` groups. They require at least
one numeric plot group; arbitrary edited SVG and nonnumeric text/highlight plot
groups have not been validated. Use `--theta-deg` for a different radial angle.
`-conf` and `-ideo` remain accepted but SVG geometry drives placement.

`-telomere` runs quarTeT `te -c plant -m 50 --noplot` on the FASTA. Bands indicate
reported left/right presence, not measured telomere lengths: each marker width
is `int(retained total genome length * -ts)` (default 0.005), clamped to its
chromosome length. This plant-specific detector is not validated for other
clades. Very small scales can round markers to zero; choose a visible scale.
Detector errors must not be interpreted as biological telomere absence. The
audited tidk 0.2.0 fails to nominate a motif on the designed one-chromosome
positive input, while a two-chromosome control passes. The wrapper fails when
no report is produced; this is a known optional dependency limitation.


The 860 bp minimal motif fixture and 30 kb one-chromosome fixture in the final
review produce no discovery candidate with tidk 0.2.0, although direct known-motif
search on the 30 kb sequence finds 100 copies at each end. The current wrapper
returns nonzero for a missing report. This qualifies the discovery workflow;
it does not imply every single-chromosome genome fails or that missing reports
mean biological absence. No discovery method/threshold was changed.
