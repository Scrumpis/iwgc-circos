# Reviewed scientific and default changes

Status: the maintainer explicitly approved the corrections below, and they are implemented in the candidate. The original before/after proposals are retained as the review record. Baseline: `fe224c81c41e4aee5905e587ed01e0c1829ca0a3` (after v0.3.0). Production behavior remained unchanged until the corresponding approval arrived.

The same Codex session inspected, implemented probes, and executed them. Numerical expectations below were derived from interval arithmetic rather than copied from the program. These are synthetic checks, not independent peer review.

## Proposed corrections

| Change | Evidence and before/after | Proposed resolution |
|---|---|---|
| GFF3 and LTR coordinates | GFF3 uses 1-based inclusive coordinates. On two 12 bp chromosomes, window 4, genes on Chr1 at 2..4, 5..8, 10..11, current normalized coverage is `0.816497,1,0.577351,0,0,0`. Counting the actual covered bases gives `sqrt(3/4),1,sqrt(2/4),0,0,0` = `0.866025,1,0.707107,0,0,0`. | Subtract one from starts when converting GFF3 / LTR `chr:start..end` to BED. Keep square-root and min-max transformations. Apply to total genes, total repeats, intact classes, and LTR overlap. |
| Gene selection | `grep gene` scans all columns: an mRNA or exon whose attributes contain `gene` is included in the gene track. | Select column 3 exactly equal to `gene`; document that annotations with only mRNA/CDS features need explicit preprocessing. |
| Reverse MUMmer query intervals | Input `1 4 10 7 4 4 99 12 12 Chr1 Chr2` becomes `Chr1 1 4 Chr2 10 10` (exit 0), losing the query span. | Repair the swap and convert one-based closed MUMmer endpoints to zero-based half-open coordinates. Retain the established orientation-free six-column output; do not claim strand preservation. |
| Link merging | Current merge condition checks only each next start against the previous end, not whether the next end precedes the previous start on the second chromosome. | Require overlap/touch on both sides and extend both lower and upper bounds. Preserve consecutive sorted greedy merging; do not introduce global clustering. |
| Top N and chromosome filtering | With totals 2,4,6 and `-top-n 2`, default ascending order keeps 2,4, despite help promising the longest links. Filtering also follows truncation. | Select longest N eligible merged links, then order those ascending or `dsc`; filter against retained karyotype IDs before selection. Preserve `-top-n 0` = all. |
| PAF helper coordinates | PAF already uses zero-based half-open intervals; helper adds one to both starts, shortening every converted span. | Emit PAF starts unchanged in six-column format. Preserve exact self-diagonal exclusion and orientation-free contract. |
| Telomere band bounds | Band size is `int(retained_genome_length * scale)` for every chromosome; on a short chromosome this may exceed its length. | Clamp each band size to the chromosome length; retain default scale and document these are presence markers, not measured telomere lengths. |

## Font default proposal

Bundled `fonts/ArialBold.ttf` identifies Monotype copyright 2006. Its internal license allows use according to the product license, not unrestricted redistribution. No redistribution grant accompanies it in this repository. Proposal: add Liberation Sans Bold under its included SIL OFL 1.1 and use it for the template and legend. Retain the historical font for now, as specified in the reviewed proposal; rights or removal remain a maintainer release decision. This changes text metrics and requires a visual review before publication. The retained Arial binary remains a redistribution blocker until rights are established or the maintainer removes it.

## Retained definitions

Do not change repeat square-root to cubing: current HEAD uses square-root, and README is inconsistent. Correct the prose. Retain global per-track min-max scaling, constant coverage/GC = 0, LTR missing = 0, identical ages = 1.5, and present ages inverted to [1,2]. Retain sliding tail cutoff of half a window. These are visualization transformations, not comparable physical densities between datasets.

## Sources

- [GFF3 specification](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md): columns 4/5 are one-based inclusive.
- [bedtools coverage](https://bedtools.readthedocs.io/en/latest/content/tools/coverage.html): fraction covered by at least one feature, not summed depths.
- [minimap2 PAF format](https://github.com/lh3/minimap2/blob/master/PAF.md): zero-based starts, half-open ends, separate strand column.

The raw initial execution is preserved in `baseline-probes.json`; commands and subsequent evidence are in `verification.md`.

## Approval and additional precision finding

The maintainer explicitly approved the documented scientific corrections and
Liberation Sans default during this audit on 2026-11-02.

An additional numerical issue was reproduced after that approval: the existing
multi-stage AWK pipelines use the six-significant-digit default print format.
`sqrt(0.9999999)` and `sqrt(1)` both serialize as `1`, so a nonconstant track can
collapse to all zero during normalization. Proposed correction: retain 17-digit
floating-point precision through transforms and compute coverage fractions from
bedtools' integer covered-base and window-length columns, retaining sufficient
precision until final output. This preserves the intended transformations,
while distinguishing nearby values. A regression will use two 10,000,000 bp
windows covered by 9,999,999 and 10,000,000 bases; expected normalized values are
exactly 0 and 1, not 0 and 0. GC should similarly use nucleotide counts to avoid
bedtools' rounded fraction column. The maintainer explicitly approved this additional precision correction during the audit.

## Final unusual-ID correction proposed

A frozen-baseline configuration probe with two chromosome records, `Chr1`
(length 12) and `scaffold_T1` (length 8), returned exit 0 but wrote
`chromosomes = Chr1`. A historical `_T[0-9]+` name heuristic excludes the second
actual chromosome from both the display list and tick statistics, even though
the records are already selected by type `chr` (telomere bands have type
`band`). Proposal: remove this redundant name exclusion and use every valid
`chr` record. This changes display/length accounting only for IDs matching that
heuristic; it leaves actual telomere bands excluded by record type. The maintainer explicitly approved including all valid chromosome records
during this audit. The redundant name exclusion is removed and tested.


## Follow-up licensing disposition

The subsequent release-validation request expressly requires removing any
unlicensed Arial binary from the candidate. `fonts/ArialBold.ttf` is deleted
in the uncommitted candidate diff; historical Git objects and original figures
are preserved. The earlier retention proposal above records the initial review,
not the final distribution policy. See `fonts/NOTICE.md` and the final GO/NO-GO
record for the current asset inventory and font-metric comparison.
