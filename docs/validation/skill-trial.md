# Fresh end-user skill trial

The requested genes + GC overview completed in `<trial>`.
The skill was explicitly supplied and loaded; this trial does not demonstrate automatic discovery.
This is an end-user usability trial, not an independent scientific review.

## Artifacts

- Labeled figure: `iwgc_circos/tmp/iwgc_circos.labeled.svg` (3000 x 3000).
- Original render: `iwgc_circos/tmp/iwgc_circos.svg` and `iwgc_circos.png` (PNG visually inspected).
- Final tracks: `iwgc_circos_data/genome.fa_gene_coverage.circos`, `genome.fa_gc.circos`, and `genome.fa_karyotype.circos`.
- Retained intermediate inputs/counts: `genome.fa_windows.bed`, `genome.fa_genes_coords.bed`, `gc_content.bed`, `genome.fa_chrs.fai`.
- Configuration: `iwgc_circos/iwgc_circos.conf`, `ticks.conf`, gap-patched `ideogram.conf`, and copied includes/fonts.
- Reproducibility: `run.sh` is the initial script, which stopped at the first runtime failure; `commands.log` captures that execution. `recovery-commands.sh` records the exact subsequent container commands. Individual step logs preserve their outputs. `identities.sha256` records image, inputs, runtime scripts and template hashes. Do not rerun these scripts over the existing outputs; substitute a fresh directory.

## What the rings mean

The labeled SVG XML contains `a.`, `b.`, `c.` and two plot groups. The actual generated config puts genes at radii 0.7500–0.9800r and GC at 0.5000–0.7300r. Thus `a.` identifies the chromosome ideogram, `b.` gene coverage (heatmap), and `c.` GC (line). Letters are ring identifiers and require this caption; they do not provide a color key.

Both Chr1 and Chr2 were retained, each 12 bp (24 bp total). No chromosome-filter option was supplied. Fixed 4 bp windows, with 4 bp advance and no sliding, yielded six zero-based half-open intervals. GFF3 starts were one-based inclusive and converted by subtracting one; only exact `gene` features contribute, so the mRNA feature is excluded. Covered bases are union coverage, not feature counts.

| Chromosome | Interval | Covered gene bases / length | Raw GC | Gene display | GC display |
|---|---|---|---|---|---|
| Chr1 | [0,4) | 3/4 | 2/4 | 0.866025 | 0.5 |
| Chr1 | [4,8) | 4/4 | 2/4 | 1 | 0.5 |
| Chr1 | [8,12) | 2/4 | 2/4 | 0.707107 | 0.5 |
| Chr2 | [0,4) | 0/4 | 4/4 | 0 | 1 |
| Chr2 | [4,8) | 0/4 | 4/4 | 0 | 1 |
| Chr2 | [8,12) | 0/4 | 0/4 (NNNN) | 0 | 0 |

Gene fractions are square-root transformed and then min-max scaled across all six windows. Here their transformed minimum is 0 and maximum 1. GC is (G+C)/actual window length, including Ns in the denominator, then min-max scaled. Here raw GC already spans 0–1, so scaled and raw fractions happen to agree. That is specific to these fixtures. In general plotted 1 does not imply 100% coverage; values across tracks/runs are not physical percentages. Constant tracks become zero. The GC line connects window values; it is a display, not per-base estimation. This tiny fixture is an execution demonstration without biological significance or general annotation validation. The minimum 1 Mb major ticks mean no meaningful internal ticks for this 12 bp example.

## Execution evidence and limitations

Software version `0.4.0-rc.1`; base HEAD `fe224c81c41e4aee5905e587ed01e0c1829ca0a3`, with pre-existing modified candidate files (full status in `commands.log`). Commit alone does not identify the executed candidate; use the recorded script/template hashes. Source and existing results were not changed by this trial. Scripts came from the supplied checkout mounted read-only, not from the image. Inputs were copied to this private writable directory before indexing. There were no downloads or external services.

SIF SHA256: `5067065f346b5233dc708bf032053a3287a38b7150ca12ca0900cb8419db7cb7`.
FASTA SHA256: `7688241fd4e02cc4e24f336ec5c45aec4b295c5fb1a07293b9eeac6c24296a2a`.
GFF3 SHA256: `30b57673cf1891293afaed51ed1479b579d4e2c08abd0e47c2498ad1c8b9a197`.

Runtime: Singularity CE 4.1.2-jammy; samtools 1.13; bedtools 2.30.0; Circos 0.69-8; GNU awk 5.1.0; container Python 3.13.2. See `dependency-versions.log`.

1. Initial restricted container startup failed with status 255: FUSE device unavailable, followed by extraction failure because no setuid installation was found. `unittest.log` preserves the error. The same local SIF ran under scoped Singularity execution escalation; no security configuration was changed.
2. Recovered prerequisite suite: 54 tests, 1 failure (`test_real_plant_telomere_detector`), 53 passed. The optional quarTeT/tidk detector reported no telomere-like repeats and produced no report. See `unittest.recovery.log`. This is an installation limitation; do not claim a fully passing installation. After inspecting that failure, I proceeded with genes/GC because neither uses telomeres, and the parent evaluator explicitly confirmed this scope decision. I did not repair or modify that dependency.
3. Smoke check returned 0 with JSON status passed, six numeric tracks, two chromosomes, rendered link and legend. This is execution/regression evidence from the provided suite, not independent scientific review. No test implementation was opened. Only the failing traceback in the execution log was inspected for troubleshooting.
4. Requested track generation, configuration, Circos rendering and letter addition each returned 0. XML inspection confirmed the three letters and two numeric groups; the original PNG was visually inspected.
5. Label addition warned that it could not solve the ideogram center and used the root center (1500,1500). It still wrote a parseable labeled SVG. Placement was not independently geometrically reviewed. No SVG-to-PNG converter exists in the supplied image, so there is no labeled PNG; the original Circos PNG exists.
6. Minor inspection recoveries: an overly specific whitespace search found no config matches (rg exit 1), then a flexible whitespace regex revealed the expected gene/GC plot order. I initially looked for the labeled SVG at the run root; that lookup failed, then the path printed by `labels.log` led to the actual nested output. Neither affected plotting.

## Instruction usability

The skill plus getting-started guide was sufficient to prepare inputs and run the four stages without reading source code, audit reports, numerical test expectations, or test implementation. Help and the allowed fixture files were inspected; runtime script bytes were hashed, not opened. The explicit source checkout/image paths made installation straightforward. Fresh directory, copying inputs before indexing, same-path read-only repo binds and retaining raw GC counts were useful instructions.

Potential improvements: give the full default labeled output path (`iwgc_circos/tmp/iwgc_circos.labeled.svg`) consistently instead of only its basename; state whether a failure confined to an optional, unused dependency blocks a minimal genes/GC workflow; explain how a user should assess the ideogram-center fallback warning. The guide already warns that the existing image lacks labeled PNG conversion and that tiny synthetic genomes have limited tick display. Automatic skill discovery was not tested.


Parent evaluator re-read the generated numeric tracks and labeled SVG: both tracks
agreed with independently enumerated base counts within 3e-6, and letters were
`a., b., c.`. This checks artifacts beyond the trial agent's report. The production
scripts subsequently changed; these trial hashes identify its earlier candidate.
