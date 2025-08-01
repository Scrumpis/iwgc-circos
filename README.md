# International Weed Genomics Consortium (IWGC) Genome Report Circos Plots
## About
This repository documents the methods used to generate the Circos plots for many of the genome reports produced by the International Genomics Consortium (IWGC). This project aims to automate the generation of commonly used Circos tracks in genomics to reduce the barrier for entry for unfamaliar with Circos and to speed-up production for those experienced with Circos.

For additional information, the [Circos](https://circos.ca/) website offers very detailed tutorials for the generation of Circos plots and tracks beyond the scope of these standard genomics tracks.  
  
## iwgc_circos_tracks.sh
**Generate track data files for any of the following:**
- Ideogram/Karyotype
  - Telomere presence
- Total gene density
- Total repeat density
- Intact transposable element density (separated by family, cat all together for total density)
- GC content
- LTR insertion age

  
### Setup
Pull container from [DockerHub](https://hub.docker.com/r/scrumpis/iwgc-circos-tracks)  
Singularity (HPC usage):
```
singularity pull iwgc-circos-tracks.sif docker://scrumpis/iwgc-circos-tracks:latest
```
Docker (Local usage):
```
docker pull scrumpis/iwgc-circos-tracks:latest
```

### Usage
```
./iwgc_circos_tracks.sh [options] <FASTA> [<GENES>] [<REPEATS>] [<INTACT>] [<LTRDATES>] <IWGC_CIRCOS_SIF> [WINDOW]
Options:
  -gene            Add gene density track (requires gene annotation GFF3)
  -repeat          Add repeat density track (requires EDTA repeat annotation: EDTA/genome.fasta.EDTA.TEanno.gff3)
  -intact          Add intact TE density track (requires EDTA annotation: EDTA/genome.fasta.mod.EDTA.intact.gff3)
  -gc              Add GC content track
  -ltr-dating      Add LTR dating track (requires EDTA annotation: EDTA/genome.mod.EDTA.raw/LTR/genome.fasta.mod.pass.list)
  -telomere        Add telomere bands to ideogram (karyotype.circos)
  -ts <value>      Telomere band size scale (default: 0.005)
  -keep-temp       Keep intermediate files
  -sliding         Use sliding windows instead of fixed
  -step <value>    Step size for sliding windows (default: WINDOW/2)
  -local           Use Docker instead of Singularity (default: false)
```
    
**Input files:**
| File Name                             | Description                                                                 |
|--------------------------------------|-----------------------------------------------------------------------------|
| `genome.fasta`                       | Required input for karyotype telomeres, and GC-content tracks. This is the only required file for the script to run. |
| `gene.annotation.gff3`               | Gene annotation file used to generate the **total gene density track**. Accepts standard GFF3 format from any annotation pipeline. |
| `genome.fasta.mod.EDTA.TEanno.gff3` | Repeat annotation file from EDTA used to create the **total repeat density track**. |
| `genome.fasta.mod.EDTA.intact.gff3` | Intact repeat annotation from EDTA used to generate **intact repeat density tracks**. |
| `genome.fasta.mod.pass.list`        | LTR insertion age from EDTA used to generate **LTR age track**. |

  
The below would produce all tracks with 500kbp sliding windows with half window size steps
```
./iwgc_circos_tracks.sh -gene -repeat -intact -gc -ltr-dating -telomere -sliding \
genome.fasta genes.gff \
repeats.fasta.mod.EDTA.TEanno.gff3 repeats.fasta.mod.EDTA.intact.gff3 repeats.fasta.mod.pass.list \
iwgc-circos-tracks_v1.sif 500000
```

### Unincluded Tracks
**Links**
Circos links tracks are not produced by this script because there is too little preprocessing required for circos and alignments are usually too custom to easily script:
- Align your regions of interest with minimap2 or similar
- If needed, convert output into a .coords file. I have used the below to convert .paf to .coords before, but other tools exist:
```
awk 'BEGIN { OFS="\t" } { printf "%s\t%d\t%d\t%s\t%d\t%d\t%s\t%.2f\n", \
$1, $3+1, $4, $6, $8+1, $9, $5, ($10/$11)*100 }' CirAr_female_haps_Chr06.paf > CirAr_female_haps_Chr06.coords
```
- Once you have a coords file, extract the columns of interest (query, qstart, qend, subject, sstart, send), collapse links which overlap in both the subject and query, and sort by qstart from least (top of file) to greatest (bottom) to overlay larger links on top of smaller ones:
```
awk 'BEGIN{OFS="\t"} {
  print $1, $2, $3, $4, $5, $6
}' CirAr_female_haps_Chr06.coords | \
sort -k1,1 -k2,2n -k4,4 -k5,5n | \
awk 'BEGIN{OFS="\t"}
{
  if ($1 != qchr || $4 != schr || $2 > qend || $5 > send) {
    if (NR > 1)
      print qchr, qstart, qend, schr, sstart, send;
    qchr = $1; qstart = $2; qend = $3;
    schr = $4; sstart = $5; send = $6;
  } else {
    if ($3 > qend) qend = $3;
    if ($6 > send) send = $6;
  }
}
END {
  print qchr, qstart, qend, schr, sstart, send;
}' | \
awk 'BEGIN{OFS="\t"} {
  qlen = $3 - $2;
  slen = $6 - $5;
  tlen = qlen + slen;
  print tlen, $0
}' | \
sort -k1,1n | \
cut -f2- > CirAr_female_haps_Chr06.coords.circos
```
- If you need to reverse the orientation of any chromosomes to untwist links, edit the ```chromosomes_reverse = ``` line of ```iwgc_circos.conf``` to indicate the chromosomes you wish to reverse (example: ```chromosomes_reverse = Chr1; Chr6; Chr9```


### Notes
- When visualizing multiple species, datasets, etc. in a single Circos plot, all files must be concatenated first
- All tracks are optional aside from the ideogram (telomeres and centromeres optional), so only required file is genome.fasta
- For automation and visualization purposes...
  - Gene desity is sqrt transformed
  - Repeat density is power 3 transformed
  - All density and LTR age tracks are normalized
  - Paried comparisons are difficult, give link to Carvense circos which should have detailed instructions of changes made since they were extensivem, if they have config files for Carvense, they can use those. Give links to all config files for all projects.
  - In special situtations where you need to edit the housekeeping.conf, such as if you have more contigs than what is typically allowed, you can either run the container interactively and update the housekeeping.conf file or you can create a housekeeping.conf file in the /circos directory and map it to the container, see Fusarium paper

### Known Warnings:
- (iwgc_circos_tracks.sh) ```File species.fasta_windows.bed has a record where naming convention (leading zero) is inconsistent with other files```
  - Seems ignorable with no effects, happens at Chr10 if previous are Chr01-09.
- (Circos) ```Use of uninitialized value in subroutine entry at /opt/circos/bin/../lib/Circos/Configuration.pm line 781```
  - Seems ignorable.
  
## Circos Plot
### Setup
Move all desired track files (.circos) into the iwgc_circos_data directory
```
iwgc_circos/
├── iwgc_circos/
    ├── tmp/
        ├── iwgc_circos.png
        ├── iwgc_circos.svg
    ├── axis.conf
    ├── background.conf
    ├── bands.conf
    ├── ideogram.conf
    ├── ideogram.label.conf
    ├── ideogram.position.conf
    ├── iwgc_circos.conf
    ├── r0r1.conf
    ├── ticks.conf
├── iwgc_circos_data/
    ├── gene_density.circos
    ├── repeat_density.circos
    ├── intact_helitrons.circos
    ├── links.circos
    ├── any.other.track.files.circos
```
  
**Manually edit files:**  
**iwgc_circos.conf (Required)**
- Change File Names:
  - Karyotype
  - Any Tracks  
- Change chromosomes following ```chromosomes =```
  - The below script can be used in the iwgc_circos_data to print chrs in correct format to stdout, then copy-paste into config file):
    ```
    grep -v '_T[0-9]\+' *karyotype.circos | awk '{print $3}' | paste -sd';' -
    ```
- Reduce ```chromosomes_units =``` for small genomes or increase for big ones if needed, try with default 1000000 (1Mbp) first

**ideogram.conf**
- Change pairwise header block to your chr exact names to put a gap for the legend (typically last chr then first chr)
```
<pairwise H2_Chr06 H1_Chr06> # Change to your chrs
spacing = 15r # spacing between chr1 and chr9 is 5x 0.1% of image
</pairwise>
```

**ticks.conf**
- If you want more or less ticks, different sizes of ticks, labels or no labels, etc.
  
**Run Circos 0.69-9:**  
*Run one level up from iwgc_circos and iwgc_circos_data*  
Docker (local use; easier refinement - quickly view png/svg, tweak .conf files, and re-run):
```
docker run --rm \
  -v "$PWD:/data" \
  scrumpis/iwgc-circos-tracks \
  circos -conf /data/iwgc_circos/iwgc_circos.conf -outputdir /data/iwgc_circos/tmp -noparanoid
```
We use --noparanoid to ignore the error about telomere bands exceeding length of chromosomes  
Singularity (HPC use):
```
singularity exec \
  --cleanenv \
  -B "$PWD":/data \
  /path/to/iwgc-circos-tracks.sif \
  circos -conf /data/iwgc_circos/iwgc_circos.conf -outputdir /data/iwgc_circos/tmp --noparanoid
```
  
## Final Touches
Use Inkspace or a similar tool capable of editing SVG images, or another method you are familiar with which will preserve quality.  
Add legend characters (a, b, c, etc.) to the central gap  
For comparisons with only two ideograms, such as the same chromosome of two haplotypes or species, you can adjust chromosome names this way  
   
## Example:
The below Circos plot illustrates [all tracks](https://github.com/Scrumpis/iwgc-circos/tree/main/examples/all-tracks), except LTR dating and links, producable by iwgc_circos_tracks.sh using the IWGC _Avena fatua_ genome assembly. Tracks can be removed and reordered using the provided iwgc_circos.conf file as a template. All values are normalized and gene density is square-root transformed and repeat density is power three transformed for automation and visualiztion purposes.
- a. Karyotype illustrating chromosomes, chromosome position, and telomere presense (distal grey boxes)
- b. Gene density, more blue is more gene rich, more yellow is less gene rich
- c. Repeat density, more red is more repeat rich, more yellow is less repeat rich
- d. Intact Gypsy LTR coverage
- e. Intact Copia LTR coverage
- f. Intact non-Gypsy or Copia LTR coverage
- g. Intact hAT-superfamily (DTA) coverage
- h. Intact CACTA-superfamily (DTC) coverage
- i. Intact Harbinger-superfamily (DTH) coverage
- j. Intact Mutator-superfamily (DTM) coverage
- k. Intact Tc1-Mariner-superfamily (DTT) coverage
- l. Intact Rolling-circle-transposon/Helitron-superfamily (Helitron) coverage
- m. GC-content
  

![image](https://github.com/Scrumpis/iwgc-circos/blob/main/examples/all-tracks/iwgc_circos/tmp/iwgc_circos.png)
