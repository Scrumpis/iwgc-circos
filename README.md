# International Weed Genomics Consortium (IWGC) Circos Plots

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
Pull container from DockerHub (can use Docker equivalent pull command too)
```
singularity pull iwgc-circos-tracks.sif docker://scrumpis/iwgc-circos-tracks:latest
```

### Usage
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
- If you need to reverse the orientation of any chromosomes to untwist links, edit the ```chromosomes_reverse = ``` line to indicate the chromosomes you wish to reverse (example: ```chromosomes_reverse = Chr1; Chr6; Chr9```


### Notes
- When visualizing multiple species, datasets, etc. in a single Circos plot, all files must be concatenated first
- All tracks are optional aside from the ideogram (telomeres and centromeres optional), so only required file is genome.fasta
- For automation and visualization purposes...
  - Gene desity is sqrt transformed
  - Repeat density is power 3 transformed
  - All density and LTR age tracks are normalized
  - Paried comparisons are difficult, give link to Carvense circos which should have detailed instructions of changes made since they were extensivem, if they have config files for Carvense, they can use those. Give links to all config files for all projects.

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
- Change pairwise header block to your chr exact names
```
<pairwise H2_Chr06 H1_Chr06> # Change chr names to your chr names to space for legend gap (typically last chr then first chr)
# spacing between chr1 and chr9 is 5x 0.1% of image
spacing = 15r 
</pairwise>
```
**Run Circos:**  
Docker (local use; easier refinement - quickly view png/svg, tweak .conf files, and re-run):
```
docker run --rm \
  -v "$PWD:/data" \
  staphb/circos \
  circos -conf /data/iwgc_circos/iwgc_circos.conf -outputdir /data/iwgc_circos/tmp --noparanoid
```
We use --noparanoid to ignore the error about telomere bands exceeding length of chromosomes  
Singularity (HPC use):
```
singularity exec \
  --cleanenv \
  -B "$PWD":/data \
  /path/to/staphb-circos.sif \
  circos -conf /data/iwgc_circos/iwgc_circos.conf -outputdir /data/iwgc_circos/tmp --noparanoid
```
  
## Final Touches
Use Inkspace or a similar tool capable of editing SVG images, or another method you are familiar with which will preserve quality.  
Add legend characters (a, b, c, etc.) to the central gap  
For comparisons with only two ideograms, such as the same chromosome of two haplotypes or species, you can adjust chromosome names this way  
   
### Examples
![image](https://github.com/user-attachments/assets/4f4ebc51-7813-4bb9-a6a4-f49da1c9b119)
  
![image](https://github.com/user-attachments/assets/e9c1d8d7-9814-48be-afea-d34ce2fcfd5a)
  
![image](https://github.com/user-attachments/assets/6dee91e6-992f-4e23-ae8b-c24a87683bf3)
  
![image](https://github.com/user-attachments/assets/eb3b2b77-ce76-4a9b-8ff9-491d5b581885)

![image](https://github.com/user-attachments/assets/b8d23985-ed7e-45ce-b741-8576e1e6f402)

