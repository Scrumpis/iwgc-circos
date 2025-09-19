# International Weed Genomics Consortium (IWGC) Genome Report Circos Plots
## About
This repository documents the methods used to generate the Circos plots for many of the genome reports produced by the [International Weed Genomics Consortium (IWGC)](https://www.weedgenomics.org/). This project aims to automate the generation of commonly used Circos tracks in genomics to reduce the barrier for entry for those unfamaliar with Circos and to speed-up production for those experienced with Circos. Useful as a template for further customization or to get a quick snapshot of patterns occurring in your genome of interest.

For additional information, the [Circos](https://circos.ca/) website offers very detailed tutorials for the generation of Circos plots and tracks beyond the scope of these standard genomics tracks.  
  
### Example:
The below Circos plot illustrates all tracks producable by iwgc-circos. Select as many of them as you would like for your output. All values are normalized, gene density is square-root transformed, and repeat density is power three transformed for automation and visualization purposes. Ideogram tick label units are the default 1Mb (50 = 50Mb). Windows are 1Mb, with sliding half steps (500Kbp). 
- a. Karyotype illustrating chromosomes, chromosome position, and telomere presense (distal grey boxes)
- b. Gene density - more blue is more gene rich, more yellow is less gene rich
- c. Repeat density - more red is more repeat rich, more yellow is less repeat rich
- d. LTR insertion age - more purple are newer insertions, less purple are older, white contain no insertions
- e. Intact Gypsy LTR coverage
- f. Intact Copia LTR coverage
- g. Intact non-Gypsy or Copia LTR coverage
- h. Intact hAT-superfamily (DTA) coverage
- i. Intact CACTA-superfamily (DTC) coverage
- j. Intact Harbinger-superfamily (DTH) coverage
- k. Intact Mutator-superfamily (DTM) coverage
- l. Intact Tc1-Mariner-superfamily (DTT) coverage
- m. Intact Rolling-circle-transposon/Helitron-superfamily (Helitron) coverage
- n. GC-content
- Syntenic links on the interior
  

![image](https://github.com/Scrumpis/iwgc-circos/blob/main/iwgc_circos/tmp/iwgc_circos.labeled.svg)

  
## Clone Repo
The scripts are built to run with this repo structure by default. It will simplify things for you to do the same.
```
git clone https://github.com/Scrumpis/iwgc-circos
```

## Setup
Pull container from [DockerHub](https://hub.docker.com/r/scrumpis/iwgc-circos-tracks). You can use [Singularity](https://apptainer.org/) or [Docker](https://docs.docker.com/desktop/) to run these scripts. Both give the same results.  
   
**Singularity** (Typically for HPC usage):
```
singularity pull iwgc-circos-tracks.sif docker://scrumpis/iwgc-circos-tracks:latest
```
**Docker** (Typically for local usage):
```
docker pull scrumpis/iwgc-circos-tracks:latest
```

## Preprocess
A genomic.fasta is the only required input file for every Circos. We include a couple of sample scripts for repeat annotation with [EDTA](https://github.com/oushujun/EDTA) and [minimap2](https://github.com/lh3/minimap2) self-alignment and .coords conversion as templates. Minimap2 is in iwgc-circos container, EDTA is not to avoid further bloating the container. A recommended EDTA container is included in the sample script. Below are descriptions of the required input files for certain tracks.

**Input files:**
| File Name                             | Description                                                                 |
|--------------------------------------|-----------------------------------------------------------------------------|
| `genome.fasta`                       | Required input for karyotype, telomeres, and GC-content tracks. This is the only required file for the base script to run. |
| `genome.gff3`                        | Gene annotation file used to generate the **total gene density track**. Accepts standard GFF3 format from any annotation pipeline. |
| `genome.fasta.mod.EDTA.TEanno.gff3`  | Repeat annotation file from EDTA used to create the **total repeat density track**. (EDTA default path: EDTA/genome.fasta.mod.EDTA.TEanno.gff3)|
| `genome.fasta.mod.EDTA.intact.gff3`  | Intact repeat annotation from EDTA used to generate **intact repeat density tracks**. (EDTA/genome.mod.EDTA.raw/LTR/genome.fasta.mod.pass.list)|
| `genome.fasta.mod.pass.list`         | LTR insertion age from EDTA used to generate **LTR age track**. (EDTA/genome.mod.EDTA.raw/LTR/genome.fasta.mod.pass.list)|
| `genome.coords`                      | .coords alignment file from minimap2 or similar used to generate **links track**. |

  
## 1. Create Track Files for Circos Plot (iwgc_circos_tracks.sh)
**Generate input data files for any of the following tracks:**
- Ideogram/Karyotype
- Telomere presence
- Total gene density
- Total repeat density
- LTR insertion age
- Intact transposable element density (separated by family)
- GC content
- Syntenic links

  
### Usage
```
  Usage: ./iwgc_circos_tracks.sh <FASTA> [options]
  
   Required:
     FASTA            Genomic FASTA file
  
   Optional:
     -gene            Add gene density track (requires gene annotation GFF3)
     -repeat          Add repeat density track (requires EDTA repeat annotation: EDTA/genome.fasta.mod.EDTA.TEanno.gff3)
     -intact          Add intact TE density track (requires EDTA intact repeat annotation: EDTA/genome.fasta.mod.EDTA.intact.gff3)
     -ltr-dating      Add LTR dating track (requires EDTA repeat annotation: EDTA/genome.mod.EDTA.raw/LTR/genome.fasta.mod.pass.list)
     -links           Add syntenic links from a .coords file (e.g., MUMmer or minimap2 output)
     -min-tlen        Minimum length of syntenic links
     -top-n           Number of top syntenic links to include (e.g., "20" for top 20 longest links)
     -link-order      Order of syntenic links - asc: smallest top to biggest bottom | dsc: big top to small bottom (default: asc)
     -gc              Add GC content track
     -telomere        Add telomere bands to ideogram (karyotype.circos)
     -ts              Telomere band size scale (default: 0.005), 0.5% of total genome size
     -window          Window size in base pairs (default: 1000000)
     -sliding         Use sliding windows instead of fixed
     -step            Step size for sliding windows (default: 0.5). The default is half window size steps
     -filter-chrs     Restrict chromosomes to those matching typical nuclear naming patterns (e.g., Chr01, Chr1, chr01B). Default: off
     -keep-temp       Keep intermediate files
     -out             Output directory for Circos track files (default: current directory)
     -h | --help      List usage options
```

     
The below commands will produce all possible track files with 1Mbp (default size) sliding windows with half window size steps (default size), telomere labels, and filter for typical chromosome names (Chr1, chr01, chr2B, etc.; not Chr00, ChrMT, N00011.1, etc.).    
  
**Singularity:**
```
singularity exec iwgc-circos-tracks.sif ./iwgc_circos_tracks.sh genome.fasta \
-gene genome.gff \
-repeat genome.fasta.mod.EDTA.TEanno.gff3 \
-intact genome.fasta.mod.EDTA.intact.gff3 \
-ltr-dating genome.fasta.mod.pass.list \
-links genome.coords \
-gc -telomere -sliding -filter-chrs \
-out iwgc_circos_data/
```
**Docker:**  
_Note: Give path to input file after ```/circos```. If local path is ```./genome.fasta```, then use ```/circos/genome.fasta```_
```
docker run --rm -v "$PWD":/circos scrumpis/iwgc-circos-tracks:latest /circos/iwgc_circos_tracks.sh /circos/data/genome.fasta \
-gene /circos/genome.gff \
-repeat /circos/genome.fasta.mod.EDTA.TEanno.gff3 \
-intact /circos/genome.fasta.mod.EDTA.intact.gff3 \
-ltr-dating /circos/genome.fasta.mod.pass.list \
-links /circos/genome.coords \
-gc -telomere -sliding -filter-chrs \
-out /circos/iwgc_circos_data/
```
  
## 2. Create Circos Plot Config Files (create_configs.sh)
This command will produce the following:
- iwgc_circos/iwgc_circos.conf file containing plots for every included track data file contained in iwgc_circos_data. So, if iwgc_circos_data/genome.fasta_gene_coverage.circos exists, gene density plot will be created.
- If the gap flag is invoked, iwgc_circos/ideogram.conf file will contain a legend gap between your last and first chromosomes in the pairwise header.
- If 8 or less of the files associated with any of the plot tracks (excludes ideogram, labels, or links) are present in iwgc_circos_data, tracks will be dynamically resized to keep the inner most tracks r0 value near 0.5 and the gap between tracks is doubled (0.02 units instead of 0.01).
- Dynamically adjusts tick amount and spacing based on genome length and chromosome amount.
  
### Setup
Move all desired track files (.circos) into the iwgc_circos_data directory if not already there. Make sure ONLY the track files you wish to visualize are present in iwgc_circos_data.  
```
iwgc-circos/
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
  
### Usage
```
  Usage: ./create_configs.sh [options]
  
   Optional:
     -template        Path, including file name, to iwgc_circos_template.config (default: ./iwgc_circos_template.config)
     -outdir          Path to output directory (default: ./iwgc_circos)
     -ideogram        Path, indluding file name, to ideogram.conf (default: $OUTDIR/iwgc_circos.conf)
     -gap             Adds gap to middle-top of plot, between last and first chromosomes (default: FALSE)
     -h | --help      List usage options
```

  
***Run in ```iwgc-circos```, one level up from ```iwgc_circos``` and ```iwgc_circos_data```***  
    
**Singularity:**
```
singularity exec iwgc-circos-tracks.sif ./create_configs.sh -gap
```
  
**Docker:**
```
docker run --rm -v "$PWD":/data -w /data \
  scrumpis/iwgc-circos-tracks:latest \
  ./create_configs.sh -gap
```
  
    
## 3. Circos Plot
The below commands will generate a Circos plot in iwgc_circos/tmp/iwgc_circos.png (and .svg).  
   
***Run in ```iwgc-circos```, one level up from ```iwgc_circos``` and ```iwgc_circos_data```***       
   
**Singularity:**
```
singularity exec \
  --cleanenv \
  -B "$PWD":/data \
  iwgc-circos-tracks.sif \
  circos -conf /data/iwgc_circos/iwgc_circos.conf -outputdir /data/iwgc_circos/tmp --noparanoid
```
  
**Docker** (easier refinement - quickly view png/svg, tweak .conf files, and re-run):
```
docker run --rm \
  -v "$PWD:/data" \
  scrumpis/iwgc-circos-tracks \
  circos -conf /data/iwgc_circos/iwgc_circos.conf -outputdir /data/iwgc_circos/tmp -noparanoid
```
  
  
## 4. Final Touches (optional)
### Add legend characters to the central gap (add_legend.sh)
Dynamically add legend characters (a., b., c., etc.) for each track present, centered vertically within tracks and horizontally within the gap (from ideogram coords).   
***Note: You may want to adjust ```angle_offset* =``` in ```iwgc_circos.conf``` if the gap looks off-center prior to adding legend characters.***  
  
```
  Usage: ./create_configs.sh [options]
  
   Optional:
     -conf            Path, including file name, to iwgc_circos.conf (default: iwgc_circos/iwgc_circos.conf)
     -svg             Path, including file name, to iwgc_circos.svg (default: iwgc_circos/tmp/iwgc_circos.svg)
     -out             Path, indluding file name, to output SVG file (default: iwgc_circos/tmp/iwgc_circos.labeled.svg)
     -font            Path, indluding file name, to legend font file (default: fonts/ArialBold.ttf)
     --dx-frac        Nudge legend characters left (+) or right (-) in fractions of total legend character size (default: 0)
     --dy-frac        Nudge legend characters up (+) or down (-) in fractions of total legend character size (default: 0)
     --theta-deg      Pivot characters around the center of the Circos plot (default: -90)
     -h | --help      List usage options
```
  
***Run in ```iwgc-circos```, one level up from ```iwgc_circos``` and ```iwgc_circos_data```***  
  
**Singularity:**
```
singularity exec iwgc-circos-tracks.sif ./add_legend.sh
```

**Docker:**
```
docker run --rm -v "$PWD":/data -w /data \
  scrumpis/iwgc-circos-tracks:latest \
  ./add_legend
```

### Manual Adjustments
If you are planning to include the Circos plot in a publication, you will likely need to make a few manual adjustments to the .config files in ```iwgc_circos/```. After making these config file adjustments, re-run Circos to generate a new plot with your updates. Below are some commonly made adjustments.  
  
**iwgc_circos.conf**
- Reduce ```chromosomes_units =``` for very small genomes or increase for very big ones. Default 1000000 (1Mbp) should cover a broad size range. You will likely have to adjust ticks.conf if you change this.
- Change r1 and r0 values within each ```<plot>``` block to adjust height/thickness of individual plot tracks, the gaps between each track, or to change the order of tracks.
- Change [plot type](https://circos.ca/documentation/tutorials/2d_tracks/), [colors of plots](https://circos.ca/documentation/tutorials/configuration/colors/images), plot background colors, reverse chromosomes, insert breaks into chromsomes, etc.
- Uncomment ```chromosome_reverse =``` and add chromosomes in a list (Chr09D; Chr01B; Chr04C) to reverse chromosomes and all associated plots.
- Syntenic links are very customizable. Change thresholds for sizes to hide or color differently, allow intrachromosomal links, add twists, flatten, etc. Defaults: 20Kbp<50Kbp blue, 50kbp<100Kbp green, 100Kbp>200Kbp yellow, 200Kbp>1Mbp orange, 1Mbp> red. 

**ideogram.label.conf**
- Increase "150p" in ```label_radius``` to move chromosome labels further from ideogram, or decrease to bring them closer.
- Adjust "chr0*" regex in ```label_format``` to match patterns you want to remove from labels. It currently removes Chr/chr/Chr0/chr0 from the beginning of all chromosomes, so Chr09B becomes 9B.

**ticks.conf**
- If you want more or less ticks, different sizes of ticks, labels or no labels, etc.
  
Use Inkspace or a similar SVG editing tool to manually edit or add any text.  
  
   
## Notes
- When visualizing multiple species, datasets, etc. in a single Circos plot, all files must be concatenated first
- For automation and visualization purposes...
  - Gene desity is sqrt transformed
  - Repeat density is power 3 transformed
  - All tracks are normalized
- In special situtations where you need to edit the housekeeping.conf, such as if you have more contigs than what is typically allowed, you can either run the container interactively and update the housekeeping.conf file or you can create a housekeeping.conf file in the /circos directory and map it to the container, see Fusarium paper

### Known Warnings:
- (iwgc_circos_tracks.sh) ```File species.fasta_windows.bed has a record where naming convention (leading zero) is inconsistent with other files```
  - Seems ignorable with no effects, happens at Chr10 or ChrChlr if previous are Chr01-09.
- (Circos) ```Use of uninitialized value in subroutine entry at /opt/circos/bin/../lib/Circos/Configuration.pm line 781```
  - Seems ignorable.
