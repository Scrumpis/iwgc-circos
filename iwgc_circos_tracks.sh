#!/bin/bash
# Quickly generate Circos plots for commonly visualized genomic features

usage() {
  echo "Usage: $0 [options] <FASTA> [options]"
  echo
  echo "Required:"
  echo "  FASTA            Genomic FASTA file"
  echo
  echo "Optional:"
  echo "  -gene            Add gene density track (requires gene annotation GFF3)"
  echo "  -repeat          Add repeat density track (requires EDTA repeat annotation: EDTA/genome.fasta.EDTA.mod.TEanno.gff3)"
  echo "  -intact          Add intact TE density track (requires EDTA intact repeat annotation: EDTA/genome.fasta.mod.EDTA.intact.gff3)"
  echo "  -ltr-dating      Add LTR dating track (requires EDTA repeat annotation: EDTA/genome.mod.EDTA.raw/LTR/genome.fasta.mod.pass.list)"
  echo "  -gc              Add GC content track"
  echo "  -telomere        Add telomere bands to ideogram (karyotype.circos)"
  echo "  -ts <value>      Telomere band size scale (default: 0.005), 0.5% of total genome size"
  echo "  -window          Window size in base pairs (default: 300000)"
  echo "  -sliding         Use sliding windows instead of fixed"
  echo "  -step <value>    Step size for sliding windows (default: 0.5). The default is half window size steps"
  echo "  -filter-chrs     Restrict chromosomes to those matching typical nuclear naming patterns (e.g., Chr01, Chr1, chr01B). Default: off"
  echo "  -keep-temp       Keep intermediate files"
  echo "  -h | --help      List usage options"
  echo
  echo "Recommended containerized usage:"
  echo "singularity exec iwgc-circos-tracks.sif ./iwgc-circos-tracks.sh <FASTA> [options]"
  echo "or"
  echo "docker run --rm -v $(pwd):/data iwgc-circos-tracks:latest ./iwgc-circos-tracks.sh <FASTA> [options]"
  echo
  echo "Note: For visualization/automation purposes, all density tracks are normalized, gene density is sqrt transformed, and repeat density is power 3 transforme
d"
  exit 1
}

INCLUDE_GENE=false
INCLUDE_REPEAT=false
INCLUDE_INTACT=false
INCLUDE_GC=false
INCLUDE_LTRDATING=false
INCLUDE_TELOMERE=false
KEEP_TEMP=false
TEMP_FILES=()
TELOMERE_SCALE=0.005
USE_SLIDING=false
STEP_SIZE=""
FILTER_CHRS=false

FASTA=$1; shift

while [[ $# -gt 0 ]]; do
  case "$1" in
    -gene) INCLUDE_GENE=true; GENES=$2; shift 2 ;;
    -repeat) INCLUDE_REPEAT=true; REPEATS=$2; shift 2 ;;
    -intact) INCLUDE_INTACT=true; INTACT=$2; shift 2 ;;
    -ltr-dating) INCLUDE_LTRDATING=true; LTRDATES=$2; shift 2 ;;
    -window) WINDOW=$2; shift 2 ;;
    -ts) TELOMERE_SCALE=$2; shift 2 ;;
    -step) STEP_SIZE=$2; shift 2 ;;
    -gc) INCLUDE_GC=true; shift ;;
    -telomere) INCLUDE_TELOMERE=true; shift ;;
    -sliding) USE_SLIDING=true; shift ;;
    -keep-temp) KEEP_TEMP=true; shift ;;
    -filter-chrs) FILTER_CHRS=true; shift ;;
    -h|--help) usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

# Default window and step sizes
WINDOW=${WINDOW:-300000}
STEP_SIZE=${STEP_SIZE:-0.5}

# Grab FASTA basename for naming output files
FASTA_BASE=$(basename "$FASTA")

# Check if required files are provided
[[ ! -f $FASTA ]] && { echo "Missing genome FASTA"; exit 1; }
[[ $INCLUDE_GENE == true && ! -f $GENES ]] && { echo "Gene track requires a GFF file"; exit 1; }
[[ $INCLUDE_REPEAT == true && ! -f $REPEATS ]] && { echo "Repeat track requires a GFF file"; exit 1; }
[[ $INCLUDE_INTACT == true && ! -f $INTACT ]] && { echo "Intact TE track requires a GFF file"; exit 1; }
[[ $INCLUDE_LTRDATING == true && ! -f $LTRDATES ]] && { echo "LTR dating track requires a coordinate file"; exit 1; }

# Create FASTA index if it doesn't exist
if [[ ! -f ${FASTA}.fai ]]; then
  samtools faidx ${FASTA}
fi

# Filter for chromosomes if requested
if [[ "$FILTER_CHRS" == true ]]; then
  grep -Ei 'Chr0?[1-9][0-9]?' ${FASTA}.fai | sort -t$'\t' -k1,1 > ${FASTA_BASE}_chrs.fai
else
  sort -t$'\t' -k1,1 ${FASTA}.fai > ${FASTA_BASE}_chrs.fai
fi

# Ensure the FASTA index contains entries
if [[ ! -s ${FASTA_BASE}_chrs.fai ]]; then
  echo "Error: No valid chromosomes or other contigs found in ${FASTA}.fai"
  exit 1
fi
TEMP_FILES+=("${FASTA_BASE}_chrs.fai")


# Generate karyotype file for Circos
cut -f1,2 ${FASTA_BASE}_chrs.fai | sort -t$'\t' | awk 'BEGIN {
  idx=0;
  for(i=2;i<=21;i++) l[idx++]="chr"i;
  for(i=2;i<=21;i++) l[idx++]="lum70chr"i;
  for(i=2;i<=21;i++) l[idx++]="lum90chr"i;
} {
  chr=$1;
  gsub(/^[Cc][Hh][Rr]0?/, "", chr);
  print "chr", "-", $1, chr, 0, $2, l[NR-1]
}' > ${FASTA_BASE}_karyotype.circos


# Add telomere bands to karyotype file if requested
if [[ $INCLUDE_TELOMERE == true ]]; then
  TELO_TMP="${FASTA_BASE}_telomere_tmp"
  mkdir -p "$TELO_TMP"
  TEMP_FILES+=("$TELO_TMP" "${FASTA_BASE}.telo.info")

  # Run quartet *inside* the temp dir
  (cd "$TELO_TMP" && quartet.py te -i ../${FASTA} -c plant -m 50 -p ../${FASTA_BASE} --noplot)

  genome_len=$(awk '{sum+=$2}END{print sum}' ${FASTA_BASE}_chrs.fai)

  if [[ "$FILTER_CHRS" == true ]]; then
    grep -v '#' ${FASTA_BASE}.telo.info | grep -Ei 'Chr0?[1-9][0-9]?'
  else
    grep -v '#' ${FASTA_BASE}.telo.info
  fi | awk '{ print $1, $2, $3, $4, "+", $6, "-" }' | \
  awk -v glen="$genome_len" '{
    t=int(glen*'$TELOMERE_SCALE');
    if($3=="both"){print $1,0,t; print $1,$2-t,$2}
    else if($3=="right"){print $1,$2-t,$2}
    else if($3=="left"){print $1,0,t}
  }' | \
  awk '{ print "band", $1, $1"_T"NR, $1"_T"NR, $2, $3, "vdgrey" }' > ${FASTA_BASE}_telomere_bands.bed

  cat ${FASTA_BASE}_telomere_bands.bed >> ${FASTA_BASE}_karyotype.circos
  TEMP_FILES+=("${FASTA_BASE}_telomere_bands.bed")
fi


# Generate sliding or fixed window file if any tracks other than karyotype are requested
if [[ $INCLUDE_GENE == true || $INCLUDE_REPEAT == true || $INCLUDE_INTACT == true || $INCLUDE_GC == true || $INCLUDE_LTRDATING == true ]]; then
  if [[ "$USE_SLIDING" == true ]]; then
    STEP_FRACTION="${STEP_SIZE:-0.5}"

# Validate step size is a float between 0 and 1
    if ! [[ "$STEP_FRACTION" =~ ^0(\.[0-9]+)?$|^1(\.0+)?$ ]]; then
        echo "Error: -step must be a decimal fraction between 0 and 1"
        exit 1
    fi

    STEP=$(awk -v frac="$STEP_FRACTION" -v win="$WINDOW" 'BEGIN { printf "%.0f", frac * win }')
    echo "Generating sliding windows with size $WINDOW and step $STEP"
    bedtools makewindows -g ${FASTA_BASE}_chrs.fai -w ${WINDOW} -s ${STEP} | \
    awk -v W=${WINDOW} '{ if (($3 - $2) >= (W / 2)) print }' > ${FASTA_BASE}_windows.bed
  else
    echo "Generating fixed windows with size $WINDOW"
    bedtools makewindows -g ${FASTA_BASE}_chrs.fai -w ${WINDOW} > ${FASTA_BASE}_windows.bed
  fi
  TEMP_FILES+=("${FASTA_BASE}_windows.bed")
fi


# Generate gene density track file if requested
if [[ $INCLUDE_GENE == true ]]; then
  grep 'gene' ${GENES} | awk '{print $1,$4,$5}' OFS='\t' | sort -k1,1 -k2,2n > ${FASTA_BASE}_genes_coords.bed
  TEMP_FILES+=("${FASTA_BASE}_genes_coords.bed")
  bedtools coverage -sorted -a ${FASTA_BASE}_windows.bed -b ${FASTA_BASE}_genes_coords.bed | \
  awk '{print $1,$2,$3,sqrt($7)}' | sort -k1,1 -k2,2n | \
  awk '{v[NR]=$4; l[NR]=$0; if(NR==1||$4<m)m=$4; if(NR==1||$4>M)M=$4} END {for(i=1;i<=NR;i++) print l[i],(M==m?0:(v[i]-m)/(M-m))}' | \
  awk '{ print $1, $2, $3, $5}' > ${FASTA_BASE}_gene_coverage.circos
fi


# Generate repeat density track file if requested
if [[ $INCLUDE_REPEAT == true ]]; then
  grep -v '#' ${REPEATS} | awk '{print $1,$4,$5}' OFS='\t' | sort -k1,1 -k2,2n > ${FASTA_BASE}_repeat_coords.bed
  TEMP_FILES+=("${FASTA_BASE}_repeat_coords.bed")
  bedtools coverage -sorted -a ${FASTA_BASE}_windows.bed -b ${FASTA_BASE}_repeat_coords.bed | \
  awk '{print $1,$2,$3,($7^3)}' | sort -k1,1 -k2,2n | \
  awk '{v[NR]=$4; l[NR]=$0; if(NR==1||$4<m)m=$4; if(NR==1||$4>M)M=$4} END {for(i=1;i<=NR;i++) print l[i],(M==m?0:(v[i]-m)/(M-m))}' | \
  awk '{ print $1, $2, $3, $5}' > ${FASTA_BASE}_repeat_coverage.circos
fi


# Generate intact TE coverage track files if requested
if [[ $INCLUDE_INTACT == true ]]; then
  awk -F'\t' -v base="${FASTA_BASE}" '
    match($9, /Classification=([^;]+)/, a) {
      split(a[1], p, "/");
      if (p[2] != "") {
        print > (base "_intactTE_" p[2] ".gff")
      }
    }' "$INTACT"

  output_files=$(ls ${FASTA_BASE}_intactTE_*.gff 2>/dev/null)
  TEMP_FILES+=($output_files)

  for gff in $output_files; do
    class=$(echo "$gff" | sed -E "s/${FASTA_BASE}_intactTE_(.*)\.gff/\1/")
    coords="${FASTA_BASE}_${class}_intactTE_coords.bed"
    awk '{print $1,$4,$5}' OFS='\t' "$gff" | sort -k1,1 -k2,2n > "$coords"
    TEMP_FILES+=("$coords")
    bedtools coverage -sorted -a "${FASTA_BASE}_windows.bed" -b "$coords" | \
    awk '{print $1,$2,$3,$7}' | sort -k1,1 -k2,2n | \
    awk '{v[NR]=$4; l[NR]=$0; if(NR==1||$4<m)m=$4; if(NR==1||$4>M)M=$4}
         END {for(i=1;i<=NR;i++) print l[i],(M==m?0:(v[i]-m)/(M-m))}' | \
    awk '{ print $1, $2, $3, $5 }' > "${FASTA_BASE}_${class}_intactTE_coverage.circos"
  done
fi


# Generate GC content track file if requested
if [[ $INCLUDE_GC == true ]]; then
  bedtools nuc -fi ${FASTA} -bed ${FASTA_BASE}_windows.bed > gc_content.bed
  TEMP_FILES+=("gc_content.bed")
  awk '{print $1,$2,$3,$5}' gc_content.bed | grep -v '#' | sort -k1,1 -k2,2n | \
  awk '{v[NR]=$4; l[NR]=$0; if(NR==1 || $4<m) m=$4; if(NR==1 || $4>M) M=$4} 
       END {for(i=1;i<=NR;i++) print l[i], (M==m ? 0 : (v[i]-m)/(M-m))}' | \
  awk '{print $1, $2, $3, $5}' > ${FASTA_BASE}_gc.circos
fi


# Generate LTR dating track file if requested
if [[ $INCLUDE_LTRDATING == true ]]; then
  awk -F'\t' 'BEGIN { OFS="\t" } !/^#/ { print $1, $12 }' ${LTRDATES} | sed -e 's/:/\t/g' -e 's/\.\./\t/g' > ${FASTA_BASE}_LTR_insertion.bed
  bedtools map -a ${FASTA_BASE}_windows.bed -b ${FASTA_BASE}_LTR_insertion.bed -c 4 -o mean > ${FASTA_BASE}_LTR_age.bed

  awk '$4 != "." {print $1,$2,$3,int($4 + 0.5)}' ${FASTA_BASE}_LTR_age.bed | sort -k1,1 -k2,2n | \
  awk '{v[NR]=$4; l[NR]=$0; if(NR==1||$4<m)m=$4; if(NR==1||$4>M)M=$4}
        END {for(i=1;i<=NR;i++) print l[i],(M==m?0:(v[i]-m)/(M-m))}' | \
  awk '{print $1, $2, $3, $5}' > ${FASTA_BASE}_LTR_age_normal.circos
  TEMP_FILES+=("${FASTA_BASE}_LTR_insertion.bed" "${FASTA_BASE}_LTR_age.bed")
fi


# Remove temporary files (default behavior)
if [[ $KEEP_TEMP == false ]]; then
  echo "Cleaning up intermediate files..."
  for file in "${TEMP_FILES[@]}"; do
    if [[ -f "$file" ]]; then
      rm -f "$file"
    elif [[ -d "$file" ]]; then
      rm -rf "$file"
    fi
  done
fi
