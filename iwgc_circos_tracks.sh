#!/bin/bash
# Generate Circos for IWGC Genome Reports with optional tracks and a single required SIF

usage() {
  echo "Usage: $0 [options] <FASTA> [<GENES>] [<REPEATS>] [<INTACT>] [<LTRDATES>] <IWGC_CIRCOS_SIF> [WINDOW]"
  echo "Options:"
  echo "  -gene            Add gene density track (requires gene annotation GFF3)"
  echo "  -repeat          Add repeat density track (requires EDTA repeat annotation: EDTA/genome.fasta.EDTA.TEanno.gff3)"
  echo "  -intact          Add intact TE density track (requires EDTA intact repeat annotation: EDTA/genome.fasta.mod.EDTA.intact.gff3)"
  echo "  -gc              Add GC content track"
  echo "  -ltr-dating      Add LTR dating track (requires EDTA repeat annotation: EDTA/genome.mod.EDTA.raw/LTR/genome.fasta.mod.pass.list)"
  echo "  -telomere        Add telomere bands to ideogram (karyotype.circos)"
  echo "  -ts <value>      Telomere band size scale (default: 0.005)"
  echo "  -keep-temp       Keep intermediate files"
  echo "  -sliding         Use sliding windows instead of fixed"
  echo "  -step <value>    Step size for sliding windows (default: WINDOW/2)"
  echo "  -local           Use Docker instead of Singularity (default: false)"
  echo
  echo "Note: For visualization purposes, all density tracks are normalized, gene desity is sqrt transformed, and repeat density is power 3 transformed"
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
USE_DOCKER=false  # Default to Singularity

# Determine how to run container commands
run_cmd() {
  if [[ "$USE_DOCKER" == true ]]; then
    docker run --rm -v "$PWD:/data" -w /data "$IWGC_CIRCOS_SIF" "$@"
  else
    singularity exec "$IWGC_CIRCOS_SIF" "$@"
  fi
}


while [[ $1 == -* ]]; do
  case "$1" in
    -gene) INCLUDE_GENE=true; shift ;;
    -repeat) INCLUDE_REPEAT=true; shift ;;
    -intact) INCLUDE_INTACT=true; shift ;;
    -gc) INCLUDE_GC=true; shift ;;
    -ltr-dating) INCLUDE_LTRDATING=true; shift ;;
    -telomere) INCLUDE_TELOMERE=true; shift ;;
    -keep-temp) KEEP_TEMP=true; shift ;;
    -ts) shift; TELOMERE_SCALE="$1"; shift ;;
    -sliding) USE_SLIDING=true; shift ;;
    -step) shift; STEP_SIZE="$1"; shift ;;
    -local) USE_DOCKER=true; shift ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

FASTA=$1; shift
[[ $INCLUDE_GENE == true ]] && { GENES=$1; shift; }
[[ $INCLUDE_REPEAT == true ]] && { REPEATS=$1; shift; }
[[ $INCLUDE_INTACT == true ]] && { INTACT=$1; shift; }
[[ $INCLUDE_LTRDATING == true ]] && { LTRDATES=$1; shift; }
IWGC_CIRCOS_SIF=$1; shift
WINDOW=${1:-300000}

[[ ! -f $FASTA ]] && { echo "Missing genome FASTA"; exit 1; }
if [[ "$USE_DOCKER" == false && ! -f $IWGC_CIRCOS_SIF ]]; then
  echo "Missing required iwgc-circos-tracks.sif for Singularity"
  exit 1
fi
if [[ "$USE_DOCKER" == true ]]; then
  if ! docker image inspect "$IWGC_CIRCOS_SIF" > /dev/null 2>&1; then
    echo "Docker image '$IWGC_CIRCOS_SIF' not found locally. Please pull or build it first."
    exit 1
  fi
fi
[[ $INCLUDE_GENE == true && ! -f $GENES ]] && { echo "Gene track requires a GFF file"; exit 1; }
[[ $INCLUDE_REPEAT == true && ! -f $REPEATS ]] && { echo "Repeat track requires a GFF file"; exit 1; }
[[ $INCLUDE_INTACT == true && ! -f $INTACT ]] && { echo "Intact TE track requires a GFF file"; exit 1; }
[[ $INCLUDE_LTRDATING == true && ! -f $LTRDATES ]] && { echo "LTR dating track requires a coordinate file"; exit 1; }

if [[ ! -f ${FASTA}.fai ]]; then
  run_cmd samtools faidx ${FASTA}
fi

grep -Ei 'Chr0?[1-9][0-9]?' ${FASTA}.fai | sort -t$'\t' -k1,1 > ${FASTA}_chrs.fai
TEMP_FILES+=("${FASTA}_chrs.fai")

cut -f1,2 ${FASTA}_chrs.fai | sort -t$'\t' | awk 'BEGIN {
  idx=0;
  for(i=2;i<=21;i++) l[idx++]="chr"i;
  for(i=2;i<=21;i++) l[idx++]="lum70chr"i;
  for(i=2;i<=21;i++) l[idx++]="lum90chr"i;
} {
  chr=$1;
  gsub(/^[Cc][Hh][Rr]0?/, "", chr);
  print "chr", "-", $1, chr, 0, $2, l[NR-1]
}' > ${FASTA}_karyotype.circos

if [[ $INCLUDE_TELOMERE == true ]]; then
  run_cmd quartet.py te -i ${FASTA} -c plant -m 50 -p ${FASTA} --noplot
  genome_len=$(awk '{sum+=$2}END{print sum}' ${FASTA}_chrs.fai)
  grep -v '#' ${FASTA}.telo.info | grep -Ei 'Chr0?[1-9][0-9]?' | awk '{ print $1, $2, $3, $4, "+", $6, "-" }' | \
  awk -v glen="$genome_len" '{
    t=int(glen*'$TELOMERE_SCALE');
    if($3=="both"){print $1,0,t; print $1,$2-t,$2}
    else if($3=="right"){print $1,$2-t,$2}
    else if($3=="left"){print $1,0,t}
  }' | \
  awk '{ print "band", $1, $1"_T"NR, $1"_T"NR, $2, $3, "vdgrey" }' > ${FASTA}_telomere_bands.bed
  cat ${FASTA}_telomere_bands.bed >> ${FASTA}_karyotype.circos
  TEMP_FILES+=("${FASTA}_telomere_bands.bed")
fi

if [[ $INCLUDE_GENE == true || $INCLUDE_REPEAT == true || $INCLUDE_INTACT == true || $INCLUDE_GC == true || $INCLUDE_LTRDATING == true ]]; then
  if [[ "$USE_SLIDING" == true ]]; then
    STEP="${STEP_SIZE:-$((WINDOW / 2))}"
    echo "Generating sliding windows with size $WINDOW and step $STEP"
    run_cmd bedtools makewindows -g ${FASTA}_chrs.fai -w ${WINDOW} -s ${STEP} | \
    awk -v W=${WINDOW} '{ if (($3 - $2) >= (W / 2)) print }' > ${FASTA}_windows.bed
  else
    echo "Generating fixed windows with size $WINDOW"
    run_cmd bedtools makewindows -g ${FASTA}_chrs.fai -w ${WINDOW} > ${FASTA}_windows.bed
  fi
  TEMP_FILES+=("${FASTA}_windows.bed")
fi

if [[ $INCLUDE_GENE == true ]]; then
  grep 'gene' ${GENES} | grep -Ei 'Chr0?[1-9][0-9]?' | awk '{print $1,$4,$5}' OFS='\t' | sort -k1,1 -k2,2n > ${GENES}_genes_coords.bed
  TEMP_FILES+=("${GENES}_genes_coords.bed")
  run_cmd bedtools coverage -sorted -a ${FASTA}_windows.bed -b ${GENES}_genes_coords.bed | \
  awk '{print $1,$2,$3,sqrt($7)}' | sort -k1,1 -k2,2n | \
  awk '{v[NR]=$4; l[NR]=$0; if(NR==1||$4<m)m=$4; if(NR==1||$4>M)M=$4} END {for(i=1;i<=NR;i++) print l[i],(M==m?0:(v[i]-m)/(M-m))}' | \
  awk '{ print $1, $2, $3, $5}' > ${FASTA}_gene_coverage.circos
fi

if [[ $INCLUDE_REPEAT == true ]]; then
  grep -v '#' ${REPEATS} | grep -Ei 'Chr0?[1-9][0-9]?' | awk '{print $1,$4,$5}' OFS='\t' | sort -k1,1 -k2,2n > ${FASTA}_repeat_coords.bed
  TEMP_FILES+=("${FASTA}_repeat_coords.bed")
  run_cmd bedtools coverage -sorted -a ${FASTA}_windows.bed -b ${FASTA}_repeat_coords.bed | \
  awk '{print $1,$2,$3,($7^3)}' | sort -k1,1 -k2,2n | \
  awk '{v[NR]=$4; l[NR]=$0; if(NR==1||$4<m)m=$4; if(NR==1||$4>M)M=$4} END {for(i=1;i<=NR;i++) print l[i],(M==m?0:(v[i]-m)/(M-m))}' | \
  awk '{ print $1, $2, $3, $5}' > ${FASTA}_repeat_coverage.circos
fi

if [[ $INCLUDE_INTACT == true ]]; then
  rm -f output_*.gff
  awk -F'\t' 'match($9, /Classification=([^;]+)/, a) {split(a[1],p,"/"); if(p[2] != "") {print > "output_"p[2]".gff"}}' $INTACT
  output_files=$(ls output_*.gff 2>/dev/null)
  TEMP_FILES+=($output_files)
  for gff in $output_files; do
    class=$(basename $gff .gff | cut -d_ -f2)
    coords="${FASTA}_${class}_intactTE_coords.bed"
    grep -Ei 'Chr0?[1-9][0-9]?' $gff | awk '{print $1,$4,$5}' OFS='\t' | sort -k1,1 -k2,2n > $coords
    TEMP_FILES+=("$coords")
    run_cmd bedtools coverage -sorted -a ${FASTA}_windows.bed -b $coords | \
    awk '{print $1,$2,$3,$7}' | sort -k1,1 -k2,2n | \
    awk '{v[NR]=$4; l[NR]=$0; if(NR==1||$4<m)m=$4; if(NR==1||$4>M)M=$4} END {for(i=1;i<=NR;i++) print l[i],(M==m?0:(v[i]-m)/(M-m))}' | awk '{ print $1, $2, $3, $5
 }' > "${FASTA}_${class}_intactTE_coverage.circos"
  done
fi

if [[ $INCLUDE_GC == true ]]; then
  run_cmd bedtools nuc -fi ${FASTA} -bed ${FASTA}_windows.bed > gc_content.bed
  TEMP_FILES+=("gc_content.bed")
  awk '{print $1,$2,$3,$5}' gc_content.bed | grep -v '#' | sort -k1,1 -k2,2n | \
  awk '{v[NR]=$4; l[NR]=$0; if(NR==1 || $4<m) m=$4; if(NR==1 || $4>M) M=$4} 
       END {for(i=1;i<=NR;i++) print l[i], (M==m ? 0 : (v[i]-m)/(M-m))}' | \
  awk '{print $1, $2, $3, $5}' > ${FASTA}_gc.circos
fi

if [[ $INCLUDE_LTRDATING == true ]]; then
  awk -F'\t' 'BEGIN { OFS="\t" } !/^#/ { print $1, $12 }' ${LTRDATES} | sed -e 's/:/\t/g' -e 's/\.\./\t/g' > ${FASTA}_LTR_insertion.bed
  run_cmd bedtools map -a ${FASTA}_windows.bed -b ${FASTA}_LTR_insertion.bed -c 4 -o mean > ${FASTA}_LTR_age.bed

  awk '$4 != "." {print $1,$2,$3,int($4 + 0.5)}' ${FASTA}_LTR_age.bed | sort -k1,1 -k2,2n | \
  awk '{v[NR]=$4; l[NR]=$0; if(NR==1||$4<m)m=$4; if(NR==1||$4>M)M=$4}
        END {for(i=1;i<=NR;i++) print l[i],(M==m?0:(v[i]-m)/(M-m))}' | \
  awk '{print $1, $2, $3, $5}' > ${FASTA}_LTR_age_normal.circos
  TEMP_FILES+=("${FASTA}_LTR_insertion.bed" "${FASTA}_LTR_age.bed")
fi

if [[ $KEEP_TEMP == false ]]; then
  echo "Cleaning up intermediate files..."
  for file in "${TEMP_FILES[@]}"; do
    [[ -f "$file" ]] && rm -f "$file"
  done
fi
