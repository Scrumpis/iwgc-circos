#!/bin/bash
set -euo pipefail
export LC_ALL=C
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# Quickly generate Circos plots for commonly visualized genomic features

usage() {
  echo "Usage: $0 <FASTA> [options]"
  echo
  echo "Required:"
  echo "  FASTA            Genomic FASTA file"
  echo
  echo "Optional:"
  echo "  -gene            Add gene density track (requires gene annotation GFF3)"
  echo "  -repeat          Add repeat density track (requires EDTA repeat annotation: EDTA/genome.fasta.EDTA.mod.TEanno.gff3)"
  echo "  -intact          Add intact TE density track (requires EDTA intact repeat annotation: EDTA/genome.fasta.mod.EDTA.intact.gff3)"
  echo "  -ltr-dating      Add LTR dating track (requires EDTA repeat annotation: EDTA/genome.mod.EDTA.raw/LTR/genome.fasta.mod.pass.list)"
  echo "  -links           Add links from 6-col BED coordinates or headerless 11/13-col MUMmer"
  echo "  -min-tlen        Minimum summed length of both syntenic intervals (bp)"
  echo "  -top-n           Number of top syntenic links to include (e.g., \"20\" for top 20 longest links)"
  echo "  -link-order      Order of syntenic links - asc: smallest top to biggest bottom | dsc: big top to small bottom (default: asc)"
  echo "  -gc              Add GC content track"
  echo "  -telomere        Add telomere bands to ideogram (karyotype.circos)"
  echo "  -ts              Telomere band size scale (default: 0.005), 0.5% of total genome size"
  echo "  -window          Window size in base pairs (default: 1000000)"
  echo "  -sliding         Use sliding windows instead of the default fixed windows"
  echo "  -step            Step size for sliding windows (default: 0.5). The default is half window size steps"
  echo "  -filter-chrs [PATTERN]"
  echo "                   Restrict chromosomes to those whose names match PATTERN."
  echo "                   If PATTERN is omitted, defaults to typical nuclear naming"
  echo "                   (e.g., Chr01, Chr1, chr01B; regex: Chr0?[1-9][0-9]?)."
  echo "                   Matching is case-insensitive; PATTERN is an extended regex"
  echo "                   applied to the chromosome name (first column). Default: off"
  echo "  -keep-temp       Keep intermediate files"
  echo "  -out             Output directory for Circos track files (default: current directory)"
  echo "  -h | --help      List usage options"
  echo
  echo "Recommended containerized usage:"
  echo "singularity exec iwgc-circos.sif ./iwgc_circos_tracks.sh <FASTA> [options]"
  echo "or"
  echo "docker run --rm -v \$(pwd):/data -w /data iwgc-circos:candidate bash ./iwgc_circos_tracks.sh <FASTA> [options]"
  echo
  echo "Note: For visualization/automation purposes, all density tracks are normalized, gene density is sqrt transformed, and repeat density is sqrt transformed"
  exit "${1:-0}"
}

INCLUDE_GENE=false
INCLUDE_REPEAT=false
INCLUDE_INTACT=false
INCLUDE_GC=false
INCLUDE_LTRDATING=false
INCLUDE_LINKS=false
LINKS_FILE=""
INCLUDE_TELOMERE=false
KEEP_TEMP=false
TEMP_FILES=()
TRACK_FILES=()
TELOMERE_SCALE=0.005
USE_SLIDING=false
STEP_SIZE=""
FILTER_CHRS=false
### New: optional pattern used when FILTER_CHRS is true
FILTER_PATTERN=""

if [[ $# -eq 0 ]]; then usage 1; fi
if [[ "$1" == -h || "$1" == --help ]]; then usage; fi
FASTA=$1; shift

while [[ $# -gt 0 ]]; do
  case "$1" in
    -gene|-repeat|-intact|-ltr-dating|-links|-min-tlen|-top-n|-link-order|-window|-ts|-step|-out)
      [[ $# -ge 2 && "$2" != -* ]] || { echo "[ERROR] $1 requires a value" >&2; exit 1; };;
  esac
  case "$1" in
    -gene) INCLUDE_GENE=true; GENES=$2; shift 2 ;;
    -repeat) INCLUDE_REPEAT=true; REPEATS=$2; shift 2 ;;
    -intact) INCLUDE_INTACT=true; INTACT=$2; shift 2 ;;
    -ltr-dating) INCLUDE_LTRDATING=true; LTRDATES=$2; shift 2 ;;
    -links) INCLUDE_LINKS=true; LINKS_FILE=$2; shift 2 ;;
    -min-tlen) MIN_TLEN=$2; shift 2 ;;
    -top-n) TOP_N=$2; shift 2 ;;
    -link-order) LINK_ORDER=$2; shift 2 ;;
    -window) WINDOW=$2; shift 2 ;;
    -ts) TELOMERE_SCALE=$2; shift 2 ;;
    -step) STEP_SIZE=$2; shift 2 ;;
    -gc) INCLUDE_GC=true; shift ;;
    -telomere) INCLUDE_TELOMERE=true; shift ;;
    -sliding) USE_SLIDING=true; shift ;;
    -keep-temp) KEEP_TEMP=true; shift ;;
    ### Updated: -filter-chrs can take an optional PATTERN argument
    -filter-chrs)
      FILTER_CHRS=true
      # If next argument exists and is not another flag, treat it as PATTERN
      if [[ -n "${2-}" && ! "$2" =~ ^- ]]; then
        FILTER_PATTERN="$2"
        shift 2
      else
        FILTER_PATTERN=""
        shift
      fi
      ;;
    -out) OUTPUT_DIR=$2; shift 2 ;;
    -h|--help) usage ;;
    *) echo "Unknown option: $1" >&2; usage 1 ;;
  esac
done

# Default window and step sizes
WINDOW=${WINDOW:-1000000}
STEP_SIZE=${STEP_SIZE:-0.5}
OUTPUT_DIR=${OUTPUT_DIR:-.}     # current directory by default
MIN_TLEN=${MIN_TLEN:-0}
TOP_N=${TOP_N:-0}
LINK_ORDER=${LINK_ORDER:-asc}   # asc|desc

# Validate before starting external tools or creating plausible outputs.
for value in "$WINDOW" "$MIN_TLEN" "$TOP_N"; do
  [[ "$value" =~ ^[0-9]+$ ]] || { echo "[ERROR] window/min-tlen/top-n must be integers" >&2; exit 1; }
done
WINDOW=$((10#$WINDOW))
(( WINDOW > 0 )) || { echo "[ERROR] window must be positive" >&2; exit 1; }
[[ "$LINK_ORDER" == asc || "$LINK_ORDER" == dsc ]] || { echo "[ERROR] link-order must be asc or dsc" >&2; exit 1; }
gawk -v x="$TELOMERE_SCALE" 'BEGIN{exit !(x ~ /^[0-9]+([.][0-9]+)?$/ && x>0 && x<=1)}' || { echo "[ERROR] ts must be in (0,1]" >&2; exit 1; }
for dependency in samtools python3 gawk; do
  command -v "$dependency" >/dev/null || { echo "[ERROR] Missing dependency: $dependency" >&2; exit 1; }
done
if $INCLUDE_GENE || $INCLUDE_REPEAT || $INCLUDE_INTACT || $INCLUDE_GC || $INCLUDE_LTRDATING; then
  command -v bedtools >/dev/null || { echo "[ERROR] Missing dependency: bedtools" >&2; exit 1; }
fi
if $INCLUDE_TELOMERE; then
  command -v quartet.py >/dev/null || { echo "[ERROR] Missing dependency: quartet.py" >&2; exit 1; }
fi
# Create output directory if it doesn't exist
if [[ ! -d $OUTPUT_DIR ]]; then
  mkdir -p "$OUTPUT_DIR"
fi

# Resolve FASTA to an absolute path so subdir cd's don't break it
if [[ "$FASTA" = /* ]]; then
  FASTA_ABS="$FASTA"
else
  FASTA_ABS="$(cd "$(dirname "$FASTA")" && pwd)/$(basename "$FASTA")"
fi

# Grab FASTA basename for naming output files
FASTA_BASE=$(basename "$FASTA")

# Check if required files are provided
[[ ! -f $FASTA ]] && { echo "Missing genome FASTA"; exit 1; }
[[ $INCLUDE_GENE == true && ! -f $GENES ]] && { echo "Gene track requires a GFF file"; exit 1; }
[[ $INCLUDE_REPEAT == true && ! -f $REPEATS ]] && { echo "Repeat track requires a GFF file"; exit 1; }
[[ $INCLUDE_INTACT == true && ! -f $INTACT ]] && { echo "Intact TE track requires a GFF file"; exit 1; }
[[ $INCLUDE_LTRDATING == true && ! -f $LTRDATES ]] && { echo "LTR dating track requires a coordinate file"; exit 1; }
[[ $INCLUDE_LINKS == true && ! -f $LINKS_FILE ]] && { echo "Links track requires a coords file"; exit 1; }

# Create FASTA index if it doesn't exist
if [[ ! -f ${FASTA}.fai ]]; then
  index_log=$(mktemp)
  if ! samtools faidx "$FASTA" 2>"$index_log"; then
    cat "$index_log" >&2; rm -f "$index_log"; exit 1
  fi
  cat "$index_log" >&2
  if grep -qi 'duplicate' "$index_log"; then
    echo "[ERROR] Duplicate FASTA sequence IDs; repair the source before indexing" >&2
    rm -f "$index_log" "${FASTA}.fai"; exit 1
  fi
  rm -f "$index_log"
fi

[[ ! "$FASTA" -nt "${FASTA}.fai" ]] || { echo "[ERROR] FASTA index is stale; rebuild with samtools faidx" >&2; exit 1; }
python3 "$SCRIPT_DIR/tools/prepare_inputs.py" index --fai "${FASTA}.fai" --retained "${FASTA}.fai"

# Filter for chromosomes if requested
### Use optional PATTERN, defaulting to Chr0?[1-9][0-9]? if empty
if [[ "$FILTER_CHRS" == true ]]; then
  if [[ -n "$FILTER_PATTERN" ]]; then
    gawk -v pat="$FILTER_PATTERN" 'BEGIN{IGNORECASE=1} $1 ~ pat' "${FASTA}.fai" \
      | sort -t$'\t' -k1,1 > "${OUTPUT_DIR}/${FASTA_BASE}_chrs.fai"
  else
    gawk 'BEGIN{IGNORECASE=1} $1 ~ /Chr0?[1-9][0-9]?/' "${FASTA}.fai" \
      | sort -t$'\t' -k1,1 > "${OUTPUT_DIR}/${FASTA_BASE}_chrs.fai"
  fi
else
  sort -t$'\t' -k1,1 "${FASTA}.fai" > "${OUTPUT_DIR}/${FASTA_BASE}_chrs.fai"
fi

# Ensure the FASTA index contains entries
if [[ ! -s "${OUTPUT_DIR}/${FASTA_BASE}_chrs.fai" ]]; then
  echo "Error: No valid chromosomes or other contigs found in ${FASTA}.fai"
  exit 1
fi
TEMP_FILES+=("${OUTPUT_DIR}/${FASTA_BASE}_chrs.fai")


# Generate karyotype file for Circos
cut -f1,2 "${OUTPUT_DIR}/${FASTA_BASE}_chrs.fai" | sort -t$'\t' | gawk 'BEGIN {
  idx=0;
  for(i=2;i<=21;i++) l[idx++]="chr"i;
  for(i=2;i<=21;i++) l[idx++]="lum70chr"i;
  for(i=2;i<=21;i++) l[idx++]="lum90chr"i;
} {
  chr=$1;
  gsub(/^[Cc][Hh][Rr]0?/, "", chr);
  print "chr", "-", $1, chr, 0, $2, l[(NR-1)%60]
}' > "${OUTPUT_DIR}/${FASTA_BASE}_karyotype.circos"


# Add telomere bands to karyotype file if requested
if [[ $INCLUDE_TELOMERE == true ]]; then
  TELO_TMP=$(mktemp -d "${OUTPUT_DIR}/${FASTA_BASE}_telomere.XXXXXX")

  TELO_INFO="${TELO_TMP}/${FASTA_BASE}.telo.info"
  TEMP_FILES+=("$TELO_TMP" "$TELO_INFO")

  # Run quartet inside the temp dir, writing .telo.info there
  # quarTeT builds shell commands internally. Fixed relative names prevent a
  # user-supplied filename from becoming shell syntax in that dependency.
  ln -s "$FASTA_ABS" "$TELO_TMP/input.fa"
  (cd "$TELO_TMP" && quartet.py te -i input.fa -c plant -m 50 -p iwgc --noplot)
  [[ -s "$TELO_TMP/iwgc.telo.info" ]] || { echo "[ERROR] quarTeT produced no telomere report; inspect its log" >&2; exit 1; }
  mv "$TELO_TMP/iwgc.telo.info" "$TELO_INFO"

  [[ -s "$TELO_INFO" ]] || { echo "[ERROR] quarTeT produced no telomere report; inspect its log" >&2; exit 1; }
  genome_len=$(gawk '{sum+=$2}END{print sum}' "${OUTPUT_DIR}/${FASTA_BASE}_chrs.fai")

  ### Apply same FILTER_PATTERN logic to telomere info
  if [[ "$FILTER_CHRS" == true ]]; then
    if [[ -n "$FILTER_PATTERN" ]]; then
      grep -v '#' "$TELO_INFO" | gawk -v pat="$FILTER_PATTERN" 'BEGIN{IGNORECASE=1} $1 ~ pat'
    else
      grep -v '#' "$TELO_INFO" | gawk 'BEGIN{IGNORECASE=1} $1 ~ /Chr0?[1-9][0-9]?/'
    fi
  else
    grep -v '#' "$TELO_INFO"
  fi | gawk '{ print $1, $2, $3, $4, "+", $6, "-" }' | \
  gawk -v glen="$genome_len" -v scale="$TELOMERE_SCALE" '{
    t=int(glen*scale); if(t>$2)t=$2;
    if($3=="both"){print $1,0,t; print $1,$2-t,$2}
    else if($3=="right"){print $1,$2-t,$2}
    else if($3=="left"){print $1,0,t}
  }' | \
  gawk '{ print "band", $1, $1"_T"NR, $1"_T"NR, $2, $3, "vdgrey" }' > "${OUTPUT_DIR}/${FASTA_BASE}_telomere_bands.bed"

  cat "${OUTPUT_DIR}/${FASTA_BASE}_telomere_bands.bed" >> "${OUTPUT_DIR}/${FASTA_BASE}_karyotype.circos"
  TEMP_FILES+=("${OUTPUT_DIR}/${FASTA_BASE}_telomere_bands.bed")
fi


# Generate sliding or fixed window file if any tracks other than karyotype are requested
if [[ "$INCLUDE_GENE" == true || "$INCLUDE_REPEAT" == true || \
      "$INCLUDE_INTACT" == true || "$INCLUDE_GC" == true || \
      "$INCLUDE_LTRDATING" == true ]]; then
  if [[ "$USE_SLIDING" == true ]]; then
    STEP_FRACTION="${STEP_SIZE:-0.5}"

# Validate step size is a float between 0 and 1
    if ! [[ "$STEP_FRACTION" =~ ^0(\.[0-9]+)?$|^1(\.0+)?$ ]]; then
        echo "Error: -step must be a decimal fraction between 0 and 1"
        exit 1
    fi

    STEP=$(gawk -v frac="$STEP_FRACTION" -v win="$WINDOW" 'BEGIN { printf "%.0f", frac * win }')
    (( STEP > 0 )) || { echo "[ERROR] sliding step rounds to zero bases" >&2; exit 1; }
    echo "Generating sliding windows with size $WINDOW and step $STEP"
    bedtools makewindows -g "${OUTPUT_DIR}/${FASTA_BASE}_chrs.fai" -w ${WINDOW} -s ${STEP} | \
    gawk -v W=${WINDOW} '{ if (($3 - $2) >= (W / 2)) print }' > "${OUTPUT_DIR}/${FASTA_BASE}_windows.bed"
  else
    echo "Generating fixed windows with size $WINDOW"
    bedtools makewindows -g "${OUTPUT_DIR}/${FASTA_BASE}_chrs.fai" -w ${WINDOW} > "${OUTPUT_DIR}/${FASTA_BASE}_windows.bed"
  fi
  [[ -s "${OUTPUT_DIR}/${FASTA_BASE}_windows.bed" ]] || { echo "[ERROR] No windows remain (check window/sliding settings)" >&2; exit 1; }
  TEMP_FILES+=("${OUTPUT_DIR}/${FASTA_BASE}_windows.bed")
fi


# Private conversions keep stale intact-class intermediates out of this run.
WORK_DIR=$(mktemp -d)
trap 'rm -rf "$WORK_DIR"' EXIT
prepare() {
  python3 "$SCRIPT_DIR/tools/prepare_inputs.py" "$1" --fai "${FASTA}.fai" \
    --retained "${OUTPUT_DIR}/${FASTA_BASE}_chrs.fai" --input "$2" "${@:3}"
}

# Generate gene density track file if requested
if [[ $INCLUDE_GENE == true ]]; then
  prepare gene "$GENES" | sort -k1,1 -k2,2n > "${OUTPUT_DIR}/${FASTA_BASE}_genes_coords.bed"
  TEMP_FILES+=("${OUTPUT_DIR}/${FASTA_BASE}_genes_coords.bed")
  bedtools coverage -sorted -a "${OUTPUT_DIR}/${FASTA_BASE}_windows.bed" -b "${OUTPUT_DIR}/${FASTA_BASE}_genes_coords.bed" | \
  gawk 'BEGIN{OFMT=CONVFMT="%.17g"} {print $1,$2,$3,sqrt($5/$6)}' | sort -k1,1 -k2,2n | \
  gawk 'BEGIN{OFMT=CONVFMT="%.17g"} {v[NR]=$4; l[NR]=$0; if(NR==1||$4<m)m=$4; if(NR==1||$4>M)M=$4} END {for(i=1;i<=NR;i++) print l[i],(M==m?0:(v[i]-m)/(M-m))}' | \
  gawk '{ print $1, $2, $3, $5}' > "${OUTPUT_DIR}/${FASTA_BASE}_gene_coverage.circos"
  TRACK_FILES+=("${OUTPUT_DIR}/${FASTA_BASE}_gene_coverage.circos")
fi


# Generate repeat density track file if requested
if [[ $INCLUDE_REPEAT == true ]]; then
  prepare repeat "$REPEATS" | sort -k1,1 -k2,2n > "${OUTPUT_DIR}/${FASTA_BASE}_repeat_coords.bed"
  TEMP_FILES+=("${OUTPUT_DIR}/${FASTA_BASE}_repeat_coords.bed")
  bedtools coverage -sorted -a "${OUTPUT_DIR}/${FASTA_BASE}_windows.bed" -b "${OUTPUT_DIR}/${FASTA_BASE}_repeat_coords.bed" | \
  gawk 'BEGIN{OFMT=CONVFMT="%.17g"} {print $1,$2,$3,sqrt($5/$6)}' | sort -k1,1 -k2,2n | \
  gawk 'BEGIN{OFMT=CONVFMT="%.17g"} {v[NR]=$4; l[NR]=$0; if(NR==1||$4<m)m=$4; if(NR==1||$4>M)M=$4} END {for(i=1;i<=NR;i++) print l[i],(M==m?0:(v[i]-m)/(M-m))}' | \
  gawk '{ print $1, $2, $3, $5}' > "${OUTPUT_DIR}/${FASTA_BASE}_repeat_coverage.circos"
  TRACK_FILES+=("${OUTPUT_DIR}/${FASTA_BASE}_repeat_coverage.circos")
fi


# Generate intact TE coverage track files if requested
if [[ $INCLUDE_INTACT == true ]]; then
  prepare intact "$INTACT" --output "$WORK_DIR"
  shopt -s nullglob
  class_files=("$WORK_DIR/"*.bed)
  shopt -u nullglob
  for coords in "${class_files[@]}"; do
    class=$(basename "$coords" .bed)
    sort -k1,1 -k2,2n "$coords" > "$coords.sorted"
    mv "$coords.sorted" "$coords"
    if [[ "$KEEP_TEMP" == true ]]; then
      cp "$coords" "${OUTPUT_DIR}/${FASTA_BASE}_${class}_intactTE_coords.bed"
    fi
    bedtools coverage -sorted -a "${OUTPUT_DIR}/${FASTA_BASE}_windows.bed" -b "$coords" | \
    gawk 'BEGIN{OFMT=CONVFMT="%.17g"} {print $1,$2,$3,$5/$6}' | sort -k1,1 -k2,2n | \
    gawk 'BEGIN{OFMT=CONVFMT="%.17g"} {v[NR]=$4; l[NR]=$0; if(NR==1||$4<m)m=$4; if(NR==1||$4>M)M=$4}
         END {for(i=1;i<=NR;i++) print l[i],(M==m?0:(v[i]-m)/(M-m))}' | \
    gawk '{ print $1, $2, $3, $5 }' > "${OUTPUT_DIR}/${FASTA_BASE}_${class}_intactTE_coverage.circos"
    TRACK_FILES+=("${OUTPUT_DIR}/${FASTA_BASE}_${class}_intactTE_coverage.circos")
  done
fi


# Generate LTR dating track file if requested
if [[ $INCLUDE_LTRDATING == true ]]; then
  prepare ltr "$LTRDATES" | sort -k1,1 -k2,2n > "${OUTPUT_DIR}/${FASTA_BASE}_LTR_insertion.bed"
  # Per-window mean age (keep all windows)
  if [[ -s "${OUTPUT_DIR}/${FASTA_BASE}_LTR_insertion.bed" ]]; then
  bedtools map -a "${OUTPUT_DIR}/${FASTA_BASE}_windows.bed" \
               -b "${OUTPUT_DIR}/${FASTA_BASE}_LTR_insertion.bed" \
               -c 4 -o mean \
  > "${OUTPUT_DIR}/${FASTA_BASE}_LTR_age.bed"
  else
    gawk 'BEGIN{OFS="\t"} {print $1,$2,$3,"."}' "${OUTPUT_DIR}/${FASTA_BASE}_windows.bed" > "${OUTPUT_DIR}/${FASTA_BASE}_LTR_age.bed"
  fi

  # Normalize ignoring zeros, invert so younger=larger, and offset to 1..2; empties=0
  gawk -v OFS='\t' '
  {
    chrom=$1; start=$2; end=$3; raw=$4;

    has = (raw != ".");            # did this window get any LTRs?
    val = has ? (raw+0) : 0;       # numeric value when present

    A[NR]=chrom; B[NR]=start; C[NR]=end; V[NR]=val; H[NR]=has;

    if (has) {
      if (!seen++) { m=val; M=val }
      if (val < m) m=val;
      if (val > M) M=val;
    }
  }
  END{
    for (i=1; i<=NR; i++) {
      if (!H[i]) {
        out = 0;                   # no insertions in this window
      } else if (M == m) {
        out = 1.5;                 # all ages identical (incl. all zeros)
      } else {
        y = 1 - ((V[i] - m) / (M - m));  # younger (smaller age) -> larger y
        out = y + 1;                      # shift to [1,2]
      }
      printf "%s\t%s\t%s\t%.6f\n", A[i], B[i], C[i], out;
    }
  }
  ' "${OUTPUT_DIR}/${FASTA_BASE}_LTR_age.bed" \
  > "${OUTPUT_DIR}/${FASTA_BASE}_LTRdating_coverage.circos"
  TRACK_FILES+=("${OUTPUT_DIR}/${FASTA_BASE}_LTRdating_coverage.circos")


  TEMP_FILES+=("${OUTPUT_DIR}/${FASTA_BASE}_LTR_insertion.bed" "${OUTPUT_DIR}/${FASTA_BASE}_LTR_age.bed")
fi


# Generate GC content track file if requested
if [[ $INCLUDE_GC == true ]]; then
  bedtools nuc -fi "$FASTA" -bed "${OUTPUT_DIR}/${FASTA_BASE}_windows.bed" > "${OUTPUT_DIR}/gc_content.bed"
  TEMP_FILES+=("${OUTPUT_DIR}/gc_content.bed")
  gawk 'BEGIN{OFMT=CONVFMT="%.17g"} !/^#/ {print $1,$2,$3,($7+$8)/$12}' "${OUTPUT_DIR}/gc_content.bed" | sort -k1,1 -k2,2n | \
  gawk 'BEGIN{OFMT=CONVFMT="%.17g"} {v[NR]=$4; l[NR]=$0; if(NR==1 || $4<m) m=$4; if(NR==1 || $4>M) M=$4}
       END {for(i=1;i<=NR;i++) print l[i], (M==m ? 0 : (v[i]-m)/(M-m))}' | \
  gawk '{print $1, $2, $3, $5}' > "${OUTPUT_DIR}/${FASTA_BASE}_gc.circos"
  TRACK_FILES+=("${OUTPUT_DIR}/${FASTA_BASE}_gc.circos")
fi

# --- Normalize links input to 6-column format if needed -----------------------
# Output schema: qchr qstart qend schr sstart send
SIXCOL_LINKS="${OUTPUT_DIR}/${FASTA_BASE}_links_input.sixcol"

if [[ $INCLUDE_LINKS == true ]]; then
  echo "Preparing links from $LINKS_FILE"

  prepare links "$LINKS_FILE" > "$SIXCOL_LINKS"

  TEMP_FILES+=("$SIXCOL_LINKS")
fi

# Generate syntenic links file if requested
if [[ $INCLUDE_LINKS == true ]]; then
  echo "Generating syntenic links (normalized) from $SIXCOL_LINKS"
  gawk 'BEGIN{OFS="\t"} { print $1,$2,$3,$4,$5,$6 }' "$SIXCOL_LINKS" \
  | sort -k1,1 -k2,2n -k4,4 -k5,5n \
  | gawk 'BEGIN{OFS="\t"}
    {
      if ($1 != qchr || $4 != schr || $2 > qend || $3 < qstart || $5 > send || $6 < sstart) {
        if (NR > 1) print qchr, qstart, qend, schr, sstart, send;
        qchr=$1; qstart=$2; qend=$3; schr=$4; sstart=$5; send=$6;
      } else {
        if ($2 < qstart) qstart=$2;
        if ($5 < sstart) sstart=$5;
        if ($3 > qend) qend=$3;
        if ($6 > send) send=$6;
      }
    }
    END { if (NR>0) print qchr, qstart, qend, schr, sstart, send; }' \
  | gawk 'BEGIN{OFS="\t"}{
      qlen=$3-$2; slen=$6-$5; tlen=qlen+slen;
      if (tlen >= '"${MIN_TLEN:-0}"') print tlen, $0
    }' \
  | sort -k1,1nr \
  | gawk -v n="$TOP_N" 'n==0 || NR<=n' \
  | { if [[ "$LINK_ORDER" == dsc ]]; then cat; else sort -k1,1n; fi; } \
  | cut -f2- > "${OUTPUT_DIR}/${FASTA_BASE}_links.circos"
fi



# Some dependency versions print errors but return 0. Require all expected
# windows and finite numeric data before reporting a complete track run.
if (( ${#TRACK_FILES[@]} )); then
  expected_windows=$(wc -l < "${OUTPUT_DIR}/${FASTA_BASE}_windows.bed")
  for track in "${TRACK_FILES[@]}"; do
    gawk -v expected="$expected_windows" '
      NF!=4 || $4 !~ /^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/ {bad=1}
      END{exit (bad || NR!=expected)}' "$track" || {
        echo "[ERROR] Incomplete or nonnumeric track: $track" >&2; exit 1;
      }
  done
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

# Exit successfully even when -keep-temp bypasses cleanup.
exit 0
