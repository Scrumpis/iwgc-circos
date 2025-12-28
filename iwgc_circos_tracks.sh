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
  echo "  -links           Add syntenic links from a .coords file (e.g., MUMmer or minimap2 output)"
  echo "  -min-tlen        Minimum length of syntenic links"
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
  echo "singularity exec iwgc-circos.sif ./iwgc-circos-tracks.sh <FASTA> [options]"
  echo "or"
  echo "docker run --rm -v \$(pwd):/data iwgc-circos-tracks:latest ./iwgc-circos-tracks.sh <FASTA> [options]"
  echo
  echo "Note: For visualization/automation purposes, all density tracks are normalized, gene density is sqrt transformed, and repeat density is power"
  echo "3 transformed"
  exit 1
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
TELOMERE_SCALE=0.005
USE_SLIDING=false
STEP_SIZE=""
FILTER_CHRS=false
### New: optional pattern used when FILTER_CHRS is true
FILTER_PATTERN=""

FASTA=$1; shift

while [[ $# -gt 0 ]]; do
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
    *) echo "Unknown option: $1"; usage ;;
  esac
done

# Default window and step sizes
WINDOW=${WINDOW:-1000000}
STEP_SIZE=${STEP_SIZE:-0.5}
OUTPUT_DIR=${OUTPUT_DIR:-.}     # current directory by default
MIN_TLEN=${MIN_TLEN:-0}
TOP_N=${TOP_N:-0}
LINK_ORDER=${LINK_ORDER:-asc}   # asc|desc

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
  samtools faidx ${FASTA}
fi

# Filter for chromosomes if requested
### Use optional PATTERN, defaulting to Chr0?[1-9][0-9]? if empty
if [[ "$FILTER_CHRS" == true ]]; then
  if [[ -n "$FILTER_PATTERN" ]]; then
    awk -v pat="$FILTER_PATTERN" 'BEGIN{IGNORECASE=1} $1 ~ pat' "${FASTA}.fai" \
      | sort -t$'\t' -k1,1 > "${OUTPUT_DIR}/${FASTA_BASE}_chrs.fai"
  else
    awk 'BEGIN{IGNORECASE=1} $1 ~ /Chr0?[1-9][0-9]?/' "${FASTA}.fai" \
      | sort -t$'\t' -k1,1 > "${OUTPUT_DIR}/${FASTA_BASE}_chrs.fai"
  fi
else
  sort -t$'\t' -k1,1 ${FASTA}.fai > "${OUTPUT_DIR}/${FASTA_BASE}_chrs.fai"
fi

# Ensure the FASTA index contains entries
if [[ ! -s "${OUTPUT_DIR}/${FASTA_BASE}_chrs.fai" ]]; then
  echo "Error: No valid chromosomes or other contigs found in ${FASTA}.fai"
  exit 1
fi
TEMP_FILES+=("${OUTPUT_DIR}/${FASTA_BASE}_chrs.fai")


# Generate karyotype file for Circos
cut -f1,2 "${OUTPUT_DIR}/${FASTA_BASE}_chrs.fai" | sort -t$'\t' | awk 'BEGIN {
  idx=0;
  for(i=2;i<=21;i++) l[idx++]="chr"i;
  for(i=2;i<=21;i++) l[idx++]="lum70chr"i;
  for(i=2;i<=21;i++) l[idx++]="lum90chr"i;
} {
  chr=$1;
  gsub(/^[Cc][Hh][Rr]0?/, "", chr);
  print "chr", "-", $1, chr, 0, $2, l[NR-1]
}' > "${OUTPUT_DIR}/${FASTA_BASE}_karyotype.circos"


# Add telomere bands to karyotype file if requested
if [[ $INCLUDE_TELOMERE == true ]]; then
  TELO_TMP="${OUTPUT_DIR}/${FASTA_BASE}_telomere_tmp"
  mkdir -p "$TELO_TMP"

  TELO_INFO="${TELO_TMP}/${FASTA_BASE}.telo.info"
  TEMP_FILES+=("$TELO_TMP" "$TELO_INFO")

  # Run quartet inside the temp dir, writing .telo.info there
  (cd "$TELO_TMP" && quartet.py te -i "${FASTA_ABS}" -c plant -m 50 -p "${FASTA_BASE}" --noplot)

  genome_len=$(awk '{sum+=$2}END{print sum}' "${OUTPUT_DIR}/${FASTA_BASE}_chrs.fai")

  ### Apply same FILTER_PATTERN logic to telomere info
  if [[ "$FILTER_CHRS" == true ]]; then
    if [[ -n "$FILTER_PATTERN" ]]; then
      grep -v '#' "$TELO_INFO" | awk -v pat="$FILTER_PATTERN" 'BEGIN{IGNORECASE=1} $1 ~ pat'
    else
      grep -v '#' "$TELO_INFO" | awk 'BEGIN{IGNORECASE=1} $1 ~ /Chr0?[1-9][0-9]?/'
    fi
  else
    grep -v '#' "$TELO_INFO"
  fi | awk '{ print $1, $2, $3, $4, "+", $6, "-" }' | \
  awk -v glen="$genome_len" -v scale="$TELOMERE_SCALE" '{
    t=int(glen*scale);
    if($3=="both"){print $1,0,t; print $1,$2-t,$2}
    else if($3=="right"){print $1,$2-t,$2}
    else if($3=="left"){print $1,0,t}
  }' | \
  awk '{ print "band", $1, $1"_T"NR, $1"_T"NR, $2, $3, "vdgrey" }' > "${OUTPUT_DIR}/${FASTA_BASE}_telomere_bands.bed"

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

    STEP=$(awk -v frac="$STEP_FRACTION" -v win="$WINDOW" 'BEGIN { printf "%.0f", frac * win }')
    echo "Generating sliding windows with size $WINDOW and step $STEP"
    bedtools makewindows -g "${OUTPUT_DIR}/${FASTA_BASE}_chrs.fai" -w ${WINDOW} -s ${STEP} | \
    awk -v W=${WINDOW} '{ if (($3 - $2) >= (W / 2)) print }' > "${OUTPUT_DIR}/${FASTA_BASE}_windows.bed"
  else
    echo "Generating fixed windows with size $WINDOW"
    bedtools makewindows -g "${OUTPUT_DIR}/${FASTA_BASE}_chrs.fai" -w ${WINDOW} > "${OUTPUT_DIR}/${FASTA_BASE}_windows.bed"
  fi
  TEMP_FILES+=("${OUTPUT_DIR}/${FASTA_BASE}_windows.bed")
fi


# Generate gene density track file if requested
if [[ $INCLUDE_GENE == true ]]; then
  grep 'gene' ${GENES} | awk '{print $1,$4,$5}' OFS='\t' | sort -k1,1 -k2,2n > "${OUTPUT_DIR}/${FASTA_BASE}_genes_coords.bed"
  TEMP_FILES+=("${OUTPUT_DIR}/${FASTA_BASE}_genes_coords.bed")
  bedtools coverage -sorted -a "${OUTPUT_DIR}/${FASTA_BASE}_windows.bed" -b "${OUTPUT_DIR}/${FASTA_BASE}_genes_coords.bed" | \
  awk '{print $1,$2,$3,sqrt($7)}' | sort -k1,1 -k2,2n | \
  awk '{v[NR]=$4; l[NR]=$0; if(NR==1||$4<m)m=$4; if(NR==1||$4>M)M=$4} END {for(i=1;i<=NR;i++) print l[i],(M==m?0:(v[i]-m)/(M-m))}' | \
  awk '{ print $1, $2, $3, $5}' > "${OUTPUT_DIR}/${FASTA_BASE}_gene_coverage.circos"
fi


# Generate repeat density track file if requested
if [[ $INCLUDE_REPEAT == true ]]; then
  grep -v '#' ${REPEATS} | awk '{print $1,$4,$5}' OFS='\t' | sort -k1,1 -k2,2n > "${OUTPUT_DIR}/${FASTA_BASE}_repeat_coords.bed"
  TEMP_FILES+=("${OUTPUT_DIR}/${FASTA_BASE}_repeat_coords.bed")
  bedtools coverage -sorted -a "${OUTPUT_DIR}/${FASTA_BASE}_windows.bed" -b "${OUTPUT_DIR}/${FASTA_BASE}_repeat_coords.bed" | \
  awk '{print $1,$2,$3,($7^3)}' | sort -k1,1 -k2,2n | \
  awk '{v[NR]=$4; l[NR]=$0; if(NR==1||$4<m)m=$4; if(NR==1||$4>M)M=$4} END {for(i=1;i<=NR;i++) print l[i],(M==m?0:(v[i]-m)/(M-m))}' | \
  awk '{ print $1, $2, $3, $5}' > "${OUTPUT_DIR}/${FASTA_BASE}_repeat_coverage.circos"
fi


# Generate intact TE coverage track files if requested
if [[ $INCLUDE_INTACT == true ]]; then
  awk -F'\t' -v base="${FASTA_BASE}" -v out="${OUTPUT_DIR}" '
    BEGIN { IGNORECASE=1 }
    match($9, /Classification=([^;]+)/, a) {
      split(a[1], p, "/");
      if (p[2] != "") {
        print > (out "/" base "_intactTE_" p[2] ".gff")
      }
    }
  ' "$INTACT"

  output_files=$(ls "${OUTPUT_DIR}/${FASTA_BASE}_intactTE_"*.gff 2>/dev/null)
  TEMP_FILES+=($output_files)

  for gff in $output_files; do
    class=$(basename "$gff" | sed -E "s/${FASTA_BASE}_intactTE_(.*)\.gff/\1/")
    coords="${OUTPUT_DIR}/${FASTA_BASE}_${class}_intactTE_coords.bed"
    awk '{print $1,$4,$5}' OFS='\t' "$gff" | sort -k1,1 -k2,2n > "$coords"
    TEMP_FILES+=("$coords")
    bedtools coverage -sorted -a "${OUTPUT_DIR}/${FASTA_BASE}_windows.bed" -b "$coords" | \
    awk '{print $1,$2,$3,$7}' | sort -k1,1 -k2,2n | \
    awk '{v[NR]=$4; l[NR]=$0; if(NR==1||$4<m)m=$4; if(NR==1||$4>M)M=$4}
         END {for(i=1;i<=NR;i++) print l[i],(M==m?0:(v[i]-m)/(M-m))}' | \
    awk '{ print $1, $2, $3, $5 }' > "${OUTPUT_DIR}/${FASTA_BASE}_${class}_intactTE_coverage.circos"
  done
fi


# Generate LTR dating track file if requested
if [[ $INCLUDE_LTRDATING == true ]]; then
  awk -F'\t' 'BEGIN{OFS="\t"} !/^#/ {print $1, $12}' "${LTRDATES}" \
| sed -e 's/:/\t/' -e 's/\.\./\t/' \
| awk -v OFS='\t' '{ if ($2 > $3) { t=$2; $2=$3; $3=t } print }' \
| sort -k1,1 -k2,2n \
> "${OUTPUT_DIR}/${FASTA_BASE}_LTR_insertion.bed"
  # Per-window mean age (keep all windows)
  bedtools map -a "${OUTPUT_DIR}/${FASTA_BASE}_windows.bed" \
               -b "${OUTPUT_DIR}/${FASTA_BASE}_LTR_insertion.bed" \
               -c 4 -o mean \
  > "${OUTPUT_DIR}/${FASTA_BASE}_LTR_age.bed"

  # Normalize ignoring zeros, invert so younger=larger, and offset to 1..2; empties=0
  awk -v OFS='\t' '
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


  TEMP_FILES+=("${OUTPUT_DIR}/${FASTA_BASE}_LTR_insertion.bed" "${OUTPUT_DIR}/${FASTA_BASE}_LTR_age.bed")
fi


# Generate GC content track file if requested
if [[ $INCLUDE_GC == true ]]; then
  bedtools nuc -fi ${FASTA} -bed "${OUTPUT_DIR}/${FASTA_BASE}_windows.bed" > "${OUTPUT_DIR}/gc_content.bed"
  TEMP_FILES+=("${OUTPUT_DIR}/gc_content.bed")
  awk '{print $1,$2,$3,$5}' "${OUTPUT_DIR}/gc_content.bed" | grep -v '#' | sort -k1,1 -k2,2n | \
  awk '{v[NR]=$4; l[NR]=$0; if(NR==1 || $4<m) m=$4; if(NR==1 || $4>M) M=$4} 
       END {for(i=1;i<=NR;i++) print l[i], (M==m ? 0 : (v[i]-m)/(M-m))}' | \
  awk '{print $1, $2, $3, $5}' > "${OUTPUT_DIR}/${FASTA_BASE}_gc.circos"
fi

# --- Normalize links input to 6-column format if needed -----------------------
# Output schema: qchr qstart qend schr sstart send
SIXCOL_LINKS="${OUTPUT_DIR}/${FASTA_BASE}_links_input.sixcol"

if [[ $INCLUDE_LINKS == true ]]; then
  echo "Preparing links from $LINKS_FILE"

  # Case A: already simple 6-col (allow comments/blank lines). Require numeric cols 2,3,5,6.
  if awk '
      BEGIN{ok=1}
      /^[[:space:]]*($|#)/{next}
      {
        if (NF!=6) {ok=0; exit}
        if (($2+0)!=$2 || ($3+0)!=$3 || ($5+0)!=$5 || ($6+0)!=$6) {ok=0; exit}
      }
      END{exit ok?0:1}
    ' "$LINKS_FILE"
  then
    # Use as-is but enforce start<=end on both sides
    awk -v OFS="\t" '
      /^[[:space:]]*($|#)/{next}
      {
        a1=$2; b1=$3; if (a1>b1){t=a1;a1=b1;b1=t}
        a2=$5; b2=$6; if (a2>b2){t=a2;a2=b2;b2=t}
        print $1,a1,b1,$4,a2,b2
      }
    ' "$LINKS_FILE" > "$SIXCOL_LINKS"

  else
    # Case B: 11-col (S1 E1 S2 E2 LEN1 LEN2 %IDY LENR LENQ REF_ID QUERY_ID), no header.
    # Require: NF>=11, first 9 fields numeric (allow decimal for %IDY), last 2 are IDs.
    awk -v OFS="\t" '
      function isn(x){ return (x ~ /^-?[0-9]+(\.[0-9]+)?$/) }  # numeric (int/float)
      /^[[:space:]]*($|#)/{next}
      {
        if (NF<11) next
        # check 1..9 numeric
        for (i=1;i<=9;i++){ gsub(/,/, "", $i); if (!isn($i)) nextline=1 }
        if (nextline){ nextline=0; next }
        s1=$1+0; e1=$2+0; s2=$3+0; e2=$4+0
        ref=$(NF-1); qry=$NF
        if (s1>e1){t=s1;s1=e1;e1=t}
        if (s2>e2){t=s2;e2=s2;s2=t}
        print ref, s1, e1, qry, s2, e2
      }
    ' "$LINKS_FILE" > "$SIXCOL_LINKS"

    if [[ ! -s "$SIXCOL_LINKS" ]]; then
      echo "Error: Could not parse $LINKS_FILE as 6-col or 11-col coords."
      exit 1
    fi
  fi

  TEMP_FILES+=("$SIXCOL_LINKS")
fi

# Generate syntenic links file if requested
if [[ $INCLUDE_LINKS == true ]]; then
  echo "Generating syntenic links (normalized) from $SIXCOL_LINKS"
  awk 'BEGIN{OFS="\t"} { print $1,$2,$3,$4,$5,$6 }' "$SIXCOL_LINKS" \
  | sort -k1,1 -k2,2n -k4,4 -k5,5n \
  | awk 'BEGIN{OFS="\t"}
    {
      if ($1 != qchr || $4 != schr || $2 > qend || $5 > send) {
        if (NR > 1) print qchr, qstart, qend, schr, sstart, send;
        qchr=$1; qstart=$2; qend=$3; schr=$4; sstart=$5; send=$6;
      } else {
        if ($3 > qend) qend=$3;
        if ($6 > send) send=$6;
      }
    }
    END { if (NR>0) print qchr, qstart, qend, schr, sstart, send; }' \
  | awk 'BEGIN{OFS="\t"}{
      qlen=$3-$2; slen=$6-$5; tlen=qlen+slen;
      if (tlen >= '"${MIN_TLEN:-0}"') print tlen, $0
    }' \
  | { if [[ "${LINK_ORDER:-asc}" == "dsc" ]]; then sort -k1,1nr; else sort -k1,1n; fi; } \
  | { if (( ${TOP_N:-0} > 0 )); then head -n "${TOP_N}"; else cat; fi; } \
  | cut -f2- \
  | {
      ### Updated: apply same optional PATTERN to link chromosomes
      if [[ "$FILTER_CHRS" == true ]]; then
        if [[ -n "$FILTER_PATTERN" ]]; then
          awk -F'\t' -v pat="$FILTER_PATTERN" 'BEGIN{IGNORECASE=1} $1 ~ pat && $4 ~ pat'
        else
          awk -F'\t' 'BEGIN{IGNORECASE=1} $1 ~ /^chr0?[1-9][0-9]?[a-z]?$/ && $4 ~ /^chr0?[1-9][0-9]?[a-z]?$/'
        fi
      else
        cat
      fi
    } > "${OUTPUT_DIR}/${FASTA_BASE}_links.circos"
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
