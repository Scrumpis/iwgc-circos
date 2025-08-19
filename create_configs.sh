#!/usr/bin/env bash
set -euo pipefail

# Defaults
TEMPLATE="iwgc_circos_template.config"
OUTDIR="iwgc_circos"
GAP=false
IDEO=""

# Parse args (single-dash flags)
while [[ $# -gt 0 ]]; do
  case "$1" in
    -template) TEMPLATE="$2"; shift 2;;
    -outdir)   OUTDIR="$2"; shift 2;;
    -ideogram) IDEO="$2"; shift 2;;   # optional override; default set after OUTFILE
    -gap)      GAP=true; shift;;
    -h|-help)
      echo "Usage: $0 [-template FILE] [-outdir DIR] [-ideogram FILE] [-gap]"
      exit 0;;
    *)
      echo "[ERROR] Unknown arg: $1" >&2; exit 1;;
  esac
done

[[ -f "$TEMPLATE" ]] || { echo "[ERROR] Template not found: $TEMPLATE" >&2; exit 1; }
mkdir -p "$OUTDIR"
OUTFILE="${OUTDIR%/}/iwgc_circos.conf"

# Default ideogram path: $OUTDIR/ideogram.conf unless -ideogram is provided
if [[ -z "${IDEO}" ]]; then
  IDEO="${OUTDIR%/}/ideogram.conf"
fi
echo "[INFO] ideogram target: $IDEO"

# === helpers ===
escape_sed_repl() { printf '%s' "$1" | sed -e 's/[&|\\]/\\&/g'; }

block_file_pattern() {
  awk 'match($0, /^[[:space:]]*file[[:space:]]*=[[:space:]]*(.*)$/, a){print a[1]}' "$1" | head -n1
}

# Return 0 if any match exists
pattern_has_match() {
  local pat="$1"
  pat="${pat%%#*}"
  pat="$(echo "$pat" | sed -E 's/^[[:space:]]+|[[:space:]]+$//g')"
  compgen -G "$pat" >/dev/null 2>&1
}

# Resolve a glob to one concrete path (first in lexicographic order). Prints "" if none.
resolve_first_match() {
  local pat="$1"
  local first
  first="$(compgen -G "$pat" | sort | head -n1 || true)"
  printf '%s' "$first"
}

rewrite_block_file_line() {
  local blk="$1" new="$2"
  awk -v new="$new" '
    BEGIN{done=0}
    /^[[:space:]]*file[[:space:]]*=/ && !done { sub(/=.*/, "= " new); print; done=1; next }
    { print }
  ' "$blk" > "$blk.tmp" && mv "$blk.tmp" "$blk"
}

is_labels_block()     { grep -qE '^[[:space:]]*type[[:space:]]*=[[:space:]]*text' "$1"; }
is_highlights_block() { grep -qE '^[[:space:]]*type[[:space:]]*=[[:space:]]*highlight' "$1"; }
block_matches_regex() { grep -qE "$2" "$1"; }

# === workspace ===
WORKDIR="$(mktemp -d)"
trap 'rm -rf "$WORKDIR"' EXIT

PRE="$WORKDIR/pre.txt"
PLOTS_BODY="$WORKDIR/plots_body.txt"
POST="$WORKDIR/post.txt"
BLOCK_DIR="$WORKDIR/blocks"
mkdir -p "$BLOCK_DIR"

# --- Evaluate ${chrs} using grep -v + awk + paste (NOT literal in output) ---
CHRS_EVAL="$(grep -v '_T[0-9]\+' iwgc_circos_data/*karyotype.circos \
  | awk 'NF>=3 {print $3}' \
  | paste -sd';' - || true)"
if [[ -z "${CHRS_EVAL// /}" ]]; then
  echo "[ERROR] Chromosome extraction returned empty. Check karyotype files." >&2
  exit 1
fi

# Substitute placeholders (${chrs}, drop ${reverse}) into a working copy
SUBBED="$WORKDIR/subbed.txt"
sed -e "s|\${chrs}|$(escape_sed_repl "$CHRS_EVAL")|g" \
    -e "s|\${reverse}||g" \
    "$TEMPLATE" > "$SUBBED"

# --- Resolve the karyotype wildcard to a concrete file (first lexicographic match) ---
kary_glob="$(sed -nE 's/^[[:space:]]*karyotype[[:space:]]*=[[:space:]]*(.*)/\1/p' "$SUBBED" | head -n1)"
if [[ -z "${kary_glob// /}" ]]; then
  echo "[ERROR] Could not find a karyotype line in template." >&2
  exit 1
fi
kary_file="$(resolve_first_match "$kary_glob")"
if [[ -z "${kary_file// /}" ]]; then
  echo "[ERROR] karyotype pattern has no matches: $kary_glob" >&2
  exit 1
fi
# Replace only the first karyotype line
awk -v newk="$kary_file" '
  BEGIN{done=0}
  /^[[:space:]]*karyotype[[:space:]]*=/ && !done { print "karyotype   = " newk; done=1; next }
  { print }
' "$SUBBED" > "$SUBBED.k" && mv "$SUBBED.k" "$SUBBED"

# --- Split into pre / plots body / post (exclude <plots> tags themselves) ---
awk -v PRE="$PRE" -v BODY="$PLOTS_BODY" -v POST="$POST" '
  BEGIN{sec="pre"}
  /<plots>/{sec="body"; next}
  /<\/plots>/{sec="post"; next}
  { if (sec=="pre") print > PRE; else if (sec=="body") print > BODY; else print > POST; }
' "$SUBBED"

# --- Extract each <plot>...</plot> block to its own file ---
awk -v dir="$BLOCK_DIR" '
  BEGIN{inblock=0; n=0}
  /<plot>/{inblock=1; n++; fn=sprintf("%s/block_%03d.txt", dir, n); print > fn; next}
  /<\/plot>/{print >> fn; inblock=0; next}
  { if(inblock){print >> fn} }
' "$PLOTS_BODY"

# --- Dynamic tracks order & regex (top -> inward) ---
DYN_KEYS=(gene_cov repeat_cov ltrdating gypsy copia unknown dta dtc dth dtm dtt helitron gc_content)
DYN_REGEX=(
  "iwgc_circos_data/\*_gene_coverage\.circos"
  "iwgc_circos_data/\*_repeat_coverage\.circos"
  "iwgc_circos_data/\*_LTRdating_coverage\.circos"
  "iwgc_circos_data/\*_Gypsy_intactTE_coverage\.circos"
  "iwgc_circos_data/\*_Copia_intactTE_coverage\.circos"
  "iwgc_circos_data/\*_unknown_intactTE_coverage\.circos"
  "iwgc_circos_data/\*_DTA_intactTE_coverage\.circos"
  "iwgc_circos_data/\*_DTC_intactTE_coverage\.circos"
  "iwgc_circos_data/\*_DTH_intactTE_coverage\.circos"
  "iwgc_circos_data/\*_DTM_intactTE_coverage\.circos"
  "iwgc_circos_data/\*_DTT_intactTE_coverage\.circos"
  "iwgc_circos_data/\*_Helitron_intactTE_coverage\.circos"
  "iwgc_circos_data/\*_gc\.circos"
)

# Geometry
UNIT=0.01
TRACK_WIDTH=$(awk -v u="$UNIT" 'BEGIN{printf "%.4f", 4*u}')   # 0.04r
TRACK_SPACING=$(awk -v u="$UNIT" 'BEGIN{printf "%.4f", 3*u}') # 0.03r
LINKS_GAP=$(awk -v u="$UNIT" 'BEGIN{printf "%.4f", 0.5*u}')   # 0.005r

# Gather blocks
mapfile -t BLOCKS < <(ls "$BLOCK_DIR"/block_*.txt 2>/dev/null | sort -V)

# Track dynamic blocks & presence
declare -A KEY_TO_BLOCK=()
declare -a PRESENT
declare -a HAS_BLOCK

for i in "${!DYN_KEYS[@]}"; do
  key="${DYN_KEYS[$i]}"; regex="${DYN_REGEX[$i]}"
  HAS_BLOCK[$i]=0; PRESENT[$i]=0
  for blk in "${BLOCKS[@]}"; do
    if block_matches_regex "$blk" "$regex"; then
      HAS_BLOCK[$i]=1
      KEY_TO_BLOCK["$key"]="$blk"
      pat="$(block_file_pattern "$blk")"
      if [[ -z "$pat" ]] || pattern_has_match "$pat"; then
        PRESENT[$i]=1
      else
        PRESENT[$i]=0
      fi
      break
    fi
  done
done

# --- Derive default r1/r0 from the template (no hard-coded arrays) ---
declare -a R1_DEFAULT R0_DEFAULT
for i in "${!DYN_KEYS[@]}"; do
  R1_DEFAULT[$i]="" ; R0_DEFAULT[$i]=""
  if [[ "${HAS_BLOCK[$i]}" -eq 1 ]]; then
    blk="${KEY_TO_BLOCK[${DYN_KEYS[$i]}]}"
    R1_DEFAULT[$i]="$(
      awk 'match($0,/^[[:space:]]*r1[[:space:]]*=[[:space:]]*([0-9.]+)r/,a){print a[1]; exit}' "$blk"
    )"
    R0_DEFAULT[$i]="$(
      awk 'match($0,/^[[:space:]]*r0[[:space:]]*=[[:space:]]*([0-9.]+)r/,a){print a[1]; exit}' "$blk"
    )"
  fi
done
TOP_R1="${R1_DEFAULT[0]:-0.97}"

# Labels/highlights: drop if file missing; otherwise resolve to concrete path now
for blk in "${BLOCKS[@]}"; do
  if is_labels_block "$blk" || is_highlights_block "$blk"; then
    pat="$(block_file_pattern "$blk")"
    if [[ -n "$pat" ]]; then
      if pattern_has_match "$pat"; then
        resolved="$(resolve_first_match "$pat")"
        rewrite_block_file_line "$blk" "$resolved"
      else
        : > "$blk"
      fi
    fi
  fi
done

# Find first missing dynamic predecessor
first_missing_idx=-1
for i in "${!DYN_KEYS[@]}"; do
  if [[ "${HAS_BLOCK[$i]}" -eq 1 && "${PRESENT[$i]}" -eq 0 ]]; then
    first_missing_idx="$i"; break
  fi
done

INNERMOST_R0=""
if [[ "$first_missing_idx" -lt 0 ]]; then
  # No missing: keep defaults; resolve file paths; inner r0 = last present default
  last_present_idx=-1
  for i in "${!DYN_KEYS[@]}"; do
    if [[ "${HAS_BLOCK[$i]}" -eq 1 && "${PRESENT[$i]}" -eq 1 ]]; then
      blk="${KEY_TO_BLOCK[${DYN_KEYS[$i]}]}"
      pat="$(block_file_pattern "$blk")"
      if [[ -n "$pat" ]]; then
        resolved="$(resolve_first_match "$pat")"
        rewrite_block_file_line "$blk" "$resolved"
      fi
      last_present_idx="$i"
    fi
  done
  if [[ "$last_present_idx" -ge 0 ]]; then
    INNERMOST_R0="${R0_DEFAULT[$last_present_idx]}"
  fi
else
  # Partial repack: before first missing -> keep defaults (resolve file); at/after -> repack + resolve file
  anchor_r1="$TOP_R1"
  last_kept_before=-1
  for i in "${!DYN_KEYS[@]}"; do
    if [[ "$i" -lt "$first_missing_idx" && "${HAS_BLOCK[$i]}" -eq 1 && "${PRESENT[$i]}" -eq 1 ]]; then
      blk="${KEY_TO_BLOCK[${DYN_KEYS[$i]}]}"
      pat="$(block_file_pattern "$blk")"
      if [[ -n "$pat" ]]; then
        resolved="$(resolve_first_match "$pat")"
        rewrite_block_file_line "$blk" "$resolved"
      fi
      last_kept_before="$i"
    elif [[ "$i" -lt "$first_missing_idx" && "${HAS_BLOCK[$i]}" -eq 1 && "${PRESENT[$i]}" -eq 0 ]]; then
      blk="${KEY_TO_BLOCK[${DYN_KEYS[$i]}]}"; : > "$blk"
    fi
  done
  if [[ "$last_kept_before" -ge 0 && -n "${R0_DEFAULT[$last_kept_before]}" ]]; then
    anchor_r1="$(awk -v r0="${R0_DEFAULT[$last_kept_before]}" -v s="$TRACK_SPACING" 'BEGIN{printf "%.4f", r0 - s}')"
  fi

  cur_r1="$anchor_r1"
  for i in "${!DYN_KEYS[@]}"; do
    key="${DYN_KEYS[$i]}"
    blk="${KEY_TO_BLOCK[$key]:-}"
    [[ -n "$blk" ]] || continue
    if [[ "$i" -lt "$first_missing_idx" ]]; then
      continue
    fi

    if [[ "${PRESENT[$i]}" -eq 0 ]]; then
      : > "$blk"
      continue
    fi

    # resolve file path first
    pat="$(block_file_pattern "$blk")"
    if [[ -n "$pat" ]]; then
      resolved="$(resolve_first_match "$pat")"
      rewrite_block_file_line "$blk" "$resolved"
    fi

    # assign packed radii
    r1="$cur_r1"
    r0=$(awk -v r1="$r1" -v w="$TRACK_WIDTH" 'BEGIN{printf "%.4f", r1 - w}')
    awk -v newr1="$r1" -v newr0="$r0" '
      BEGIN{r1done=0; r0done=0}
      /^[[:space:]]*r1[[:space:]]*=[[:space:]]*[0-9.]+r/ && !r1done { sub(/[0-9.]+r/, sprintf("%.4fr", newr1)); r1done=1 }
      /^[[:space:]]*r0[[:space:]]*=[[:space:]]*[0-9.]+r/ && !r0done { sub(/[0-9.]+r/, sprintf("%.4fr", newr0)); r0done=1 }
      {print}
    ' "$blk" > "$blk.tmp" && mv "$blk.tmp" "$blk"

    cur_r1=$(awk -v r0="$r0" -v s="$TRACK_SPACING" 'BEGIN{printf "%.4f", r0 - s}')
    INNERMOST_R0="$r0"
  done
fi

# --- Rebuild <plots> region ---
PLOTS_REBUILT="$WORKDIR/plots_rebuilt.txt"
{
  echo "<plots>"
  echo
  for blk in "${BLOCKS[@]}"; do
    if [[ -s "$blk" ]]; then
      cat "$blk"
      echo; echo
    fi
  done
  echo "</plots>"
} > "$PLOTS_REBUILT"

# --- Stitch: PRE + plots + POST ---
cat "$PRE" "$PLOTS_REBUILT" "$POST" > "$WORKDIR/stitched.conf"

# === Resolve <links> file and adjust radius (or drop links if no file) ===
links_pat="$(awk '
  BEGIN{inlinks=0}
  /<links>/{inlinks=1; next}
  /<\/links>/{inlinks=0}
  inlinks && /^[[:space:]]*file[[:space:]]*=/ {
    sub(/^[[:space:]]*file[[:space:]]*=[[:space:]]*/, "", $0)
    print; exit
  }
' "$WORKDIR/stitched.conf" || true)"

if [[ -n "${links_pat// /}" ]]; then
  links_file="$(resolve_first_match "$links_pat")"
  if [[ -z "${links_file// /}" ]]; then
    # No links file -> drop entire <links>...</links> block
    awk '
      BEGIN{inlinks=0}
      /<links>/{inlinks=1; next}
      /<\/links>/{inlinks=0; next}
      { if(!inlinks) print }
    ' "$WORKDIR/stitched.conf" > "$WORKDIR/nolinks.conf"
    mv "$WORKDIR/nolinks.conf" "$WORKDIR/stitched.conf"
    cp "$WORKDIR/stitched.conf" "$OUTFILE"
  else
    # Rewrite file= and adjust radius if INNERMOST_R0 known
    awk -v newfile="$links_file" -v have_r0="${INNERMOST_R0:-}" -v gap="$LINKS_GAP" '
      BEGIN{inlinks=0; file_done=0; radius_done=0}
      /<links>/{inlinks=1}
      /<\/links>/{inlinks=0}
      {
        if(inlinks==1 && file_done==0 && /^[[:space:]]*file[[:space:]]*=/){
          sub(/=.*/, "= " newfile); file_done=1; print; next
        }
        if(inlinks==1 && have_r0 != "" && radius_done==0 && /^[[:space:]]*radius[[:space:]]*=[[:space:]]*[0-9.]+r/){
          newrad = sprintf("%.4f", have_r0 - gap)
          sub(/[0-9.]+r/, newrad "r"); radius_done=1; print; next
        }
        print
      }
    ' "$WORKDIR/stitched.conf" > "$OUTFILE"
  fi
else
  # No <links> block (or no file line) -> write as-is
  cp "$WORKDIR/stitched.conf" "$OUTFILE"
fi

echo "[OK] Wrote $OUTFILE"

# === Optional: patch ideogram.conf with pairwise, if -gap was provided ===
if $GAP; then
  echo "[INFO] -gap set: patching pairwise into: $IDEO"
  if [[ ! -f "$IDEO" ]]; then
    echo "[WARN] ideogram not found at: $IDEO — skipping pairwise patch."
  else
    echo "[INFO] using karyotype: $kary_file"

    # Robustly derive FIRST (first 'chr' row col3) and LAST (last 'chr' row col3)
    FIRST_ID="$(awk 'BEGIN{IGNORECASE=1} {sub(/\r$/,"")} $1=="chr" && $0!~/^#/ && NF>=3 {print $3; exit}' "$kary_file" || true)"
    LAST_ID="$(awk  'BEGIN{IGNORECASE=1} {sub(/\r$/,"")} $1=="chr" && $0!~/^#/ && NF>=3 {id=$3} END{if(id!="") print id}' "$kary_file" || true)"

    if [[ -z "${FIRST_ID:-}" || -z "${LAST_ID:-}" ]]; then
      echo "[ERROR] Could not derive chromosome IDs from: $kary_file" >&2
      echo "[ERROR] (looked for lines with first column == 'chr' and >=3 fields)" >&2
    else
      echo "[INFO] pairwise will be: <pairwise $LAST_ID $FIRST_ID>"
    fi

    pairwise_block="$(
      cat <<EOF
<pairwise ${LAST_ID:-chrZ} ${FIRST_ID:-chrA}>
# spacing between edge chromosomes; computed from karyotype ($kary_file)
spacing = 15r
</pairwise>
EOF
    )"

    # Idempotent patch: replace or insert inside <spacing>. If no <spacing>, append a minimal one.
    awk -v block="$pairwise_block" '
      { sub(/\r$/, "", $0) }                               # normalize CRLF
      BEGIN{ in_spacing=0; in_pair=0; replaced=0; saw_spacing=0 }
      /<spacing[^>]*>/ { in_spacing=1; saw_spacing=1 }
      in_spacing && /<pairwise[[:space:]][^>]*>/ { in_pair=1 }
      in_pair && /<\/pairwise>/ {
        if (!replaced) { print block; replaced=1 }
        in_pair=0; next
      }
      in_pair { next }
      in_spacing && /<\/spacing>/ {
        if (!replaced) { print block; replaced=1 }
        in_spacing=0; print; next
      }
      { print }
      END{
        if (!saw_spacing) {
          print ""
          print "<spacing>"
          print block
          print "</spacing>"
        }
      }
    ' "$IDEO" > "$IDEO.tmp" && mv "$IDEO.tmp" "$IDEO"

    if grep -qE '^<pairwise[[:space:]]' "$IDEO"; then
      echo "[OK] Patched $IDEO with pairwise: ${LAST_ID:-NA} ${FIRST_ID:-NA}"
    else
      echo "[WARN] Pairwise block not found after patch. Dumping <spacing> for debug:"
      awk '{
        sub(/\r$/,"",$0)
        if ($0 ~ /<spacing[^>]*>/) {show=1}
        if (show) print
        if ($0 ~ /<\/spacing>/) {show=0}
      }' "$IDEO"
    fi
  fi
fi
