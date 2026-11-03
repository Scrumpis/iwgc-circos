#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C
PNG_REQUIRED=false

# ---- Inputs ----
CONF="iwgc_circos/iwgc_circos.conf"
SVG="iwgc_circos/tmp/iwgc_circos.svg"
IDEO="iwgc_circos/ideogram.conf"
OUT="iwgc_circos/tmp/iwgc_circos.labeled.svg"
PNG_OUT=""
FONT_TTF="fonts/LiberationSans-Bold.ttf"

# Placement & sizing
DX_FRAC="-0.15"        # SVG dx = -DX_FRAC * font-size
DY_FRAC="0"        # down nudge as fraction of font-size
THETA_DEG=""          # if empty, default to -90° (top)

usage(){ cat <<EOF
Usage: $0 [-conf FILE] [-svg FILE] [-ideo FILE] [-out FILE] [-png FILE] [-font FILE]
          [--dx-frac F] [--dy-frac F] [--theta-deg DEG]
EOF
exit 0; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    -conf|-svg|-ideo|-out|-png|-font|--dx-frac|--dy-frac|--theta-deg)
      [[ $# -ge 2 ]] || { echo "[ERROR] $1 requires a value" >&2; exit 1; };;
  esac
  case "$1" in
    -conf) CONF="$2"; shift 2;;
    -svg)  SVG="$2"; shift 2;;
    -ideo) IDEO="$2"; shift 2;;
    -out)  OUT="$2"; shift 2;;
    -png)  PNG_REQUIRED=true; PNG_OUT="$2"; shift 2;;
    -font) FONT_TTF="$2"; shift 2;;
    --dx-frac) DX_FRAC="$2"; shift 2;;
    --dy-frac) DY_FRAC="$2"; shift 2;;
    --theta-deg) THETA_DEG="$2"; shift 2;;
    -h|--help) usage;;
    *) echo "[ERROR] Unknown arg: $1" >&2; exit 1;;
  esac
done
for number in "$DX_FRAC" "$DY_FRAC" "${THETA_DEG:--90}"; do
  [[ "$number" =~ ^[-+]?[0-9]+([.][0-9]+)?$ ]] || { echo "[ERROR] label offsets/angle must be numeric" >&2; exit 1; }
done
[[ "$SVG" != "$OUT" ]] || { echo "[ERROR] SVG input and output must differ" >&2; exit 1; }
[[ -f "$SVG"  ]] || { echo "[ERROR] SVG not found: $SVG" >&2; exit 1; }
python3 - "$SVG" "$OUT" <<'PYXML'
import sys
from pathlib import Path
import xml.etree.ElementTree as ET
if Path(sys.argv[1]).resolve() == Path(sys.argv[2]).resolve():
    raise SystemExit('[ERROR] SVG input and output resolve to the same file')
root=ET.parse(sys.argv[1]).getroot()
if root.tag != '{http://www.w3.org/2000/svg}svg':
    raise SystemExit('[ERROR] Input must be an SVG document')
if any(x.get('id') == 'iwgc-gap-labels' for x in root.iter()):
    raise SystemExit('[ERROR] Input SVG already contains gap labels; use the original')
PYXML
[[ -f "$FONT_TTF" ]] || echo "[WARN] Font not found at $FONT_TTF — falling back to Liberation Sans/sans-serif."
[[ -n "${PNG_OUT:-}" ]] || PNG_OUT="${OUT%.svg}.png"

WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT
numonly(){ sed -E 's/[^0-9.+-eE]//g'; }
one_line(){ tr -d '\r' | gawk 'NR==1{print; exit}'; }

# ---------- Solve true center from two ideogram chords ----------
center_from_ideogram_chords(){
  local svg="$1"
  gawk '
    BEGIN{RS="</g>"; ORS=""}
    {
      blk=$0
      if (blk !~ /<g[^>]*id="ideograms"/) next
      n=0
      # use two paths; from each, take M x,y and the FIRST A ... end x,y
      while (match(blk, /<path[^>]*d="([^"]+)"/, p)) {
        d=p[1]
        if (match(d, /[mM][[:space:]]*([-0-9.+eE]+)[ ,]([-0-9.+eE]+)/, m) &&
            match(d, /[aA][^A-Za-z]*[-0-9.+eE]+[ ,][-0-9.+eE]+[^A-Za-z]*[, ][01][ ,]*,[01][ ,]*([-0-9.+eE]+)[ ,]([-0-9.+eE]+)/, a)) {
          x1=m[1]+0; y1=m[2]+0; x2=a[1]+0; y2=a[2]+0
          chords[++n]=x1" "y1" "x2" "y2
          if (n==2) break
        }
        blk=substr(blk, RSTART+RLENGTH)
      }
      if (n<2) { print "NaN NaN"; next }
      split(chords[1], c1, " "); split(chords[2], c2, " ")
      mx1=(c1[1]+c1[3])/2; my1=(c1[2]+c1[4])/2; nx1=-(c1[4]-c1[2]); ny1=(c1[3]-c1[1])
      mx2=(c2[1]+c2[3])/2; my2=(c2[2]+c2[4])/2; nx2=-(c2[4]-c2[2]); ny2=(c2[3]-c2[1])
      det = nx1*(-ny2) - ny1*(-nx2)
      if (det==0) { print "NaN NaN"; next }
      dx=mx2-mx1; dy=my2-my1
      t = (dx*(-ny2) - dy*(-nx2)) / det
      cx = mx1 + t*nx1; cy = my1 + t*ny1
      printf "%.6f %.6f", cx, cy
      next
    }
  ' "$svg"
}

# ---------- Read ideogram radii (px) from A rx,ry in <g id="ideograms"> ----------
ideo_mid_thick_from_svg(){
  local svg="$1"
  gawk '
    BEGIN{RS="</g>"; ORS=""; have=0; minr=1e18; maxr=0}
    {
      blk=$0
      if (blk !~ /<g[^>]*id="ideograms"/) next
      while (match(blk, /[aA][[:space:]]*([-0-9.+eE]+)[ ,]([-0-9.+eE]+)/, a)) {
        r=a[1]+0; have=1
        if (r<minr) minr=r; if (r>maxr) maxr=r
        blk=substr(blk, RSTART+RLENGTH)
      }
    }
    END{
      if (!have){ print "NaN NaN"; exit }
      mid=(minr+maxr)/2; th=maxr-minr
      printf "%.6f %.6f", mid, th
    }
  ' "$svg"
}

# ---------- NEW: Read all plot ring mids/thickness (px) from <g id="plotN"> ----------
extract_plot_mids_from_svg(){
  local svg="$1"
  gawk '
    function upd_r(line,  a) {
      while (match(line, /[aA][[:space:]]*([-0-9.+eE]+)[ ,]([-0-9.+eE]+)/, a)) {
        r=a[1]+0
        if (r<minr) minr=r
        if (r>maxr) maxr=r
        have=1
        line=substr(line, RSTART+RLENGTH)
      }
    }
    BEGIN{
      FS="\n"; RS="\n"
      inplot=0; depth=0; idx=-1
      have=0; minr=1e18; maxr=0
    }
    {
      line=$0

      if (!inplot) {
        # handle id="plot0" OR id="plot_0"
        if (match(line, /<g[^>]*id="plot_?([0-9]+)"/, m)) {
          inplot=1; idx=m[1]+0
          have=0; minr=1e18; maxr=0
          # count opens/closes on this same line
          opens = gsub(/<g[ >]/, "", line)
          closes = gsub(/<\/g>/, "", line)
          depth = opens - closes
          upd_r(line)
        }
        next
      }

      # inside a plotN group
      upd_r(line)
      opens = gsub(/<g[ >]/, "", line)
      closes = gsub(/<\/g>/, "", line)
      depth += opens - closes

      if (depth<=0) {
        if (have) {
          mid=(minr+maxr)/2; th=maxr-minr
          printf "%d\t%.6f\t%.6f\n", idx, mid, th
        }
        # reset for next plot
        inplot=0; depth=0; idx=-1
        have=0; minr=1e18; maxr=0
      }
    }
    END{
      if (inplot && have) {
        mid=(minr+maxr)/2; th=maxr-minr
        printf "%d\t%.6f\t%.6f\n", idx, mid, th
      }
    }
  ' "$svg" | sort -n -k1,1
}

# ---------- Get gap endpoints (first ideogram start, last ideogram end) ----------
gap_endpoints_from_svg(){
  local svg="$1"
  gawk '
    BEGIN{ingrp=0; have_first=0; x2=""; y2=""}
    /<g[^>]*id="ideograms"/ {ingrp=1}
    ingrp && /<\/g>/ {ingrp=0}
    ingrp && match($0, /<path[^>]*d="([^"]+)"/, p){
      d = p[1]
      # first path start "M x,y"
      if (!have_first && match(d, /[mM][[:space:]]*([-0-9.+eE]+)[ ,]([-0-9.+eE]+)/, m)) {
        x1=m[1]+0; y1=m[2]+0; have_first=1
      }
      # update to the end x,y of the FIRST arc "A rx,ry rot laf,sf x,y"
      if (match(d, /[aA][^A-Za-z]*[-0-9.+eE]+[ ,][-0-9.+eE]+[^A-Za-z]*[, ][01][ ,]*,[01][ ,]*([-0-9.+eE]+)[ ,]([-0-9.+eE]+)/, a)) {
        x2=a[1]+0; y2=a[2]+0
      }
    }
    END{
      if (have_first && x2!="") {
        printf "%.6f %.6f %.6f %.6f", x1, y1, x2, y2
      } else {
        print "NaN NaN NaN NaN"
      }
    }
  ' "$svg"
}

alpha_label(){ local n="$1" s="" base=26 a=97; n=$((n+1)); while (( n>0 )); do local r=$(( (n-1)%base )); s="$(printf \\$(printf "%03o" $((a+r))))${s}"; n=$(( (n-1)/base )); done; printf "%s." "$s"; }

# ---------- PNG ----------
svg_to_png(){ local in="$1" out="$2"
  if command -v rsvg-convert >/dev/null 2>&1; then rsvg-convert -o "$out" "$in"
  elif command -v inkscape >/dev/null 2>&1; then inkscape "$in" --export-type=png --export-filename="$out" >/dev/null 2>&1
  elif command -v magick >/dev/null 2>&1; then magick -density 300 "$in" -background none -flatten "$out"
  elif command -v convert >/dev/null 2>&1; then convert -density 300 "$in" -background none -flatten "$out"
  else echo "[WARN] No SVG->PNG converter found (install rsvg-convert, inkscape, or ImageMagick)." >&2; return 1; fi
}

# ===================== main =====================
# 1) True center
read -r CX CY <<<"$(center_from_ideogram_chords "$SVG" || true)"
if [[ -z "${CX:-}" || "$CX" == "NaN" ]]; then
  # fall back to root box center if needed
  root="$(grep -m1 -oE '<svg[^>]*>' "$SVG" || true)"
  W="$(printf '%s' "$root" | grep -oE 'width="[^"]*"'  -m1 | cut -d'"' -f2 | numonly)"
  H="$(printf '%s' "$root" | grep -oE 'height="[^"]*"' -m1 | cut -d'"' -f2 | numonly)"
  CX="$(gawk -v w="$W" 'BEGIN{printf "%.6f", (w==""?3000:w)/2}')"
  CY="$(gawk -v h="$H" 'BEGIN{printf "%.6f", (h==""?3000:h)/2}')"
  echo "[WARN] Could not solve ideogram center; using root center ($CX,$CY)"
else
  echo "[INFO] Center from ideograms: ($CX,$CY)"
fi

# 2) Ideogram mid & thickness (px)
read -r IDEO_MID_PX IDEO_TH_PX <<<"$(ideo_mid_thick_from_svg "$SVG" || true)"
if [[ -z "${IDEO_MID_PX:-}" || "$IDEO_MID_PX" == "NaN" ]]; then
  echo "[ERROR] Could not read ideogram radii from SVG."; exit 1
fi
echo "[INFO] Ideogram mid/thickness (px): $IDEO_MID_PX / $IDEO_TH_PX"

# 3) Plot mids/thickness (px) from SVG
PLOTS_TXT="$WORK/plots_px.txt"
extract_plot_mids_from_svg "$SVG" > "$PLOTS_TXT"
if [[ ! -s "$PLOTS_TXT" ]]; then
  echo "[ERROR] No <g id=\"plotN\"> radii found in SVG."; exit 1
fi
echo "[INFO] Plot mids/thickness:"
cat "$PLOTS_TXT"

# 4) Compose label radii (px): a. (ideogram) + each plotN in order
ALL_MIDS_PX="$WORK/mids_px.txt"
ALL_H_PX="$WORK/heights_px.txt"

# ideogram mid goes first
printf "%.6f\n" "$IDEO_MID_PX" > "$ALL_MIDS_PX"
# size of a.: use OUTERMOST plot thickness so it matches others
OUTER_TH_PX="$(gawk 'NR==1{print $3; exit}' "$PLOTS_TXT")"
printf "%.6f\n" "$OUTER_TH_PX" > "$ALL_H_PX"

# then each plot mid & its own thickness
gawk '{printf "%.6f\n", $2}' "$PLOTS_TXT" >> "$ALL_MIDS_PX"
gawk '{printf "%.6f\n", $3}' "$PLOTS_TXT" >> "$ALL_H_PX"

# 5) Angle
if [[ -n "$THETA_DEG" ]]; then
  GAP_THETA="$(gawk -v d="$THETA_DEG" 'BEGIN{pi=atan2(0,-1); printf "%.10f", d*pi/180.0}')"
  echo "[INFO] Gap angle override (deg): $THETA_DEG -> $GAP_THETA"
else
  GAP_THETA="$(gawk 'BEGIN{pi=atan2(0,-1); printf "%.10f", -90*pi/180.0}')"
  echo "[INFO] Gap angle defaulting to -90°: $GAP_THETA"
fi

# 5.5) NEW: compute dynamic gap X from ideogram first-start and last-end
read -r X1 Y1 X2 Y2 <<<"$(gap_endpoints_from_svg "$SVG" || true)"
if [[ -z "${X1:-}" || "$X1" == "NaN" || -z "${X2:-}" || "$X2" == "NaN" ]]; then
  echo "[WARN] Could not parse ideogram gap endpoints; falling back to radial X."
  XGAP_MID=""
else
  XGAP_MID="$(gawk -v a="$X1" -v b="$X2" 'BEGIN{printf "%.6f", (a+b)/2.0}')"
  echo "[INFO] Gap endpoints: first-start=($X1,$Y1) last-end=($X2,$Y2); Xmid=$XGAP_MID"
fi

# 6) Build and inject labels
INJ="$WORK/inject.svgfrag"
{
  echo '<style type="text/css"><![CDATA['
  if [[ -f "$FONT_TTF" ]]; then
    FONT_DATA=$(base64 < "$FONT_TTF" | tr -d '\n')
    echo "  @font-face { font-family: 'GapLabelBold'; src: url('data:font/ttf;base64,$FONT_DATA') format('truetype'); font-weight: bold; }"
  fi
  echo '  .iwgc-gap-label { font-family: "GapLabelBold", "Liberation Sans", sans-serif; text-anchor: middle; dominant-baseline: middle; }'
  echo ']]></style>'
  echo '<g id="iwgc-gap-labels" aria-label="interior track labels" role="group">'

  idx=0
  paste "$ALL_MIDS_PX" "$ALL_H_PX" | while read -r mid_px h_px; do
    fsz="$(gawk -v h="$h_px" 'BEGIN{v=h*0.6; if(v<10)v=10; if(v>48)v=48; printf "%.2f", v}')"

    # y stays dynamic per ring (unchanged)
    ypix="$(gawk -v cy="$CY" -v r="$mid_px" -v th="$GAP_THETA" 'BEGIN{printf "%.6f", cy + r*sin(th)}')"

    # x is dynamic gap-midpoint (if available), else previous radial fallback
    if [[ -n "${XGAP_MID:-}" && -z "$THETA_DEG" ]]; then
      xpix="$XGAP_MID"
    else
      xpix="$(gawk -v cx="$CX" -v r="$mid_px" -v th="$GAP_THETA" 'BEGIN{printf "%.6f", cx + r*cos(th)}')"
    fi

    dxpx="$(gawk -v f="$fsz" -v frac="$DX_FRAC" 'BEGIN{printf "%.3f", -frac*f}')"
    dypx="$(gawk -v f="$fsz" -v frac="$DY_FRAC" 'BEGIN{printf "%.3f",  frac*f}')"
    printf '  <text class="iwgc-gap-label" x="%.6f" y="%.6f" dx="%.3f" dy="%.3f" font-size="%s">%s</text>\n' \
      "$xpix" "$ypix" "$dxpx" "$dypx" "$fsz" "$(alpha_label "$idx")"
    idx=$((idx+1))
  done

  echo '</g>'
} > "$INJ"

gawk -v injfile="$INJ" '
  BEGIN{ while ((getline l < injfile) > 0) inj = inj l "\n"; close(injfile) }
  /<\/svg>/ { print inj; print; next }
  { print }
' "$SVG" > "$OUT"

echo "[OK] Wrote labeled SVG -> $OUT"

# Optional PNG
if svg_to_png "$OUT" "$PNG_OUT"; then
  echo "[OK] Wrote PNG -> $PNG_OUT"
else
  echo "[WARN] PNG export unavailable or failed." >&2
  if $PNG_REQUIRED; then exit 1; fi
fi
