#!/usr/bin/env bash
# Execute this INSIDE the newly built full dependency image, from the source checkout.
# Usage: bash tools/check_container.sh /fresh/output/directory
set -euxo pipefail
repo=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
out=${1:?Usage: check_container.sh /fresh/output/directory}
[[ "$out" == /* ]] || { echo "Output path must be absolute" >&2; exit 1; }
mkdir -p "$out"
[[ -z "$(ls -A "$out")" ]] || { echo 'Output directory must be empty' >&2; exit 1; }
cd "$repo"
exec > >(tee "$out/checks.log") 2>&1
trap 'rc=$?; printf "check_container_exit_code=%s\n" "$rc"' EXIT
for tool in samtools bedtools circos gawk python3 quartet.py tidk rsvg-convert; do
  command -v "$tool" >/dev/null || { echo "Missing full-image dependency: $tool" >&2; exit 1; }
done
test ! -e fonts/ArialBold.ttf
test -s fonts/LiberationSans-Bold.ttf
test -s fonts/LICENSE.LiberationSans
cmp fonts/LiberationSans-Bold.ttf /usr/local/share/fonts/iwgc/LiberationSans-Bold.ttf
cmp fonts/LICENSE.LiberationSans /usr/share/doc/iwgc-circos-font/LICENSE.LiberationSans
cmp patches/circos-zero-width-svg.patch /usr/share/doc/iwgc-circos/circos-zero-width-svg.patch
cmp patches/NOTICE.md /usr/share/doc/iwgc-circos/NOTICE.md
cmp patches/LICENSE.Circos /usr/share/doc/iwgc-circos/LICENSE.Circos
for script in iwgc_circos_tracks.sh create_configs.sh add_legend.sh preprocess/minimap2_self_alignment.sh tools/check_container.sh; do
  bash -n "$script"
done
python3 tools/validate_skill.py
python3 -m unittest discover -s tests -v
python3 tests/independent_check.py --output "$out/independent"
python3 tests/smoke.py --output "$out/smoke"
# Require the converter result as well as the structural SVG checks in smoke.py.
python3 - "$out/smoke/iwgc_circos/tmp/labeled.png" <<'PY'
import struct, sys
from pathlib import Path
p=Path(sys.argv[1]); data=p.read_bytes()
assert data[:8] == b'\x89PNG\r\n\x1a\n'
assert struct.unpack('>II',data[16:24]) == (3000,3000)
PY
# The observed old-image SVG/native histogram mismatch must not pass publication.
python3 tests/check_render_gaps.py "$out/smoke/iwgc_circos/tmp/iwgc_circos.png" "$out/smoke/iwgc_circos/tmp/labeled.png"
# Execute the human guide's minimal genes/GC example, with exactly its options.
example="$out/example"
mkdir -p "$example/iwgc_circos/tmp" "$example/iwgc_circos_data" "$example/fonts"
cp iwgc_circos_template.config "$example/"
cp iwgc_circos/*.conf "$example/iwgc_circos/"
cp fonts/LiberationSans-Bold.ttf fonts/LICENSE.LiberationSans "$example/fonts/"
cp tests/fixtures/genome.fa tests/fixtures/genes.gff3 "$example/"
(
  cd "$example"
  bash "$repo/iwgc_circos_tracks.sh" genome.fa -gene genes.gff3 -gc -window 4 -keep-temp -out iwgc_circos_data
  bash "$repo/create_configs.sh" -gap
  circos -conf iwgc_circos/iwgc_circos.conf -outputdir iwgc_circos/tmp
  bash "$repo/add_legend.sh"
)
python3 - "$example" <<'PY'
from pathlib import Path
import sys, math, xml.etree.ElementTree as ET
p=Path(sys.argv[1]); ns='{http://www.w3.org/2000/svg}'
for suffix,want in [('gene_coverage',[math.sqrt(.75),1,math.sqrt(.5),0,0,0]),('gc',[.5,.5,.5,1,1,0])]:
 got=[float(line.split()[3]) for line in (p/f'iwgc_circos_data/genome.fa_{suffix}.circos').read_text().splitlines()]
 assert len(got)==6 and all(abs(a-b)<1e-12 for a,b in zip(got,want)), (suffix,got)
r=ET.parse(p/'iwgc_circos/tmp/iwgc_circos.labeled.svg').getroot()
assert [e.text for e in r.iter(ns+'text') if e.get('class')=='iwgc-gap-label']==['a.','b.','c.']
PY
{
  python3 --version
  samtools --version | sed -n '1,2p'
  bedtools --version
  circos -version
  tidk --version
  quartet.py -h
  rsvg-convert --version
  dpkg-query -W
  if command -v conda >/dev/null; then conda list --explicit; fi
} > "$out/dependencies.txt"
echo 'Full-image suite, independent arithmetic, rendered smoke, labeled PNG, and human example: PASS'
