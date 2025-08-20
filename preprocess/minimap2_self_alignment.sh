#!/usr/bin/env bash
# Usage:
#   ./minimap2_self_alignment.sh -f genome.fa [-p output_prefix] [-t threads]
# Removes rows where (qID==sID && qStart==sStart && qEnd==sEnd). - This needs to be tested before publishing

set -euo pipefail

# === Defaults ===
fasta=""
prefix=""
threads=8

# === Parse flags ===
while getopts "f:p:t:" opt; do
  case "$opt" in
    f) fasta="$OPTARG" ;;
    p) prefix="$OPTARG" ;;
    t) threads="$OPTARG" ;;
    *) echo "❌ Invalid option: -$OPTARG"; exit 1 ;;
  esac
done

# === Check required input ===
if [[ -z "$fasta" ]]; then
  echo "❌ Error: -f <fasta> is required"
  exit 1
fi

# === Auto-determine prefix if not provided ===
if [[ -z "$prefix" ]]; then
  prefix=$(basename "$fasta" .fa)
  prefix="${prefix%.fasta}"
fi

# === Outputs ===
paf="${prefix}.paf"
coords="${prefix}.coords"

echo "🔁 Running minimap2 self-alignment"
echo "  ➤ FASTA:   $fasta"
echo "  ➤ PREFIX:  $prefix"
echo "  ➤ THREADS: $threads"

date

minimap2 -x asm5 -t "$threads" "$fasta" "$fasta" > "$paf"

echo "📄 Converting PAF → coords and filtering exact self-diagonal hits..."
awk 'BEGIN { OFS="\t" }
{
  qstart = $3 + 1;  qend = $4;
  sstart = $8 + 1;  send = $9;

  # Skip exact self hits: same contig and same interval
  if ($1 == $6 && qstart == sstart && qend == send) next;

  print $1, qstart, qend, $6, sstart, send
}' "$paf" > "$coords"

echo "✅ Done."
echo " - $paf"
echo " - $coords"

date
