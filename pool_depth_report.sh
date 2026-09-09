#!/usr/bin/env bash
set -euo pipefail

SITES_POS="sites.pos.txt"        # CHROM<TAB>POS
BAMDIR="dedup_rg"
OUT="pool_depth_report.txt"

: > "$OUT"   # truncate/create

for BAM in "$BAMDIR"/*-PL.dedup.bam; do
  NAME=$(basename "$BAM" .dedup.bam)

  {
    echo "=================================================="
    echo "POOL: $NAME"
    echo "=================================================="

    echo "--- summary (raw depth, no MAPQ/BAQ filter) ---"
    samtools depth -a -b "$SITES_POS" "$BAM" \
      | awk '{print $3}' | sort -n \
      | awk '{a[NR]=$1}
         END{
           print "n covered:", NR;
           print "median:  ", a[int(NR/2)];
           print "90th pct:", a[int(NR*0.9)];
           print "max:     ", a[NR];
           s=0; for(i=1;i<=NR;i++)s+=a[i]; print "mean:    ", s/NR}'

    echo ""
    echo "--- distribution (filtered: -q 20 -Q 30) ---"
    samtools depth -a -b "$SITES_POS" -q 20 -Q 30 "$BAM" \
      | awk '{d=$3;
          if(d==0)b0++;
          else if(d<5)b1++;
          else if(d<10)b2++;
          else if(d<20)b3++;
          else if(d<50)b4++;
          else if(d<100)b5++;
          else b6++}
        END{
          print "0:      "b0+0;
          print "1-4:    "b1+0;
          print "5-9:    "b2+0;
          print "10-19:  "b3+0;
          print "20-49:  "b4+0;
          print "50-99:  "b5+0;
          print ">=100:  "b6+0}'

    echo ""
  } >> "$OUT"

  echo "done: $NAME"   # progress to terminal
done

echo ""
echo "report written to $OUT"