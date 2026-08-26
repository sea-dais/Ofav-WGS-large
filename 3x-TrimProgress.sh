# Trimming job progress
ls -lhS trimmed/
ls -l pylauncher_tmp3437752/success* | wc

# how much was trimmed
grep -A2 "passing filters" pylauncher_tmp3437752/out0     # or wherever your cutadapt logs are
grep "Total basepairs" pylauncher_tmp3437752/out*
grep -E "too short|Pairs written|basepairs|Quality-trimmed" pylauncher_tmp3437752/out0
grep "too short" pylauncher_tmp3437752/out* | grep -oE "\([0-9.]+%\)"
## most samples around 20% or less discarded for too short
grep -l "too short.*31\.[0-9]%" pylauncher_tmp3437752/out*
# both adults have 30% discarded reads that were too short 

cd pylauncher_tmp3437752
for f in out*; do
  n=${f#out}
  sample=$(grep -oE 'trimmed/[^ ]*_R1.trim.gz' "$f" | sed 's|trimmed/||; s|_R1.trim.gz||')
  status=$([ -f "success$n" ] && echo OK || echo ???)
  echo "$n  $sample  $status"
done

cd pylauncher_tmp3437752
for f in out*; do
  n=${f#out}
  [ -f "success$n" ] && continue          # skip OK tasks
  # sample base name, e.g. F3-h-T3_S29_L001
  sample=$(grep -oE 'trimmed/[^ ]*_R1.trim.gz' "$f" | sed 's|trimmed/||; s|_R1.trim.gz||')
  echo "=== task $n: $sample ==="
  ls -lh ../files/${sample}_R1_001.fastq.gz ../trimmed/${sample}_R1.trim.gz 2>/dev/null
done