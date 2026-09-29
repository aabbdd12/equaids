#!/bin/bash
# bf_parallel.sh -- the brute-force influence function in K parallel Stata
# processes (one block of households each), then the comparison with the
# analytic standard errors.
#   bash bf_parallel.sh spec tag n K
#     spec  the specification do-file (spec.do, see bf_worker.do)
#     tag   a name for the output files (raw/bf_<tag>_<k>.mmat)
#     n     the number of households
#     K     the number of processes
# STATA: the path of the Stata executable (environment variable; default
# StataNow 19 SE).
# Each block runs in its own Stata process at BelowNormal priority; the
# wrappers and logs go to run/ (one log per block, the combination in
# run/bfc_<tag>.log).
SPEC="$1"; TAG="$2"; N="$3"; K="$4"
STATA="${STATA:-C:\Program Files\StataNow19\StataSE-64.exe}"
BF="$(cd "$(dirname "$0")" && pwd)"
RUN="$BF/run"
mkdir -p "$RUN" "$BF/raw" "$BF/out"
BFM="$(cygpath -m "$BF")"

launch() {  # launch <wrapper name>: runs run/<name>.do, waits for it
  powershell -NoProfile -Command "\$p = Start-Process -FilePath '$STATA' -ArgumentList '/e do \"$(cygpath -w "$RUN/$1.do")\"' -WorkingDirectory '$(cygpath -w "$RUN")' -PassThru; \$p.PriorityClass = 'BelowNormal'; \$p.WaitForExit()" > /dev/null
}

s0=$(date +%s)
size=$(( (N + K - 1) / K ))
for k in $(seq 1 "$K"); do
  first=$(( (k - 1) * size + 1 ))
  last=$(( k * size )); [ "$last" -gt "$N" ] && last=$N
  printf 'cd "%s"\ndo bf_worker.do %s %d %d %s %d\n' "$BFM" "$SPEC" "$first" "$last" "$TAG" "$k" > "$RUN/bfw_${TAG}_$k.do"
  rm -f "$RUN/bfw_${TAG}_$k.log" "$BF/raw/bf_${TAG}_$k.mmat"
  launch "bfw_${TAG}_$k" &
done
wait
s1=$(date +%s)
for k in $(seq 1 "$K"); do
  grep -h "^block $k:" "$RUN/bfw_${TAG}_$k.log" | tail -1
  [ -f "$BF/raw/bf_${TAG}_$k.mmat" ] || echo "block $k: NO OUTPUT (see run/bfw_${TAG}_$k.log)"
done
printf 'cd "%s"\ndo bf_combine.do %s %s %d\n' "$BFM" "$SPEC" "$TAG" "$K" > "$RUN/bfc_$TAG.do"
rm -f "$RUN/bfc_$TAG.log"
launch "bfc_$TAG"
grep -h "^re-runs not converged" "$RUN/bfc_$TAG.log"
sed -n '/^brute force: /,/^end of do-file/p' "$RUN/bfc_$TAG.log" | grep -v -E '^end of do-file|^[.:] |^> |^$'
echo "wall time: $((s1 - s0)) s for the $K blocks"
