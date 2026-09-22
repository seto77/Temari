#!/usr/bin/env bash
# reproduce/run.sh — recompute the Temari benchmark and compare it with the published reference CSV.
#
# Target: dataset F v7.0.0 (DOI 10.5281/zenodo.22643468), E0 = 200 keV, 10 channels x 17 s-nodes = 170 values
# (docs/notes/benchmark_spec_2026-08-19.md). The engine is tools/benchmark_recompute.jl; this script only
# prepares the inputs, runs it from an empty working directory (= an empty SCF cache) and classifies the outcome.
#
#   bash reproduce/run.sh [OUTDIR]
#
# Environment (all optional):
#   JULIA               the Julia command, may contain arguments (default "julia"; e.g. "julia +1.11.9")
#   THREADS             Julia threads (default "auto")
#   CHANNELS            space-separated subset of the 10 channels, for a quick setup check (the result is then
#                       labelled SUBSET and is not a benchmark result)
#   TEMARI_NUCLEAR_DATA directory holding the two IAEA tables (default refs/data). Missing tables are downloaded.
#
# Exit status: 0 = every requested channel within the gate (max|dF| <= 1e-5)
#              1 = at least one channel outside the gate
#              2 = inputs or environment (Julia missing, table download failed, table or CSV sha256 differs)
#              3 = the engine stopped without a verdict (a defect, not a result)
set -u

# `pwd -W` gives C:/... under Git Bash (MSYS); Julia on Windows cannot read /c/... passed through the environment
abspath_of() { (cd "$1" && (pwd -W 2>/dev/null || pwd)); }
ROOT=$(abspath_of "$(dirname "$0")/..")
read -r -a JCMD <<< "${JULIA:-julia}"
THREADS=${THREADS:-auto}
ALL_CHANNELS="K_Z6 K_Z14 K_Z26 K_Z50 L1_Z20 L2_Z20 L3_Z79 L3_Z86 M1_Z30 M5_Z86"
CHANNELS=${CHANNELS:-$ALL_CHANNELS}
EXPECT_VERSION=7.0.0
CSV_SHA256=510c1cb7aca0ef0d2df65df5761bd1c87fbe5aedd5f352eef26b921f946b49a1   # tables/F_200keV_preview.csv
DATA=${TEMARI_NUCLEAR_DATA:-$ROOT/refs/data}
# name|url|approved sha256 (the same values as APPROVED_TABLE_SHA256 in src/l1_nucleus_radius.jl)
TABLES=(
  "charge_radii_IAEA_AngeliMarinova.csv|https://www-nds.iaea.org/radii/charge_radii.csv|683d028dae93235376ddf7a345f736561ee7c0dba070091912c034bdf740eee6"
  "iaea_livechart_ground_states.csv|https://nds.iaea.org/relnsd/v1/data?fields=ground_states&nuclides=all|8aee5dc431af1e35fcb49746387b83e927b3c300e7787defbda621a08212c795"
)

die2() { echo "reproduce: $*" >&2; exit 2; }
sha256_of() {
  if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | cut -d' ' -f1
  elif command -v shasum >/dev/null 2>&1; then shasum -a 256 "$1" | cut -d' ' -f1
  else die2 "neither sha256sum nor shasum is available"; fi
}

if [ $# -ge 1 ]; then OUT=$1; mkdir -p "$OUT" || die2 "cannot create $OUT"
else OUT=$(mktemp -d "${TMPDIR:-/tmp}/temari_reproduce.XXXXXX") || die2 "mktemp failed"; fi
OUT=$(abspath_of "$OUT")
[ -e "$OUT/result.jsonl" ] && die2 "$OUT/result.jsonl already exists (use a new OUTDIR)"
WORK="$OUT/work"
mkdir "$WORK" || die2 "$WORK already exists (use a new OUTDIR; the working directory must start empty)"

command -v "${JCMD[0]}" >/dev/null 2>&1 || die2 "Julia not found: ${JCMD[*]} (set JULIA)"
JVER=$("${JCMD[@]}" --startup-file=no -e 'print(VERSION)' 2>/dev/null) || die2 "cannot run ${JCMD[*]}"
case "$JVER" in
  1.11.9|1.12.*) ;;
  *) echo "reproduce: warning: Julia $JVER — the reference was produced with 1.11.9 (1.12.6 gives identical bits here)" >&2 ;;
esac

mkdir -p "$DATA" || die2 "cannot create $DATA"
for t in "${TABLES[@]}"; do
  IFS='|' read -r name url want <<< "$t"
  f="$DATA/$name"
  if [ ! -f "$f" ]; then
    echo "reproduce: downloading $name from IAEA"
    curl -fsSL --retry 2 -m 120 -A "temari-reproduce" -o "$f.part" "$url" || { rm -f "$f.part"; die2 "download failed: $url"; }
    mv "$f.part" "$f"
  fi
  got=$(sha256_of "$f")
  [ "$got" = "$want" ] || die2 "$name has sha256 $got, approved $want.
  IAEA may have updated the table, or the file is not the one Temari was built with. The finite-nucleus
  radii of dataset F v7.0.0 come from the approved bytes, so the benchmark cannot be run against other ones."
done

n_req=$(echo $CHANNELS | wc -w | tr -d " ")
label=BENCHMARK
[ "$CHANNELS" = "$ALL_CHANNELS" ] || label="SUBSET ($n_req of 10 channels; not a benchmark result)"
COMMIT=$(git -C "$ROOT" rev-parse HEAD 2>/dev/null || echo unknown)
{
  echo "label       $label"
  echo "started_utc $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "julia       $JVER  threads=$THREADS"
  echo "commit      $COMMIT"
  echo "csv_sha256  $CSV_SHA256"
  for t in "${TABLES[@]}"; do IFS='|' read -r name url want <<< "$t"; echo "table       $name $want"; done
  echo "channels    $CHANNELS"
} > "$OUT/summary.txt"

t0=$(date +%s)
# shellcheck disable=SC2086
( cd "$WORK" && TEMARI_NUCLEAR_DATA="$DATA" "${JCMD[@]}" --startup-file=no --project="$ROOT" -t "$THREADS" --gcthreads=1 \
    "$ROOT/tools/benchmark_recompute.jl" $CHANNELS \
    --expect-version "$EXPECT_VERSION" --expect-csv-sha256 "$CSV_SHA256" --out "$OUT/result.jsonl" ) \
    > "$OUT/log.txt" 2> "$OUT/err.txt"
rc=$?
secs=$(( $(date +%s) - t0 ))

result=$(grep -E '^RESULT: [0-9]+/[0-9]+ PASS' "$OUT/log.txt" | tail -1)
n_rows=0; [ -f "$OUT/result.jsonl" ] && n_rows=$(grep -c '' "$OUT/result.jsonl")
n_scf=$(grep -c '^\[SCF' "$OUT/log.txt")
{
  echo "seconds     $secs"
  echo "engine_rc   $rc"
  echo "result      ${result:-none}"
  echo "jsonl_rows  $n_rows"
  echo "scf_solved  $n_scf   (SCF runs in this process: > 0 shows the cache started empty)"
} >> "$OUT/summary.txt"

verdict=3
if [ "$rc" -eq 2 ]; then
  verdict=2
elif [ -n "$result" ]; then
  npass=$(sed -E 's/^RESULT: ([0-9]+)\/([0-9]+).*/\1/' <<< "$result")
  ntot=$(sed -E 's/^RESULT: ([0-9]+)\/([0-9]+).*/\2/' <<< "$result")
  if [ "$ntot" -ne "$n_req" ] || [ "$n_rows" -ne "$n_req" ]; then
    verdict=3   # the count does not match what was asked for
  elif [ "$rc" -eq 0 ] && [ "$npass" -eq "$n_req" ]; then
    verdict=0
  elif [ "$rc" -eq 1 ] && [ "$npass" -lt "$n_req" ]; then
    verdict=1
  fi
fi
echo "verdict     $verdict" >> "$OUT/summary.txt"

cat "$OUT/summary.txt"
case $verdict in
  0) echo "reproduce: PASS — $label, every channel within 1e-5 of the reference CSV. Output: $OUT" ;;
  1) echo "reproduce: FAIL — see docs/notes/benchmark_spec_2026-08-19.md section 5.3 for what each kind of mismatch means. Output: $OUT" ;;
  2) echo "reproduce: input or environment problem (see $OUT/err.txt)"; tail -5 "$OUT/err.txt" ;;
  3) echo "reproduce: the engine stopped without a verdict (see $OUT/err.txt)"; tail -5 "$OUT/err.txt" ;;
esac
exit $verdict
