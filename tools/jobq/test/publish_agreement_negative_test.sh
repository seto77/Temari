#!/bin/bash
# publish_agreement_negative_test.sh — Python / checker を起動できない衝突を一致扱いしない。
# 本物の共有・worker は使わず、scratch の 1 票と最小 gen_production stub だけを走らせる。
set -u

here=$(cd "$(dirname "$0")" && pwd); jobq_dir=$(cd "$here/.." && pwd)
worker=${JOBQ_WORKER_UNDER_TEST:-$jobq_dir/worker.sh}
JULIA=${JOBQ_JULIA_CHANNEL:-+1.11.9}
scratch=$(mktemp -d "${TMPDIR:-/tmp}/jobq-publish-negative.XXXXXX") || exit 1
trap 'rm -rf "$scratch"' EXIT
ROOT=$scratch/root; SPOOL=$ROOT/spool; LOCAL=$scratch/local
WID=publish-negative; CAMPAIGN=publish_negative; BASE=${CAMPAIGN}_000001.e001; BASE2=${CAMPAIGN}_000001.e002
CODE_SHA=aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa; CODE16=${CODE_SHA:0:16}
npass=0; nfail=0
check() { local d=$1; shift; if "$@" >/dev/null 2>&1; then npass=$((npass+1)); printf 'PASS  %s\n' "$d"; else nfail=$((nfail+1)); printf 'FAIL  %s\n' "$d"; fi; }
nfiles() { find "$1" -maxdepth 1 -type f -name "$2" 2>/dev/null | wc -l | tr -d ' '; }

mkdir -p "$SPOOL/queue" "$SPOOL/running" "$SPOOL/results/$CAMPAIGN" "$SPOOL/done" "$SPOOL/failed" \
         "$SPOOL/control" "$SPOOL/hosts" "$SPOOL/campaigns" "$LOCAL/setup" "$LOCAL/code/$CODE16/src" \
         "$LOCAL/code/$CODE16/tools" "$LOCAL/logs" "$LOCAL/state" "$LOCAL/work"
cat > "$LOCAL/worker.conf" <<CONF
JOBQ_ROOT=$ROOT
JOBQ_SPOOL=$SPOOL
JOBQ_LOCAL=$LOCAL
WORKER_ID=$WID
SLOTS=1
THREADS=1
PYTHON=/definitely/not/a/python.exe
STALL_SECONDS=60
MAX_ATTEMPTS=1
HEARTBEAT_INTERVAL=1
RETRY_BACKOFF=1
DEGRADED_SLEEP=1
CONF
cp "$jobq_dir/../agreement_check.py" "$LOCAL/code/$CODE16/tools/agreement_check.py"
cat > "$LOCAL/code/$CODE16/src/gen_production.jl" <<'JULIA_STUB'
out = ARGS[findfirst(==("--out"), ARGS) + 1]
mkpath(out)
write(joinpath(out, "F_M5_Z30.json"), """{"dataset_version":"6.0.0","generator_source_fingerprint":"0123456789abcdef","spec_sha256":"749fadc500000000000000000000000000000000000000000000000000000000","F":[1.0,0.5000000000000001]}\n""")
println("gen_production: 1/1 チャネル (lane 0/1, tags=M5, HIGH, スレッド 1)")
println("完了: 1 計算 / 0 skip (既存)")
JULIA_STUB
cat > "$SPOOL/queue/$BASE.json" <<TICKET
{"schema":1,"campaign":"$CAMPAIGN","jobseq":1,"claim_epoch":1,"task":"temari.gen_production","code_sha256":"$CODE_SHA","code_commit":"","args":{"tags":["M5"],"lane":0,"lane_count":1,"profile":"v6_high","expected_dataset_version":"6.0.0"},"created_utc":"2026-08-23T00:00:00Z","issued_by":"negative-test"}
TICKET
cat > "$SPOOL/queue/$BASE2.json" <<TICKET
{"schema":1,"campaign":"$CAMPAIGN","jobseq":1,"claim_epoch":2,"task":"temari.gen_production","code_sha256":"$CODE_SHA","code_commit":"","args":{"tags":["M5"],"lane":0,"lane_count":1,"profile":"v6_high","expected_dataset_version":"6.0.0"},"created_utc":"2026-08-23T00:00:01Z","issued_by":"negative-test"}
TICKET
published="$SPOOL/results/$CAMPAIGN/F_M5_Z30.json"
printf '%s\n' '{"dataset_version":"6.0.0","generator_source_fingerprint":"0123456789abcdef","spec_sha256":"749fadc500000000000000000000000000000000000000000000000000000000","F":[1.0,0.5]}' > "$published"
published_sha=$(sha256sum "$published" | cut -c1-64)
printf '{"schema":1,"result_sha256":"%s","hostname":"prior-host","cpu":"prior-cpu","worker_id":"prior","owner":"prior-s0-b1","code_sha256":"%s"}\n' \
       "$published_sha" "$CODE_SHA" > "$published.manifest.json"

env JOBQ_LOCAL="$LOCAL" JOBQ_ROOT="$ROOT" JOBQ_SPOOL="$SPOOL" JOBQ_MAX_IDLE_LOOPS=1 \
    JOBQ_QUEUECTL="$jobq_dir/queuectl.jl" JOBQ_PIN="$jobq_dir/PIN.json" JOBQ_WATCH_INTERVAL=1 \
    bash "$worker" 0 > "$scratch/worker.log" 2>&1
wrc=$?
check "worker は判定不能を処理して終了する" test "$wrc" -eq 0
check "判定不能の票は DONE にならない" test "$(nfiles "$SPOOL/done/$CAMPAIGN" '*.json')" -eq 0
check "判定不能の e001/e002 はともに FAIL receipt になる" test "$(nfiles "$SPOOL/failed/$CAMPAIGN" "${CAMPAIGN}_000001.e00[12].*.json")" -eq 2
receipt=$(find "$SPOOL/failed/$CAMPAIGN" -maxdepth 1 -type f -name "$BASE.*.json" | head -1)
check "FAIL reason は could not judge と明記する" grep -q 'could not judge' "$receipt"
check "候補は dup でなく unjudged に隔離する" test "$(nfiles "$SPOOL/failed/$CAMPAIGN/unjudged" 'F_M5_Z30.json.*')" -ge 1
records_are_distinct() {
  local refs n ref
  refs=$(grep -h -o 'failed/[^" ]*/unjudged/[^" ]*\.agreement\.json' "$SPOOL/failed/$CAMPAIGN"/${CAMPAIGN}_000001.e00[12].*.json 2>/dev/null | sort -u)
  n=$(printf '%s\n' "$refs" | sed '/^$/d' | wc -l | tr -d ' ')
  [ "$n" -eq 2 ] || return 1
  [ "$(nfiles "$SPOOL/failed/$CAMPAIGN/unjudged" '*.agreement.json')" -eq 2 ] || return 1
  while IFS= read -r ref; do [ -z "$ref" ] || { [ -f "$SPOOL/$ref" ] && grep -q '"verdict": "unjudged"' "$SPOOL/$ref"; } || return 1; done <<EOF
$refs
EOF
}
check "e001/e002 receipt は上書きされない別々の機械可読 record を指す" records_are_distinct
check "測って不一致の dup 証拠には混ぜない" test "$(nfiles "$SPOOL/failed/$CAMPAIGN/dup" '*')" -eq 0
check "先客の成果物は上書きしない" test "$(sha256sum "$published" | cut -c1-64)" = "$published_sha"

printf 'publish_agreement_negative_test: PASS %d / FAIL %d\n' "$npass" "$nfail"
[ "$nfail" -eq 0 ]
