#!/bin/bash
# publish_agreement_mutation_test.sh — fail-closed の 1 行を無効化すると負のテストが落ちることを実演。
set -u

here=$(cd "$(dirname "$0")" && pwd); jobq_dir=$(cd "$here/.." && pwd)
src=$jobq_dir/worker.sh
scratch=$(mktemp -d "${TMPDIR:-/tmp}/jobq-publish-mutant.XXXXXX") || exit 1
trap 'rm -rf "$scratch"' EXIT
mutant=$scratch/worker-mutant.sh; out=$scratch/negative.out
anchor='return 4 # AGREEMENT_UNJUDGED_FAIL_CLOSED'
n=$(grep -cF "$anchor" "$src" 2>/dev/null || true)
if [ "$n" -ne 1 ]; then
  printf 'FAIL  mutation anchor count is %s, expected exactly 1\n' "$n"
  exit 1
fi
sed 's/return 4 # AGREEMENT_UNJUDGED_FAIL_CLOSED/return 0 # MUTANT_ACCEPTS_UNJUDGED/' "$src" > "$mutant"
if ! grep -qF 'return 0 # MUTANT_ACCEPTS_UNJUDGED' "$mutant"; then
  printf 'FAIL  mutant marker is absent after replacement\n'
  exit 1
fi
if JOBQ_WORKER_UNDER_TEST="$mutant" bash "$here/publish_agreement_negative_test.sh" > "$out" 2>&1; then
  printf 'FAIL  negative test passed against a mutant that accepts an unjudged collision\n'
  tail -n 20 "$out"
  exit 1
fi
if ! grep -q '^FAIL  ' "$out"; then
  printf 'FAIL  mutant run failed for an unrelated reason (no failed assertion)\n'
  tail -n 20 "$out"
  exit 1
fi
printf 'PASS  exactly one fail-closed line was mutated, marker verified, and the negative test rejected the mutant\n'
