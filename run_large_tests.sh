#!/usr/bin/env bash

err_exit() {
 echo -e "Error: $1" >&2
 exit 1
}

prog="./gclib-test"
if [[ ! -x $prog ]]; then
 make -j4 gclib-test || err_exit "Build failed."
fi

mkdir -p tests/wrk || err_exit "cannot create tests/wrk"
rm -f tests/wrk/large_genome.fa tests/wrk/large_genome.fa.fai tests/wrk/large_annotation.gtf tests/wrk/large_ranges.tsv tests/wrk/large_*.tsv

echo "=== Generating large-offset fixtures"
"$prog" generate-large-fixtures tests/wrk/large_genome.fa tests/wrk/large_annotation.gtf tests/wrk/large_ranges.tsv \
  || err_exit "could not generate large-offset fixtures"

run_test() {
 local label="$1"
 local outfile="$2"
 local expfile="$3"
 shift 3
 echo "=== Running $label"
 "$prog" "$@" > "$outfile" || err_exit "command failed for $label"
 if diff -u "$expfile" "$outfile" ; then
   echo "  OK."
 else
   err_exit "output mismatch for $label"
 fi
}

run_test "large fasta-index" tests/wrk/large_fasta_index.tsv tests/expected/large_fasta_index.tsv \
  fasta-index tests/wrk/large_genome.fa
run_test "large fasta-ranges" tests/wrk/large_fasta_ranges.tsv tests/expected/large_fasta_ranges.tsv \
  fasta-ranges tests/wrk/large_genome.fa tests/wrk/large_ranges.tsv
run_test "large transcript-seq" tests/wrk/large_annotation.spliced.tsv tests/expected/large_annotation.spliced.tsv \
  transcript-seq tests/wrk/large_annotation.gtf tests/wrk/large_genome.fa

echo "All large-offset tests passed successfully!"
