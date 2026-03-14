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
rm -f tests/wrk/genome.fa tests/wrk/genome.fa.fai tests/wrk/*.tsv
cp tests/data/genome.fa tests/wrk/genome.fa || err_exit "cannot copy genome.fa"

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

run_test "fasta-index" tests/wrk/fasta_index.tsv tests/expected/fasta_index.tsv \
  fasta-index tests/wrk/genome.fa
run_test "fasta-ranges" tests/wrk/fasta_ranges.tsv tests/expected/fasta_ranges.tsv \
  fasta-ranges tests/wrk/genome.fa tests/data/fasta_ranges.tsv
run_test "annotation.gff summary" tests/wrk/annotation_gff.summary.tsv \
  tests/expected/annotation_gff.summary.tsv gff-summary tests/data/annotation.gff
run_test "annotation.gff.gz summary" tests/wrk/annotation_gff_gz.summary.tsv \
  tests/expected/annotation_gff.summary.tsv gff-summary tests/data/annotation.gff.gz
run_test "annotation.gtf summary" tests/wrk/annotation_gtf.summary.tsv \
  tests/expected/annotation_gtf.summary.tsv gff-summary tests/data/annotation.gtf
run_test "annotation.gtf.gz summary" tests/wrk/annotation_gtf_gz.summary.tsv \
  tests/expected/annotation_gtf.summary.tsv gff-summary tests/data/annotation.gtf.gz
run_test "transcripts.gtf summary" tests/wrk/transcripts_gtf.summary.tsv \
  tests/expected/transcripts_gtf.summary.tsv gff-summary tests/data/transcripts.gtf
run_test "transcripts.gtf.gz summary" tests/wrk/transcripts_gtf_gz.summary.tsv \
  tests/expected/transcripts_gtf.summary.tsv gff-summary tests/data/transcripts.gtf.gz
run_test "annotation.gff transcript-seq" tests/wrk/annotation_gff.spliced.tsv \
  tests/expected/annotation_gff.spliced.tsv transcript-seq tests/data/annotation.gff tests/wrk/genome.fa
run_test "annotation.gff.gz transcript-seq" tests/wrk/annotation_gff_gz.spliced.tsv \
  tests/expected/annotation_gff.spliced.tsv transcript-seq tests/data/annotation.gff.gz tests/wrk/genome.fa
run_test "annotation.gtf transcript-seq" tests/wrk/annotation_gtf.spliced.tsv \
  tests/expected/annotation_gtf.spliced.tsv transcript-seq tests/data/annotation.gtf tests/wrk/genome.fa
run_test "annotation.gtf.gz transcript-seq" tests/wrk/annotation_gtf_gz.spliced.tsv \
  tests/expected/annotation_gtf.spliced.tsv transcript-seq tests/data/annotation.gtf.gz tests/wrk/genome.fa

echo "All tests passed successfully!"
