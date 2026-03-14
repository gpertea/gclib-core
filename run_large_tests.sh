#!/usr/bin/env bash

err_exit() {
 echo -e "Error: $1" >&2
 exit 1
}

prog="./gclib-test"
shift_bases=5368709120
srcdir="tests/data"
expdir="tests/expected"
wrkdir="tests/wrk"
if [[ ! -x $prog ]]; then
 make -j4 gclib-test || err_exit "Build failed."
fi

mkdir -p "$wrkdir" || err_exit "cannot create $wrkdir"
rm -f "$wrkdir"/genome_large.fa "$wrkdir"/genome_large.fa.fai \
  "$wrkdir"/annotation_large.gff "$wrkdir"/annotation_large.gff.gz \
  "$wrkdir"/annotation_large.gtf "$wrkdir"/annotation_large.gtf.gz \
  "$wrkdir"/transcripts_large.gtf "$wrkdir"/transcripts_large.gtf.gz \
  "$wrkdir"/fasta_ranges_large.tsv "$wrkdir"/large_*.tsv \
  "$wrkdir"/large_genome.fa "$wrkdir"/large_genome.fa.fai \
  "$wrkdir"/large_annotation.gtf "$wrkdir"/large_ranges.tsv

shift_feature_file() {
 local infile="$1"
 local outfile="$2"
 SHIFT_BASES="$shift_bases" perl -F'\t' -lane '
  BEGIN { $shift=$ENV{SHIFT_BASES}+0; }
  if ($_ ne "" && $F[0] !~ /^#/ && @F >= 5) {
    $F[3]+=$shift;
    $F[4]+=$shift;
    $_=join("\t", @F);
  }
  print $_;
 ' "$infile" > "$outfile" || err_exit "could not create $outfile"
}

shift_ranges_file() {
 local infile="$1"
 local outfile="$2"
 SHIFT_BASES="$shift_bases" perl -F'\t' -lane '
  BEGIN { $shift=$ENV{SHIFT_BASES}+0; }
  if ($_ ne "" && $F[0] !~ /^#/ && @F >= 5) {
    $F[2]+=$shift;
    $F[3]+=$shift;
    $_=join("\t", @F);
  }
  print $_;
 ' "$infile" > "$outfile" || err_exit "could not create $outfile"
}

shift_index_expected() {
 local infile="$1"
 local outfile="$2"
 SHIFT_BASES="$shift_bases" perl -F'\t' -lane '
  BEGIN { $shift=$ENV{SHIFT_BASES}+0; }
  $F[1]+=$shift;
  print join("\t", @F);
 ' "$infile" > "$outfile" || err_exit "could not create $outfile"
}

shift_ranges_expected() {
 local infile="$1"
 local outfile="$2"
 SHIFT_BASES="$shift_bases" perl -F'\t' -lane '
  BEGIN { $shift=$ENV{SHIFT_BASES}+0; }
  $F[2]+=$shift;
  $F[3]+=$shift;
  print join("\t", @F);
 ' "$infile" > "$outfile" || err_exit "could not create $outfile"
}

shift_summary_expected() {
 local infile="$1"
 local outfile="$2"
 SHIFT_BASES="$shift_bases" perl -F'\t' -lane '
  BEGIN {
    $shift=$ENV{SHIFT_BASES}+0;
    sub shift_ranges {
      my ($s, $shift)=@_;
      return $s if !defined($s) || $s eq "." || $s eq "";
      my @ranges=split(/,/, $s);
      for my $r (@ranges) {
        my ($start, $end)=split(/-/, $r, 2);
        $r=($start+$shift)."-".($end+$shift);
      }
      return join(",", @ranges);
    }
  }
  $F[2]+=$shift;
  $F[3]+=$shift;
  $F[7]=shift_ranges($F[7], $shift);
  $F[9]=shift_ranges($F[9], $shift);
  print join("\t", @F);
 ' "$infile" > "$outfile" || err_exit "could not create $outfile"
}

generate_large_fasta() {
 local infile="$1"
 local outfile="$2"
 local seqname=""
 local line_len=""
 local full_lines=0
 local remainder=0
 local nline=""
 seqname=$(sed -n '1{s/^>//;p;q;}' "$infile")
 line_len=$(awk 'NR==2 { print length($0); exit }' "$infile")
 [[ -n $seqname ]] || err_exit "could not read FASTA header from $infile"
 [[ -n $line_len && $line_len -gt 0 ]] || err_exit "could not determine FASTA line length from $infile"
 full_lines=$((shift_bases / line_len))
 remainder=$((shift_bases % line_len))
 nline=$(printf '%*s' "$line_len" '' | tr ' ' 'N')
  {
   printf '>%s\n' "$seqname"
   {
     yes "$nline" | head -n "$full_lines" | tr -d '\n'
     if (( remainder > 0 )); then
       printf '%*s' "$remainder" '' | tr ' ' 'N'
     fi
     tail -n +2 "$infile" | tr -d '\n'
   } | fold -w "$line_len"
 } > "$outfile" || err_exit "could not create $outfile"
}

echo "=== Building large-offset fixtures from committed test data"
generate_large_fasta "$srcdir/genome.fa" "$wrkdir/genome_large.fa"
shift_feature_file "$srcdir/annotation.gff" "$wrkdir/annotation_large.gff"
shift_feature_file "$srcdir/annotation.gtf" "$wrkdir/annotation_large.gtf"
shift_feature_file "$srcdir/transcripts.gtf" "$wrkdir/transcripts_large.gtf"
shift_ranges_file "$srcdir/fasta_ranges.tsv" "$wrkdir/fasta_ranges_large.tsv"
gzip -c "$wrkdir/annotation_large.gff" > "$wrkdir/annotation_large.gff.gz" \
  || err_exit "could not create $wrkdir/annotation_large.gff.gz"
gzip -c "$wrkdir/annotation_large.gtf" > "$wrkdir/annotation_large.gtf.gz" \
  || err_exit "could not create $wrkdir/annotation_large.gtf.gz"
gzip -c "$wrkdir/transcripts_large.gtf" > "$wrkdir/transcripts_large.gtf.gz" \
  || err_exit "could not create $wrkdir/transcripts_large.gtf.gz"

shift_index_expected "$expdir/fasta_index.tsv" "$wrkdir/large_fasta_index.expected.tsv"
shift_ranges_expected "$expdir/fasta_ranges.tsv" "$wrkdir/large_fasta_ranges.expected.tsv"
shift_summary_expected "$expdir/annotation_gff.summary.tsv" "$wrkdir/annotation_large_gff.summary.expected.tsv"
shift_summary_expected "$expdir/annotation_gtf.summary.tsv" "$wrkdir/annotation_large_gtf.summary.expected.tsv"
shift_summary_expected "$expdir/transcripts_gtf.summary.tsv" "$wrkdir/transcripts_large_gtf.summary.expected.tsv"

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

run_test "large fasta-index" "$wrkdir/large_fasta_index.tsv" "$wrkdir/large_fasta_index.expected.tsv" \
  fasta-index "$wrkdir/genome_large.fa"
run_test "large fasta-ranges" "$wrkdir/large_fasta_ranges.tsv" "$wrkdir/large_fasta_ranges.expected.tsv" \
  fasta-ranges "$wrkdir/genome_large.fa" "$wrkdir/fasta_ranges_large.tsv"
run_test "large annotation.gff summary" "$wrkdir/annotation_large_gff.summary.tsv" \
  "$wrkdir/annotation_large_gff.summary.expected.tsv" gff-summary "$wrkdir/annotation_large.gff"
run_test "large annotation.gff.gz summary" "$wrkdir/annotation_large_gff_gz.summary.tsv" \
  "$wrkdir/annotation_large_gff.summary.expected.tsv" gff-summary "$wrkdir/annotation_large.gff.gz"
run_test "large annotation.gtf summary" "$wrkdir/annotation_large_gtf.summary.tsv" \
  "$wrkdir/annotation_large_gtf.summary.expected.tsv" gff-summary "$wrkdir/annotation_large.gtf"
run_test "large annotation.gtf.gz summary" "$wrkdir/annotation_large_gtf_gz.summary.tsv" \
  "$wrkdir/annotation_large_gtf.summary.expected.tsv" gff-summary "$wrkdir/annotation_large.gtf.gz"
run_test "large transcripts.gtf summary" "$wrkdir/transcripts_large_gtf.summary.tsv" \
  "$wrkdir/transcripts_large_gtf.summary.expected.tsv" gff-summary "$wrkdir/transcripts_large.gtf"
run_test "large transcripts.gtf.gz summary" "$wrkdir/transcripts_large_gtf_gz.summary.tsv" \
  "$wrkdir/transcripts_large_gtf.summary.expected.tsv" gff-summary "$wrkdir/transcripts_large.gtf.gz"
run_test "large annotation.gff transcript-seq" "$wrkdir/annotation_large_gff.spliced.tsv" \
  "$expdir/annotation_gff.spliced.tsv" transcript-seq "$wrkdir/annotation_large.gff" "$wrkdir/genome_large.fa"
run_test "large annotation.gff.gz transcript-seq" "$wrkdir/annotation_large_gff_gz.spliced.tsv" \
  "$expdir/annotation_gff.spliced.tsv" transcript-seq "$wrkdir/annotation_large.gff.gz" "$wrkdir/genome_large.fa"
run_test "large annotation.gtf transcript-seq" "$wrkdir/annotation_large_gtf.spliced.tsv" \
  "$expdir/annotation_gtf.spliced.tsv" transcript-seq "$wrkdir/annotation_large.gtf" "$wrkdir/genome_large.fa"
run_test "large annotation.gtf.gz transcript-seq" "$wrkdir/annotation_large_gtf_gz.spliced.tsv" \
  "$expdir/annotation_gtf.spliced.tsv" transcript-seq "$wrkdir/annotation_large.gtf.gz" "$wrkdir/genome_large.fa"

echo "All large-offset tests passed successfully!"
