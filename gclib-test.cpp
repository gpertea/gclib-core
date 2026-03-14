#include "GBase.h"
#include "GStr.h"
#include "GVec.hh"
#include "GList.hh"
#include "GFastaIndex.h"
#include "GFaSeqGet.h"
#include "gff.h"
#include "gdna.h"

#define USAGE "Usage:\n\
  gclib-test generate-large-fixtures <fasta> <annotation.gtf> <ranges.tsv>\n\
  gclib-test fasta-index <fasta>\n\
  gclib-test fasta-ranges <fasta> <ranges.tsv>\n\
  gclib-test gff-summary <annotation.gff|annotation.gtf>\n\
  gclib-test transcript-seq <annotation.gff|annotation.gtf> <fasta>\n\
"

static const char* LARGE_SEQNAME = "chrHuge";
static const int64_t LARGE_LINE_LEN = 1 << 20;
static const int64_t LARGE_FULL_LINES = 4096;
static const int64_t LARGE_BASE_COORD = LARGE_LINE_LEN * LARGE_FULL_LINES;
static const char* LARGE_TAIL =
    "ACGTACGTACGTACGTACGTACGTACGTACGT"
    "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
    "AGCTTAGCAGCTTAGCAGCTTAGCAGCTTAGC"
    "CCGGTTAACCGGTTAACCGGTTAACCGGTTAA"
    "TGCATGCATGCATGCATGCATGCATGCATGCA"
    "AACCAACCAACCAACCAACCAACCAACCAACC";

static void writeLargeFixtureFasta(const char* fasta) {
 FILE* fw=fopen(fasta, "wb");
 if (fw==NULL) GError("Error: could not create large FASTA %s\n", fasta);
 char* wbuf=NULL;
 GMALLOC(wbuf, 8*1024*1024);
 setvbuf(fw, wbuf, _IOFBF, 8*1024*1024);
 fprintf(fw, ">%s\n", LARGE_SEQNAME);

 char* linebuf=NULL;
 GMALLOC(linebuf, LARGE_LINE_LEN+1);
 memset(linebuf, 'N', LARGE_LINE_LEN);
 linebuf[LARGE_LINE_LEN]='\n';
 for (int64_t i=0;i<LARGE_FULL_LINES;++i) {
   if (fwrite(linebuf, 1, LARGE_LINE_LEN+1, fw)!=(size_t)(LARGE_LINE_LEN+1))
     GError("Error writing large FASTA body to %s\n", fasta);
 }
 const size_t tail_len=strlen(LARGE_TAIL);
 if (fwrite(LARGE_TAIL, 1, tail_len, fw)!=tail_len)
   GError("Error writing large FASTA tail to %s\n", fasta);
 fputc('\n', fw);
 fclose(fw);
 GFREE(linebuf);
 GFREE(wbuf);
}

static void writeLargeFixtureGtf(const char* gtf) {
 FILE* fw=fopen(gtf, "wb");
 if (fw==NULL) GError("Error: could not create large GTF %s\n", gtf);
 const int64_t plus_exon1_start=LARGE_BASE_COORD+1;
 const int64_t plus_exon1_end=LARGE_BASE_COORD+16;
 const int64_t plus_exon2_start=LARGE_BASE_COORD+101;
 const int64_t plus_exon2_end=LARGE_BASE_COORD+120;
 const int64_t minus_exon1_start=LARGE_BASE_COORD+33;
 const int64_t minus_exon1_end=LARGE_BASE_COORD+48;
 const int64_t minus_exon2_start=LARGE_BASE_COORD+145;
 const int64_t minus_exon2_end=LARGE_BASE_COORD+160;

 fprintf(fw, "%s\ttest\ttranscript\t%" PRId64 "\t%" PRId64 "\t.\t+\t.\tgene_id \"gene_plus\"; transcript_id \"tx_plus\";\n",
         LARGE_SEQNAME, plus_exon1_start, plus_exon2_end);
 fprintf(fw, "%s\ttest\texon\t%" PRId64 "\t%" PRId64 "\t.\t+\t.\tgene_id \"gene_plus\"; transcript_id \"tx_plus\";\n",
         LARGE_SEQNAME, plus_exon1_start, plus_exon1_end);
 fprintf(fw, "%s\ttest\texon\t%" PRId64 "\t%" PRId64 "\t.\t+\t.\tgene_id \"gene_plus\"; transcript_id \"tx_plus\";\n",
         LARGE_SEQNAME, plus_exon2_start, plus_exon2_end);
 fprintf(fw, "%s\ttest\tCDS\t%" PRId64 "\t%" PRId64 "\t.\t+\t0\tgene_id \"gene_plus\"; transcript_id \"tx_plus\";\n",
         LARGE_SEQNAME, plus_exon1_start, plus_exon1_end);
 fprintf(fw, "%s\ttest\tCDS\t%" PRId64 "\t%" PRId64 "\t.\t+\t0\tgene_id \"gene_plus\"; transcript_id \"tx_plus\";\n",
         LARGE_SEQNAME, plus_exon2_start, plus_exon2_end);

 fprintf(fw, "%s\ttest\ttranscript\t%" PRId64 "\t%" PRId64 "\t.\t-\t.\tgene_id \"gene_minus\"; transcript_id \"tx_minus\";\n",
         LARGE_SEQNAME, minus_exon1_start, minus_exon2_end);
 fprintf(fw, "%s\ttest\texon\t%" PRId64 "\t%" PRId64 "\t.\t-\t.\tgene_id \"gene_minus\"; transcript_id \"tx_minus\";\n",
         LARGE_SEQNAME, minus_exon1_start, minus_exon1_end);
 fprintf(fw, "%s\ttest\texon\t%" PRId64 "\t%" PRId64 "\t.\t-\t.\tgene_id \"gene_minus\"; transcript_id \"tx_minus\";\n",
         LARGE_SEQNAME, minus_exon2_start, minus_exon2_end);
 fprintf(fw, "%s\ttest\tCDS\t%" PRId64 "\t%" PRId64 "\t.\t-\t0\tgene_id \"gene_minus\"; transcript_id \"tx_minus\";\n",
         LARGE_SEQNAME, minus_exon1_start, minus_exon1_end);
 fprintf(fw, "%s\ttest\tCDS\t%" PRId64 "\t%" PRId64 "\t.\t-\t0\tgene_id \"gene_minus\"; transcript_id \"tx_minus\";\n",
         LARGE_SEQNAME, minus_exon2_start, minus_exon2_end);
 fclose(fw);
}

static void writeLargeFixtureRanges(const char* rangesf) {
 FILE* fw=fopen(rangesf, "wb");
 if (fw==NULL) GError("Error: could not create large ranges file %s\n", rangesf);
 fprintf(fw, "huge_e1\t%s\t%" PRId64 "\t%" PRId64 "\t+\n",
         LARGE_SEQNAME, LARGE_BASE_COORD+1, LARGE_BASE_COORD+16);
 fprintf(fw, "huge_e2\t%s\t%" PRId64 "\t%" PRId64 "\t+\n",
         LARGE_SEQNAME, LARGE_BASE_COORD+101, LARGE_BASE_COORD+120);
 fprintf(fw, "huge_rev\t%s\t%" PRId64 "\t%" PRId64 "\t-\n",
         LARGE_SEQNAME, LARGE_BASE_COORD+145, LARGE_BASE_COORD+160);
 fclose(fw);
}

static int cmdGenerateLargeFixtures(const char* fasta, const char* gtf, const char* rangesf) {
 writeLargeFixtureFasta(fasta);
 writeLargeFixtureGtf(gtf);
 writeLargeFixtureRanges(rangesf);
 return 0;
}

static const char* nullStr(const char* s) {
 return (s==NULL || s[0]==0) ? "." : s;
}

static const char* getGeneId(GffObj* gffrec) {
 char* s=gffrec->getGeneID();
 if (s!=NULL && s[0]!=0) return s;
 s=gffrec->getAttr("gene_id");
 if (s!=NULL && s[0]!=0) return s;
 s=gffrec->getAttr("geneID");
 if (s!=NULL && s[0]!=0) return s;
 return ".";
}

static const char* getGeneName(GffObj* gffrec) {
 char* s=gffrec->getGeneName();
 if (s!=NULL && s[0]!=0) return s;
 s=gffrec->getAttr("gene_name");
 if (s!=NULL && s[0]!=0) return s;
 s=gffrec->getAttr("gene");
 if (s!=NULL && s[0]!=0) return s;
 s=gffrec->getAttr("Name");
 if (s!=NULL && s[0]!=0) return s;
 return ".";
}

static void addSeg(GStr& s, int64_t start, int64_t end) {
 if (!s.is_empty()) s+=',';
 s+=start;
 s+='-';
 s+=end;
}

static const char* getExonRanges(GffObj* gffrec, GStr& s) {
 s.clear();
 for (int i=0;i<gffrec->exons.Count();++i) {
   addSeg(s, gffrec->exons[i]->start, gffrec->exons[i]->end);
 }
 return s.is_empty() ? "." : s.chars();
}

static const char* getCDSRanges(GffObj* gffrec, GVec<GffExon>& cds, GStr& s) {
 s.clear();
 cds.Clear();
 if (!gffrec->hasCDS()) return ".";
 gffrec->getCDSegs(cds);
 for (int i=0;i<cds.Count();++i) {
   addSeg(s, cds[i].start, cds[i].end);
 }
 return s.is_empty() ? "." : s.chars();
}

static void printFastaIndex(GFastaIndex& faidx) {
 GList<GFastaRec> reclist(true,false,true);
 faidx.records.startIterate();
 GFastaRec* rec=NULL;
 while ((rec=faidx.records.NextData())!=NULL) {
   reclist.Add(rec);
 }
 for (int i=0;i<reclist.Count();++i) {
#ifdef _WIN32
   fprintf(stdout, "%s\t%" PRId64 "\t%ld\t%d\t%d\n", reclist[i]->seqname,
           reclist[i]->seqlen, (long)reclist[i]->fpos,
           reclist[i]->line_len, reclist[i]->line_blen);
#else
   fprintf(stdout, "%s\t%" PRId64 "\t%lld\t%d\t%d\n", reclist[i]->seqname,
           reclist[i]->seqlen, (long long)reclist[i]->fpos,
           reclist[i]->line_len, reclist[i]->line_blen);
#endif
 }
}

static int cmdFastaIndex(const char* fasta) {
 GStr fainame(fasta);
 fainame+=".fai";
 GFastaIndex faidx(fasta, fainame.chars());
 if (!faidx.hasIndex()) {
   if (faidx.buildIndex()==0)
     GError("Error: no FASTA records found in %s\n", fasta);
   if (faidx.storeIndex(fainame.chars())==0)
     GError("Error: could not write FASTA index %s\n", fainame.chars());
 }
 GFastaIndex loaded(fasta, fainame.chars());
 if (!loaded.hasIndex())
   GError("Error: could not reload FASTA index %s\n", fainame.chars());
 printFastaIndex(loaded);
 return 0;
}

static int cmdFastaRanges(const char* fasta, const char* rangesf) {
 GFastaDb gfasta(fasta, true);
 GLineReader lr(rangesf);
 char* line=NULL;
 while ((line=lr.nextLine())!=NULL) {
   if (line[0]==0 || line[0]=='#') continue;
   GDynArray<char*> fields;
   strsplit(line, fields, '\t');
   if (fields.Count()<5)
     GError("Error parsing ranges line %d in %s\n", lr.readcount(), rangesf);
   int64_t start=0;
   int64_t end=0;
   char* p=fields[2];
   if (!parseUInt(p, start) || *p!=0)
     GError("Error parsing start coordinate at line %d in %s\n", lr.readcount(), rangesf);
   p=fields[3];
   if (!parseUInt(p, end) || *p!=0)
     GError("Error parsing end coordinate at line %d in %s\n", lr.readcount(), rangesf);
   bool rev=(fields[4][0]=='-');
   GFaSeqGet* faseq=gfasta.fetch(fields[1]);
   if (faseq==NULL)
     GError("Error: could not fetch sequence %s from %s\n", fields[1], fasta);
   char* seq=faseq->copyRange(start, end, rev, false);
   if (seq==NULL)
     GError("Error: could not fetch range %s %" PRId64 "-%" PRId64 "\n", fields[1], start, end);
   fprintf(stdout, "%s\t%s\t%" PRId64 "\t%" PRId64 "\t%c\t%d\t%s\n",
           fields[0], fields[1], start, end,
           rev ? '-' : '+', (int)strlen(seq), seq);
   GFREE(seq);
 }
 return 0;
}

static int cmdGffSummary(const char* fname) {
 GffReader reader(fname, true, true);
 reader.readAll();
 GVec<GffExon> cds;
 GStr exons;
 GStr cdslist;
 for (int i=0;i<reader.gflst.Count();++i) {
   GffObj* gffrec=reader.gflst[i];
   if (gffrec==NULL || !gffrec->isValidTranscript()) continue;
   getExonRanges(gffrec, exons);
   getCDSRanges(gffrec, cds, cdslist);
   fprintf(stdout, "%s\t%s\t%" PRId64 "\t%" PRId64 "\t%c\t%s\t%" PRId64 "\t%s\t%" PRId64 "\t%s\t%s\t%s\n",
           nullStr(gffrec->getID()),
           nullStr(gffrec->getGSeqName()),
           gffrec->start, gffrec->end, gffrec->strand,
           nullStr(gffrec->getFeatureName()),
           gffrec->exons.Count(),
           exons.is_empty() ? "." : exons.chars(),
           cds.Count(),
           cdslist.is_empty() ? "." : cdslist.chars(),
           getGeneId(gffrec), getGeneName(gffrec));
 }
 return 0;
}

static int cmdTranscriptSeq(const char* fname, const char* fasta) {
 GffReader reader(fname, true, true);
 GFastaDb gfasta(fasta, true);
 reader.readAll();
 for (int i=0;i<reader.gflst.Count();++i) {
   GffObj* gffrec=reader.gflst[i];
   if (gffrec==NULL || !gffrec->isValidTranscript()) continue;
   GFaSeqGet* faseq=gfasta.fetch(gffrec->getGSeqName());
   if (faseq==NULL)
     GError("Error: could not fetch sequence %s from %s\n", gffrec->getGSeqName(), fasta);
   int splen=0;
   char* spliced=gffrec->getSpliced(faseq, false, &splen);
   if (spliced==NULL)
     GError("Error: could not get spliced sequence for %s\n", gffrec->getID());
   char* cdsseq=NULL;
   int cdslen=0;
   if (gffrec->hasCDS()) {
     cdsseq=gffrec->getSpliced(faseq, true, &cdslen);
     if (cdsseq==NULL)
       GError("Error: could not get CDS sequence for %s\n", gffrec->getID());
     fprintf(stdout, "%s\t%s\t%c\t%" PRId64 "\t%d\t%s\t%d\t%s\n",
             nullStr(gffrec->getID()), nullStr(gffrec->getGSeqName()),
             gffrec->strand, gffrec->exons.Count(),
             splen, spliced, cdslen, cdsseq);
   }
   else {
     fprintf(stdout, "%s\t%s\t%c\t%" PRId64 "\t%d\t%s\t.\t.\n",
             nullStr(gffrec->getID()), nullStr(gffrec->getGSeqName()),
             gffrec->strand, gffrec->exons.Count(),
             splen, spliced);
   }
   GFREE(spliced);
   GFREE(cdsseq);
 }
 return 0;
}

int main(int argc, char* argv[]) {
 gDnaInit();
 if (argc<3) {
   fprintf(stderr, "%s", USAGE);
   return 1;
 }
 if (strcmp(argv[1], "generate-large-fixtures")==0) {
   if (argc!=5) {
     fprintf(stderr, "%s", USAGE);
     return 1;
   }
   return cmdGenerateLargeFixtures(argv[2], argv[3], argv[4]);
 }
 if (strcmp(argv[1], "fasta-index")==0) {
   if (argc!=3) {
     fprintf(stderr, "%s", USAGE);
     return 1;
   }
   return cmdFastaIndex(argv[2]);
 }
 if (strcmp(argv[1], "fasta-ranges")==0) {
   if (argc!=4) {
     fprintf(stderr, "%s", USAGE);
     return 1;
   }
   return cmdFastaRanges(argv[2], argv[3]);
 }
 if (strcmp(argv[1], "gff-summary")==0) {
   if (argc!=3) {
     fprintf(stderr, "%s", USAGE);
     return 1;
   }
   return cmdGffSummary(argv[2]);
 }
 if (strcmp(argv[1], "transcript-seq")==0) {
   if (argc!=4) {
     fprintf(stderr, "%s", USAGE);
     return 1;
   }
   return cmdTranscriptSeq(argv[2], argv[3]);
 }
 fprintf(stderr, "%s", USAGE);
 return 1;
}
