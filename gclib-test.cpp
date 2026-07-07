#include "GBase.h"
#include "GStr.h"
#include "GVec.hh"
#include "GList.hh"
#include "GFastaIndex.h"
#include "GFaSeqGet.h"
#include "gff.h"
#include "gdna.h"

#define USAGE "Usage:\n\
  gclib-test fasta-index <fasta>\n\
  gclib-test fasta-ranges <fasta> <ranges.tsv>\n\
  gclib-test bgzip <in.fa|-> <out.fa.gz>\n\
  gclib-test gff-summary <annotation.gff|annotation.gtf>\n\
  gclib-test transcript-seq <annotation.gff|annotation.gtf> <fasta>\n\
"

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
   fprintf(stdout, "%s\t%" PRId64 "\t%lld\t%d\t%d\n", reclist[i]->seqname,
           reclist[i]->seqlen, (long long)reclist[i]->fpos,
           reclist[i]->line_len, reclist[i]->line_blen);
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

// write one .fai record line (5 columns: name, seqlen, uncompressed offset,
// line bases, line bytes) — identical layout to a plain-FASTA .fai
static void writeFaiRec(FILE* fai, const char* name, int64_t seqlen,
                        uint64_t ofs, int llen, int lblen) {
 fprintf(fai, "%s\t%" PRId64 "\t%" PRIu64 "\t%d\t%d\n",
         name, seqlen, ofs, llen, lblen);
}

// Compress a FASTA (file or stdin) into a BGZF file <out>, emitting the sibling
// <out>.gzi (block index) and <out>.fai (record index with UNCOMPRESSED offsets)
// in a single streaming pass — no external bgzip/samtools required, and the huge
// input is never materialized uncompressed on disk.
static int cmdBgzip(const char* infa, const char* outgz) {
 FILE* fi = (strcmp(infa, "-")==0) ? stdin : fopen(infa, "rb");
 if (fi==NULL) GError("Error: cannot open input FASTA %s\n", infa);
 GBgzf bg;
 if (!bg.openWrite(outgz)) GError("Error: cannot create BGZF file %s\n", outgz);
 GStr faiName(outgz); faiName+=".fai";
 GStr gziName(outgz); gziName+=".gzi";
 FILE* fai = fopen(faiName.chars(), "w");
 if (fai==NULL) GError("Error: cannot create FASTA index %s\n", faiName.chars());

 // streaming FASTA scanner (tracks UNCOMPRESSED byte offsets)
 bool lineStart=true, curIsHeader=false, collectingName=false, haveRec=false;
 bool firstSeqLine=false;
 GStr curName;
 int64_t curSeqlen=0;
 uint64_t curSeqofs=0, uoff=0;
 int curLlen=0, curLblen=0, curLineChars=0, curLineBytes=0;

 const size_t CHUNK=1<<20;
 char* buf=NULL; GMALLOC(buf, CHUNK);
 size_t n=0;
 while ((n=fread(buf, 1, CHUNK, fi))>0) {
   bg.write(buf, (int64_t)n);
   for (size_t i=0;i<n;i++) {
     char c=buf[i];
     if (lineStart) {
       lineStart=false;
       if (c=='>') {
         if (haveRec)
           writeFaiRec(fai, curName.chars(), curSeqlen, curSeqofs, curLlen, curLblen);
         curIsHeader=true; collectingName=true; curName.clear();
         curSeqlen=0; firstSeqLine=true; curLlen=0; curLblen=0;
         haveRec=true;
         uoff++;
         continue; //consume the '>' (not part of the name)
       }
       curIsHeader=false; curLineChars=0; curLineBytes=0;
     }
     if (curIsHeader) {
       if (collectingName) {
         if (c=='\n'||c=='\r'||c==' '||c=='\t') collectingName=false;
         else curName+=c;
       }
       if (c=='\n') { curSeqofs=uoff+1; lineStart=true; curIsHeader=false; }
     } else {
       if (c=='\n') {
         curLineBytes++;
         if (firstSeqLine) { curLlen=curLineChars; curLblen=curLineBytes; firstSeqLine=false; }
         curSeqlen+=curLineChars;
         lineStart=true;
       } else if (c=='\r') {
         curLineBytes++;
       } else {
         curLineChars++; curLineBytes++;
       }
     }
     uoff++;
   }
 }
 if (haveRec) //flush the last record
   writeFaiRec(fai, curName.chars(), curSeqlen, curSeqofs, curLlen, curLblen);
 GFREE(buf);
 if (fi!=stdin) fclose(fi);
 fclose(fai);
 bg.dumpIndex(gziName.chars());
 bg.close();
 if (!haveRec) GError("Error: no FASTA records found in %s\n", infa);
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
   fprintf(stdout, "%s\t%s\t%" PRId64 "\t%" PRId64 "\t%c\t%zu\t%s\n",
           fields[0], fields[1], start, end,
           rev ? '-' : '+', strlen(seq), seq);
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
   int64_t splen=0;
   char* spliced=gffrec->getSpliced(faseq, false, &splen);
   if (spliced==NULL)
     GError("Error: could not get spliced sequence for %s\n", gffrec->getID());
   char* cdsseq=NULL;
   int64_t cdslen=0;
   if (gffrec->hasCDS()) {
     cdsseq=gffrec->getSpliced(faseq, true, &cdslen);
     if (cdsseq==NULL)
       GError("Error: could not get CDS sequence for %s\n", gffrec->getID());
     fprintf(stdout, "%s\t%s\t%c\t%" PRId64 "\t%" PRId64 "\t%s\t%" PRId64 "\t%s\n",
             nullStr(gffrec->getID()), nullStr(gffrec->getGSeqName()),
             gffrec->strand, gffrec->exons.Count(),
             splen, spliced, cdslen, cdsseq);
   }
   else {
     fprintf(stdout, "%s\t%s\t%c\t%" PRId64 "\t%" PRId64 "\t%s\t.\t.\n",
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
 if (strcmp(argv[1], "bgzip")==0) {
   if (argc!=4) {
     fprintf(stderr, "%s", USAGE);
     return 1;
   }
   return cmdBgzip(argv[2], argv[3]);
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
