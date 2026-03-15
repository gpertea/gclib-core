#ifndef GFASEQGET_H
#define GFASEQGET_H
#include "GFastaIndex.h"

#define MAX_FASUBSEQ 0x20000000
//max 512MB sequence data held in memory at a time

class GSubSeq {
 public:
  int64_t sqstart; //1-based coord of subseq start on sequence
  int64_t sqlen;   //length of subseq loaded
  char* sq; //actual subsequence data will be stored here
                // (with end-of-line characters removed)

  /*char* xseq; //the exposed pointer to the last requested subsequence start
  off_t xstart; //the coordinate start for the last requested subseq
  off_t xlen; //the last requested subseq len*/
  GSubSeq() {
     sqstart=0;
     sqlen=0;
     sq=NULL;
     /* xseq=NULL;
     xstart=0;
     xlen=0;*/
     }
  void forget() { //forget about pointer data, so we can reuse it
  	sq=NULL;
  	sqstart=0;
  	sqlen=0;
  }
  ~GSubSeq() {
     GFREE(sq);
     }
  // genomic, 1-based coordinates:
  void setup(int64_t sstart, int64_t slen, int64_t sovl=0, int64_t qfrom=0, int64_t qto=0, int64_t maxseqlen=0);
    //check for overlap with previous window and realloc/extend appropriately
    //returns offset from seq that corresponds to sstart
    // the window will keep extending until MAX_FASUBSEQ is reached
};

//
class GFaSeqGet {
  char* fname; //file name where the sequence resides
  FILE* fh;
  off_t fseqstart; //file offset where the sequence actually starts
  int64_t seq_len; //total sequence length, if known (when created from GFastaIndex)
  int64_t line_len; //length of each line of text
  int64_t line_blen; //binary length of each line
                 // = line_len + number of EOL character(s)
  GSubSeq* lastsub;
  void initialParse(off_t fofs=0, bool checkall=true);
  const char* loadsubseq(int64_t cstart, int64_t& clen);
  void finit(const char* fn, off_t fofs, bool validate);
 public:
  //GStr seqname; //current sequence name
  char* seqname;
  GFaSeqGet(): fname(NULL), fh(NULL), fseqstart(0), seq_len(0),
		  line_len(0), line_blen(0), lastsub(NULL), seqname(NULL) {
  }

  GFaSeqGet(const char* fn, off_t fofs, bool validate=false):fname(NULL), fh(NULL),
		    fseqstart(0), seq_len(0), line_len(0), line_blen(0),
			lastsub(NULL), seqname(NULL) {
     finit(fn,fofs,validate);
  }

  GFaSeqGet(const char* fn, bool validate=false):fname(NULL), fh(NULL),
		    fseqstart(0), seq_len(0), line_len(0), line_blen(0),
			lastsub(NULL), seqname(NULL) {
     finit(fn,0,validate);
  }

  GFaSeqGet(const char* faname, int64_t seqlen, off_t fseqofs, int l_len, int l_blen);
  //constructor from GFastaIndex record

  GFaSeqGet(FILE* f, off_t fofs=0, bool validate=false);

  ~GFaSeqGet() {
    if (fname!=NULL) {
       GFREE(fname);
       fclose(fh);
    }
    GFREE(seqname);
    delete lastsub;
  }

  const char* seq(int64_t cstart=1, int64_t clen=0) {
		  int64_t cend = clen==0 ? 0 : cstart+clen-1;
		  return getRange(cstart, cend);
	  }

  const char* subseq(int64_t cstart, int64_t& clen);
  const char* getRange(int64_t cstart=1, int64_t cend=0) {
	      if (cend==0) cend=(seq_len>0)?seq_len : MAX_FASUBSEQ;
	      if (cstart>cend) { Gswap(cstart, cend); }
	      int64_t clen=cend-cstart+1;
	      //int rdlen=clen;
	      return subseq(cstart, clen);
	  }

  //caller is responsible for deallocating the return string
  char* copyRange(int64_t cstart, int64_t cend, bool revCmpl=false, bool upCase=false);

  //uncached, read and return allocated buffer
  //caller is responsible for deallocating the return string
  char* fetchSeq(int64_t* retlen=NULL) {
  	int64_t clen=(seq_len>0) ? seq_len : MAX_FASUBSEQ;
  	if (lastsub) { delete lastsub; lastsub=NULL; }
  	subseq(1, clen);
  	if (retlen) *retlen=clen;
  	char* r=lastsub->sq;
  	lastsub->forget();
  	if (clen>0) {
  	   r[clen]=0;
  	}
  	else {
  		r=NULL;
  	}
  	return r;
  }

  void loadall(int64_t max_len=0) {
    //TODO: better read the whole sequence differently here - line by line
    //so when EOF or another '>' line is found, the reading stops!
    int64_t clen=(seq_len>0) ? seq_len : ((max_len>0) ? max_len : MAX_FASUBSEQ);
    subseq(1, clen);
    }
  void load(int64_t cstart, int64_t cend) {
     //cache as much as possible
      if (seq_len>0 && cend>seq_len) cend=seq_len; //correct a bad request
      int64_t clen=cend-cstart+1;
      subseq(cstart, clen);
     }
  int64_t getsublen() { return lastsub!=NULL ? lastsub->sqlen : 0 ; }
  int64_t getseqlen() { return seq_len; } //known when loaded with GFastaIndex
  off_t getseqofs() { return fseqstart; }
  int64_t getLineLen() { return line_len; }
  int64_t getLineBLen() { return line_blen; }
  //reads a subsequence starting at genomic coordinate cstart (1-based)
 };

//multi-fasta sequence handling
class GFastaDb {
 public:
  char* fastaPath;
  GFastaIndex* faIdx; //could be a cdb .cidx file
  //int last_fetchid;
  const char* last_seqname;
  GFaSeqGet* faseq;
  //GCdbYank* gcdb;
  GFastaDb(const char* fpath=NULL, bool forceIndexFile=true):fastaPath(NULL), faIdx(NULL), last_seqname(NULL),
		  faseq(NULL) {
     //gcdb=NULL;
     init(fpath, forceIndexFile);
  }

  void init(const char* fpath, bool writeIndexFile=true) {
     if (fpath==NULL || fpath[0]==0) return;
     //last_fetchid=-1;
     last_seqname=NULL;
     if (!fileExists(fpath))
       GError("Error: file/directory %s does not exist!\n",fpath);
     fastaPath=Gstrdup(fpath);
     //GStr gseqpath(fpath);
     if (fileExists(fastaPath)>1) { //exists and it's not a directory
            char* fainame=Gstrdup(fastaPath,4);
            size_t fainamelen=strlen(fainame);
            //int fainame_len=strlen(fainame);
            if (trimSuffix(fastaPath, ".fai")) {
               //.fai index file given directly
               if (!fileExists(fastaPath))
                  GError("Error: cannot find fasta file for index %s !\n", fastaPath);
            }
            else { //append .fai as needed
            	 strcpy(fainame+fainamelen, ".fai");
            	 fainamelen+=4;
            }
            //GMessage("creating GFastaIndex with fastaPath=%s, fainame=%s\n", fastaPath, fainame.chars());
            faIdx=new GFastaIndex(fastaPath, fainame);
            char* fainamecwd=fainame; //will hold just the file name without the path
            char* plast=strrchr(fainamecwd, '/');
            if (plast!=NULL) {
              fainamecwd=plast+1; //point to the file name only
            }
            if (!faIdx->hasIndex()) { //could not load index file .fai
               //try current directory (Warning: might not be the correct index for that file!)
               if (plast==NULL) {
                 if (fileExists(fainamecwd)>1) {
                    faIdx->loadIndex(fainamecwd);
                  }
               }
            } //tried to load index
            if (!faIdx->hasIndex()) { //no index file to be loaded, build the index
                 //if (forceIndexFile)
                 //   GMessage("No fasta index found for %s. Rebuilding, please wait..\n",fastaPath);
                 faIdx->buildIndex(); //build index in memory only
                 if (faIdx->getCount()==0) GError("Error: no fasta records found!\n");
                 if (writeIndexFile) {
                     //GMessage("Fasta index rebuilt.\n");
                     FILE* fcreate=fopen(fainame, "w");
                     char* idxfname=fainame;
                     if (fcreate==NULL) {
                        GMessage("Warning: cannot create fasta index file %s! (permissions?)\n", fainame);
                        if (fainame!=fainamecwd) {
                        	//try cwd
                        	idxfname=fainamecwd;
                        	GMessage("   Attempting to create the index in the current directory..\n");
                        	if ((fcreate=fopen(fainamecwd, "w"))==NULL)
                        		GError("Error: cannot create fasta index file %s!\n", fainamecwd);
                        }
                     }
                     if (fcreate!=NULL) {
                    	 if (faIdx->storeIndex(fcreate)<faIdx->getCount())
                              GMessage("Warning: error writing the index file %s!\n", idxfname);
                    	 else GMessage("FASTA index file %s created.\n", idxfname);
                     }
                 } //file storage of index requested
            } //creating FASTA index
    GFREE(fainame);
    } //multi-fasta file
  }

  GFaSeqGet* fetchFirst(const char* fname, bool checkFasta=false) {
	 faseq=new GFaSeqGet(fname, checkFasta);
	 faseq->loadall();
	 //last_fetchid=gseq_id;
	 GFREE(last_seqname);
	 last_seqname=Gstrdup(faseq->seqname);
	 return faseq;
  }

 char* getFastaFile(const char* gseqname) {
	if (fastaPath==NULL) return NULL;
		size_t gnl=strlen(gseqname);
		if (gnl>(size_t)INT_MAX-8)
		  GError("Error: FASTA sequence name too long: %s\n", gseqname);
		char* s=Gstrdup(fastaPath, (int)(gnl+8));
		size_t slen=strlen(s);
	if (s[slen-1]!='/') {//CHPATHSEP ?
	  s[slen]='/';
	  slen++;
	  s[slen]='\0';
	}
	//s.append(gseqname);
    strcpy(s+slen, gseqname);
    slen+=gnl;
	if (!fileExists(s)) {
		//s.append(".fa")
		strcpy(s+slen, ".fa");
		slen+=3;
	}
	if (!fileExists(s)) { strcpy(s+slen, "sta"); slen+=3; }
	if (fileExists(s)) return Gstrdup(s);
	 else {
	   GMessage("Warning: cannot find genomic sequence file %s/%s{.fa,.fasta}\n",fastaPath, s);
	   return NULL;
	 }
	GFREE(s);
  }

 GFaSeqGet* fetch(const char* gseqname) {
    if (fastaPath==NULL) return NULL;
    if (last_seqname!=NULL && (strcmp(gseqname, last_seqname)==0)
    		&& faseq!=NULL) return faseq;
    delete faseq;
    faseq=NULL;
    //last_fetchid=-1;
    GFREE(last_seqname);
    last_seqname=NULL;
    //char* gseqname=GffObj::names->gseqs.getName(gseq_id);
    if (faIdx!=NULL) { //fastaPath was the multi-fasta file name and it must have an index
        GFastaRec* farec=faIdx->getRecord(gseqname);
        if (farec!=NULL) {
             faseq=new GFaSeqGet(fastaPath,farec->seqlen, farec->fpos,
                               farec->line_len, farec->line_blen);
             if (farec->seqlen>0 && farec->seqlen<=MAX_FASUBSEQ) {
               faseq->loadall(); //just cache the whole sequence when it fits in the cache window
             }
             //last_fetchid=gseq_id;

             last_seqname=Gstrdup(gseqname);
        }
        else {
          GMessage("Warning: couldn't find fasta record for '%s'!\n",gseqname);
          return NULL;
        }
    }
    else { //directory with FASTA files named as gseqname
        char* sfile=getFastaFile(gseqname);
        if (sfile!=NULL) {
      	   faseq=new GFaSeqGet(sfile);
           faseq->loadall();
           //last_fetchid=gseq_id;
           GFREE(sfile);
           }
    } //one fasta file per contig

    //else GMessage("Warning: fasta index not available, cannot retrieve sequence %s\n",
    //		gseqname);
    return faseq;
  }

   ~GFastaDb() {
     GFREE(fastaPath);
     GFREE(last_seqname);
     //delete gcdb;
     delete faIdx;
     delete faseq;
     }
};


GFaSeqGet* fastaSeqGet(GFastaDb& gfasta, const char* seqid);

#endif
