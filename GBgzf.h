#ifndef GBGZF_H
#define GBGZF_H
// Compact, self-contained BGZF (blocked gzip) reader/writer for random access
// to bgzip-compressed FASTA files. Modeled on the htslib BGZF format so files
// are interoperable with samtools' `bgzip`/`faidx` (same block framing and the
// same binary `.gzi` block index), but implemented over plain stdio + zlib with
// no htslib dependency. Single-threaded, 64-bit offsets.
//
// A BGZF file is a sequence of standard gzip members, each carrying an extra
// "BC" subfield giving the total block size, followed by a 28-byte empty EOF
// marker block. The `.gzi` companion maps uncompressed offsets to physical
// (compressed) offsets so a byte position in the decompressed stream can be
// located without scanning; the `.fai` for a bgzf FASTA holds those same
// uncompressed offsets, so GFastaIndex is reused unchanged.

#include "GBase.h"
#include "GVec.hh"
#include <zlib.h>

class GBgzf {
  struct GBgzBlockIdx { uint64_t caddr; uint64_t uaddr; }; //block start offsets

  FILE* fp;
  bool is_write;
  int wlevel;
  bool zs_init;
  z_stream zs;

  // read state
  unsigned char* cblock; //raw block bytes (header stripped payload), <=BGZF_MAX_BLOCK
  unsigned char* ublock; //current decompressed block, <=BGZF_MAX_BLOCK
  int ublock_len;        //valid decompressed bytes in ublock
  int ublock_off;        //read cursor within ublock
  uint64_t cur_caddr;    //compressed offset of current block start
  uint64_t cur_uaddr;    //uncompressed offset of current block start
  int pushback;          //single-byte pushback for getc/ungetc (-1 if none)

  // write state
  unsigned char* wbuf;   //uncompressed accumulation buffer (<=BGZF_BLOCK)
  int wbuf_len;
  uint64_t wtotal_u;     //total uncompressed bytes written so far
  uint64_t wtotal_c;     //total compressed bytes written so far

  GVec<GBgzBlockIdx> idx; //block index (idx[0] is always {0,0})
  bool have_idx;

  int fill_block();      //read+inflate one block at current fp pos;
                         // returns decompressed len (0 for empty block), -1 at EOF
  bool flush_block();    //deflate+write wbuf as one bgzf block, record index entry
 public:
  GBgzf();
  ~GBgzf();

  static bool isBgzf(const char* fname); //sniff for gzip magic + BGZF "BC" subfield

  // -- reading
  bool openRead(const char* fname);
  bool loadIndex(const char* gziname); //load the binary .gzi block index
  bool hasIndex() { return have_idx; }
  int64_t useek(uint64_t uoffset);     //seek by uncompressed offset (needs index)
  int64_t read(void* buf, int64_t len);//read decompressed bytes across blocks
  int getc();                          //one decompressed byte, or -1 at EOF
  void ungetc(int c);

  // -- writing
  bool openWrite(const char* fname, int level=-1);
  int64_t write(const void* buf, int64_t len);
  bool dumpIndex(const char* gziname); //write the .gzi (call before close())
  bool close();
  bool isWrite() { return is_write; }
};

#endif
