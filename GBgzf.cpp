#include "GBgzf.h"

// BGZF constants (see the SAM/BGZF spec)
#define BGZF_BLOCK      0xff00  // max uncompressed bytes packed per block (65280)
#define BGZF_MAX_BLOCK  0x10000 // max total/decompressed block size (65536)
#define BGZF_HDR_LEN    18      // fixed gzip header + 6-byte BC extra field
#define BGZF_FTR_LEN    8       // CRC32 + ISIZE trailer

// little-endian pack/unpack helpers (BGZF is always little-endian)
static inline uint16_t le_u16(const unsigned char* p) {
  return (uint16_t)(p[0] | (p[1] << 8));
}
static inline uint32_t le_u32(const unsigned char* p) {
  return (uint32_t)p[0] | ((uint32_t)p[1] << 8) |
         ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24);
}
static inline void put_u16(unsigned char* p, uint16_t v) {
  p[0] = (unsigned char)(v & 0xff);
  p[1] = (unsigned char)((v >> 8) & 0xff);
}
static inline void put_u32(unsigned char* p, uint32_t v) {
  p[0] = (unsigned char)(v & 0xff);
  p[1] = (unsigned char)((v >> 8) & 0xff);
  p[2] = (unsigned char)((v >> 16) & 0xff);
  p[3] = (unsigned char)((v >> 24) & 0xff);
}

// fixed 18-byte BGZF block header template (BSIZE filled in at bytes 16-17)
static const unsigned char bgzf_hdr[BGZF_HDR_LEN] = {
  31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, 66, 67, 2, 0, 0, 0
};
// 28-byte empty EOF marker block
static const unsigned char bgzf_eof[28] = {
  31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, 66, 67, 2, 0, 27, 0,
  3, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

GBgzf::GBgzf(): fp(NULL), is_write(false), wlevel(Z_DEFAULT_COMPRESSION),
   zs_init(false), cblock(NULL), ublock(NULL), ublock_len(0), ublock_off(0),
   cur_caddr(0), cur_uaddr(0), pushback(-1), wbuf(NULL), wbuf_len(0),
   wtotal_u(0), wtotal_c(0), idx(), have_idx(false) {
}

GBgzf::~GBgzf() {
  close();
  GFREE(cblock);
  GFREE(ublock);
  GFREE(wbuf);
}

bool GBgzf::isBgzf(const char* fname) {
  FILE* f = fopen(fname, "rb");
  if (f == NULL) return false;
  unsigned char h[BGZF_HDR_LEN];
  size_t n = fread(h, 1, BGZF_HDR_LEN, f);
  fclose(f);
  if (n < BGZF_HDR_LEN) return false;
  // gzip magic + deflate + FEXTRA, plus the "BC" (66,67) extra subfield id
  return (h[0] == 31 && h[1] == 139 && h[2] == 8 && (h[3] & 4) != 0 &&
          h[12] == 66 && h[13] == 67);
}

bool GBgzf::openRead(const char* fname) {
  fp = fopen(fname, "rb");
  if (fp == NULL) return false;
  is_write = false;
  GMALLOC(cblock, BGZF_MAX_BLOCK);
  GMALLOC(ublock, BGZF_MAX_BLOCK);
  zs.zalloc = Z_NULL; zs.zfree = Z_NULL; zs.opaque = Z_NULL;
  if (inflateInit2(&zs, -15) != Z_OK)
    GError("Error (GBgzf): inflateInit2 failed\n");
  zs_init = true;
  cur_caddr = 0;
  cur_uaddr = 0;
  ublock_len = 0;
  ublock_off = 0;
  pushback = -1;
  return true;
}

// read+inflate a single physical block at the current file position.
// returns decompressed length (>=0), or -1 at physical end of file.
int GBgzf::fill_block() {
  cur_caddr = (uint64_t)ftello(fp);
  unsigned char hdr[BGZF_HDR_LEN];
  size_t n = fread(hdr, 1, BGZF_HDR_LEN, fp);
  if (n == 0) { ublock_len = 0; ublock_off = 0; return -1; }
  if (n < BGZF_HDR_LEN)
    GError("Error (GBgzf): truncated block header\n");
  if (hdr[0] != 31 || hdr[1] != 139 || hdr[2] != 8 || (hdr[3] & 4) == 0 ||
      hdr[12] != 66 || hdr[13] != 67)
    GError("Error (GBgzf): not a valid BGZF block\n");
  int xlen = le_u16(hdr + 10);
  int bsize = le_u16(hdr + 16) + 1; //total block size
  if (xlen != 6) {
    // uncommon: extra field carries subfields beyond BC; skip the remainder
    if (fseeko(fp, (off_t)(xlen - 6), SEEK_CUR) != 0)
      GError("Error (GBgzf): cannot skip extra field\n");
  }
  int cpayload = bsize - BGZF_HDR_LEN; //CDATA + CRC32 + ISIZE
  if (cpayload <= 0 || cpayload > BGZF_MAX_BLOCK)
    GError("Error (GBgzf): invalid block size %d\n", bsize);
  size_t got = fread(cblock, 1, (size_t)cpayload, fp);
  if ((int)got != cpayload)
    GError("Error (GBgzf): truncated block payload\n");
  int cdata_len = cpayload - BGZF_FTR_LEN;
  uint32_t isize = le_u32(cblock + cdata_len + 4);
  if (isize > (uint32_t)BGZF_MAX_BLOCK)
    GError("Error (GBgzf): block ISIZE too large\n");
  if (isize == 0) { ublock_len = 0; ublock_off = 0; return 0; } //empty block
  inflateReset(&zs);
  zs.next_in = cblock;
  zs.avail_in = (uInt)cdata_len;
  zs.next_out = ublock;
  zs.avail_out = (uInt)BGZF_MAX_BLOCK;
  int zr = inflate(&zs, Z_FINISH);
  if (zr != Z_STREAM_END || zs.total_out != isize)
    GError("Error (GBgzf): block inflate failed (zlib %d)\n", zr);
  ublock_len = (int)isize;
  ublock_off = 0;
  return ublock_len;
}

int64_t GBgzf::useek(uint64_t uoffset) {
  if (!have_idx)
    GError("Error (GBgzf): useek requires a loaded .gzi index\n");
  pushback = -1;
  // fast path: target already within the currently loaded block
  if (ublock_len > 0 && uoffset >= cur_uaddr &&
      uoffset < cur_uaddr + (uint64_t)ublock_len) {
    ublock_off = (int)(uoffset - cur_uaddr);
    return (int64_t)uoffset;
  }
  // binary search for the last block whose uaddr <= uoffset
  int64_t lo = 0, hi = idx.Count() - 1, best = 0;
  while (lo <= hi) {
    int64_t mid = (lo + hi) / 2;
    if (idx[mid].uaddr <= uoffset) { best = mid; lo = mid + 1; }
    else hi = mid - 1;
  }
  if (fseeko(fp, (off_t)idx[best].caddr, SEEK_SET) != 0)
    GError("Error (GBgzf): seek failed\n");
  fill_block();
  cur_uaddr = idx[best].uaddr;
  int off = (int)(uoffset - cur_uaddr);
  if (off > ublock_len) off = ublock_len; //clamp (request past block end)
  ublock_off = off;
  return (int64_t)uoffset;
}

int64_t GBgzf::read(void* buf, int64_t len) {
  if (len <= 0) return 0;
  unsigned char* out = (unsigned char*)buf;
  int64_t total = 0;
  if (pushback >= 0) { out[total++] = (unsigned char)pushback; pushback = -1; }
  while (total < len) {
    if (ublock_off >= ublock_len) { //current block exhausted
      cur_uaddr += (uint64_t)ublock_len;
      int r = fill_block();
      if (r < 0) break; //physical EOF
      continue;         //re-evaluate (handles empty blocks too)
    }
    int64_t avail = ublock_len - ublock_off;
    int64_t want = len - total;
    int64_t cp = (want < avail) ? want : avail;
    memcpy(out + total, ublock + ublock_off, (size_t)cp);
    ublock_off += (int)cp;
    total += cp;
  }
  return total;
}

int GBgzf::getc() {
  if (pushback >= 0) { int c = pushback; pushback = -1; return c; }
  unsigned char b;
  if (read(&b, 1) == 1) return (int)b;
  return -1;
}

void GBgzf::ungetc(int c) {
  pushback = c;
}

bool GBgzf::loadIndex(const char* gziname) {
  FILE* f = fopen(gziname, "rb");
  if (f == NULL) return false;
  unsigned char cnt8[8];
  if (fread(cnt8, 1, 8, f) != 8) { fclose(f); return false; }
  uint64_t n = (uint64_t)le_u32(cnt8) | ((uint64_t)le_u32(cnt8 + 4) << 32);
  idx.Clear();
  GBgzf::GBgzBlockIdx first; first.caddr = 0; first.uaddr = 0;
  idx.Add(first); //implicit first block at {0,0}
  for (uint64_t i = 0; i < n; i++) {
    unsigned char rec[16];
    if (fread(rec, 1, 16, f) != 16) { fclose(f); return false; }
    GBgzf::GBgzBlockIdx e;
    e.caddr = (uint64_t)le_u32(rec)     | ((uint64_t)le_u32(rec + 4) << 32);
    e.uaddr = (uint64_t)le_u32(rec + 8) | ((uint64_t)le_u32(rec + 12) << 32);
    idx.Add(e);
  }
  fclose(f);
  have_idx = true;
  return true;
}

bool GBgzf::openWrite(const char* fname, int level) {
  fp = fopen(fname, "wb");
  if (fp == NULL) return false;
  is_write = true;
  wlevel = (level < 0) ? Z_DEFAULT_COMPRESSION : level;
  GMALLOC(wbuf, BGZF_BLOCK);
  wbuf_len = 0;
  wtotal_u = 0;
  wtotal_c = 0;
  zs.zalloc = Z_NULL; zs.zfree = Z_NULL; zs.opaque = Z_NULL;
  if (deflateInit2(&zs, wlevel, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY) != Z_OK)
    GError("Error (GBgzf): deflateInit2 failed\n");
  zs_init = true;
  idx.Clear();
  GBgzf::GBgzBlockIdx first; first.caddr = 0; first.uaddr = 0;
  idx.Add(first); //first block starts at {0,0}
  return true;
}

// deflate the accumulated wbuf as a single BGZF block and write it out
bool GBgzf::flush_block() {
  if (wbuf_len == 0) return true;
  unsigned char cbuf[BGZF_MAX_BLOCK];
  deflateReset(&zs);
  zs.next_in = wbuf;
  zs.avail_in = (uInt)wbuf_len;
  zs.next_out = cbuf;
  zs.avail_out = (uInt)BGZF_MAX_BLOCK;
  int zr = deflate(&zs, Z_FINISH);
  if (zr != Z_STREAM_END)
    GError("Error (GBgzf): deflate failed (zlib %d)\n", zr);
  int cdata_len = (int)zs.total_out;
  int bsize = BGZF_HDR_LEN + cdata_len + BGZF_FTR_LEN;
  if (bsize > BGZF_MAX_BLOCK)
    GError("Error (GBgzf): compressed block exceeds max size\n");
  unsigned char hdr[BGZF_HDR_LEN];
  memcpy(hdr, bgzf_hdr, BGZF_HDR_LEN);
  put_u16(hdr + 16, (uint16_t)(bsize - 1)); //BSIZE
  unsigned char ftr[BGZF_FTR_LEN];
  put_u32(ftr, (uint32_t)crc32(crc32(0L, Z_NULL, 0), wbuf, (uInt)wbuf_len));
  put_u32(ftr + 4, (uint32_t)wbuf_len); //ISIZE
  if (fwrite(hdr, 1, BGZF_HDR_LEN, fp) != BGZF_HDR_LEN ||
      fwrite(cbuf, 1, (size_t)cdata_len, fp) != (size_t)cdata_len ||
      fwrite(ftr, 1, BGZF_FTR_LEN, fp) != BGZF_FTR_LEN)
    GError("Error (GBgzf): write failed\n");
  wtotal_c += (uint64_t)bsize;
  wtotal_u += (uint64_t)wbuf_len;
  wbuf_len = 0;
  // record the start offsets of the NEXT block
  GBgzf::GBgzBlockIdx e; e.caddr = wtotal_c; e.uaddr = wtotal_u;
  idx.Add(e);
  return true;
}

int64_t GBgzf::write(const void* buf, int64_t len) {
  const unsigned char* in = (const unsigned char*)buf;
  int64_t done = 0;
  while (done < len) {
    int space = BGZF_BLOCK - wbuf_len;
    int64_t want = len - done;
    int cp = (want < space) ? (int)want : space;
    memcpy(wbuf + wbuf_len, in + done, (size_t)cp);
    wbuf_len += cp;
    done += cp;
    if (wbuf_len == BGZF_BLOCK) flush_block();
  }
  return done;
}

bool GBgzf::dumpIndex(const char* gziname) {
  // flush any pending data first so the index is complete
  if (is_write && wbuf_len > 0) flush_block();
  FILE* f = fopen(gziname, "wb");
  if (f == NULL) return false;
  // idx holds a trailing {wtotal_c,wtotal_u} entry for the (not-yet-written)
  // block after the last data block; and idx[0] is the implicit {0,0}. The .gzi
  // stores only the real block-start offsets after the first: idx[1 .. count-2].
  int64_t nblocks = idx.Count() - 2; //exclude idx[0]={0,0} and the trailing sentinel
  if (nblocks < 0) nblocks = 0;
  unsigned char cnt8[8];
  put_u32(cnt8, (uint32_t)((uint64_t)nblocks & 0xffffffffULL));
  put_u32(cnt8 + 4, (uint32_t)(((uint64_t)nblocks >> 32) & 0xffffffffULL));
  if (fwrite(cnt8, 1, 8, f) != 8) { fclose(f); return false; }
  for (int64_t i = 1; i <= nblocks; i++) {
    unsigned char rec[16];
    put_u32(rec, (uint32_t)(idx[i].caddr & 0xffffffffULL));
    put_u32(rec + 4, (uint32_t)((idx[i].caddr >> 32) & 0xffffffffULL));
    put_u32(rec + 8, (uint32_t)(idx[i].uaddr & 0xffffffffULL));
    put_u32(rec + 12, (uint32_t)((idx[i].uaddr >> 32) & 0xffffffffULL));
    if (fwrite(rec, 1, 16, f) != 16) { fclose(f); return false; }
  }
  fclose(f);
  return true;
}

bool GBgzf::close() {
  if (fp == NULL) return true;
  if (is_write) {
    if (wbuf_len > 0) flush_block();
    fwrite(bgzf_eof, 1, sizeof(bgzf_eof), fp); //EOF marker block
  }
  if (zs_init) {
    if (is_write) deflateEnd(&zs); else inflateEnd(&zs);
    zs_init = false;
  }
  fclose(fp);
  fp = NULL;
  return true;
}
