# gclib-core

`gclib-core` is the shared low-level subset of the `gclib` codebase used by
`gffcompare` and other tools.

This repository is the single source of truth for the common headers, template
classes, and core implementation files that were previously copied into
downstream projects.

## Consumer layout

- `gffcompare` consumes this repository directly as its `gclib/` submodule.
- The full `gclib` repository consumes this repository as its `core/`
  submodule and keeps non-core code in the parent repository.

## Maintenance rule

Changes to shared code must be made here first. Downstream repositories should
only update their submodule pointer and must not keep edited copies of these
files.

## Project overview

This repository contains the small, dependency-light C++ core that downstream
tools use for:

- generic container and utility code
- FASTA indexing and random-access sequence retrieval
- genomic interval handling
- GFF/GTF/BED parsing and transcript models
- nucleotide and codon helpers

The code is written in an older `gclib` style that favors lightweight custom
data structures over the C++ standard library. This is intentional and should
be preserved unless there is a strong reason to change it.

## Main components

### Base utilities

- `GBase.h` / `GBase.cpp`
  Core utility layer: memory helpers, error handling, string helpers, parsing
  helpers, file helpers, and small foundational types like `GSeg`.
- `GArgs.*`
  Minimal argument parsing utilities used by tools built on top of the library.
- `GStr.*`
  Custom mutable string class used throughout the codebase instead of
  `std::string`.

### Core containers

- `GVec.hh`
  Dynamic array and pointer-vector templates (`GVec`, `GPVec`).
- `GList.hh`
  Sorted/unsorted list and array-style wrappers built on the vector layer.
- `GHashMap.hh`
  Hash tables and sets used for name lookup and record indexing.
- `GBitVec.h`
  Compact bit-vector implementation used by overlap and splice-structure logic.

These containers are the project’s normal abstraction layer. Prefer using them
consistently instead of introducing STL containers into the shared core.

### Sequence and interval handling

- `GFastaIndex.*`
  FASTA `.fai` indexing: stores sequence length, file offsets, and line layout.
- `GFaSeqGet.*`
  FASTA random access: retrieves genomic subsequences efficiently from indexed
  FASTA files.
- `GIntervalTree.hh`
  Interval tree support for genomic interval queries.
- `gdna.*`
  DNA helpers such as reverse-complement and IUPAC base handling.
- `codons.*`
  Codon translation helpers.

### Annotation model and parsing

- `gff.h` / `gff.cpp`
  The main annotation layer:
  - line parsing for GFF3, GTF, and BED-derived transcript records
  - feature/name registries
  - exon/CDS/transcript object model (`GffObj`, `GffLine`, `GffReader`)
  - transcript comparison, overlap classification, spliced sequence extraction,
    CDS projection, and summary/printing helpers

This is the largest and most behavior-rich unit in the repository.

## Coordinate and size conventions

Genomic coordinates and coordinate-derived lengths in the core code are now
expected to use `int64_t`. Collection sizes and indexes in the core containers
were also widened to 64-bit.

When touching coordinate-sensitive code, preserve 64-bit arithmetic end-to-end,
especially in:

- `GSeg` and interval math
- FASTA offsets and sequence lengths
- exon/CDS/transcript spans
- container counts and indexes

## Coding style

The codebase follows a pragmatic C-with-classes style:

- prefer simple headers and inline templates
- avoid exceptions and RTTI
- avoid heavy STL usage in shared core code
- prefer existing `gclib` classes (`GVec`, `GList`, `GHash`, `GStr`, `GSeg`)
  over introducing new abstractions
- keep dependencies minimal and portability high
- preserve existing naming and formatting patterns when editing older files

Memory management is explicit in much of the code. Ownership conventions matter,
especially for pointer containers like `GPVec` and `GList`.

## Repo role

This repository is the canonical place for shared core changes. If a bug fix or
API change affects code mirrored into downstream tools, it should be made here
first and then propagated by submodule update.
