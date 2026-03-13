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
