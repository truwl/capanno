Changelog
==========

<!--
Newest changes should be on top.

This document is user facing. Please word the changes in such a way
that users understand how the changes affect the new version.
-->

version 1.4.0
---------------------------
+ Replaced Stringtie merge with gffcompare for merging of stringtie
  assembly results.
+ Downgrade stringtie to version 1.3.6.
+ Make it possible to not run Stringtie.
+ When multiple TPM/FPKM values are returned for a single gene by
  stringtie, they will now be added together in the multi-sample
  expression tables. Previously only the last value encountered would be
  used.
+ Updated default docker image for collect-columns (now uses
  version 1.0.0 instead of 0.2.0).


version 1.3.0
---------------------------
+ Updated default docker images.
+ Tasks were updated to contain the `time_minutes` runtime attribute and
  associated `timeMinutes` input, describing the maximum time the task will
  take to run.

version 1.2.0
---------------------------
+ Add proper copyright headers to WDL files. So the free software license
  is clear to end users who wish to adapt and modify.
+ Added inputs overview to the docs.
+ Added parameter_meta.
+ Added wdl-aid to linting.
+ Updated default htseq image to version 0.11.2.
+ Add miniwdl to linting.

version 1.1.0
---------------------------
+ Update tasks so they pass the correct memory requirements to the
  execution engine. Memory requirements are set on a per-task (not
  per-core) basis.

version 1.0.0
---------------------------
+ Updated documentation.
