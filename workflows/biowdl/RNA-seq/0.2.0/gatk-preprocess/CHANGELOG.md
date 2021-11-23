Changelog
==========

<!--

Newest changes should be on top.

This document is user facing. Please word the changes in such a way
that users understand how the changes affect the new version.
-->

version 3.0.0
-----------------
+ Increase the time needed for baserecalibration
+ Update default docker images.
+ The scatter regions has been removed. Regions for scattering need to be 
  provided externally. This is more efficient for pipelines that work on 
  multiple samples.
+ Tasks were updated to contain the `time_minutes` runtime attribute and
  associated `timeMinutes` input, describing the maximum time the task will
  take to run.
+ Only scatter when the number of scatters turns out to be 1.


version 2.0.0
-----------------
+ Add a scatterSizeMillions parameter to make it easier to set larger scatter 
  sizes.
+ Add proper copyright headers to WDL files. So the free software license
  is clear to end users who wish to adapt and modify.
+ Remove redundant orderedscatters task. Ordering is now handled by the 
  scatterRegions task.
+ Remove option to recalibrate or not. The bam file is always recalibrated.
+ Correctly apply base recalibration after the reads have been spliced.
+ Remove structs from input and output.
+ Added inputs overview to the docs.
+ Added parameter_meta.
+ Added wdl-aid to linting.
+ Added miniwdl to linting.

version 1.1.0
---------------------------
+ Update tasks so they pass the correct memory requirements to the 
  execution engine. Memory requirements are set on a per-task (not
  per-core) basis.

version 1.0.0
---------------------------
+ fixed the md5sum key in gatheredBam (md5 -> md5sum)
+ Updated documentation to reflect latest changes
