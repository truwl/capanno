Changelog
==========

<!--

Newest changes should be on top.

This document is user facing. Please word the changes in such a way
that users understand how the changes affect the new version.
-->

version 2.1.0
---------------------------
+ Use the latest versions of samtools and picard by default.
+ Add a `reports` output for use in downstream pipelines.
+ Move `collectAlignmentSummaryMetrics` & `meanQualityByCycle` to top level
  workflow so they can be set to false for TALON-WDL pipeline.
+ Tasks were updated to contain the `time_minutes` runtime attribute and
  associated `timeMinutes` input, describing the maximum time the task will
  take to run.

version 2.0.0
---------------------------
+ Add proper copyright headers to WDL files. So the free software license
  is clear to end users who wish to adapt and modify.
+ Remove complex structs from the input.
+ Added inputs overview to the docs.
+ Added parameter_meta.
+ Added wdl-aid to linting.
+ Added miniwdl to linting.

version 1.2.0
---------------------------
+ Use the latest version of the tasks repository.

version 1.1.0
---------------------------
+ Update tasks so they pass the correct memory requirements to the 
  execution engine. Memory requirements are set on a per-task (not
  per-core) basis.


version 1.0.0
---------------------------
+ Update documentation to reflect latest changes
+ Picard: Use version 2.20.5 of the biocontainer as this includes the R dependency
