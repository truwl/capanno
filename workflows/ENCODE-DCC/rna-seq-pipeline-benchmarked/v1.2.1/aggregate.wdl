version 1.0

task melt {
  input {
      String job_id
      String workflow_instance_identifier
      String workflow_identifier
	  String rep
      File qcfile
      File Rscript_aggregate
  }
  output {
    File talltable = "truwlbenchmarks{~rep}.txt"
  }
  command <<<
    Rscript ~{Rscript_aggregate} ~{job_id} ~{workflow_instance_identifier} ~{workflow_identifier} ~{rep} ~{qcfile} truwlbenchmarks{~rep}.txt
  >>>
  runtime {
    docker: "rocker/tidyverse:4.1.0"
    memory: "1 GB"
    cpu: 1
  }
}
