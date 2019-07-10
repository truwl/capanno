from collections import OrderedDict
from itertools import chain


ENCODE_atac_seq = OrderedDict([
    ('ST_43baaf.f7', 'cwl-scripts/ENCODE-DCC/atac-seq-pipeline/v1.1/encode_ataqc/encode_ataqc.cwl'),
    ('ST_49b9f4.f7', 'cwl-scripts/ENCODE-DCC/atac-seq-pipeline/v1.1/encode_bam2ta/encode_bam2ta.cwl'),
    ('ST_cff563.f7', 'cwl-scripts/ENCODE-DCC/atac-seq-pipeline/v1.1/encode_bowtie2/encode_bowtie2.cwl'),
    ('ST_5f34e1.f7', 'cwl-scripts/ENCODE-DCC/atac-seq-pipeline/v1.1/encode_filter/encode_filter.cwl'),
    ('ST_8646fc.f7', 'cwl-scripts/ENCODE-DCC/atac-seq-pipeline/v1.1/encode_idr/encode_idr.cwl'),
    ('ST_da7194.f7', 'cwl-scripts/ENCODE-DCC/atac-seq-pipeline/v1.1/encode_macs2_atac/encode_macs2_atac.cwl'),
    ('ST_29b258.f7', 'cwl-scripts/ENCODE-DCC/atac-seq-pipeline/v1.1/encode_naive_overlap/encode_naive_overlap.cwl'),
    ('ST_78d274.f7', 'cwl-scripts/ENCODE-DCC/atac-seq-pipeline/v1.1/encode_pool_ta/encode_pool_ta.cwl'),
    ('ST_3b3a64.f7', 'cwl-scripts/ENCODE-DCC/atac-seq-pipeline/v1.1/encode_qc_report/encode_qc_report.cwl'),
    ('ST_ae51d1.f7', 'cwl-scripts/ENCODE-DCC/atac-seq-pipeline/v1.1/encode_reproducibility_qc/encode_reproducibility_qc.cwl'),
    ('ST_14c72b.f7', 'cwl-scripts/ENCODE-DCC/atac-seq-pipeline/v1.1/encode_spr/encode_spr.cwl'),
    ('ST_263ee1.f7', 'cwl-scripts/ENCODE-DCC/atac-seq-pipeline/v1.1/encode_trim_adapter/encode_trim_adapter.cwl'),
    ('ST_23a390.f7', 'cwl-scripts/ENCODE-DCC/atac-seq-pipeline/v1.1/encode_xcor/encode_xcor.cwl'),
])



GA4GH_Workflow_Execution_Challenge_helloworld_checker = OrderedDict([
    ('ST_99b1ff.e4', 'cwl-scripts/GA4GH_Workflow_Execution_Challenge/helloworld/1.0.2/hello_world/hello_world.cwl'),
    ('ST_d03368.e4', 'cwl-scripts/GA4GH_Workflow_Execution_Challenge/helloworld-checker/1.1.2/helloworld_check/helloworld_check.cwl')
])

GA4GH_Workflow_Execution_Challenge_md5sum_checker = OrderedDict([
    ('ST_493ada.e4', 'cwl-scripts/GA4GH_Workflow_Execution_Challenge/md5sum-checker/master__bfd7e5c/md5sum-checker/md5sum-checker.cwl'),
    ('ST_71882e.e4', 'cwl-scripts/GA4GH_Workflow_Execution_Challenge/md5sum/master__ca3823d/my_md5sum/my_md5sum.cwl')
])

GA4GH_Workflow_Execution_Challenge = OrderedDict(chain(
    GA4GH_Workflow_Execution_Challenge_helloworld_checker.items(),
    GA4GH_Workflow_Execution_Challenge_md5sum_checker.items(),
))

everything = OrderedDict(chain(
    ENCODE_atac_seq.items(),
    GA4GH_Workflow_Execution_Challenge.items()
))
