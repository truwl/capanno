from collections import OrderedDict
from itertools import chain


ENCODE_atac_seq = OrderedDict([
    ('ST_43baaf.f7', 'cwl-scripts/ENCODE_DCC/atac-seq-pipeline/v1.1/encode_ataqc.cwl'),
    ('ST_49b9f4.f7', 'cwl-scripts/ENCODE_DCC/atac-seq-pipeline/v1.1/encode_bam2ta.cwl'),
    ('ST_cff563.f7', 'cwl-scripts/ENCODE_DCC/atac-seq-pipeline/v1.1/encode_bowtie2.cwl'),
    ('ST_5f34e1.f7', 'cwl-scripts/ENCODE_DCC/atac-seq-pipeline/v1.1/encode_filter.cwl'),
    ('ST_8646fc.f7', 'cwl-scripts/ENCODE_DCC/atac-seq-pipeline/v1.1/encode_idr.cwl'),
    ('ST_da7194.f7', 'cwl-scripts/ENCODE_DCC/atac-seq-pipeline/v1.1/encode_macs2_atac.cwl'),
    ('ST_29b258.f7', 'cwl-scripts/ENCODE_DCC/atac-seq-pipeline/v1.1/encode_naive_overlap.cwl'),
    ('ST_78d274.f7', 'cwl-scripts/ENCODE_DCC/atac-seq-pipeline/v1.1/encode_pool_ta.cwl'),
    ('ST_3b3a64.f7', 'cwl-scripts/ENCODE_DCC/atac-seq-pipeline/v1.1/encode_qc_report.cwl'),
    ('ST_ae51d1.f7', 'cwl-scripts/ENCODE_DCC/atac-seq-pipeline/v1.1/encode_reproducibility_qc.cwl'),
    ('ST_14c72b.f7', 'cwl-scripts/ENCODE_DCC/atac-seq-pipeline/v1.1/encode_spr.cwl'),
    ('ST_263ee1.f7', 'cwl-scripts/ENCODE_DCC/atac-seq-pipeline/v1.1/encode_trim_adapter.cwl'),
    ('ST_23a390.f7', 'cwl-scripts/ENCODE_DCC/atac-seq-pipeline/v1.1/encode_xcor.cwl'),
])

# ENCODE_DCC = None

GA4GH_Workflow_Execution_Challenge_hello_world_with_checker = OrderedDict([
    ('ST_99b1ff.e4', 'cwl-scripts/GA4GH_Workflow_Execution_Challenge/hello_world_with_checker/1.0/bash-hello_world.cwl'),
    ('ST_d03368.e4', 'cwl-scripts/GA4GH_Workflow_Execution_Challenge/hello_world_with_checker/1.0/python-helloworld_check.cwl')
])

GA4GH_Workflow_Execution_Challenge_md5sum_checker = OrderedDict([
    ('ST_493ada.e4', 'cwl-scripts/GA4GH_Workflow_Execution_Challenge/md5_sum_checker/master_2018_3_27/check_md5sum.cwl'),
    ('ST_71882e.e4', 'cwl-scripts/GA4GH_Workflow_Execution_Challenge/md5_sum_checker/master_2018_3_27/my_md5sum.cwl')
])

GA4GH_Workflow_Execution_Challenge = OrderedDict(chain(
    GA4GH_Workflow_Execution_Challenge_hello_world_with_checker.items(),
    GA4GH_Workflow_Execution_Challenge_md5sum_checker.items(),
))

everything = OrderedDict(chain(
    ENCODE_atac_seq.items(),
    GA4GH_Workflow_Execution_Challenge.items()
))
