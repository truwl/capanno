# Used for creating and updating db items.

from collections import OrderedDict
from itertools import chain

gnu_tools = OrderedDict([
    ('TL_d077f2.47', 'cwl-tools/cat/8.25/cat/cat.cwl'),
    ('TL_cbb11e.47', 'cwl-tools/echo/8.25/echo/echo.cwl'),
    ('TL_9579db.89', 'cwl-tools/gawk/4.1/gawk/gawk.cwl'),
    ('TL_c8f1ee.47', 'cwl-tools/md5sum/8.25/common/md5sum-metadata.yaml'),
    ('TL_c8f1ee_d4.47', 'cwl-tools/md5sum/8.25/md5sum/md5sum.cwl'),
    ('TL_c8f1ee_0b.47', 'cwl-tools/md5sum/8.25/md5sum_check/md5sum-check.cwl'),
    ('TL_cadc8c.47', 'cwl-tools/sort/8.25/sort/sort.cwl'),
    ('TL_e7d707.47', 'cwl-tools/tr/8.25/tr/tr.cwl'),
])


samtools = OrderedDict([
    ('TL_ec2a8d.0b', 'cwl-tools/samtools/1.3/common/samtools-metadata.yaml'),
    ('TL_ec2a8d_a6.0b', 'cwl-tools/samtools/1.3/samtools_flagstat/samtools-flagstat.cwl'),
    ('TL_ec2a8d_6a.0b', 'cwl-tools/samtools/1.3/samtools_index/samtools-index.cwl'),
    ('TL_ec2a8d_ca.0b', 'cwl-tools/samtools/1.3/samtools_sort/samtools-sort.cwl'),
    ('TL_ec2a8d_46.0b', 'cwl-tools/samtools/1.3/samtools_view/samtools-view.cwl'),
])

STAR = OrderedDict([
    ('TL_8ab263.82', 'cwl-tools/STAR/2.5/common/STAR-metadata.yaml'),
    ('TL_8ab263_a4.82', 'cwl-tools/STAR/2.5/STAR_alignReads/STAR-alignReads.cwl'),
    ('TL_8ab263_95.82', 'cwl-tools/STAR/2.5/STAR_genomeGenerate/STAR-genomeGenerate.cwl'),
    ('TL_8ab263_d6.82', 'cwl-tools/STAR/2.5/STAR_inputAlignmentsFromBAM/STAR-inputAlignmentsFromBAM.cwl'),

])

trimmomatic = OrderedDict([
    ('TL_a2e03a.af', 'cwl-tools/trimmomatic/0.38/common/trimmomatic-metadata.yaml'),
    ('TL_a2e03a_53.af', 'cwl-tools/trimmomatic/0.38/trimmomatic_SE/trimmomatic-SE.cwl'),
])

everything = OrderedDict(chain(
    gnu_tools.items(),
    samtools.items(),
    STAR.items(),
    trimmomatic.items(),
))
