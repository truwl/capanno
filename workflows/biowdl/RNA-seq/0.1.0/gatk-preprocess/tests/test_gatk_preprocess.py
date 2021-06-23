from pathlib import Path

import pysam
import pytest


@pytest.mark.workflow("recalibrated_bam")
def test_applybqsr_used(workflow_dir):
    bam_path = workflow_dir / Path("test-output") / Path("test.bam")
    bam_file = pysam.AlignmentFile(str(bam_path), "rb")
    programs = [program.get('ID') for program in bam_file.header.get('PG')]
    assert "GATK ApplyBQSR" in programs
    assert "GATK SplitNCigarReads" not in programs


@pytest.mark.workflow("split_n_cigar", "no_scatter")
def test_splitncigar_used(workflow_dir):
    bam_path = workflow_dir / Path("test-output") / Path("test.bam")
    bam_file = pysam.AlignmentFile(str(bam_path), "rb")
    programs = [program.get('ID') for program in bam_file.header.get('PG')]
    assert "GATK ApplyBQSR" in programs
    assert "GATK SplitNCigarReads" in programs
