version 1.0

# Copyright (c) 2018 Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
import "tasks/bedtools.wdl" as bedtools
import "tasks/chunked-scatter.wdl" as chunkedscatter

workflow CalculateRegions {
    input {
        File referenceFasta
        File referenceFastaDict
        File referenceFastaFai
        File? XNonParRegions
        File? YNonParRegions
        File? regions
        Int scatterSizeMillions = 1000
        # scatterSize is on number of bases. The human genome has 3 000 000 000 bases.
        # 1 billion gives approximately 3 scatters per sample.
        Int? scatterSize
        Map[String, String] dockerImages = {
            "bedtools": "quay.io/biocontainers/bedtools:2.23.0--hdbcaa40_3",
            "chunked-scatter": "quay.io/biocontainers/chunked-scatter:1.0.0--py_0"
        }
    }
    meta {allowNestedInputs: true}

    Boolean knownParRegions = defined(XNonParRegions) && defined(YNonParRegions)

    if (knownParRegions) {
        call bedtools.MergeBedFiles as mergeBeds {
            input:
                bedFiles = select_all([XNonParRegions, YNonParRegions]),
                dockerImage = dockerImages["bedtools"]
        }
        # We define the 'normal' regions by creating a regions file that covers
        # everything except the XNonParRegions and the YNonParRegions.
        call bedtools.Complement as inverseBed {
            input:
                inputBed = mergeBeds.mergedBed,
                faidx = referenceFastaFai,
                outputBed = "autosomal_regions.bed",
                dockerImage = dockerImages["bedtools"]
        }


        if (defined(regions)) {
            call bedtools.Intersect as intersectAutosomalRegions {
                input:
                    regionsA = inverseBed.complementBed,
                    regionsB = select_first([regions]),
                    faidx = referenceFastaFai,
                    outputBed = "intersected_autosomal_regions.bed",
                    dockerImage = dockerImages["bedtools"]
            }

            call bedtools.Intersect as intersectX {
                input:
                    regionsA = select_first([XNonParRegions]),
                    regionsB = select_first([regions]),
                    faidx = referenceFastaFai,
                    outputBed = "intersected_x_non_par_regions.bed",
                    dockerImage = dockerImages["bedtools"]
            }

            call bedtools.Intersect as intersectY {
                input:
                    regionsA = select_first([YNonParRegions]),
                    regionsB = select_first([regions]),
                    faidx = referenceFastaFai,
                    outputBed = "intersected_y_non_par_regions.bed",
                    dockerImage = dockerImages["bedtools"]
            }
        }
    }

    # When there are non-PAR regions and there are regions of interest, use the intersect of the autosomal regions and the regions of interest.
    # When there are non-PAR regions and there are no specified regions of interest, use the autosomal regions.
    # When there are no non-PAR regions, use the optional regions parameter.
    File? calculatedAutosomalRegions = if knownParRegions
                      then select_first([intersectAutosomalRegions.intersectedBed, inverseBed.complementBed])
                      else regions

    call chunkedscatter.ScatterRegions as scatterAutosomalRegions {
        input:
            inputFile = select_first([calculatedAutosomalRegions, referenceFastaFai]),
            scatterSize = scatterSize,
            scatterSizeMillions = scatterSizeMillions,
            dockerImage = dockerImages["chunked-scatter"]
    }

    output {
        File? Xregions = if defined(regions) then intersectX.intersectedBed else XNonParRegions
        File? Yregions = if defined(regions) then intersectY.intersectedBed else YNonParRegions
        File? autosomalRegions = calculatedAutosomalRegions
        Array[File] autosomalRegionScatters = scatterAutosomalRegions.scatters
    }

    parameter_meta {
        regions: {description: "A bed file describing the regions to operate on. Will be used to intersect the other regions.", category: "common"}
        XNonParRegions: {description: "Bed file with the non-PAR regions of X", category: "common"}
        YNonParRegions: {description: "Bed file with the non-PAR regions of Y", category: "common"}
        referenceFasta: { description: "The reference fasta file", category: "required" }
        referenceFastaFai: { description: "Fasta index (.fai) file of the reference", category: "required" }
        referenceFastaDict: { description: "Sequence dictionary (.dict) file of the reference", category: "required" }
        dockerImages: { description: "specify which docker images should be used for running this pipeline", category: "advanced" }
        scatterSize: {description: "The size of the scattered regions in bases. Scattering is used to speed up certain processes. The genome will be seperated into multiple chunks (scatters) which will be processed in their own job, allowing for parallel processing. Higher values will result in a lower number of jobs. The optimal value here will depend on the available resources.",
              category: "advanced"}
        scatterSizeMillions:{ description: "Same as scatterSize, but is multiplied by 1000000 to get scatterSize. This allows for setting larger values more easily",
                              category: "advanced"}
    
    }
}
