ENCODE Transcription Factor and Histone ChIP-Seq processing pipeline.

This ChIP-Seq pipeline is based off the ENCODE (phase-3) transcription factor and histone ChIP-seq pipeline
specifications (by Anshul Kundaje) in this [google doc](https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit#).

---
#### Job Costs
Costs for this workflow can vary widely. Primary determinants of cost include the number of samples (multiple can be specified in a single job), input file sizes, number of controls, and what options are chosen (workflow can only do alignment step, xcor, etc.). 

Examples:

 |  Public example     |   Job cost  | Notes |
| ------------------ | ---------------- |-------------------| 
|  [ENCSR135FYV Control ChIP-seq alignonly](https://bit.ly/3H6ze1L) | $ 12 | This is a control sample that only does alignment. (No peak calling) |
| [Re-run of ENCODE ChIP-seq analysis ENCSR777VNA](https://bit.ly/3FwnAgA) | $ 45 | This is a SE mouse sample with 2 samples and 2 controls. |
|[ENCODE ChIP-seq analysis ENCSR859FDL](https://bit.ly/3mqwnJq) | $ 96 | This is a PE human cell line sample with **2 biological replicates** and 2 controls.


