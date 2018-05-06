# pcawgAdditiveVariance

This repository consist of code relevant for additive variance analysis performed on PCAWG mutations.

Following dependencies are required to run this workflow.

FunSeq2 (http://funseq2.gersteinlab.org/)

Python

Matlab

GCTA (http://cnsgenomics.com/software/gcta/#Overview)

This workflow consist of two components.

1) pre-processing step 
    In this step we generate summary file for each cancer cohort.
   
    *******************
    
    Two input files are needed for this step. 
    a) PCAWG driver mutation list and b)FunSeq2 output file(in BED format) 
    
    Usage:
    generateSummaryInfo.py -d <driverFile> -I <funSeqOutFile> -O <outSummaryFile>
    generateSummaryInfo.py (-h | --help)


2) post-processing steps

   Pipeline is run by calling additive_variance_demo.m

   *******************

   Full pipeline requires the following inputs:

    - in bedFiles folder:
      cohortName.null.bed
      cohortName.obs.bed

   - in summaryFiles folder:
     cohortName.null.summary.txt
     cohortName.obs.summary.txt

   and generates the output cohortName.txt in the results folder.

  A gcta executable is required (Linux version is included, Windows and Mac versions are available from cnsgenomics.com/software/gcta), which should be placed in the gctaFiles folder.

*******************

Current pipeline is set up to call only final stage, which summarizes the results from precomputed intermediate outputs.
Outputs are computed from the Breast-AdenoCa cohort with randomized samples used for both null and obs conditions.
Results text file shows calculated additive variance for each funseq threshold, which is ~0 (1e-6), along with associated p-values (0.5 indicates that no significant genetic variance was found).
Values of -1 in the results file for funseq thresholds 5 and 6 indicate that insufficient data was found at these thresholds.
