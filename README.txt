Demo of additive variance code

*******************

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

*******************

