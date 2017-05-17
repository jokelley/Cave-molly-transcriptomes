# Cave-molly-transcriptomes
Scripts for "Complexities of gene expression patterns in natural populations of an extremophile fish (Poecilia mexicana, Poeciliidae)"

Passow et al 2017, Molecular Ecology

ANALYSIS PIPELINE
 - Trim_map_expression_script.txt

EXTRACT COUNTS
 - extractcounts.pl

Used within the analysis pipeline

DIFFERENTIAL EXPRESSION ANALYSIS
 - DE_tissue_final.R

COUNTS MATRIX
 - tissue.txt.gz 

Transcript abundance was measured using the program eXpress (http://bio.math.berkeley.edu/eXpress/overview.html). 
For more information see methods section in Passow et al. 2017 and linux output (Trim_map_expression_scripts.txt)
To extract counts information we used the perl script extractcounts.pl


GO-RESULTS
 - GO-results folder

The folder contains files for section 3. Gene Function [Physiological pathway] - level within the R script DE_tissue_final.R
 
Gene ontology terms were obtained using BLAST2GO PRO (https://www.blast2go.com) with the follow protocol: 

(1) For the BLAST2GO input, we used the human swissprot accessions (as .txt files)

(2) We then uploaded the .txt files individually into Blast2Go PRO and ran the function "run annotation" to generate Gene ontology (GO) IDs.

(3) We then used the function Validate Annotation to remove redundant GO terms from the dataset. This was to ensure that only the most specific annotations for a given sequence were saved. 

(4) We then implemented the function "Remove First Level Annotations" to remove the three main top-level GO terms (molecular function, cellular component and biological process). 

(5) The annotations were then exported as .txt files from BLAST2GO and saved with a _GO.txt extension. 

BLAST OUTPUT
 - TableS4.csv

For all transcripts the BLAST output was saved to TableS4.csv. For BLAST output parameters, see methods in Passow et al. 2017 
