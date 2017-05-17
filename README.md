# Cave-molly-transcriptomes
Files for "Complexities of gene expression patterns in natural populations of an extremophile fish (Poecilia mexicana, Poeciliidae)"

Passow et al 2017, Molecular Ecology

## Analysis pipeline
 - Trim_map_expression_script.txt

Used within the analysis pipeline


## Differential expression analysis 

 - DE_tissue_final.R

## Counts matrix 
 - [tissue.txt.gz](https://github.com/jokelley/Cave-molly-transcriptomes/blob/master/tissue.txt.gz)

Transcript abundance was measured using the program eXpress (http://bio.math.berkeley.edu/eXpress/overview.html). 
For more information see methods section in Passow et al. 2017 and linux output (Trim_map_expression_scripts.txt)
To extract counts information we used the perl script [extractcounts.pl](https://github.com/jokelley/Cave-molly-transcriptomes/blob/master/extractcounts.pl)


## GO results 
 - GO-results folder

The folder contains files for section "3. Physiological pathways (gene function) shared between contrasting environmental conditions" within the R script DE_tissue_final.R
 
Gene ontology terms were obtained using BLAST2GO PRO (https://www.blast2go.com) with the follow protocol: 

(1) For the BLAST2GO input, we used the human swissprot accessions (as .txt files)

(2) We then uploaded the .txt files individually into BLAST2GO PRO and ran the function "run annotation" to generate Gene ontology (GO) IDs.

(3) We then used the function Validate Annotation to remove redundant GO terms from the dataset. This was to ensure that only the most specific annotations for a given sequence were saved. 

(4) We then implemented the function "Remove First Level Annotations" to remove the three main top-level GO terms (molecular function, cellular component and biological process). 

(5) The annotations were then exported as .txt files from BLAST2GO and saved with a _GO.txt extension. 

## BLAST output 
 - TableS4.csv

For all transcripts the BLAST output was saved to TableS4.csv. For BLAST output parameters, see methods in Passow et al. 2017 
