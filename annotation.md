# Denovo repeat library
A comprehensive denovo repeat library was prepared for the genome for repeat characterization.

Workflow:

## A. Repeat library preparation
1. [Repeatmoduler](https://github.com/upendrabhattarai/Earwig_genome_project/blob/main/Denovo_repeat_library/Repeatmoduler.md)
2. [LTRharvest & LTRdigest](https://github.com/upendrabhattarai/Earwig_genome_project/blob/main/Denovo_repeat_library/LTRharvest%26LTRdigest.md)
3. TransposonPSI
4. Sine database
5. Concatenating, filtering, and classifying repeats
6. RepeatClassifier
7. Repeat masking the genome
8. RepeatMasker

   
## B. Genome annotation
Maker2 pipeline was used for genome annotation. mRNA-seq data is denovo assembled using Trinity. Other relavant publicly available datasets were downloaded and used as input.

1. Processing mRNA-seq data
2. GeneMark-ES
3. Braker
4. Configuring and running Maker2 pipeline
