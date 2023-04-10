# Analysis Write-ups

A repository containing work summaries for various projects I have done.  Most will be statistics or bioinformatics projects completed for others.  I may also store summaries of personal projects here in the future.
 
Started by [@jessicarowell](https://github.com/jessicarowell)

## TOC
* [Baitcap Analysis](#baitcap-analysis)
* [Translate Transcript 2 Genomic Coordinates](#translate-transcript-2-genomic-coordinates)
* [Machine Learning Example](#machine-learning-example)

## Baitcap Analysis

This is an analysis done to account for a change in wet lab methodology that happened in this project.
A detailed description is included in the html file.

Topic(s): R software, laboratory, basic statistics
File(s): baitcap_analysis_2021.html 

## Translate Transcript 2 Genomic Coordinates

This script takes two input files that result from a mapping a set of transcripts to a reference genome.  It parses the CIGAR string and maps the transcript positions to their corresponding genomic positions.  There are packages to do this, but I did it manually to demonstrate an example of a complete script with logging.  It isn't built for efficiency since it parses the string multiple times per transcript.  

Topic(s): Python, genomics, scripting
Usage: `python3 translate-t2g.py -h`

## Machine Learning Example

This R Markdown report explores the *ERBB2*/HER2 oncogene in relation to HER2 status among patients with breast cancer. The goal was to analyze HER2 in TCGA breast cancer data and assess how RNA expression, copy number, and clinical status are correlated and whether RNA or DNA is predictive of HER2 clinical status.

Topic(s): R software, gene expression, machine learning
File(s): HER2_breast_cancer_2022.html, HER2_breast_cancer_2022.pdf, HER2_breast_cancer_2022.Rmd
