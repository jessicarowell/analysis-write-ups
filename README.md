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

## Translate Transcript 2 Genomic Coordinates

This script takes two input files that result from a mapping a set of transcripts to a reference genome.  It parses the CIGAR string and maps the transcript positions to their corresponding genomic positions.  There are packages to do this, but I did it manually to demonstrate an example of a complete script with logging.  It isn't built for efficiency since it parses the string multiple times per transcript.  

Usage: `python3 translate-t2g.py -h`

## Machine Learning Example

Example available upon request.
