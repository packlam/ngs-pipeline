# ngs-pipeline

This repository contains some of the Python code that I wrote while I was working as a bioinformatics analyst for a major healthcare provider.

* **coverage.py** makes use of the [Pandas](https://pandas.pydata.org/) library to read in data regarding individual samples and calculate the [coverage](https://en.wikipedia.org/wiki/Coverage_(genetics)) level at specific areas of interest defined in a [BED file](https://www.ensembl.org/info/website/upload/bed.html)
* **generate_results.py** reads in raw data from a [VCF file](https://en.wikipedia.org/wiki/Variant_Call_Format) and [Alamut Batch](https://www.interactive-biosoftware.com/alamut-batch/) before merging the required data columns, carrying out a number of different filtering steps, and producing a human-readable tab-separated data file at the end
* **run_pipeline.py** is a bioinformatics pipeline script that takes raw DNA sequencing data and runs it through various third party software, as well as the above Python scripts, in order to generate all of the files required for analysing and interpreting the genomic data for each sample
