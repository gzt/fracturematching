# fracturematching
Repository of Data and Code for Fracture Matching Using Matrix Variate Classification

This repository contains the code and data to reproduce the analyses and figures in the paper.

* The R package version used for the paper is contained in `MixMatrix_0.1.0.tar.gz`. The other code is 
written presuming this package has been unpacked in this directory.
* The underlying data is contained in the files with filenames starting with "correlation" in the `/data/` folder.
The three repetitions on one knife dataset are labelled `Correlationsset1.csv`, `Correlationsset2.csv`, and `Correlationsset3.csv`, 
The two steel bar data sets are labelled with `S-`. The data fields are:
	- `5-10` - correlations in the 5-10 $mm^{-1}$ range 
	- `10-20` - correlations in the 10-20 $mm^{-1}$ range
	- `knife` - the labels of the base:tip parts.
	- `img` - which image in the sequence (from 1 to 9)
	- `match` - whether it is a `match` or `nonmatch`
	- `set` - the label of the data set
* The file `Nature-Figures-Compilation.Rmd` is an R Markdown file which will read the data and output to 
produce most of the figures in the paper. This includes some additional diagnostic plots. It depends on
results from the CV analysis and the overlap analysis -- results from those are included in `/data/` and
`/oddsconsec/` but code to reproduce those analyses is included.
* The file `CV-Analysis.R` is an R script to reproduce the cross validation analysis and 
related figures. The results have been placed in the `/data/` folder.
* The files `OVERLAP-analysis-and-results.Rmd`, `OVERLAP-analysis-and-results-13579.Rmd`, 
and `OVERLAP-analysis-and-results-159.Rmd` perform the overlap analysis for the 9 image, 5 image (50\% overlap),
and 3 image (0\% overlap) settings for one setting of the $df$ parameter. This parameter can be edited 
to produce the entire analysis from the paper. The results have been placed in the `/oddsconsec/` folder.
