# PrEST standard evaluation script

This repository contains the evaluation script for the PrEST standards for protein inference calibration.

Usage:
```
python generate_prest_report.py <results_file> <vial>
```
* The results file is a file starting with a header line, followed by lines: protein1 \<tab\> q-value1. An example results file is included in the repository (`vialA.tsv`)
* Vial is one of 'A', 'B' or 'AB'.
* Protein groups can be specified as a protein identifier of the form "protein1a,protein1b" (without quotes). The employed null hypothesis is that at least one of the proteins in the protein group is present.

