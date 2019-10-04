# PeakZCL

PeakZCL is a pipeline that implements the Zero-Crossing Lines method based on wavelets for analyzing NGS signals.

## Installation

This script need no installation, but has the following requirements:
    -R version 3.5.3 or later
    -R packages:
	data.table
	signal
	optparse
	fields
	wavethresh
	Rwave
	wavelets
	wmtsa
	multtest
	VennDiagram
	IRanges
	ChIPpeakAnno


### This tool analyzes signals with coverage at nucleotide level (SEE option -d of the tool genomeCoverageBed).



##  BASIC Usage
	./PeakZCL.sh -f /path/to/CHIP/directory/

## Before running your analysis you must to copy the file "chr_names.txt" to the CHIP directory. In this file you can select, (writing it or deleting the rest) the chromosome/s that you need to analyze. 

For more information type:
	./PeakZCL.sh

## License
MIT License

Copyright (c) 2019 José González Gomariz

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

