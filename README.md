AcroPoly, version 1.5, Ehsan Motazedi (ehsan.motazedi@gmail.com), last modified: August 10 2018

Introduction
============

**AcroPoly** \(executable *AcroPoly.py*\) is an expectation-maximisation \(EM\) based algorithm to obtain multi-allelic haplotype markers for \(a collection of\) diploid or polyploid individuals using next generation sequencing data. Haplotype markers are defined for windows consisting of *w* SNPs, so that each possible haplotype of the SNP alleles corresponds to a marker allele. The method considers a Poisson process for the number of reads that support each haplotype allele in a window and applies EM to estimate the rates of these Poisson processes \(while adjusting for the total sequencing depth\). Dosage scores are assigned in this way to each haplotype allele. For RNA sequence data, the estimated rates reflect the expression level of each haplotype while for DNA sequence data the estimates are normalized to obtain a fractional dosage \(between 0 and 1\) for each haplotype.

Given the reads mapped to a specific genomic or transcript region and the variants in the reads, the method scans the region with a sliding window of size, *w* \(specified by -*l*, --*length* option\) that is each time shifted *s* SNPs (specified by -*s*, --*shift* option\) to the right (starting from the 5' to 3' end) until the whole region is scanned. The haplotype alleles and their dosages are reported in a data table format for each region, so that they can be easily used for genetic association analysis. Alternatively, the best deterministic phasing estimates (obtained by rounding the dosage scores) can be reported for each haplotyping window using -*d*, --*deterministic* option. 

To analyze alignment files that contain several contigs, a wrapper *AcroPoly_wrapped.pl* can be used, which extracts the desired contig from the bam and vcf files and in addition can filter markers based upon their genotype missing rate in the population (see ./AcroPoly_wrapped.pl -h).

Citation:
=====================

To cite AcroPoly, please refer to *Estimating haplotype frequencies in polyploid populations from sequence data by expectation maximisation*. In Motazedi, E., 2019. Haplotype estimation in polyploids using DNA sequence data (Doctoral dissertation, Wageningen University). 

Input parameters:
=====================

For the input parameters, see the help of the software:

./AcroPoly.py --help, -h

***Example:***

./AcroPoly.py TetraMerged.bam TetraMerged.vcf -iter 2000 -mbq 19 --mmq 13 -maxIS 500 -mincov 4 -eps 1e-6 -d -v

Copyright and License:
=====================
Copyright (c) 2018, Ehsan Motazedi, Wageningen UR.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files \(the "Software"\), to deal in the Software without restriction, including without limitation the rights to use, copy, share and modify the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. This permission is only valid for academic and personal use, i.e. not for commercial purposes.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
