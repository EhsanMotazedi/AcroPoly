#!/usr/bin/env perl
#Written by Ehsan Motazedi, Wageningen University & Research
#10-10-2018
use strict;
use warnings;

my $outfile="/mnt/scratch/motaz001/Alstroemria/NewVCF_Alstroemeria.vcf";
open(my $fh, "<", "/mnt/scratch/motaz001/Alstroemria/Merged6000_All-samples.vcf") or die "Can't open < Merged6000_All-samples.vcf: $!";

open(my $fh2, ">", $outfile) or die "Can't open $outfile: $!";

while (my $line = <$fh>){
	chomp $line;
	if ((substr $line, 0, 1) eq "#" and (substr $line, 1, 1) ne "#"){
		my @l = split /\t/, $line;
		print $fh2 (join "\t", @l[0,1,3,4,9..$#l])."\n"
	} elsif ((substr $line, 0, 1) ne "#") {
		my @recs=split /\t/, $line;
		my %alleles = map {$_ => 1} (split /,/, $recs[4]); #filter out multi-allelic variants
		my $n_alt_alleles = keys %alleles;
		if ($n_alt_alleles==1){
			print $fh2 (join "\t", @recs[0,1,3,4])."\t";
			for my $i (9..$#recs){
				my @geno = split /\/|\|/, ((split /:/, $recs[$i])[0]);
				my %alleles = map {$_ => 1} @geno; #set missing values as missing
				my $dosage;
				if (exists($alleles{"."})) {
					$dosage = "."
				} else { 
					$dosage = eval join "+", @geno # get the dosages for non-missing bi-allelic variants 
				}
				if ($i lt $#recs) {
					print $fh2 "$dosage\t"
				} else {
					print $fh2 "$dosage\n"
				}
			}
		}
	}
}

close $fh;
close $fh2;
