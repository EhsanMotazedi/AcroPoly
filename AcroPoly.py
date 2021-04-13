#!/usr/bin/env python
# Written by Ehsan Motazedi, Wageningen University & Research, 26-04-2018
# Last modified: 10-08-2018

import argparse
import bamprocess_regression
import os
import re
import sys
import tempfile

from os.path import isfile, join

def check_file():
	""" Checks if the given file path exists."""
	class Check_File(argparse.Action):
		def __call__(self, parser, args, values, option_string=None):
			if values is None:
				setattr(args, self.dest, '')
			elif False in [os.path.exists(_val) for _val in values]:
				raise argparse.ArgumentError(self, "no file called '{}' was found!".format(_val))
			else:
				setattr(args, self.dest, values)
	return Check_File

def check_positive(value):
	""" Checks if the passed argument is a positive integer."""
	try:
		ivalue = int(float(value))
	except ValueError:
		raise argparse.ArgumentTypeError("{} is not a positive integer value!".format(value))
	except TypeError:
		raise argparse.ArgumentTypeError("{} is not a positive integer value!".format(value))
	else:
		if ivalue <= 0:
			raise argparse.ArgumentTypeError("{} is not a positive integer value!".format(value))
		else:
			pass
	return ivalue

def check_positive_float(value):
	""" Checks if the passed argument is a positive real number."""
	fvalue=None
	try:
		fvalue = float(value)
	except ValueError:
		raise argparse.ArgumentTypeError("{} is not a positive real value!".format(value))
	except TypeError:
		raise argparse.ArgumentTypeError("{} is not a positive real value!".format(value))
	else:
		if fvalue <= 0:
			raise argparse.ArgumentTypeError("{} is not a positive real value!".format(value))
		else:
			pass
	return fvalue

if __name__ == "__main__":
	parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='AcroPoly method to estimate multi-allelic haplotype dosage scores from a set of reads (individual \
	or multi-sample bamfile) and a variant file (an individual or multi-sample VCF file) accross a region. AcroPoly scans the SNPs on the target region with a given \
	window size and window shift (move) size and calculates the dosage scores of multi-allelic haplotypes for each window from the sequence reads using the \
	Expectation Maximisation (EM) algorithm. The results will be reported separately for each window and optionally, the most likely phasing can be reported for each window \
	instead of the dosage scores for the possible haplotypes.')
	parser.add_argument('bam', metavar='bam file', action=check_file(), type=str, nargs=1, help='the bam file containing single or \
	paired-end reads of one or a collection of individuals.')
	parser.add_argument('vcf', metavar='vcf file', action=check_file(), type=str, nargs=1, help='the vcf file containing the variants \
	in the reads compared to some reference sequence.')
	parser.add_argument('output', metavar='output file', action='store', type=str, nargs='?', help='the base name of the output files to report haplotype markers and \
	dosages, or estimated phasings if -d, --deterministic is set. All windows will be written to the standard output if no output name is specified.')
	parser.add_argument('-l','--length', dest='windowsize', action='store', type=int, nargs=1, metavar="int+",
	help='the length of the windows (in number of SNPs) for defining multi-allelic haplotype markers.',
	default=3)
	parser.add_argument('-s','--shift', dest='shiftsize', action='store', type=check_positive, nargs=1, metavar="int+",
	help='the step size (in number of SNPs) to each time shift the current haplotype window to the right (from 5''\' to 3''\' end) for scanning the target region.',
	default=1)
	parser.add_argument('-cov','--cov','-mincov','--mincov', dest='min_cov', action='store', type=check_positive, nargs=1, metavar="int+",
	help='the minimum sequencing coverage required to estimate dosages for the possible haplotype alleles in each window.',
	default=4)
	parser.add_argument('-iter','--iter','-itermax','--itermax', dest='itermax', action='store', type=check_positive, nargs=1, metavar="int+",
	help='the maximum number of iterations allowed for the EM algorithm.',
	default=1000)
	parser.add_argument('-rel', '-relative', '--rel', '--relative', dest = 'use_relative', action='store_true',
	help='use relative convergence tolerance for the EM algorithm.',
	default=False)
	process = parser.add_mutually_exclusive_group(required=False)
	process.add_argument('-expression', '--expression', dest = 'expression', action='store_true',
	help='set if the reads belong to expression, i.e. NOT genomic, data. The reported dosages will reflect the normalized count of each allele for \
	expression data. NOT to be used with -d, --deterministic!', default=False)
	process.add_argument('-d', '--deterministic', dest = 'deterministic', action='store_true',
	help='report the most likely phasing estimate for each window in a conventional deterministic format as in SDhaP, HapTree, .... Not meaningful with expression data and \
	therefore incompatible with -expression, --expression.', default=False)
	parser.add_argument('-eps','--eps', dest='eps', action='store', type=check_positive_float, nargs=1, metavar="Real+",
	help='the absolute tolerance of iteration convergence for the EM algorithm.',
	default=1e-06)
	parser.add_argument('-epsrel','--epsrel', dest='eps_rel', action='store', type=check_positive_float, nargs=1, metavar="Real+",
	help='the relative tolerance of iteration convergence for the EM algorithm, used only when relative convergence is chosen by -rel, --relative.',
	default=1e-04)
	parser.add_argument('-mmq','--mmq', dest='mmq', action='store', type=check_positive, nargs=1, metavar="int+",
	help='the minimum mapping quality of a read to be considred for getting haplotype allele dosages.', default=20)
	parser.add_argument('-mbq','--mbq', dest='mbq', action='store', type=check_positive, nargs=1, metavar="int+",
	help='the minimum base-calling quality to consider a base within a read.',
	default=13)
	parser.add_argument('-maxIS','--maxIS', dest='maxIS', action='store', type=check_positive, nargs=1, metavar="int+",
	help='maximum insert-size for a paired-end read to be considered as a single fragment for phasing.',
	default=3000)
	parser.add_argument('-v', '--verbose', dest = 'verbose', action='store_true', help='report details of the EM runs.', 
        default=False)
	try:
	    if not len(sys.argv)>1:
		    parser.parse_args('-h'.split())
	    elif ('-h' in sys.argv) or ('--help' in sys.argv):
		    parser.parse_args('-h'.split())
	    else:
		    args = vars(parser.parse_args()) 
	except SystemExit:
	    raise
	for key, val in args.items(): #in Unix: cat ./AcroPoly.py |grep nargs=1|awk 'BEGIN{FS="dest"}{print $2}'|awk 'BEGIN{FS=","}{print $1}'|paste -s -d" "|sed "s/=/,/g"
		if key in ('bam', 'vcf', 'windowsize' ,'shiftsize' ,'min_cov' ,'max_haps' ,'min_dos' ,'itermax' ,'eps' ,'eps_rel' ,'mmq' ,'mbq' ,'maxIS'):
			try:
				exec(key + '=val[0]')
			except TypeError:
				exec(key + '=val')
		else:
			exec(key + '=val')
	if windowsize<2:
		parser.error('the window size (-l, --length) must be at least 2!')
	if verbose and use_relative:
		garbage = sys.stderr.write("Relative convergence tolerance applied to EM with epsilon={0:.2e}...\n".format(eps_rel))
	elif verbose:
		garbage = sys.stderr.write("Absolute convergence tolerance applied to EM with epsilon={0:.2e}...\n".format(eps))
	if expression:
		garbage = sys.stderr.write("The reads are assumed to be transcriptomic!\n".format(eps))
	else:
		garbage = sys.stderr.write("The reads are assumed to be genomic!\n".format(eps))
	if deterministic:
		sys.stderr.write("WARNING: deterministic haplotypes will be reported by rounding fractional dosages of each window!\n") 
	if not output:
		tmp=True
		output = tempfile.mkdtemp()+os.sep+'tmp'
		base_dir = (os.sep).join(output.split(os.sep)[0:-1])
	else:
		tmp=False	
	qoffset = 33
	if deterministic: # report a phasing estimate for each scanning window
		outhaps = bamprocess_regression.get_haplos(bam, vcf, windowsize, shiftsize, eps, eps_rel, itermax, use_relative, maxIS, mbq, mmq, qoffset, min_cov, expression, deterministic, verbose)
		if tmp:
			for _sample, _haps  in outhaps.items():
				sys.stdout.write('**************'+_sample+':\n'+'\n'.join(_haps)+'\n')
			os.rmdir(base_dir)
		else:
			for _sample, _haps  in outhaps.items():
				_file = output + '_' + _sample
				with open(_file, 'w') as _fh:
					_fh.write('\n'.join(_haps)+'\n')
	else: #report the estimated dosages for each possible haplotype within each scanning window
		bamprocess_regression.scan_region(bam, vcf, windowsize, shiftsize, output, eps, eps_rel, itermax, use_relative,  maxIS, mbq, mmq, qoffset, min_cov, expression, verbose)
		if tmp:
			pattern = re.compile(base_dir+os.sep+'tmp_window_[0-9]*_from_position_[0-9]*_to_position_[0-9]*.dat')
			for _entity in os.listdir(base_dir):
				if isfile(join(base_dir, _entity)) and pattern.match(join(base_dir, _entity)):
					winnum = _entity.split('tmp_window_')[-1].split('.dat')[0]
					with open(join(base_dir, _entity)) as _fh:
						sys.stdout.write("Alleles and dosages in window {0:s}:\n".format(' '.join(winnum.split('_'))))
						for _line in _fh:
							sys.stdout.write(_line)
					os.unlink(join(base_dir, _entity))
			os.rmdir(base_dir)
	try:
            os.close(sys.stdout.fileno())
            #sys.stdout.close()
	except:
	    pass
	try:
	    #sys.stderr.close()
            os.close(sys.stderr.fileno())
	except:
	    pass
