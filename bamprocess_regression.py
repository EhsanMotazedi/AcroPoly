# Functions used by AcroPoly algorithm to process the bam file, scan the genomic region, get EM dosages and report them in the desired format.
# Written by Ehsan Motazedi, Wageningen UR, 21-01-2018.
# Last mofified: 08-09-2018.

import copy
import itertools
import math
import multiprocessing as mpthread
import numpy as np
import pysam
import re
import subprocess
import sys
import threading
from collections import OrderedDict
from cStringIO import StringIO
from genotypes import getAllelesPop, getGenotypesPop
from haplotypes import Haplotypes
from math import log
from reads import Read
#from scipy.misc import factorial as fact

MAXCORES = mpthread.cpu_count()-1 # Max number of available cores
NCORES = 8 # desired number of cores to use
logfactorial_table = {n: sum(log(_n) for _n in range(1,n+1)) for n in range(0,1001)}

Global_Sema = None # Global Semaphore not initialized yet 

class BlockException(Exception):
	def __init__(self, value):
		super(BlockException, self).__init__(value)
		self.args = (value,)
	def __str__(self):
		return "{}".format(':'.join(str(_arg) for _arg in self.args))
	def __repr__(self):
		return self.args

class Capturing(list):
    """class defined to temporarily redirect stdout and capture what is written in a list. To be used in 'with' statement block."""
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self
    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        del self._stringio    # free up some memory
	self._stdout.flush()
        sys.stdout = self._stdout

def GetSemaLock(useglobal=True):
	""" Return the semaphore and lock objects needed for multi threading."""
	global NCORES
	global MAXCORES
	global Global_Sema
	if min(NCORES, MAXCORES)>1:
		if useglobal: # use a global semaphore
			if Global_Sema is None:
				Global_Sema = mpthread.BoundedSemaphore(min(NCORES, MAXCORES))
		else: # use a local semaphore at each call to BranchPop
			sema = mpthread.BoundedSemaphore(min(NCORES, MAXCORES))
		lock = mpthread.Lock()
		Parallel = True
	else:   # If parallel execution is not possible (as just one core could be used), concurrent execution is performed using threads.
		if useglobal: 
			if Global_Sema is None:
				Global_Sema = threading.BoundedSemaphore(NCORES) 
		else:
			sema = threading.BoundedSemaphore(NCORES)
		lock = threading.Lock()
		Parallel = False
	if useglobal:
		return Global_Sema, lock , Parallel 
	else:
		return sema, lock , Parallel 

def thread_func(sema, lock, q, func, *args, **kwargs):
	""" thread to call func with args and put the return value in q. """
	_locked = False
	try:
		b = func(*args, **kwargs)
		_locked = lock.acquire()
		q.append(b)
	except:
		raise
	finally:
		if _locked:
			lock.release()	
		sema.release()

class Infix:
	def __init__(self, function):
		self.function = function
	def __ror__(self, other):
		return Infix(lambda x, self=self, other=other: self.function(other, x))
	def __or__(self, other):
		return self.function(other)
	def __rlshift__(self, other):
		return Infix(lambda x, self=self, other=other: self.function(other, x))
	def __rshift__(self, other):
		return self.function(other)
	def __call__(self, value1, value2):
		return self.function(value1, value2)

mulx = Infix(lambda X, Y: [[sum(a*b for a,b in zip(X_row,Y_col)) for Y_col in zip(*Y)] for X_row in X]) # matrix multiplication |x|

def adjust_seq(read):
	""" Adjust the read sequence according to the mapped positions, i.e. get read of the insertions and clipped bases."""
	cig = list(_x for _x in re.split('[0-9]{1,}', read.cigarstring) if _x)
	cign = list(_x for _x in re.split('[^0-9]', read.cigarstring) if _x)
	cig = list(cig[_n]*int(cign[_n]) for _n in range(0, len(cig)))
	cig = ''.join(_x for _x in cig if 'D' not in _x) # deleted nucleotides from the reference are not present in the read sequence
	adj_seq=[]
	for _n, _x in enumerate(read.seq):
		if cig[_n]=='M': # Allow only match/mismatch of the read nucleotides, i.e. no clipping, no insertion.
			adj_seq.append(_x)
	return ''.join(adj_seq)

def frag_gen(varpos, allelelst, genolst, coordinates, nucleotides, qscores):
	""" Generate SNP-fragments from (paired-end) reads.""" 
	if (not coordinates) or (not nucleotides):
		return Read(), {}
	var_codes = []
	var_num = []
	var_q = []
	for _n, _x in enumerate(varpos):
		if _x < coordinates[0]:
			continue
		if _x > coordinates[-1]:
			break
		if _x in coordinates:
			if set(genolst[_n].GetGenes()).intersection(set(['.','-'])): # Throw away missing genotypes or genotypes with one or more missing allele(s)
				continue
			#if len(set(genolst[_n].GetGenes()))<2: # Do not include in the SNP-fragments belonging to a population member its homozygous alleles
			#	continue
			try:
				var_codes.append(allelelst[_n][2][nucleotides[coordinates.index(_x)]])
				var_num.append(allelelst[_n][0])
				var_q.append(qscores[coordinates.index(_x)])
			except KeyError: # if the called nucleotide is wrong, i.e. does not exist in VCF alleles
				pass
	return Read({_x:str(_y) for _x, _y in zip(var_num, var_codes)}), {_x:str(_y) for _x, _y in zip(var_num, var_q)} # return the reads {SNP number: allele} and the associated quality scores {SNP number: Qscore}

class InputError(Exception):
	""" Handle invalid input specifications.""" 
	def __init__(self, msg, *args):
		super(InputError, self).__init__(args)
		self.msg = msg
	def __str__(self):
		return "InputError: {}\n".format(self.msg)
	def __repr__(self):
		return (self.msg,)+self.args+('\n',)

def get_expectation(haplodict, M):
	"""calculate the expected count of each possible haplotype within a window, using the current haplotype rate
	estimates and the compatibility matrix M. Part of the EM-algorithm to get the expectation of the unknown counts
	in  the likelihood function for the reads: P(R|k_1,...,k_N,u_1,...,u_N) ~ Pois(k_1|lu_1)...Pois(k_N|lu_N) (with 
	l the length of the haplotyes in SNPs). Having M, the matrix showing copatiblity of reads (each read one row)
	with the haplotype (each hapoltype one column, values of course 0 and 1), it is easy to calculate k_i assuming 
	a multinomial model for the assignment of each read. U is the vector of haplotype specific Poisson rates. Expected number of 
	each haplotype h_j,i.e. k_j, will be: sum over the reads (U_j.I(r_i|H_j) divided by the sum of U_k for h_k 
	compatible with r_i)."""
	c = len(M) # total number of reads overlapping with the current region for the current sample
	N = len(haplodict) # number of possible haplotypes for the region
	count_dict  = OrderedDict()
	haplo_freq = [_u for _u in haplodict.values()]
	frag_prob = np.transpose(np.dot(np.diag(haplo_freq), np.transpose(M)))
	#frag_prob_norm = M |mulx| [[_f] for _f in haplo_freq] # sum of U_k for h_k compatible with r_i
	frag_prob_norm = np.dot(M, [[_f] for _f in haplo_freq]) # sum of U_k for h_k compatible with r_i
	frag_prob = [[frag_prob[_i][_j] / float(frag_prob_norm[_i][0]) for _j in range(0, N)] for _i in range (0, c)]
	expectations = np.dot([[1 for _i in range(0,c)]], frag_prob)
	for _j, _h in enumerate(haplodict.keys()):
		count_dict[_h] = expectations[0][_j]
	return count_dict

def get_frags(bamfile, vcffile, maxIS=3000, mbq=13, mmq=20, qoffset=33):
	""" 
Make the SNP-fragment list from a multi-sample bam file and its corresponding VCF file  with input options: mmq, mbq, maxIS and qoffset.
Modified from the original version for read-based regression. Specifically, read length does NOT need to be at least 2 here, and homozygous sites of 
individual are also included.
mmq : minimum read mapping quality to consider a read for phasing, default 20\n
qvoffset <33/64> : quality value offset, 33/64 depending on how quality values were encoded, default is 33\n
mbq  : minimum base quality to consider a base for haplotype fragment, default 13\n
maxIS : maximum insert-size for a paired-end read to be considered as a single fragment for phasing, default 3000.\n 
	"""
	try:
		all_reads = pysam.Samfile(bamfile, 'rb')
	except IOError:
		raise InputError('The input BAM file was not found!')
	ReadHeader = subprocess.Popen(["samtools","view","-H", bamfile], shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
	header, err_header = ReadHeader.communicate()
	if ReadHeader.returncode!=0:
		raise InputError('Failed to read the header from the bam file! Original error message:\n'+err_header)
	if isinstance(header, bytes):
		header=bytes.decode(header)
	else:
		pass
	RGIDs, SMnames = [], []
	for _headerline in header.splitlines(): # pasre the header of the bam file to extract the ID and SM fields of each Read Group
		if _headerline[0:3]=='@RG':
			RGID_added, SM_added = False, False
			for _n, _field in enumerate(_headerline.split()):
				if 'ID' in _field:
					if not RGID_added:
						RGIDs.append(''.join(_headerline.split()[_n].split(':')[1:])) # add read group ID
						RGID_added = True
					else:
						raise InputError('Double ID fields detected in @RG header line!')
				elif 'SM' in _field:
					if not SM_added:
						SMnames.append(''.join(_headerline.split()[_n].split(':')[1:])) # add the sample name associated with the read group ID
						SM_added = True
					else:
						raise InputError('Double SM fields detected in @RG header line!')
			if SM_added and RGID_added:
				pass
			elif SM_added:
				raise InputError('ID field missing in @RG header line!')
			elif RGID_added:
				raise InputError('SM field missing in @RG header line!')
			else:
				raise InputError('ID and SM fields missing in @RG header line!')
	if len(RGIDs)!=len(set(RGIDs)):
		raise InputError('Duplicate read group IDs detected in the bam header!')
	GroupedReadsWithID = [[] for _id in RGIDs] # separate reads belonging to each Read Group
	for _read in all_reads:
		GroupedReadsWithID[RGIDs.index(dict(_read.get_tags())['RG'])].append(_read)
	GroupedReads, GroupedSM = [], [] # combine the reads with different RGID but the same SM as they are assumed to belong to the same sample
	for _SMname, _ReadGroup in zip(SMnames, GroupedReadsWithID):
		if _SMname not in GroupedSM:
			GroupedReads.append(_ReadGroup)
			GroupedSM.append(_SMname)
		else:
			GroupedReads[GroupedSM.index(_SMname)]+=_ReadGroup
	del GroupedReadsWithID
	try:
		genolst = getGenotypesPop(vcffile, GroupedSM)
		allelelst = getAllelesPop(vcffile, GroupedSM)
	except IOError:
		raise InputError('The VCF file was not found!')
	except:
		raise
	Frag_lsts, Q_lsts = [], []
	varpos = list(_x[1] for _x in allelelst)
	for _group, reads in enumerate(GroupedReads):
		frag_lst, q_lst = [], []
		reads = sorted(reads, key= lambda _x: (_x.qname, _x.flag & 0x900)) # sort the alignments using their names, with the primary alignments being placed first.
		_rnum = 0
		rNUM = len(reads)
		while _rnum < rNUM: # scan through the alignments to find the pairs/singles 
			break_mate = False
			is_proper_pair = False
			read = copy.deepcopy(reads[_rnum])
			if read.is_unmapped or read.is_duplicate or (read.flag & 0x900): # throw away unmapped reads, duplicates and secondary/supplementary alignments
				_rnum+=1
				continue
			try:
				if read.qname == reads[_rnum+1].qname: # means the read is paired to a mate or has multiple/supplemenatry alignments
					if reads[_rnum+1].is_unmapped or reads[_rnum+1].is_duplicate or (reads[_rnum+1].flag & 0x900): # if the next read is unmapped, a duplicate or not primarym: skip it
						pass
					else:
						is_proper_pair = True	# means the read is paired to a proper mate
						mate = copy.deepcopy(reads[_rnum+1])
					_rnum+=2
				else: # means the read is single
					_rnum+=1
			except IndexError: # could occur for the last read in the alignments' list
				_rnum+=1
			if is_proper_pair:
				if (max(mate.positions+read.positions)-min(mate.positions+read.positions)+1)>maxIS: # Check the maximum insert-size to consider the mates as a single fragment
					break_mate = True
				if read.mapping_quality >= mmq:
					try:
						coordinates, nucleotides, quals = list(zip(*[(_x, _y, _z) for _x, _y, _z in zip(read.positions, adjust_seq(read).upper(), list(ord(_x)-qoffset for _x in read.qual)) if _z>=mbq])) 
					except ValueError as e:
						if e.args[0][0:len("need more than 0 values to unpack")]=="need more than 0 values to unpack" or e.args[0][0:len("not enough values to unpack")]=="not enough values to unpack":
							coordinates, nucleotides, quals= [(), (), ()] 
						else:
							raise
				else:
					coordinates, nucleotides, quals = [(), (), ()] 
				if mate.mapping_quality >= mmq:
					try: 
						coordinates_mate, nucleotides_mate, quals_mate = list(zip(*[(_x, _y, _z) for _x, _y, _z in zip(mate.positions, adjust_seq(mate).upper(), list(ord(_x)-qoffset for _x in mate.qual)) if _z>=mbq]))
					except ValueError as e:
						if e.args[0][0:len("need more than 0 values to unpack")]=="need more than 0 values to unpack" or e.args[0][0:len("not enough values to unpack")]=="not enough values to unpack":
							coordinates_mate, nucleotides_mate, quals_mate = [(), (), ()] 
						else:
							raise
				else:
					coordinates_mate, nucleotides_mate, quals_mate = [(), (), ()] 
				if break_mate:
					pass
				else: # merge the sub-reads if the insert-size is less than maxIS
					try: 
						coordinates, nucleotides, quals = list(zip(*sorted(zip(coordinates+coordinates_mate, nucleotides + nucleotides_mate, quals+quals_mate), key = lambda x: x[0])))
					except ValueError as e:
						if e.args[0][0:len("need more than 0 values to unpack")]=="need more than 0 values to unpack" or e.args[0][0:len("not enough values to unpack")]=="not enough values to unpack":
							coordinates, nucleotides, quals = [(), (), ()] 
						else:
							raise
			else: 
				break_mate = True
				if read.mapping_quality >= mmq:
					try:
						coordinates, nucleotides, quals = list(zip(*[(_x, _y, _z) for _x, _y, _z in zip(read.positions, adjust_seq(read).upper(), list(ord(_x)-qoffset for _x in read.qual)) if _z>=mbq]))
					except ValueError as e:
						if e.args[0][0:len("need more than 0 values to unpack")]=="need more than 0 values to unpack" or e.args[0][0:len("not enough values to unpack")]=="not enough values to unpack":
							coordinates, nucleotides, quals = [(), (), ()] 
						else:
							raise
				else:
					coordinates, nucleotides, quals = [(), (), ()] 
				coordinates_mate, nucleotides_mate, quals_mate = [(), (), ()] 
			if break_mate:
				pass
			else:
				unique_q = []
				unique_c = []
				unique_n = []
				for _n, _c in enumerate(coordinates): # remove the duplicates from overlapping positions 
					try:
						if unique_c[-1]!=_c:
							unique_c.append(_c)
							unique_n.append(nucleotides[_n])
							unique_q.append(quals[_n])
						elif unique_n[-1]==nucleotides[_n]:
							unique_q[-1] = min(126-qoffset, unique_q[-1]+quals[_n])
						else: # if the called nucleotides differ at overlapping sites, use the one with the highest Phred score and adjust the Phred score.
							if quals[_n]>unique_q[-1]:
								_new_q_score = round(-10*log(1-10**(-unique_q[-1]/10)*(1-10**(-quals[_n]/10)), 10), 5) # Q=-10log(p,10)
								if _new_q_score >= mbq:
									unique_n[-1] = nucleotides[_n]
									unique_q[-1] = _new_q_score
								else:
									del(unique_c[-1], unique_n[-1], unique_q[-1])
							else:
								_new_q_score = round(-10*log(1-(1-10**(-unique_q[-1]/10))*10**(-quals[_n]/10), 10), 5)
								if _new_q_score >= mbq:
									unique_q[-1] = _new_q_score
								else:
									del(unique_c[-1], unique_n[-1], unique_q[-1])
					except IndexError:
						unique_c.append(_c)
						unique_n.append(nucleotides[_n])
						unique_q.append(quals[_n])
				coordinates, nucleotides, quals = [unique_c, unique_n, unique_q]
			coordinates = list(_x+1 for _x in coordinates) # Convert the zero-based BAM coordinates to 1-based, as the coordinates are 1-based in the VCF (like the SAM format).
			new_frag, new_q = frag_gen(varpos, allelelst, genolst[_group], coordinates, nucleotides, quals)
			frag_lst.append(new_frag)
			q_lst.append(new_q)
			if break_mate:
				coordinates_mate = list(_x+1 for _x in coordinates_mate)
				new_frag_mate, new_q_mate = frag_gen(varpos, allelelst, genolst[_group], coordinates_mate, nucleotides_mate, quals_mate)
				frag_lst.append(new_frag_mate)
				q_lst.append(new_q_mate)
		try:
			frag_lst, q_lst = [_lst for _lst in zip(*[(_x, _y) for _x, _y in zip(frag_lst, q_lst) if not _x.isNULL()])]
		except ValueError as e:
			if e.args[0][0:len("need more than 0 val:ques to unpack")]=="need more than 0 values to unpack" or e.args[0][0:len("not enough values to unpack")]=="not enough values to unpack":
				frag_lst, q_lst = [], []
		Frag_lsts.append(frag_lst)
		Q_lsts.append(q_lst)
	return Frag_lsts, Q_lsts, GroupedSM

def get_haplos(bamfile, vcffile, windowsize = 3, shiftsize = 1, eps=1e-06, eps_rel=1e-04, maxit=5000, use_relative=False, maxIS=500, mbq=10, mmq=9, qoffset=33, min_cov = 4, expression=True, deterministic=False, verbose=False):
	"""scan the SNPs given by vcffile, using a window of size winodwsize and window shifts of size shiftsize, and estimate haplotype frequencies for each window. 
        Finally, bridge the windows by starting from the leftmost window and stacking consecutive windows upon each other."""
	output = []
	sema, lock , Parallel = GetSemaLock(True)
	threads = []
	if Parallel:
		manager = mpthread.Manager()
	Frag_lsts, Q_lsts, Samples = get_frags(bamfile, vcffile, maxIS, mbq, mmq, qoffset)
	if not Samples:
		raise InputError('No sample name included in the bam header!')
	alleleSet = [str(_allele) for _allele in set(_x for _y in getAllelesPop(vcffile, Samples) for _x in _y[2].values())] # obtain the set of alleles observed at all variation sites
	varpos = [] # obtain the variant positions from the vcf file
	ploidy = [0 for _sample in Samples]
	contigs = []
	with open(vcffile, 'rU') as _vcfh:
		for _line in _vcfh:
			if '#'!=_line.strip()[0]:
				varpos.append(_line.strip().split()[1])
				contigs.append(_line.strip().split()[0])
				if not all(ploidy): # determine the ploidy levels of each sample
					_genotype_fields = [_x.split(':')[0] for _x in _line.strip().split()[9:]]
					if not all(_x=='.' for _x in _genotype_fields):
						for _n, _x in enumerate(_genotype_fields):
							if not ploidy[_n]:
								ploidy[_n] = len(_x.split('/')) if _x!='.' else ploidy[_n]
	if len(set(contigs))>1:
		raise BlockException('All of the variants must be located on the same contig for phasing! Check the VCF file!')
	contig = contigs[0]
	if not all(ploidy):
		raise BlockException('Ploidy levels could not be detected for some of the samples! Check the VCF file!')
	l = len(varpos)
	if Parallel:
		_window_haps = [manager.list() for _start in range(0, l-windowsize+shiftsize, shiftsize)] # shiftsize must be <= windowsize-2 for informative overlapping between the windows
	else:
		_window_haps = [[] for _start in range(0, l-windowsize+shiftsize, shiftsize)]   # ordinary lists suffice for multiple threads belonging to the same process
	_index = -1
	_start = -1*shiftsize
	while _start<(l-windowsize):
		_index+=1 
		_start+=shiftsize
		if Parallel: 
			t = mpthread.Process(target=thread_func, name = "Haplotype window {0:d}".format(_index+1),
			args = (sema, lock, _window_haps[_index], make_data_for_window, 
			varpos[_start:_start+windowsize],Frag_lsts, Q_lsts, Samples, vcffile, None, eps, eps_rel, maxit, use_relative),  
			kwargs = {'queue':True, 'verbose':verbose, 'maxIS':maxIS, 'mmq':mmq, 'mbq':mbq, 'qoffset':qoffset, 'expression':expression, 'varpos':varpos,
			'alleles':alleleSet, 'min_cov':min_cov, 'Parallel':Parallel})
		else:
			t = threading.Thread(target=thread_func, name = "Haplotype window {0:d}".format(_index+1), 
			args = (sema, lock, _window_haps[_index], make_data_for_window, 
			varpos[_start:_start+windowsize],Frag_lsts, Q_lsts, Samples, vcffile, None, eps, eps_rel, maxit, use_relative),  
			kwargs = {'queue':True, 'verbose':verbose, 'maxIS':maxIS, 'mmq':mmq, 'mbq':mbq, 'qoffset':qoffset, 'expression':expression, 'varpos':varpos,
			'alleles':alleleSet, 'min_cov':min_cov, 'Parallel':Parallel})
		sema.acquire()
		threads.append(t)
		t.start()
	for _thread in threads:
		_thread.join()
	if Parallel:
		window_haps = [list(_lprox)[0] if len(_lprox) else [] for _lprox in _window_haps]  # empty dictionaries reported when coverage is too low
	else:
		window_haps = [_l[0] if len(_l) else [] for _l in _window_haps] # empty dictionaries reported when coverage is too low
	Haplotype_Blocks = dict()
	for _sample in Samples:
		_window_dics = {}
		for _n, _win in enumerate(window_haps):
			_window_dics[_n] = _win[_sample]
                Haplotype_Blocks[_sample] = bridge_windows(_window_dics, shiftsize, windowsize, alleleSet, deterministic) # now combine the atomic windows one by one
	_sample_num = -1
	for _sample in Samples:
		_sample_num+=1
		if not deterministic:
			with Capturing() as _output:
				sys.stdout.write('Haplotype estimates for {}:\n'.format(_sample))
				for _block_num in range(0, len(Haplotype_Blocks[_sample])):
					sys.stdout.write('--->Block {0:d} on {1:s}:\n'.format(_block_num+1, contig))
					_n = 0
					for _H in Haplotype_Blocks[_sample][_block_num]:
						_n+=1
						sys.stdout.write('\tH_{0:d}: {1:s}\n\tStart: {3:d} ({5:s} bp), Stop: {4:d} ({6:s} bp), Prob H_{0:d}: {2:.5f}\n'.format(_n, ''.join(_H.GetVS()[0]), _H.GetRL(), _H.GetStart()+1, _H.GetStop()+1, varpos[_H.GetStart()], varpos[_H.GetStop()]))
				sys.stdout.write('*******************\n')
				sys.stdout.flush()
		else:
			with Capturing() as _output:
				for _block_num in range(0, len(Haplotype_Blocks[_sample])):
					columns = [[], []]+ [[] for _x in range(0, ploidy[_sample_num])]
					sys.stdout.write('BLOCK_{0:d}_from_{1:s}bp_to_{2:s}bp\t{5:d}_SNPs\t{3:d}\t{4:d}\n'.format(_block_num+1, varpos[Haplotype_Blocks[_sample][_block_num][0].GetStart()],  varpos[Haplotype_Blocks[_sample][_block_num][0].GetStop()], Haplotype_Blocks[_sample][_block_num][0].GetStart()+1, Haplotype_Blocks[_sample][_block_num][0].GetStop()+1, len(Haplotype_Blocks[_sample][_block_num][0].GetVS()[0]))) # write the header
					#columns[0] = ['.' for _pos in range(0, len(Haplotype_Blocks[_sample][_block_num][0].GetVS()[0]))] # contig column in the individual solution
					columns[0] = [contig for _pos in range(0, len(Haplotype_Blocks[_sample][_block_num][0].GetVS()[0]))] # contig column in the individual solution
					columns[1] = [varpos[_pos] for _pos in range(Haplotype_Blocks[_sample][_block_num][0].GetStart(), Haplotype_Blocks[_sample][_block_num][0].GetStop()+1)] # variant position column in the individual solution
					for _h in range(2, ploidy[_sample_num]+2):
						columns[_h] = ['-' for _pos in range(Haplotype_Blocks[_sample][_block_num][0].GetStart(), Haplotype_Blocks[_sample][_block_num][0].GetStop()+1)] # initially undetermined haplotypes
					_h = 2
					try:
						for _H in Haplotype_Blocks[_sample][_block_num]:
							for _hh in range(0, int(round(_H.GetRL()*ploidy[_sample_num]))): # convert fraction dosages to integer dosages
								try:
									columns[_h+_hh] = _H.GetVS()[0] # write a homologue d times with d being its integer dosage
								except IndexError:
									sys.stderr.write("WARNING: Dubious rounding of fractional dosages to integer! Maybe better to use fractional dosages instead of the deterministic!\n")
									break
							_h += int(round(_H.GetRL()*ploidy[_sample_num]))
					except ValueError:
						sys.stderr.write("WARNING: Undefined haplotype probabilities encountered in block {0:d}, sample {1:s}!\n".format(_block_num, _sample))
					if any(_h==['.' for _pos in range(0, len(Haplotype_Blocks[_sample][_block_num][0].GetVS()[0]))] for _h in columns[2:]):
						sys.stderr.write("WARNING: Some of the haplotypes could not be reliably determined and are therefore missing! Maybe better to use fractional dosages instead of the deterministic!\n")
					for _row in range(0, len(columns[0])):
						sys.stdout.write('\t'.join(str(_col[_row]) for _col in columns)+'\n')
			output.append(_output)
	if not deterministic:
		return output # return the log report with global haplotype markers and their dosages
	else:
		return {_sample:_output for _sample, _output in zip(Samples, output)} # return a dictionary with sample names as keys and determined haplotypes as values

def bridge_windows(hap_window_dics, shiftsize=1, winsize=3, alleleset=('0','1'), deterministic=False):
        """combine haplotype windows, starting from the leftmost window, using the probabilities of the possible haplotypes within each window.
	If deterministic haplotypes are required, report a phasing estimate for each haplotype block (window) obtained by calculating the expected dosages from the probabilities."""
	Haplo_Blocks = [] # homologues and probs estimated for each haploblock
	starts = [] #start SNP position for each haploblock 
	stops = [] #stop SNP position for each haploblock
	current_found = False
	_start_win_num = -1
	while not current_found:
		_start_win_num+=1
		try:
			if not any(math.isnan(_prob) for _prob in hap_window_dics[_start_win_num].values()):
				current_window = hap_window_dics[_start_win_num]
				current_found = True
		except IndexError:
			sys.stderr.write("ERROR: No haplotype windows with proper probabilites given to bridge!\n")
			return Haplo_Blocks
	next_window_number = _start_win_num
	current_added = False # has the current window been added to the haplotype blocks?
	len_hap = len(hap_window_dics.keys())
	while next_window_number < (len_hap-1): # for test, just stack windows
		next_window_number += 1
		next_found = False
		while (not next_found) and (next_window_number < len_hap):
			if not (any(math.isnan(_prob) for _prob in hap_window_dics[next_window_number].values())):
				#if (stops==[]) or (next_window_number*shiftsize > stops[-1]): # this results in non-overlapping windows
				if 1:
					next_found = True
					current_window = hap_window_dics[next_window_number] # start a new block with the next window
					starts.append(next_window_number*shiftsize) # start position of the new haplotype block stored
					stops.append(next_window_number*shiftsize+min(len(next(iter(current_window))), winsize)-1) # stop position of the current haplotype block stored
					Haplo_Blocks.append(current_window)
				else:
					next_window_number+=1
			else:
				next_window_number+=1
	Haplo_Blocks_new = []
	for _n, _block in enumerate(Haplo_Blocks):
		Haplo_Blocks_new.append([])
		for _H in _block.keys():
#			print(starts[_n], stops[_n], _H)
			Haplo_Blocks_new[-1].append(Haplotypes(starts[_n], stops[_n], _block[_H], 0, None, None, _H))
	return Haplo_Blocks_new

def scan_region(bamfile, vcffile, windowsize = 3, shiftsize=1, outprefix = 'scan', eps=1e-06, eps_rel=1e-04, maxit=5000, use_relative=False, maxIS=500, mbq=10, mmq=9, qoffset=33, min_cov = 4, expression = True, verbose=False):
	"""scan the SNPs given by vcffile, using a window of size winodwsize. The haplotypes (expression) dosages are estimated for the SNPs in each window from the bam file and 
	written to a dat files for that window. Windows will slide one or more nucleotides at each scanning step according to the shiftsize (l-windowsize+1 windows in total with 
	shiftsize = 1 and l//windowsize+1 windows with shiftsize = windowsize)."""
	output = []
	sema, lock , Parallel = GetSemaLock(True)
	threads = []
	if Parallel:
		manager = mpthread.Manager()
	Frag_lsts, Q_lsts, Samples = get_frags(bamfile, vcffile, maxIS, mbq, mmq, qoffset)
	if not Samples:
		raise InputError('No sample name included in the bam header!')
	alleleSet = [str(_allele) for _allele in set(_x for _y in getAllelesPop(vcffile, Samples) for _x in _y[2].values())] # obtain the set of alleles observed at all variation sites
	varpos = [] # obtain the variant positions from the vcf file
	with open(vcffile, 'rU') as _vcfh:
		for _line in _vcfh:
			if '#'!=_line.strip()[0]:
				varpos.append(_line.strip().split()[1])
	l = len(varpos)
	if Parallel:
		garbage = [manager.list() for _start in range(0, l-windowsize+shiftsize, shiftsize)]
	else:
		garbage = [[] for _start in range(0, l-windowsize+shiftsize, shiftsize)]
	_index = -1
	_start = -1*shiftsize
	while _start<(l-windowsize):
		_index+=1 
		_start+=shiftsize
		current_varpos = varpos[_start:_start+windowsize]
		if Parallel: 
			t = mpthread.Process(target=thread_func, name = "Haplotype window {0:d}".format(_index+1),
			args = (sema, lock, garbage[_index], make_data_for_window, 
			current_varpos,Frag_lsts, Q_lsts, Samples, vcffile, outprefix+'_window_'+str(_index+1)+'_from_position_'+str(current_varpos[0])+'_to_position_'+str(current_varpos[-1])+'.dat', eps, eps_rel, maxit, use_relative), 
			kwargs = {'maxIS':maxIS, 'mbq':mbq, 'mmq':mmq, 'qoffset':qoffset, 'alleles':alleleSet, 'min_cov':min_cov, 'expression':expression, 
			'varpos':varpos, 'verbose':verbose, 'Parallel':Parallel})
		else:
			t = threading.Thread(target=thread_func, name = "Haplotype window {0:d}".format(_index+1),
			args = (sema, lock, garbage[_index], make_data_for_window, 
			current_varpos,Frag_lsts, Q_lsts, Samples, vcffile, outprefix+'_window_'+str(_index+1)+'_from_position_'+str(current_varpos[0])+'_to_position_'+str(current_varpos[-1])+'.dat', eps, eps_rel, maxit, use_relative), 
			kwargs = {'maxIS':maxIS, 'mbq':mbq, 'mmq':mmq, 'qoffset':qoffset, 'alleles':alleleSet, 'min_cov':min_cov, 'expression':expression, 
			'varpos':varpos, 'verbose':verbose, 'Parallel':Parallel})
		sema.acquire()
		threads.append(t)
		t.start()
	for _thread in threads:
		_thread.join()

def make_data_for_window(window,Frag_lsts, Q_lsts, Samples, vcffile, outfile = None, eps = 1e-06, eps_rel = 1e-04, maxit=1000, use_relative= False, maxIS=500, mbq=10, mmq=9, qoffset=33, alleles = ('0','1'), min_cov = 4, expression = True, varpos = None, queue = False, verbose=False, Parallel=False):
	"""Gets SNP-fragments and then generates a data frame: sample names are the rows and the dosages of
	each windowed haplotype are stored in the columns. Window gives SNP positions to be included in the windowed haplotype 
	as a list. Assuming bi-allelic SNPs (0/1), 2**len(window) haplotypes will be possible the dosage of each is determined
	from the SNP-fragments. The results will be written to the given outfile or to standard output."""
	if not outfile:
		outhandel = sys.stdout
	else:
		outhandel = open(outfile, 'w')
	if not varpos:
		varpos = []
		with open(vcffile, 'rU') as _vcfh:
			for _line in _vcfh:
				if '#'!=_line.strip()[0]:
					varpos.append(_line.strip().split()[1])
	window = {varpos.index(str(_base)): _winpos for _winpos, _base in enumerate(window)}
	all_haplos = [_winhaplo for _winhaplo in itertools.product(*[alleles for _winpos in range(0,len(window))])] # all of the haplotypes possible in the window
	if not queue:
		garbage = outhandel.write("Sample"+'\t'+'\t'.join('H_'+''.join(_haplo) for _haplo in all_haplos)+'\n')  #write the header for the data frame file
	l = len(window) # length of the haplotypes (just equal to the size of the SNP scanning window)
	N = len(all_haplos)  # number of the haplotypes
	norm = 1.
	if queue:
		Outqueue = OrderedDict()
	else:
		Outqueue = None
	for _sn, _sample in enumerate(Samples): # initialize the EM-algorithm and set once and for all the compatibility matrix
		to_delete = [] # delete SNP-fragments that do not overlap with the current window
		coverage = len(Frag_lsts[_sn]) # total number of reads for sample _sn
		Matrix = [[0 for _j in range(0, N)] for _i in range(0,coverage)] # compatibility matrix for the reads (_i) vs haplotypes (_j) for sample _sn
		sample_haplo_dict = OrderedDict({_h:1./N for _h in all_haplos}) # Ordered dictionary to store the dosages of the window haplotypes for each sample
		for _i in range(0, len(Frag_lsts[_sn])):
			_SNP_frag = Frag_lsts[_sn][_i].GetDict() # get the dictionary corresponding to the SNP-fragment (throwing '-' alleles away)
			_overlap  =  set(window.keys()) & set(_SNP_frag.keys()) # set of SNP positions common between the SNP-fragment and the current window
			if not _overlap: # if no overlap with the window
				to_delete.append(_i)
			else:
				_j = -1
				for _h in sample_haplo_dict.keys(): # for each read, detect window haplotypes compatible with that read
					_j+=1
					_read_belongs_2_hap = True
					for _base in _overlap:
						if _SNP_frag[_base]!=_h[window[_base]]:
							_read_belongs_2_hap = False
							break
					if _read_belongs_2_hap:
						Matrix[_i][_j] = 1
		Frag_lsts[_sn] = np.delete(Frag_lsts[_sn], to_delete, axis= None).tolist()
		Matrix = np.delete(Matrix, to_delete, axis=0)
		sample_haplo_dict = EM(sample_haplo_dict, Matrix, l, eps, eps_rel, maxit, min_cov, use_relative=use_relative, verbose=verbose) #, 1e-12, use_relative = True)
		if expression:
			norm = 1.
		else: # normalize the haplotype specific expression rates to obtain haplotype dosages if the data is genomic, i.e. NOT RNA-seq expression data
			norm = sum(sample_haplo_dict[_haplo] for _haplo in all_haplos)
		if queue:
			output = OrderedDict()
			for _haplo in all_haplos:
				output[_haplo] = sample_haplo_dict[_haplo]/norm
			Outqueue[_sample] = output
		else:
			garbage = outhandel.write(_sample+'\t'+'\t'.join("{0:.3f}".format(sample_haplo_dict[_haplo]/norm) for _haplo in all_haplos)+'\n')
	if Parallel or outfile: # close stdout descriptor in case several processes are run or close the outfile handle in case output is given 
		try:
			os.close(outhandel.fileno())
		except:
			pass
	return Outqueue

def EM(haplo_dict, M, l, eps = 1e-06, eps_rel = 1e-04, itermax = 1000, min_cov = 4, use_relative=False, verbose = True):
	""" EM algorithm to obtain haplotype frequencies, using the NGS reads that cover a region and the initial Poisson rates
	for the hidden haplotypes (given in sample_haplo_dict). Absolute convergence test is used with eps value unless use_relative
	is set True in which case relative (less stric with log kernel) convergence is tested with eps_rel.
	Q(R|U) = E[log(poisson(k_1|lu_1)*poisson(k_2|lu_2)*...*poisson(k_N|lu_N))] =
	E[log(exp(-l*u_1)*(l*u_1)**k_1/k_1! * exp(-l*u_2)*(l*u_2)**k_2/k_2! *...* exp(-l*u_N)*(l*u_N)**k_N/k_N!)] =
	E[-l*u_1 + k_1*(log(l)+log(u_1)) - log(k_1!)+
	-l*u_2 + k_2*(log(l)+log(u_2)) - log(k_2!)+
	...
	-l*u_N + k_N*(log(l)+log(u_N)) - log(k_N!)] =
	-l*(u_1+u_2+...+u_N) + c*log(l) + sum(E[k_i|U]*log(u_i) for i in 1,2,...,N) - sum(log(E[k_i|U]!) for i in 1,2,...N)
	=> Q(R|U) = -l*(u_1+u_2+...+u_N) + sum(E[k_i|U]*log(u_i) for i in 1,2,...,N) - sum(log(E[k_i|U]!) for i in 1,2,...N) => u_i = E[k_i|U]/l (Maximization)
	E(k_i|U) obtained from the expectation step."""
	process_name = mpthread.current_process().name
	thread_name = threading.current_thread().name
	converged = False # convergence not reached
	iteration = 1 # current iteration
	c = len(M) # total read count for the region
	if c < min_cov:
		sys.stderr.write('[{0:s}:{1:s}] ERROR: Coverage was too low to estimate haplotypes!\n'.format(process_name, thread_name))
		for _key in haplo_dict.keys():
			haplo_dict[_key] = float('NaN')
		return haplo_dict
	while (not converged and iteration<=itermax):
		count_haplo_dict = get_expectation(haplo_dict, M) # Expectation (E) step
		log_count = np.log([_u for _u in haplo_dict.values()])
		valid_haps = [_n for _n in range(0, len(haplo_dict.values())) if not (math.isnan(log_count[_n]) or math.isinf(log_count[_n]))]
		current_kernel = -l*sum(haplo_dict.values()[_n] for _n in valid_haps) + np.dot([count_haplo_dict.values()[_n] for _n in valid_haps], [log_count[_n] for _n in valid_haps]) - sum(logfact(count_haplo_dict.values()[_n]) for _n in valid_haps) #current kernel value
		new_haplo_dict = OrderedDict() 
		for _h in haplo_dict.keys(): #maximization (M) step
			new_haplo_dict[_h] = count_haplo_dict[_h]/float(l)
		new_count_haplo_dict = get_expectation(new_haplo_dict, M)
		new_log_count = np.log([_u for _u in new_haplo_dict.values()])
		new_valid_haps = [_n for _n in range(0, len(new_haplo_dict.values())) if not (math.isnan(new_log_count[_n]) or math.isinf(new_log_count[_n]))]
		new_kernel = -l*sum(new_haplo_dict.values()[_n] for _n in new_valid_haps) + np.dot([new_count_haplo_dict.values()[_n] for _n in new_valid_haps], [new_log_count[_n] for _n in new_valid_haps]) - sum(logfact(new_count_haplo_dict.values()[_n]) for _n in new_valid_haps) #update kernel value
		diff = new_kernel - current_kernel
		if diff<0 and verbose:
                        sys.stderr.write("[{1:s}:{2:s}] WARNING: likelihood decreased at EM iteration {0:d}!\n".format(iteration, process_name, thread_name))
		if not use_relative:
                        converged = abs(diff) < eps
		else:
			converged = float(diff)/abs(max(new_kernel, current_kernel)) < eps_rel
		haplo_dict = OrderedDict()
		for _key in new_haplo_dict.keys():
			haplo_dict[_key] = new_haplo_dict[_key]
		iteration+=1
	if verbose and not converged:
		sys.stderr.write("[{1:s}:{2:s}] WARNING: convergence NOT acheived after {0:d} iterations!\n".format(itermax, process_name, thread_name))
        elif verbose:
		sys.stderr.write("[{1:s}:{2:s}] Convergence achieved at {0:d}'th iteration!\n".format(iteration-1, process_name, thread_name))
        if verbose and not use_relative:
                sys.stderr.write("[{1:s}:{2:s}] Absolute difference in log likelihood at termination: {0:5.4e}\n".format(diff, process_name, thread_name))
	elif verbose:
                sys.stderr.write("[{1:s}:{2:s}] Relative difference in log likelihood at termination: {0:5.4e}\n".format(float(diff)/abs(max(new_kernel, current_kernel)), process_name, thread_name))
	return haplo_dict 

def logfact(n):
	"""calculate the log factorial function in a way to avoid overflow."""
	global logfactorial_table
	if n<1000.5:
		return logfactorial_table[round(n)]
	else:
		return (n-1/2)*log(n)-n+(1/2)*log(2*math.pi)+1/(12*n) # approximate by Gamma function
