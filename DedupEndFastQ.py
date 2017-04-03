#! python

## DedupEndFastQ.py
## Removes duplicate-ended reads from fastq.gz file
## Inputs:
##  - One .fastq.gz file with sequencing reads, 
##  OR - Two paired .fastq.gz files with read pairs in corresponding positions
## - n5 and n3, the number of nts at each end of the read to use for deduplication
##
## Outputs:
##  - a deduplicated file for each out. 
## 
## This is not designed to be memory-efficient.
## 
## Edward Wallace ewjwallace@gmail.com, 2017

import sys, os, csv, gzip, shutil, argparse
from itertools import islice

if __name__=="__main__" :    
    # test1sing: python DedupEndFastQ.py -r1 data/Sample_4reads_R1.fastq.gz -n5 2 -n3 2 -o TestSingleDedup4reads
    # test1sing2: python DedupEndFastQ.py -r1 data/Sample_4reads_R1.fastq.gz -n5 2 -n3 2 -t Slartibartfast
    # test2sing: python DedupEndFastQ.py -r1 data/Sample_init10000_R1.fastq.gz -o TestSingleDedup10000
    # test1pair: python DedupEndFastQ.py -r1 data/Sample_4reads_R1.fastq.gz -r2 data/Sample_4reads_R2.fastq.gz  -n5 2 -n3 2 -o TestPairDedup4reads
    # test1pair2: python DedupEndFastQ.py -r1 data/Sample_4reads_R1.fastq.gz -r2 data/Sample_4reads_R2.fastq.gz  -n5 2 -n3 2 -t Slartibartfast
    # test2pair: python DedupEndFastQ.py -r1 data/Sample_init10000_R1.fastq.gz -r2 data/Sample_init10000_R2.fastq.gz -o TestPairDedup10000
        
    # define input options
    parser = argparse.ArgumentParser(description="Demultiplex reads from fastq.gz by inline barcodes")
    parser.add_argument("-r1", "--read1",dest="r1_fn",nargs='?',help="read 1 filename in fastq.gz format")
    parser.add_argument("-r2", "--read2",dest="r2_fn",default=False,nargs='?',help="read 2 pair filename in fastq.gz format")
    parser.add_argument("-n5", "--nnts5end",dest="n5",type=int,default=20,nargs='?',help="number of nts to use at 5' end")
    parser.add_argument("-n3", "--nnts3end",dest="n3",type=int,default=20,nargs='?',help="number of nts to use at 3' end")
    parser.add_argument("-t", "--target", dest="target",nargs='?',default=False,help="target output filename stem")
    parser.add_argument("-o", "--outdir", dest="outdir",nargs='?',default=False,help="output directory")
    # it would be nice to have mismatches, but not necessary
    # parser.add_argument("-m", "--mismatches", dest="mismatches", default=0, type=int, help="number of mismatches permitted for duplicate")
    options = parser.parse_args()
    
    print "Deduplicating reads for file:\n" + options.r1_fn
    print "counting {} nt from 5' end and {} nt from 3' end".format(options.n5,options.n3)
    # print "allowed mismatches = {}".format(options.mismatches)
    
    # check read 1 fastq is present, if so open it
    if not os.path.isfile(options.r1_fn):
        raise IOError("# Error: read 1 file {} does not exist".format(options.r1_fn))
    r1_fgf = gzip.open(options.r1_fn, 'rt')
    
    # check read 2 fastq is supplied and present, if so open it
    is_paired_end = bool(options.r2_fn)
    if is_paired_end :
        if not os.path.isfile(options.r2_fn):
            raise IOError("# Error: read 2 file {} does not exist".format(options.r2_fn))
        r2_fgf = gzip.open(options.r2_fn, 'rt')
    
    # import pdb; pdb.set_trace()
    
    if options.target :
        if is_paired_end :
            r1outname = options.target + "_R1.fastq.gz"
            r2outname = options.target + "_R2.fastq.gz"
        else : 
            r1outname = options.target + ".fastq.gz"
    else : 
        r1outname = os.path.basename(options.r1_fn)
        if is_paired_end :
            r2outname = os.path.basename(options.r2_fn)
    
    if options.outdir : 
        # make output directory and then file handle for each read end
        if os.path.dirname(options.r1_fn) == options.outdir :
            raise IOError("# Error: output directory should be different from input")        
        try:
            os.mkdir(options.outdir)
        except Exception:
            raise IOError("# Error: output directory {} cannot be created".format(options.outdir))
        r1outname = options.outdir + "/" + r1outname
        if is_paired_end :
            r2outname = options.outdir  + "/" + r2outname
    # make one output file
    r1_out = gzip.open(r1outname,"wt")
    if is_paired_end :
        # make an output file for each paired read
        r2_out = gzip.open(r2outname,"wt")
    
    ntotreads = 0 # counts reads
    ndups = 0 # counts duplicates
    shortseqset = set() # holds the short sequences used for deduplication
    
    ## This loop is the heart of the program
    while True:
        # get fastq record/read (4 lines)
        fqrec1 = list(islice(r1_fgf, 4))
        if not fqrec1 :
            break
        if is_paired_end :
            fqrec2 = list(islice(r2_fgf, 4))
            shortseq = fqrec1[1][:options.n5] + fqrec2[1][:-1][-options.n3:]
        else :
            shortseq = fqrec1[1][:options.n5] + fqrec1[1][:-1][-options.n3:]
        
        # check if new shortseq is already there
        if shortseq in shortseqset :
            ndups += 1
            # import pdb; pdb.set_trace()
        else :
            shortseqset.add(shortseq)
            r1_out.writelines( fqrec1 )
            if is_paired_end :
                r2_out.writelines( fqrec2 )
        
        # count number of processed reads, output every millionth
        ntotreads += 1
        if (ntotreads % 1000000) == 0:
            print "{} reads read".format(ntotreads)
    
    print "removed {} duplicates from {} reads".format(ndups,ntotreads)
    
    # close output handles and fasta file
    r1_out.close()
    r1_fgf.close()
    if is_paired_end :
        r2_out.close()
        r2_fgf.close()

    print "Done"
