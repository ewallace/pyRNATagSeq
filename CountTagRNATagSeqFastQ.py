#! python

## CountTagRNATagSeq2FastQ.py
## Counts initial 9-nt barcode (TagRead) for fastq file
## Inputs:
##  - .fastq or .fastq.gz file with sequencing reads
##  - output file name (optional)
##  - length of TagRead (optional)
##
## Outputs:
##  - a tab delimited text file containing counts with each 5'-end Tag
## 
## Edward Wallace ewjwallace@gmail.com, 2017

import sys, os, csv, gzip, shutil, argparse
from itertools import islice, product

def init_fastq_record(fqr,n=9):
    """return initial n letters from fastq record"""
    return fqr[1][:n]

if __name__=="__main__" :    
    # test1: python CountTagRNATagSeqFastQ.py -o Sample_4reads_TagCounts.txt data/Sample_4reads_R1.fastq.gz
    # test2: python CountTagRNATagSeqFastQ.py -o Sample_10000_TagCounts.txt data/Sample_init10000_R1.fastq.gz
    # check reads add up by:
    #    import pandas as pd; cts = pd.read_csv("Sample_10000_TagCounts.txt",delimiter="\t"); cts.Count.sum() == 10000;
    
    # define input options
    parser = argparse.ArgumentParser(description="Demultiplex reads from fastq.gz by inline barcodes")
    # Required arguments
    parser.add_argument(dest="fastq_fn",nargs='?',help="reads filename in fastq.gz format")
    # Optional arguments
    parser.add_argument("-o", "--outfile", dest="out_fn",nargs='?',default="TagCodesCounted.txt",help="output file for Counts")
    parser.add_argument("-tl", "--taglength", dest="taglength",nargs='?',default=9,help="length of TagRead")
    options = parser.parse_args()
    
    print "Counting TagReads for file:\n" + options.fastq_fn
    
    if not os.path.isfile(options.fastq_fn):
	 	raise IOError("# Error: file {} does not exist".format(options.samplesheet_fn))
    
    # setup input fastq
    if not os.path.isfile(options.fastq_fn):
	 	raise IOError("# Error: file {} does not exist".format(options.fastq_fn))
	
	# Make list of all potential tags
    # TagIter = product(["T","C","A","G"], repeat = options.taglength)
    Tags = [''.join(Tag) for Tag in product(["T","C","A","G"], repeat = options.taglength)]
    TagCounts = dict((Tag , 0 ) for Tag in Tags)
    # any keys with an "N" will be mapped to Bad Tag, all N's. 
    BadTag = "".join(["N" for i in range(options.taglength)])
    TagCounts[BadTag] = 0
    
    # open fastq_file
    with gzip.open(options.fastq_fn, 'rt') as fqf :
        ntotreads = 0
        while True:
            # get fastq record/read (4 lines)
            fqrec = list(islice(fqf, 4))
            if not fqrec :
                break
            # count number of processed reads
            ntotreads += 1
            if (ntotreads % 1000000) == 0:
                print "{} reads processed".format(ntotreads)
            # TagReads.append( init_fastq_record(fqrec,n=options.taglength) )
            # Increment the counter
            try : 
                TagCounts[ init_fastq_record( fqrec, n=options.taglength ) ] += 1
            except KeyError :
                TagCounts[BadTag] += 1
    
    print "All {} reads processed and counted".format(ntotreads)
    
    # This next line is very slow
    # TagCounts = [[x,TagReads.count(x)] for x in set(TagReads)]
    
    print "Writing output".format(ntotreads)

    with open(options.out_fn, "wb") as out_f:
       writer = csv.writer(out_f,dialect="excel-tab")
       writer.writerow( ["TagRead","Count"] )
       writer.writerows( [ [Tag,Count] for Tag, Count in TagCounts.items() ])
    
    # missing: function call/comments in output
    
    print "Done"


