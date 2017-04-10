#! python

## CutDedupTrimQuant.py

## For paired-end RNATagSeq reads with random barcodes on 5' end of read 2.
## After read assignment to libraries: 
##   - cut adaptors using flexbar
##   - deduplicate reads with DedupEndFastQ.py
##   - trim random barcodes from reads and filter for length
##   - quantify at transcript level with kallisto (optional)
##   - align reads to genome with hisat2 (optional)
##   - sort and index hisat2 bamoutput, and make bedgraph for genome browser
## 
## This script leaves many files in the intermediate tempdir; check and delete after use

# hind = /homes/genomes/s.cerevisiae/sacCer3/hisat_indexes/sacCer3_hisat2
# kind = /homes/genomes/s.cerevisiae/sacCer3/kallisto_indexes/sacCer3_orf_coding

# example: python pyRNATagSeq/CutDedupTrimQuant.py -o testCutDedupTrimQuant1 -f SplitCatSelect/EWP01_T01_Lall
# example: python pyRNATagSeq/CutDedupTrimQuant.py -o testCutDedupTrimQuant1 -f testSplitCatSelect/EWP01_T01_Lall_init100000
# example: python pyRNATagSeq/CutDedupTrimQuant.py -o testCutDedupTrimQuant1 -f testSplitCatSelect/EWP01_T01_Lall_init100000 -kind /homes/genomes/s.cerevisiae/sacCer3/kallisto_indexes/sacCer3_orf_coding -hind /homes/genomes/s.cerevisiae/sacCer3/hisat_indexes/sacCer3_hisat2 
import sys, os, csv, gzip, shutil, argparse, warnings



if __name__=="__main__" :  
    parser = argparse.ArgumentParser(description="Demultiplex reads from fastq.gz by inline barcodes")
    parser.add_argument("-f", "--filestem",dest="filestem",nargs='?',help="filestem for paired-end reads")
    parser.add_argument("-td", "--tempdir", dest="tempdir",nargs='?',default="tmp",help="temp directory name")
    parser.add_argument("-kind", "--kallistoindex", dest="kallistoindex",nargs='?',default=False,help="kallisto index")
    parser.add_argument("-hind", "--hisat2index", dest="hisat2index",nargs='?',default=False,help="hisat 2 index")
    parser.add_argument("-o", "--outdir", dest="outdir",nargs='?',default=False,help="output directory")
    options = parser.parse_args()
    
    filestembit = os.path.basename(options.filestem)
    ftarget_cut = options.tempdir + "/" + filestembit + "_cut"
    ftarget_dedup = options.tempdir + "/" + filestembit + "_dedup"
    ftarget_trim = options.tempdir + "/" + filestembit + "_trim"
    ftarget_sam = options.tempdir + "/" + filestembit + "_hisat2" 
    ftarget_bam = options.outdir + "/" + filestembit + "_hisat2" 
    ftarget_kallisto = options.outdir + "/" + filestembit 
    
    # import pdb; pdb.set_trace()
    
    try:
        os.mkdir(options.outdir)
    except Exception:
        warnings.warn("Warning: output directory {} cannot be created".format(options.outdir))
    
    try:
        os.mkdir(options.tempdir)
    except Exception:
        warnings.warn("Warning: temp directory {} cannot be created".format(options.tempdir))
        
    # cut adapters
    
    os.system("flexbar --adapter-trim-end RIGHT -at 2 -ao 5 -l ALL -n 8 -r " + \
    options.filestem + "_R1.fastq.gz -p " + options.filestem + "_R2.fastq.gz " + \
    " -t " + ftarget_cut
    )
    # deduplex
    cut_target = options.tempdir + "/" + filestembit
    os.system("python pyRNATagSeq/DedupEndFastQ.py -r1 " + \
    ftarget_cut + "_1.fastq -r2 "  + \
    ftarget_cut + "_2.fastq -t " + ftarget_dedup
    )
    # trim
    os.system("flexbar --pre-trim-right 5 --min-read-length 40  -n 8  -r " + \
    ftarget_dedup + "_R1.fastq.gz -p " + ftarget_dedup + "_R2.fastq.gz " + 
    " -t " + ftarget_trim
    )
    
    # compress everything in tempdir to reduce footprint?
    
    if (options.kallistoindex) : 
        # kallisto quant
        os.system("kallisto quant --rf-stranded -t 4 -i " + options.kallistoindex + \
        " -o " + ftarget_kallisto + " " + \
        ftarget_trim + "_1.fastq " + ftarget_trim + "_2.fastq ")
    
    
    if (options.hisat2index) : 
        # HISAT2 quant
        os.system("hisat2 -k 5 -p 4 --rna-strandness F --max-intronlen 2000 \
        --non-deterministic -x " + options.hisat2index + 
        " -1 " + ftarget_trim + "_1.fastq -2 " + ftarget_trim + "_2.fastq " + \
        "-S " + ftarget_sam + ".sam")
        ## convert sam (text) output to bam file (compressed binary)
        os.system("samtools view -b " + ftarget_sam + ".sam | " + \
        ## sort bam file on genome and write
        " samtools sort -@ 3 -O bam -o " + ftarget_bam + ".bam ")
        ## index bamfile
        os.system("samtools index " + ftarget_bam + ".bam ")
        
        ## calculate genome coverage for plus strand, compress it
        os.system("bedtools genomecov -ibam " + ftarget_bam + \
        ".bam -bga -strand + | gzip > " + ftarget_bam + \
        "_plus.bedgraph.gz")
        ## calculate genome coverage for minus strand, compress it
        os.system("bedtools genomecov -ibam " + ftarget_bam + \
        ".bam -bga -strand - | gzip > " + ftarget_bam + \
        "_minus.bedgraph.gz")
