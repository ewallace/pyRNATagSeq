# pyRNATagSeq
For counting and splitting RNATagSeq reads.

Edward Wallace, ewjwallace@gmail.com, 2017.
https://github.com/ewallace
Please let me know if you use this, and report bugs!

RNA sequencing (RNASeq) library preparation using the RNATagSeq method of Shishkin, et al. (2015). *Simultaneous generation of many RNA-seq libraries in a single reaction.* Nature Methods, 12(4), 323–325. http://doi.org/10.1038/nmeth.3313 (Broad institute).
Edward Wallace's mild updates to the library prep protocol are at: http://drummondlab.org/protocols/protocol/rnatagseq-library-prep

I wrote this because fastx_barcode_splitter, from fastx-toolkit, was prohibitively slow.

## Contents

CountTagRNATagSeqFastQ.py counts the reads in a fastq.gz with each _possible_ TagRead. Useful for debugging and counting mismatches.

SplitRNATagSeqFastQ.py splits (single or paired end) fastq.gz reads by RNATag in 5' end of read 1, depending on user-supplied sheet of SampleIDs and TagReads

data/Sample_4reads_R1.fastq.gz - Artificial sample with 4 read 1s for debugging

data/Sample_4reads_R2.fastq.gz - 4 read 2s corresponding to Sample_4reads_R1.fastq.gz

data/Sample_init10000_R1.fastq.gz - initial 10000 read 1s from a paired-end S. cerevisiae dataset

data/Sample_init10000_R2.fastq.gz - 10000 read 2s corresponding to Sample_init10000_R1.fastq.gz

data/TagSeqBarcodedOligos2015.txt - TagSeq barcoded oligos used in Shishkin et al. Edit this with your meaningful sample ids!


## Examples

Count:
python CountTagRNATagSeqFastQ.py -o Sample_4reads_TagCounts.txt data/Sample_4reads_R1.fastq.gz
python CountTagRNATagSeqFastQ.py -o Sample_10000_TagCounts.txt data/Sample_init10000_R1.fastq.gz

Split single-end:
python SplitRNATagSeqFastQ.py -r1 data/Sample_4reads_R1.fastq.gz -ss data/TagSeqBarcodedOligos2015.txt -o TestSingleSplit4reads
python SplitRNATagSeqFastQ.py -r1 data/Sample_init10000_R1.fastq.gz -ss data/TagSeqBarcodedOligos2015.txt -o TestSingleSplit10000

Split paired-end:
python SplitRNATagSeqFastQ.py -r1 data/Sample_4reads_R1.fastq.gz -r2 data/Sample_4reads_R2.fastq.gz -ss data/TagSeqBarcodedOligos2015.txt -o TestPairSplit4reads
python SplitRNATagSeqFastQ.py -r1 data/Sample_init10000_R1.fastq.gz -r2 data/Sample_init10000_R2.fastq.gz -ss data/TagSeqBarcodedOligos2015.txt -o TestPairSplit10000