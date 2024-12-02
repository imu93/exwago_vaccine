# Scripts and files associated with the vaccine manuscript
---
## To align trimm adpters from sRNA-seq libraries we used reaper from the Kraken package (Davis et al., 2013):
```
$srcPath/reaper -i $file -geom no-bc -3pa $adp3 -3p-global 12/2/1 -3p-prefix 8/2/1 -3p-head-to-tail 1 -nnn-check 3/5 -format-clean @%I%n%C%n+%n%Q%n -polya  5 -qqq-check 35/10 -tri 35""
```
	1.	-i $file
Specifies the input file containing raw sequencing reads. $file is a placeholder for the actual file name.
	2.	-geom no-bc
Indicates that there are no barcodes in the reads. This disables barcode parsing and processing.
	3.	-3pa $adp3
Specifies the 3’ adapter sequence ($adp3) to be trimmed from the reads. $adp3 is a placeholder for the actual adapter sequence.
	4.	-3p-global 12/2/1
Defines the scoring scheme for global alignment during 3’ adapter trimming:
	•	12: Match score.
	•	2: Mismatch penalty.
	•	1: Gap penalty.
	5.	-3p-prefix 8/2/1
Sets the scoring scheme for 3’ adapter trimming using a prefix alignment:
	•	8: Match score.
	•	2: Mismatch penalty.
	•	1: Gap penalty.
	6.	-3p-head-to-tail 1
Specifies whether to trim adapters in head-to-tail alignment mode:
	•	1: Enables head-to-tail alignment (adapter sequences concatenated).
	7.	-nnn-check 3/5
Filters reads based on the occurrence of consecutive N (ambiguous bases):
	•	3: Maximum allowable consecutive N.
	•	5: Fractional threshold of ambiguous bases to discard a read.
	8.	-format-clean @%I%n%C%n+%n%Q%n
Specifies the output format for cleaned reads, matching the standard FASTQ format:
	•	@%I: Read ID.
	•	%C: Read sequence.
	•	+%n: Separator line.
	•	%Q: Quality scores.
	9.	-polya 5
Trims poly-A tails from reads if the tail is at least 5 bases long.
	10.	-qqq-check 35/10
Filters reads based on quality:
	•	35: Quality threshold.
	•	10: Minimum percentage of bases in the read that must meet the quality threshold.
	11.	-tri 35
Trims low-quality bases from the ends of reads until the remaining sequence has an average quality score of at least 35.
 
## sRNA-seq alignmet
We used bowtie to align reads in a range of 15 to 50 nt using pullseq and bowtie1 
```
for file in *.fastq.gz; do outFile=$(echo $file | sed 's|.fastq.gz|.ps_15_50.fastq|'); pullseq -i $file -m 15 -a 50 > $outFile; pigz $outFile; done
```
```
for file in *.ps_15_50.fastq.gz; do
   outFile=$(echo $file | sed 's/.fastq.gz//')
   echo "bowtie -S -q -v 1 -M 1 --all -t --threads $thr --best --strata -x $genome $file | samtools view -F 4 -b > ${outFile}.M1.us.bam"
   bowtie -S -q -v 1 -M 1 --all -t --threads $thr --best --strata -x $genome $file | samtools view -F 4 -b > ${outFile}.M1.us.bam
   samtools sort -o $outFile.M1.bam -@ $thr $outFile.M1.us.bam
   rm -rf $outFile.M1.us.bam
done
```

The resulting bam files were transformed into GenomicRanges using the `01.bam2Rds.R` script which produce GRanges objects with the following structure:
seqnames    ranges strand |       len                     Id flag              seq
Acey_s0001_scaf     71-86      - |        16 VH00144:497:AACVH52M..	16 GAAGAAACTGGGTCTT

Where teh colums are
1.seqnames: Chr/Scafold
2.range: start and end of alignment
3.strand: strand of alignment
4.len: read length
5.ID: read ID
6.flag: alignment flag
7.seq: read sequence

Using this files we fist produce matrices of raw counts with the `02.get_fn_mtx.R` and then we transformed this counts to cpms using the cpm function from edgeR () with the `03.fn_long_cpm.R` script 

