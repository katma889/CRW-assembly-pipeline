# CRW-assembly-pipeline
## we sequenced the clover root weevil( Sitona obsoletus) using the data from 4 flow cell runs with 4 different individual and the data from all 4 runs were combined and basecalled with the guppy 5 version.
script for Guppy 5 
```
#!/bin/bash -e

#SBATCH --job-name=guppy_crw                 #name of the job
#SBATCH --account=uoo02772              #my project number in nesi
#SBATCH --time=10:00:00                 #wall time
#SBATCH --partition=gpu                 #guppy runs faster in gpu partition in nesi, than other partition
#SBATCH --gres=gpu:1                    #some configuration for gpu partition, that i don't understand, suggested by nesi support
#SBATCH --mem=6G                                # memory 6gb
#SBATCH --ntasks=4                              #ntask set to 4
#SBATCH --cpus-per-task=1               #cpu per task set to 1
#SBATCH --output=%x-%j.out              #%x gives job name and %j gives job number, this is slurm output file
#SBATCH --error=%x-%j.err               #similar slurm error file
#SBATCH --mail-type=ALL
#SBATCH --mail-user=katma889@student.otago.ac.nz

module load ont-guppy-gpu/5.0.7
guppy_basecaller -i ../ -s . --flowcell FLO-MIN106 --kit SQK-LSK109 --num_callers 4 -x auto --recursive --trim_barcodes --disable_qscore_filtering

```
##  Along we the merged fastqc files from guppy we also got the sequencing summary.txt file which we process further with pycoQC-2.5.2 which hives the html file in output makes the viewing data for quality easier.
Script for pycoqc
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --partition=large
#SBATCH --job-name pycoqc
#SBATCH --mem=50G
#SBATCH --time=01:00:00
#SBATCH --account=uoo02772
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

export PATH="/nesi/nobackup/uoo02752/nematode/bin/miniconda3/bin:$PATH"

pycoQC -f ../sequencing_summary.txt -o pycoQC_output.html

```
###After viwing the quality of our output data via html file we further proceed to remove the reads mapping to the lambda phage genome from our fastq files using Nanolyse. This is because we used DNA CS while running our sample in the minion flow cells. The script for nanolyse id given below;

Script for Nanolyse
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --partition=large
#SBATCH --job-name nanolyse.job
#SBATCH --mem=50G
#SBATCH --time=08:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

export PATH="/nesi/nobackup/uoo02752/nematode/bin/miniconda3/bin:$PATH"

cat ../crw.ont.merged.fastq | NanoLyse --reference ./dna_cs.fasta | gzip > crw_nanopore_filtered.fastq.gz

```
Similarly the dna_cs. fasta we used in our experiment is given below;
```
>DNA_CS
GCCATCAGATTGTGTTTGTTAGTCGCTTTTTTTTTTTGGAATTTTTTTTTTGGAATTTTTTTTTTGCGCTAACAACCTCCTGCCGTTTTGCCCGTGCATATCGGTCACGAACAAATCTGATTACTAAACACAGTAGCCTGGATTTGTTCTATCAGTAATCGACCTTATTCCTAATTAAATAGAGCAAATCCCCTTATTGGGGGTAAGACATGAAGATGCCAGAAAAACATGACCTGTTGGCCGCCATTCTCGCGGCAAAGGAACAAGGCATCGGGGCAATCCTTGCGTTTGCAATGGCGTACCTTCGCGGCAGATATAATGGCGGTGCGTTTACAAAAACAGTAATCGACGCAACGATGTGCGCCATTATCGCCTAGTTCATTCGTGACCTTCTCGACTTCGCCGGACTAAGTAGCAATCTCGCTTATATAACGAGCGTGTTTATCGGCTACATCGGTACTGACTCGATTGGTTCGCTTATCAAACGCTTCGCTGCTAAAAAAGCCGGAGTAGAAGATGGTAGAAATCAATAATCAACGTAAGGCGTTCCTCGATATGCTGGCGTGGTCGGAGGGAACTGATAACGGACGTCAGAAAACCAGAAATCATGGTTATGACGTCATTGTAGGCGGAGAGCTATTTACTGATTACTCCGATCACCCTCGCAAACTTGTCACGCTAAACCCAAAACTCAAATCAACAGGCGCCGGACGCTACCAGCTTCTTTCCCGTTGGTGGGATGCCTACCGCAAGCAGCTTGGCCTGAAAGACTTCTCTCCGAAAAGTCAGGACGCTGTGGCATTGCAGCAGATTAAGGAGCGTGGCGCTTTACCTATGATTGATCGTGGTGATATCCGTCAGGCAATCGACCGTTGCAGCAATATCTGGGCTTCACTGCCGGGCGCTGGTTATGGTCAGTTCGAGCATAAGGCTGACAGCCTGATTGCAAAATTCAAAGAAGCGGGCGGAACGGTCAGAGAGATTGATGTATGAGCAGAGTCACCGCGATTATCTCCGCTCTGGTTATCTGCATCATCGTCTGCCTGTCATGGGCTGTTAATCATTACCGTGATAACGCCATTACCTACAAAGCCCAGCGCGACAAAAATGCCAGAGAACTGAAGCTGGCGAACGCGGCAATTACTGACATGCAGATGCGTCAGCGTGATGTTGCTGCGCTCGATGCAAAATACACGAAGGAGTTAGCTGATGCTAAAGCTGAAAATGATGCTCTGCGTGATGATGTTGCCGCTGGTCGTCGTCGGTTGCACATCAAAGCAGTCTGTCAGTCAGTGCGTGAAGCCACCACCGCCTCCGGCGTGGATAATGCAGCCTCCCCCCGACTGGCAGACACCGCTGAACGGGATTATTTCACCCTCAGAGAGAGGCTGATCACTATGCAAAAACAACTGGAAGGAACCCAGAAGTATATTAATGAGCAGTGCAGATAGAGTTGCCCATATCGATGGGCAACTCATGCAATTATTGTGAGCAATACACACGCGCTTCCAGCGGAGTATAAATGCCTAAAGTAATAAAACCGAGCAATCCATTTACGAATGTTTGCTGGGTTTCTGTTTTAACAACATTTTCTGCGCCGCCACAAATTTTGGCTGCATCGACAGTTTTCTTCTGCCCAATTCCAGAAACGAAGAAATGATGGGTGATGGTTTCCTTTGGTGCTACTGCTGCCGGTTTGTTTTGAACAGTAAACGTCTGTTGAGCACATCCTGTAATAAGCAGGGCCAGCGCAGTAGCGAGTAGCATTTTTTTCATGGTGTTATTCCCGATGCTTTTTGAAGTTCGCAGAATCGTATGTGTAGAAAATTAAACAAACCCTAAACAATGAGTTGAAATTTCATATTGTTAATATTTATTAATGTATGTCAGGTGCGATGAATCGTCATTGTATTCCCGGATTAACTATGTCCACAGCCCTGACGGGGAACTTCTCTGCGGGAGTGTCCGGGAATAATTAAAACGATGCACACAGGGTTTAGCGCGTACACGTATTGCATTATGCCAACGCCCCGGTGCTGACACGGAAGAAACCGGACGTTATGATTTAGCGTGGAAAGATTTGTGTAGTGTTCTGAATGCTCTCAGTAAATAGTAATGAATTATCAAAGGTATAGTAATATCTTTTATGTTCATGGATATTTGTAACCCATCGGAAAACTCCTGCTTTAGCAAGATTTTCCCTGTATTGCTGAAATGTGATTTCTCTTGATTTCAACCTATCATAGGACGTTTCTATAAGATGCGTGTTTCTTGAGAATTTAACATTTACAACCTTTTTAAGTCCTTTTATTAACACGGTGTTATCGTTTTCTAACACGATGTGAATATTATCTGTGGCTAGATAGTAAATATAATGTGAGACGTTGTGACGTTTTAGTTCAGAATAAAACAATTCACAGTCTAAATCTTTTCGCACTTGATCGAATATTTCTTTAAAAATGGCAACCTGAGCCATTGGTAAAACCTTCCATGTGATACGAGGGCGCGTAGTTTGCATTATCGTTTTTATCGTTTCAATCTGGTCTGACCTCCTTGTGTTTTGTTGATGATTTATGTCAAATATTAGGAATGTTTTCACTTAATAGTATTGGTTGCGTAACAAAGTGCGGTCCTGCTGGCATTCTGGAGGGAAATACAACCGACAGATGTATGTAAGGCCAACGTGCTCAAATCTTCATACAGAAAGATTTGAAGTAATATTTTAACCGCTAGATGAAGAGCAAGCGCATGGAGCGACAAAATGAATAAAGAACAATCTGCTGATGATCCCTCCGTGGATCTGATTCGTGTAAAAAATATGCTTAATAGCACCATTTCTATGAGTTACCCTGATGTTGTAATTGCATGTATAGAACATAAGGTGTCTCTGGAAGCATTCAGAGCAATTGAGGCAGCGTTGGTGAAGCACGATAATAATATGAAGGATTATTCCCTGGTGGTTGACTGATCACCATAACTGCTAATCATTCAAACTATTTAGTCTGTGACAGAGCCAACACGCAGTCTGTCACTGTCAGGAAAGTGGTAAAACTGCAACTCAATTACTGCAATGCCCTCGTAATTAAGTGAATTTACAATATCGTCCTGTTCGGAGGGAAGAACGCGGGATGTTCATTCTTCATCACTTTTAATTGATGTATATGCTCTCTTTTCTGACGTTAGTCTCCGACGGCAGGCTTCAATGACCCAGGCTGAGAAATTCCCGGACCCTTTTTGCTCAAGAGCGATGTTAATTTGTTCAATCATTTGGTTAGGAAAGCGGATGTTGCGGGTTGTTGTTCTGCGGGTTCTGTTCTTCGTTGACATGAGGTTGCCCCGTATTCAGTGTCGCTGATTTGTATTGTCTGAAGTTGTTTTTACGTTAAGTTGATGCAGATCAATTAATACGATACCTGCGTCATAATTGATTATTTGACGTGGTTTGATGGCCTCCACGCACGTTGTGATATGTAGATGATAATCATTATCACTTTACGGGTCCTTTCCGGTGAAAAAAAAGGTACCAAAAAAAACATCGTCGTGAGTAGTGAACCGTAAGC
```
###Then we used porechop to remove the adapters from our reads. The scrpit for porechop is given below ;
Script for Porechop
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --partition=large
#SBATCH --job-name porechop
#SBATCH --mem=50G
#SBATCH --time=04:00:00
#SBATCH --account=uoo02772
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load Porechop/0.2.4-gimkl-2020a-Python-3.8.2
porechop -i ../crw_nanopore_filtered.fastq.gz -o crw_ont_nanolyse_porechop.fastq.gz --threads 10
```
## We used FLye to assemble the long read data from Oxford Minion. The script for flye assembly algorithm is given below;
Script for Flye 
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --partition=hugemem
#SBATCH --job-name flye.crwV3
#SBATCH --mem=150G
#SBATCH --time=72:00:00
#SBATCH --account=uoo02772
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=katma889@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load Flye/2.9-gimkl-2020a-Python-3.8.2

flye --nano-hq ../crw_ont_nanolyse_porechop_nanofilt.fastq.gz -o ./Flye -t 10 -i 3 

```
This Flye output also includes a main file assembly.fasta which is further used for running the Quast. Quast is usually used to evaluate the assembly quality even without a reference genome. The script of Quast is given below;
Script for Quast
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 6
#SBATCH --partition=large
#SBATCH --job-name quast.crw
#SBATCH --mem=30G
#SBATCH --time=09:00:00
#SBATCH --account=uoo02772
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=katma889@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load QUAST
#=========> running quast for assembly quality

quast.py -t 10 --eukaryote --large --conserved-genes-finding --k-mer-stats \
assembly.fasta \
-o quastqless 
```
## By running the above script it yielded a Quast folder wwhich yielded a main file called report.txt.
The above mentioned report.txt yielded us  the Complete BUSCO percent of 91.42 and partial BUSCO percent of 5.61 with the number of contigs equals to 82815.
## Then we used Purgehaplotigs to remove the haplotigs from our assembly. It helps us to to identify and reassign the duplicate contigs to improve our assembly. The script that we ran is given below;
Script for Purgehaplotigs
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --job-name purgehap.crw
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --account=uoo02772
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=katma889@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load SAMtools/1.12-GCC-9.2.0
module load minimap2/2.20-GCC-9.2.0
module load BEDTools/2.29.2-GCC-9.2.0

#minimap2 -t 10 -ax map-ont assembly.fasta crw_ont_nanolyse_porechop_nanofilt.fastq.gz \
#--secondary=no | samtools sort -m 5G -o aligned.bam -T tmp.ali

export PATH="/nesi/nobackup/uoo02772/bin/miniconda3/bin:$PATH"
#purge_haplotigs hist -b aligned.bam -g assembly.fasta -t 10

#purge_haplotigs cov -i aligned.bam.gencov -l 0 -m 20 -h 199 -o coverage_stats.csv

purge_haplotigs purge -g assembly.fasta -c coverage_stats.csv -b aligned.bam

#awk '{print $1",s,"}' gapclosed.fasta.pilon3.fasta.fai > cov_stat.csv
#purge_haplotigs purge -g gapclosed.fasta.pilon3.fasta -c cov_stat.csv -b aligned.bam
```
### This yielded us the file called curated.fasta which we further ran quast on it. This purge haplotigs bring down the contigs number to 51390 from 82815. However, the complete BUSCO percent was sligthly reduded to 90.10 and little increase on partial BUSCO to 6.27. Therefore we further used the RagTag algorithm  a toolset for automating assembly scaffolding and patching our long read assembly. The script for RAgTag is given below;
Script for RagTag
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --job-name ragtag.crw
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --account=uoo02772
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=katma889@student.otago.ac.nz
#SBATCH --hint=nomultithread

export PATH="/nesi/nobackup/uoo02752/nematode/bin/miniconda3/bin:$PATH"

ragtag.py scaffold curated.haplotigs.fasta curated.fasta
```
By running above script we got ragtag.scaffold.fasta as our main output and we can check the stats file ( for example ragtag.scaffold.stats) to see  the number of scaffold it removed. Our result is given below;
```
placed_sequences        placed_bp	unplaced_sequences	unplaced_bp     gap_bp  gap_sequences
15186   629971398	36204   511046512	377200  3772
```
## We further ran quast in the output file from ragtag that is ragtag.scaffold.fasta and it further reduced to number of contigs to 47618 with the complete BUSCO percent and partial BUSCO percent to 90.10 and 5.94 respectively.
## Then we used lrscaff to further scaffold the assembly from ragtag using long reads that is crw_ont_nanolyse_porechop_nanofilt.fasta in our case. The script for lrscaff is given below;
Script for LRScaff
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --job-name lrscaf.crw
#SBATCH --mem=80G
#SBATCH --time=72:00:00
#SBATCH --account=uoo02772
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=katma889@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load minimap2

minimap2 -t 10 ragtag.scaffold.fasta crw_ont_nanolyse_porechop_nanofilt.fasta > ./aln.mm
export PATH="/nesi/nobackup/uoo02752/bin/lrscaf/target/:$PATH"

java -Xms80g -Xmx80g -jar /nesi/nobackup/uoo02752/bin/lrscaf/target/LRScaf-1.1.11.jar --contig ragtag.scaffold.fasta --alignedFile aln.mm -t mm -p 10 --output ./scaffolds1

minimap2 -t 10 ./scaffolds1/scaffolds.fasta crw_ont_nanolyse_porechop_nanofilt.fasta > ./scaffolds1/aln.mm
export PATH="/nesi/nobackup/uoo02752/bin/lrscaf/target/:$PATH"

java -Xms80g -Xmx80g -jar /nesi/nobackup/uoo02752/bin/lrscaf/target/LRScaf-1.1.11.jar --contig ./scaffolds1/scaffolds.fasta --alignedFile ./scaffolds1/aln.mm -t mm -p 10 --output ./scaffolds1/scaffolds2

minimap2 -t 10 ./scaffolds1/scaffolds2/scaffolds.fasta crw_ont_nanolyse_porechop_nanofilt.fasta > ./scaffolds1/scaffolds2/aln.mm
export PATH="/nesi/nobackup/uoo02752/bin/lrscaf/target/:$PATH"

java -Xms80g -Xmx80g -jar /nesi/nobackup/uoo02752/bin/lrscaf/target/LRScaf-1.1.11.jar --contig ./scaffolds1/scaffolds2/scaffolds.fasta --alignedFile ./scaffolds1/scaffolds2/aln.mm -t mm -p 10 --output ./scaffolds1/scaffolds2/scaffolds3

minimap2 -t 10 ./scaffolds1/scaffolds2/scaffolds3/scaffolds.fasta crw_ont_nanolyse_porechop_nanofilt.fasta > ./scaffolds1/scaffolds2/scaffolds3/aln.mm
export PATH="/nesi/nobackup/uoo02752/bin/lrscaf/target/:$PATH"

java -Xms80g -Xmx80g -jar /nesi/nobackup/uoo02752/bin/lrscaf/target/LRScaf-1.1.11.jar --contig ./scaffolds1/scaffolds2/scaffolds3/scaffolds.fasta --alignedFile ./scaffolds1/scaffolds2/scaffolds3/aln.mm -t mm -p 10 --output ./scaffolds1/scaffolds2/scaffolds3/scaffolds4

minimap2 -t 10 ./scaffolds1/scaffolds2/scaffolds3/scaffolds4/scaffolds.fasta crw_ont_nanolyse_porechop_nanofilt.fasta > ./scaffolds1/scaffolds2/scaffolds3/scaffolds4/aln.mm
export PATH="/nesi/nobackup/uoo02752/bin/lrscaf/target/:$PATH"

java -Xms80g -Xmx80g -jar /nesi/nobackup/uoo02752/bin/lrscaf/target/LRScaf-1.1.11.jar --contig ./scaffolds1/scaffolds2/scaffolds3/scaffolds4/scaffolds.fasta --alignedFile ./scaffolds1/scaffolds2/scaffolds3/scaffolds4/aln.mm -t mm -p 10 --output ./scaffolds1/scaffolds2/scaffolds3/scaffolds4/scaffolds5
```
### we checked the output from the lrscaff from quast and scaffolds5_scaffolds resulted in the reduction in the number of contigs to 26788 with a complete BUSCO percent of 94.72 and partial BUSCO 2.64 percent respectively. 

## we further scaffold our genome usimng long DNA sequences from ONT (crw_ont_nanolyse_porechop_nanofilt.fasta) using RAILS v1.5.1 and Cobbler v0.6.1. Here RAIlS is known as all-in -one scaffolder and gap-filler whereas Cobbler is tool that patch gaps automatically. The script used for this is given below;
Script for RAILS and Cobbler
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --job-name rails.cob.crw
#SBATCH --mem=50G
#SBATCH --time=72:00:00
#SBATCH --account=uoo02772
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=katma889@student.otago.ac.nz
#SBATCH --hint=nomultithread


module load Perl/5.30.1-GCC-9.2.0
module load minimap2
module load SAMtools/1.13-GCC-9.2.0
export PATH="/nesi/nobackup/uoo02752/nematode/bin/RAILS/bin:$PATH"

sh runRAILSminimapSTREAM.sh scaffolds.fasta crw_ont_nanolyse_porechop_nanofilt.fasta 250 0.80 500 2 ont \
/scale_wlg_persistent/filesets/opt_nesi/CS400_centos7_bdw/SAMtools/1.13-GCC-9.2.0/bin/samtools 10

```
#### By running the script given above first we will get the output of cobbler as crw_ont_nanolyse_porechop_nanofilt.fasta_vs_scaffolds.fasta_250_0.80_gapsFill.fa in my case which was further used by RAILS to get the final output as crw_ont_nanolyse_porechop_nanofilt.fasta_vs_scaffolds.fasta_250_0.80_rails.scaffolds.fa at the end. The summary of my RAILS log file is given below;
```
/RAILS Summary ---------------
Number of merges induced : 47
Average closed gap length (bp) : 2916.23
Closed gap length st.dev +/- : 6705.50
Total bases added : 137063
Largest gap resolved (bp) : 38027
Shortest gap resolved (bp) : 6
---------------------------------------------
*0 bp gaps are not counted towards the average
crw_ont_nanolyse_porechop_nanofilt.fasta_vs_scaffolds.fasta_250_0.80_rails.log 

```
We then ran quast for the output file from the rails. The QUAST report is given below;
```
Assembly                    crw_ont_nanolyse_porechop_nanofilt.fasta_vs_scaffolds.fasta_250_0.80_rails.scaffolds
# contigs (>= 0 bp)         26741                                                                               
# contigs (>= 1000 bp)      23273                                                                               
# contigs (>= 5000 bp)      15510                                                                               
# contigs (>= 10000 bp)     13401                                                                               
# contigs (>= 25000 bp)     10776                                                                               
# contigs (>= 50000 bp)     8265                                                                                
Total length (>= 0 bp)      1507291159                                                                          
Total length (>= 1000 bp)   1504938917                                                                          
Total length (>= 5000 bp)   1486560638                                                                          
Total length (>= 10000 bp)  1471422085                                                                          
Total length (>= 25000 bp)  1427853847                                                                          
Total length (>= 50000 bp)  1335778686                                                                          
# contigs                   17613                                                                               
Largest contig              4214430                                                                             
Total length                1494738492                                                                          
GC (%)                      31.82                                                                               
N50                         178346                                                                              
N75                         98144                                                                               
L50                         2465                                                                                
L75                         5281                                                                                
# N's per 100 kbp           8276.39                                                                             
Complete BUSCO (%)          94.06                                                                               
Partial BUSCO (%)           2.97

```
## Then we further used the LRScaff to further boost the contiguity of our assemly using ONT long filetered reads.The script for LRScaff is given below;
script for LRScaff
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 16
#SBATCH --job-name lr-gapEW
#SBATCH --mem=40G
##SBATCH --time=00:15:00
#SBATCH --time=50:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bhaup057@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load BWA/0.7.17-gimkl-2017a
export PATH=/nesi/nobackup/uoo02752/bin/LR_Gapcloser/src/:$PATH

sh LR_Gapcloser.sh -i crw.scaffolds.fasta -l crw_ont_nanolyse_porechop_nanofilt.fasta -s n -t 16 -r 10

```
By running above script we got iteration -1 to iteration-10 folders with the common filename gapclosed.fasta in each of them. Then we ran quast (iteration10) to check the quality of the assembly, LRScaff reduced gaps (N's per 100 kbp) to 3002.56 from 8276.39   using crw_ont_nanolyse_porechop_nanofilt.fasta for gap filling and increased the partial BUSCO to 3.30 from 2.97 while complete busco and number of contigs the same as previous.

## we further scaffold the output (gapclosed.fasta) by our supernova assembly from 10X (scaffold crw.10x.all.pseudo.fasta) ragtag. 
Script for ragtag
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --job-name ragtag.10x
#SBATCH --mem=50G
#SBATCH --time=05:00:00
#SBATCH --account=uoo02772
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=katma889@student.otago.ac.nz
#SBATCH --hint=nomultithread

export PATH="/nesi/nobackup/uoo02752/nematode/bin/miniconda3/bin:$PATH"

ragtag.py scaffold crw.10x.all.pseudo.fasta gapclosed.fasta

```
### By running above script we got ragtag.scaffold.fasta as output file which we further used for scafollding agian using ragtag. We renamed this output as assembly.fasta and also changed the sequence header to make it compatible to use in scaffolding using ragtag. The script for ragtag for second times is same except using assembly.fasta instead of gapclosed.fasta. By running the ragtag script this time we used ragtag.scaffold.fasta. 

## We further used ARBitR for further merging and scaffolding our current genome assembly from ragtag. As it takes alignment file in the bam/sam format (for example possorted_bam.bam in our case)  with 10X Chromium barcodes when provided with genome fasta file used for mapping, it will sort and merge the provided contigs into scaffolds. Therefore, first we created a reference data ragtag.scaffold.fasta and alignment in the folder id CRW using longranger to be used by the ARBitR. 

Script for longranger

```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --job-name longR_CRW
#SBATCH --mem=50G
#SBATCH --time=72:00:00
#SBATCH --account=uoo02772
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=katma889@student.otago.ac.nz
#SBATCH --hint=nomultithread

export PATH=/nesi/project/uoo02752/bin/longranger-2.2.2:$PATH

#longranger mkref ragtag.scaffold.fasta

longranger align --id=CRW \
--fastqs=/nesi/nobackup/uoo02772/crw/10x/1.raw.hiseq.novaseq \
--reference=/nesi/nobackup/uoo02772/crw/2.nanopore/1.CRW_nanopore_rawdata/guppy.5/nanolyse/porechop/nanoqc/nanofilt/flye/Flye/purgehaplotigs/ragtag_output/lrscaff/scaffolds1/scaffolds2/scaffolds3/scaffolds4/scaffolds5/rails.cobbler/lrgapcloser/Output/sn.10x.ragtag/ragtag_output/ragtag.2/ragtag_output/arbitr/refdata-ragtag.scaffold

```
## Then we ran ARBitR on the aligned scaffold fasta file using possorted_bam.bam.

Script for ARBitR

```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --job-name arbitr.crw
#SBATCH --mem=5G
#SBATCH --time=05:00:00
#SBATCH --account=uoo02772
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=katma889@student.otago.ac.nz
#SBATCH --hint=nomultithread

export PATH="/nesi/nobackup/uoo02752/nematode/bin/ARBitR/src:$PATH"
export PATH="/nesi/nobackup/uoo02752/nematode/bin/miniconda3/bin:$PATH"

arbitr.py -i ragtag.scaffold.fasta -o output.arbitr.scaffolds ./CRW/outs/possorted_bam.bam

```
## We further ran arks on the output from the above script.
Script for arks

```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 16
##SBATCH --qos=debug
#SBATCH --job-name arks.crw
#SBATCH --mem=50G
#SBATCH --time=72:00:00
##SBATCH --time=00:15:00
#SBATCH --account=uoo02772
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=katma889@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load BWA/0.7.17-GCC-9.2.0
module load SAMtools/1.13-GCC-9.2.0
module load BEDTools/2.29.2-GCC-9.2.0
module load LINKS/1.8.7-GCC-9.2.0
export PATH=/nesi/nobackup/uoo02752/nematode/nematode_nanopore/0.all_fast5/gupppy.5/pycoqc/nanolyse/porechop/nanoqc/flye/M.neg_flye/purgehaplotigs/lrscaff/scaffolds1/scaffolds2/scaffolds3/scaffolds4/scaffolds5/lrgapcloser/Output/rails.cobler/ragtag/ragtag_output/arbitr/arbitr.default/arbitr.2/arcs.links/arks-1.0.4/Examples:$PATH
export PATH=/nesi/nobackup/uoo02752/nematode/nematode_nanopore/0.all_fast5/gupppy.5/pycoqc/nanolyse/porechop/nanoqc/flye/M.neg_flye/purgehaplotigs/lrscaff/scaffolds1/scaffolds2/scaffolds3/scaffolds4/scaffolds5/lrgapcloser/Output/rails.cobler/ragtag/ragtag_output/arbitr/arbitr.default/arbitr.2/arcs.links/arks-1.0.4/Arks:$PATH

arks-make arks draft=output.arbitr.scaffolds reads=barcoded threads=16

```
### Then we ran arks version 1.1.0 to scaffold the genome assemblies produced by arbirr using our 10X Chromium Genomics reads. The script is given below;  

Scripts for arks

```
(base) [katma889@mahuika01 arks]$ less arks.sl

#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 16
##SBATCH --qos=debug
#SBATCH --job-name arks.crw
#SBATCH --mem=50G
#SBATCH --time=72:00:00
##SBATCH --time=00:15:00
#SBATCH --account=uoo02772
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=katma889@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load BWA/0.7.17-GCC-9.2.0
module load SAMtools/1.13-GCC-9.2.0
module load BEDTools/2.29.2-GCC-9.2.0
module load LINKS/1.8.7-GCC-9.2.0
export PATH=/nesi/nobackup/uoo02752/nematode/nematode_nanopore/0.all_fast5/gupppy.5/pycoqc/nanolyse/porechop/nanoqc/flye/M.neg_flye/purgehaplotigs/lrscaff/scaffolds1/scaffolds2/scaffolds3/scaffolds4/scaffolds5/lrgapcloser/Output/rails.cobler/ragtag/ragtag_output/arbitr/arbitr.default/arbitr.2/arcs.links/arks-1.0.4/Examples:$PATH
export PATH=/nesi/nobackup/uoo02752/nematode/nematode_nanopore/0.all_fast5/gupppy.5/pycoqc/nanolyse/porechop/nanoqc/flye/M.neg_flye/purgehaplotigs/lrscaff/scaffolds1/scaffolds2/scaffolds3/scaffolds4/scaffolds5/lrgapcloser/Output/rails.cobler/ragtag/ragtag_output/arbitr/arbitr.default/arbitr.2/arcs.links/arks-1.0.4/Arks:$PATH

arks-make arks draft=output.arbitr.scaffolds reads=barcoded threads=16

```
