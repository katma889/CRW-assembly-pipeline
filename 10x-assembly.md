## Assembling linked reads from 10X Genomics

For Sitona obsoletus whole genome project, we combinated the linked read data from 10x genomics obtained from sequencing in Novaseq and Hiseq.
As the first sequencing resulted only about 36X coverage of the eastimated genome size of this weevil, therefore the same libraries (leftover) were resequenced in HiSEQ to obtain additional data.
Combining both data resulted in a data of over 60X times coverage of the estimated genome of the weevil.   
The script for 'Supernova v.2.1.1' which I used for assembling linked reads is given below:

```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=1
#SBATCH --job-name Supernova_crw
#SBATCH --mem=400G
#SBATCH --time=168:00:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=katma889@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load Supernova/2.1.1

supernova run --id=crw_10xSN \
              --fastqs=/path/to/linked-read/fastq/files 
             
```
Upon completing the assembly which resulted in fasta file was produced to 'pseudohap'style with the script below:

```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task=2
#SBATCH --ntasks=1
#SBATCH --job-name Supernova_crw_fasta
#SBATCH --mem=1G
#SBATCH --time=00:10:00
#SBATCH --account=uoo02752
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=katma889@student.otago.ac.nz
#SBATCH --hint=nomultithread

module load Supernova/2.1.1

supernova mkoutput \
        --style=pseudohap \
        --asmdir=path/to/working-directory-of-supernova/crw_10xSN/outs/assembly \      # this is a path to the output directory `assembly` created by supernova
        --outprefix=./crw_10xSN
```
