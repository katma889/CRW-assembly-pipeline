# CRW-assembly-pipeline
## we sequenced the clover root weevil( Sitona obsoletus) using the data from 4 flow cell runs with 4 different individual and the data from all 4 runs were combined and basecalled with the guppy 5 version.
script
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
