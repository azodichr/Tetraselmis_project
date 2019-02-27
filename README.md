# Tetraselmis Genome Assembly & Annotation
**Table of Contents**

1. [Genome Assembly](#genome-assembly)
2. [Genome Polishing](#genome-polishing)
3. [Genome Annotation](#genome-annotation)
4. [BUSCO Analysis](#busco-analysis)
5. [Expression level time course](#expression-level-time-course)


# Genome Assembly
Completed by Nick Panchy


# Genome Polishing

**Plan: Polish genome 3 times with Pilon using 10xGenomics short reads***

**Back up 10xGenomics data**
Raw data: /data/scratch/azodi/TetraselmisProject/10xGenomics_RawData/fastq_path/HTCLVBBXX/Tetrastriata/Tetrastriata_S1_L006_*_001.fastq.gz
Barcodes Trimmed: /data/scratch/azodi/TetraselmisProject/10xGenomics_BarcodeTrimmed/barcoded.fastq.gz
Genome saved on calculon: /data/scratch/panchy/TetraselmisAssembly/Tetra664_04212017/unitigging/5.3-concensus_dg24_db24_gsize380000000_el499/tetra.contigs.fasta

## Round 1 

### Align the transcripts to the genome using Bowtie

```
bowtie2-build tetra.contigs.fasta tetra.contigs.fasta.bt
sbatch job_bowtie_round1 
```
Example bowtie jobs:
```
bowtie2 -x tetra.contigs.fasta.bt --phred33 -S tetra.contigs_10x_1.sam -q -U barcoded.fastq.gz
python ~shius/codes/qsub_hpc.py -f submit -u azodichr -c job_bowtie_round2.sh -w 100:40:00 -m 100 -p 8 -J bowtie2 -mo bowtie2 -wd /mnt/scratch/azodichr/10xgenom_2232_20180405/01_bam/
```

Example sort job:
```
#PBS -q main
#PBS -l walltime=12:00:00,mem=100gb
#PBS -N sort_samtools
#PBS -d /mnt/scratch/azodichr/10xgenom_2232_20180405/01_bam/

module load samtools
cd /mnt/scratch/azodichr/10xgenom_2232_20180405/01_bam/
samtools view -Su tetra.contigs.10x_2.sam | samtools sort - tetra.contigs.10x_2.sam.sorted
samtools index tetra.contigs.10x_2.sam.sorted.bam
```

### Polish genome using Pilon
submit: job_pilon_1.sh

Example pilon job:
```
#!/bin/bash -login

#PBS -q main
#PBS -l walltime=3:59:00,mem=180gb
#PBS -N pilon
#PBS -d /mnt/scratch/azodichr/10xgenom_2232_20180405/01_bam/

module load Java/1.7.0_51
cd /mnt/scratch/azodichr/10xgenom_2232_20180405/01_bam/
java -Xmx160G -jar pilon-1.21.jar --genome 'tetra.contigs.fasta' --bam 'tetra.contigs_10x_1.sam.sorted.bam' --diploid --outdir 'Tet_10x_pilon_round1' --output tetra.contigs.10x.pilon_1.fa --threads 32 --debug
```

## Round 2 

```
bowtie2-build tetra.contigs.10x.pilon_1.fa tetra.contigs.10x.pilon_1.fa.bt
sbatch job_bowtie_round2.sh
sbatch job_samtools_2.sh
sbatch job_pilon_1.sh
```

## Round 3 
```
bowtie2-build tetra.contigs.10x.pilon_2.fa.fasta tetra.contigs.10x.pilon_2.fa.fasta.bt
sbatch bowtie3.sb
sbatch bowtie_sort3.sb
# Note: 77.41% overall alignment rate*
sbatch job_pilon_3.sh
```

## Get Stats from 3 rounds of pilon 

Pull summary of corrected small events out of run file
```
grep 'Corrected ' pilon.o60971579 > pilon_stats_smallevent_1
grep 'Corrected ' pilon2.o61967548 > pilon_stats_smallevent_2
grep 'Corrected ' slurm-4765115.out > pilon_stats_smallevent_3
```

Rename in R
```
r <- read.csv('pilon_stats_smallevent_3', header=F, sep=' ')
r1 <- r[c('V2','V5','V9','V11','V15')]
names(r1) <- c('SNPs','SIn','SIn_bp','SDel','SDel_bp')
print(colSums(r1))
```

Pull Number of gaps closed and breaks fixed.
```
grep 'ClosedGap' pilon.o60971579 | wc -l 
grep 'PartialFill' pilon.o60971579 | wc -l 
grep 'BreakFix' pilon.o60971579 | wc -l 

grep 'ClosedGap' pilon2.o61967548 | wc -l 
grep 'PartialFill' pilon2.o61967548 | wc -l 
grep 'BreakFix' pilon2.o61967548 | wc -l 

grep 'ClosedGap' slurm-4765115.out | wc -l 
grep 'PartialFill' slurm-4765115.out | wc -l 
grep 'BreakFix' slurm-4765115.out | wc -l 
```



# Genome Annotation

This pipeline was developed based on the [Augustus retraining protocol](http://augustus.gobics.de/binaries/retraining.html) and from [this GitHub post](https://www.biostars.org/p/261203/)
## 1. Initial MAKER Analysis
*See Expression level time course section for details on QC and trimming RNA-Seq data

Needed input: Genome (*/mnt/home/azodichr/02_Tetraselmis/01_FinalDrafts/tetra.contigs.10x.pilon_3.fasta*), [EST Evidence](#est-evidence), and [protein homology evidence](#protein-homology-evidence)




### EST Evidence
wkdir: /mnt/gs18/scratch/users/azodichr/10xgenom_2232_20180405/02_Assembly/

**Set up alignment using HISAT2**

```module load GCC/6.4.0-2.28  
module load OpenMPI/2.1.2
module load GMAP-GSNAP/2018-05-11
gmap_build -d Tet_Pilon3_genome -D . -k 15 /mnt/gs18/scratch/users/azodichr/10xgenom_2232_20180405/01_bam/Tet_10x_pilon_round3/tetra.contigs.10x.pilon_3.fa.fasta

module load GCC/4.7.2
module load hisat2
hisat2-build /mnt/gs18/scratch/users/azodichr/10xgenom_2232_20180405/01_bam/Tet_10x_pilon_round3/tetra.contigs.10x.pilon_3.fa.fasta tetra.contigs.10xpilon3

module load GCC/6.4.0-2.28  
module load OpenMPI/2.1.2
module load GMAP-GSNAP/2018-05-11
gmap -d tetra.contigs.10xpilon3 -D . -B 5 -A

ln -s ~/02_Tetraselmis/02_TimeCourse_Transcriptomics/01_DataProcessing/02_Trim/h*E .
```

**Align the trimmed RNA-Seq reads to the genome**

```
module load GCC/4.7.2
module load hisat2
hisat2 -x tetra.contigs.10xpilon3 --phred33 -q -S Tet_transcriptome.sam -1 h0_e1_R1_PE,h0_e2_R1_PE,h0_e3_R1_PE,h12_e1_R1_PE,h12_e2_R1_PE,h12_e3_R1_PE,h168_e1_R1_PE,h168_e2_R1_PE,h168_e3_R1_PE,h1_e1_R1_PE,h1_e2_R1_PE,h1_e3_R1_PE,h24_e1_R1_PE,h24_e2_R1_PE,h24_e3_R1_PE,h2_e1_R1_PE,h2_e2_R1_PE,h2_e3_R1_PE,h48_e1_R1_PE,h48_e2_R1_PE,h48_e3_R1_PE,h4_e1_R1_PE,h4_e2_R1_PE,h4_e3_R1_PE,h6_e1_R1_PE,h6_e2_R1_PE,h6_e3_R1_PE,h96_e1_R1_PE,h96_e2_R1_PE,h96_e3_R1_PE -2 h0_e1_R2_PE,h0_e2_R2_PE,h0_e3_R2_PE,h12_e1_R2_PE,h12_e2_R2_PE,h12_e3_R2_PE,h168_e1_R2_PE,h168_e2_R2_PE,h168_e3_R2_PE,h1_e1_R2_PE,h1_e2_R2_PE,h1_e3_R2_PE,h24_e1_R2_PE,h24_e2_R2_PE,h24_e3_R2_PE,h2_e1_R2_PE,h2_e2_R2_PE,h2_e3_R2_PE,h48_e1_R2_PE,h48_e2_R2_PE,h48_e3_R2_PE,h4_e1_R2_PE,h4_e2_R2_PE,h4_e3_R2_PE,h6_e1_R2_PE,h6_e2_R2_PE,h6_e3_R2_PE,h96_e1_R2_PE,h96_e2_R2_PE,h96_e3_R2_PE -U h0_e1_R1_SE,h0_e1_R2_SE,h0_e2_R1_SE,h0_e2_R2_SE,h0_e3_R1_SE,h0_e3_R2_SE,h12_e1_R1_SE,h12_e1_R2_SE,h12_e2_R1_SE,h12_e2_R2_SE,h12_e3_R1_SE,h12_e3_R2_SE,h168_e1_R1_SE,h168_e1_R2_SE,h168_e2_R1_SE,h168_e2_R2_SE,h168_e3_R1_SE,h168_e3_R2_SE,h1_e1_R1_SE,h1_e1_R2_SE,h1_e2_R1_SE,h1_e2_R2_SE,h1_e3_R1_SE,h1_e3_R2_SE,h24_e1_R1_SE,h24_e1_R2_SE,h24_e2_R1_SE,h24_e2_R2_SE,h24_e3_R1_SE,h24_e3_R2_SE,h2_e1_R1_SE,h2_e1_R2_SE,h2_e2_R1_SE,h2_e2_R2_SE,h2_e3_R1_SE,h2_e3_R2_SE,h48_e1_R1_SE,h48_e1_R2_SE,h48_e2_R1_SE,h48_e2_R2_SE,h48_e3_R1_SE,h48_e3_R2_SE,h4_e1_R1_SE,h4_e1_R2_SE,h4_e2_R1_SE,h4_e2_R2_SE,h4_e3_R1_SE,h4_e3_R2_SE,h6_e1_R1_SE,h6_e1_R2_SE,h6_e2_R1_SE,h6_e2_R2_SE,h6_e3_R1_SE,h6_e3_R2_SE,h96_e1_R1_SE,h96_e1_R2_SE,h96_e2_R1_SE,h96_e2_R2_SE,h96_e3_R1_SE,h96_e3_R2_SE

sbatch job_hisat2.sb

mv slurm-5404049.out tetra.contigs.10xpilon3_hisat2.stats
sbatch job_samtools.sb
```
*NOTE: 99.27% overall alignment rate !!!!!!!!!!!*


**Run Trinity Guided Transcriptome assembly**

```
Trinity --seqType fq --max_memory 100G --output /mnt/gs18/scratch/users/azodichr/10xgenom_2232_20180405/02_Assembly/trinity_181228 --CPU 10 --left h0_e1_R1_PE,h0_e2_R1_PE,h0_e3_R1_PE,h12_e1_R1_PE,h12_e2_R1_PE,h12_e3_R1_PE,h168_e1_R1_PE,h168_e2_R1_PE,h168_e3_R1_PE,h1_e1_R1_PE,h1_e2_R1_PE,h1_e3_R1_PE,h24_e1_R1_PE,h24_e2_R1_PE,h24_e3_R1_PE,h2_e1_R1_PE,h2_e2_R1_PE,h2_e3_R1_PE,h48_e1_R1_PE,h48_e2_R1_PE,h48_e3_R1_PE,h4_e1_R1_PE,h4_e2_R1_PE,h4_e3_R1_PE,h6_e1_R1_PE,h6_e2_R1_PE,h6_e3_R1_PE,h96_e1_R1_PE,h96_e2_R1_PE,h96_e3_R1_PE --right h0_e1_R2_PE,h0_e2_R2_PE,h0_e3_R2_PE,h12_e1_R2_PE,h12_e2_R2_PE,h12_e3_R2_PE,h168_e1_R2_PE,h168_e2_R2_PE,h168_e3_R2_PE,h1_e1_R2_PE,h1_e2_R2_PE,h1_e3_R2_PE,h24_e1_R2_PE,h24_e2_R2_PE,h24_e3_R2_PE,h2_e1_R2_PE,h2_e2_R2_PE,h2_e3_R2_PE,h48_e1_R2_PE,h48_e2_R2_PE,h48_e3_R2_PE,h4_e1_R2_PE,h4_e2_R2_PE,h4_e3_R2_PE,h6_e1_R2_PE,h6_e2_R2_PE,h6_e3_R2_PE,h96_e1_R2_PE,h96_e2_R2_PE,h96_e3_R2_PE --KMER_SIZE 27 --min_kmer_cov 2 --full_cleanup

sbatch job_trinity.sb  # Took ~18 hours

module purge
module load intel/2017b Trinity/2.6.6
/opt/software/Trinity/2.6.6/util/TrinityStats.pl trinity_181228.Trinity.fasta > trinity_181228.Trinity.fasta.stats
```

### Protein Homology Evidence

**CEG (Core Eukaryotic Genes)**
http://korflab.ucdavis.edu/Datasets/genome_completeness/index.html#SCT2
```
wget http://korflab.ucdavis.edu/Datasets/genome_completeness/core/248.prots.fa.gz
mv 248.prots.fa.gz CEG.fa
```

**GreenCut Genes**
GreenCut list from Ben (/mnt/home/azodichr/02_Tetraselmis/06_MAKER/GreenCut2JBC2011.xls)
Arabidopsis: pull Arabidopsis IDs and get longest peptide: 
```python ~shius/codes/FastaManager.py -f getseq2 -fasta ~/Sequences/Arabidopsis/TAIR10_pep_20101214_updated.txt.longest.mod.fa -name greencut_at```

Chlammy: Pull JGI 3 names and convert them to JGI v5.5
http://pathways.mcdb.ucla.edu/algal/id_conversion.html

Select all non plastid CEG and GreenCut genes:
```
python ~shius/codes/FastaManager.py -f getseq2 -name greencut_noPlastid_cd -fasta ~/Sequences/Chlammy/Creinhardtii_281_v5.5.protein_primaryTranscriptOnly.fa.mod.fa
cat CEG.fa greencut_noPlastid_cd.fa > annot_proteins.fa
```

** Add chlorella repiprocal best matches to Chlammy CEG and Greencut (non-plastid) genes **
*Note - also tried to add ostreococcus genes, but only one hit reciprocal best match hit, so just using chlorella*
Downloaded from JGI
```
wd: /mnt/home/azodichr/02_Tetraselmis/06_MAKER/process_training_proteins
module purge
module load BLAST/2.2.26-Linux_x86_64

formatdb -i Chlorella_NC64A.best_proteins.fasta -p T
formatdb -i Ostta4221_3_GeneCatalog_proteins_20161028.aa.fasta -p T
formatdb -i greencut_noPlastid_cd.fa -p T

blastall -p blastp -i greencut_noPlastid_cd.fa -d Ostta4221_3_GeneCatalog_proteins_20161028.aa.fasta -e 0.001 -m 8 -o Os_greencut_hits.txt
blastall -p blastp -i greencut_noPlastid_cd.fa -d Chlorella_NC64A.best_proteins.fasta -e 0.001 -m 8 -o Cnc64A_greencut_hits.txt

python ~/GitHub/Utilities/ParseBlast.py -f get_reciprocal -blast Os_greencut_hits.txt
python ~/GitHub/Utilities/ParseBlast.py -f get_reciprocal -blast Cnc64A_greencut_hits.txt

cat Cnc64A_greencut_hits.txt.recip ../CEG.fa greencut_noPlastid_cd.fa > proteins_CEG_GCcr_GCc64.fa
```
**1813 proteins to use as protein homology evidence**


### Run MAKER
*maker_opts.ctl*
est2genome=1
protein2genome=1

wd: 06_MAKER/09_Pilon10x3_prot_3species/


## 2. Train gene prediction models (Augustus and SNAP)
### SNAP
Export 'confident' gene models from MAKER and rename to something meaningful
```
module purge
module load GCC/7.2.0-2.29 MPICH/3.2.1
module load maker/2.31.9
maker2zff -x 0.25 -l 50 -d ../09_Pilon10x3_prot_3species/tetra.contigs.10x.pilon_3.fa.maker.output/tetra.contigs.10x.pilon_3.fa_master_datastore_index.log
for i in *; do mv $i $(echo $i | sed 's/genome/tetra10xpilonx3.zff.length50_aed0.25/'); done
```

Gather stats, validate, and collect training sequences and annotations plus 1kb surrounding for training
```
fathom tetra10xpilonx3.zff.length50_aed0.25.ann tetra10xpilonx3.zff.length50_aed0.25.dna -gene-stats > gene-stats.log
fathom tetra10xpilonx3.zff.length50_aed0.25.ann tetra10xpilonx3.zff.length50_aed0.25.dna -validate > validate.log
fathom tetra10xpilonx3.zff.length50_aed0.25.ann tetra10xpilonx3.zff.length50_aed0.25.dna -categorize 1000 > categorize.log
fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log
```

Create the training parameters and assemble the HMM
```
mkdir params
cd params
forge ../export.ann ../export.dna > ../forge.log
cd ../
hmm-assembler.pl tetra10xpilonx3.zff.length50_aed0.25 params > tetra10xpilonx3.zff.length50_aed0.25.hmm
```


### Augustus

#### Refine and format set of annotated proteins to train and test on
A. Remove genes from ab initio MAKER round 1 that have >90% amino acid similarity (over 90% of the length of the gene) with other genes (cause redundant genes can cause overfitting).
```
grep -P "maker\tgene" tetra.contigs.10x.pilon_3.fa.all.gff > tetra.contigs.10x.pilon_3.fa.all.gff_genesOnly
python ~/GitHub/Utilities/FastaManager.py -f gff_to_coord -gff tetra.contigs.10x.pilon_3.fa.all.gff_genesOnly
python ~/GitHub/Utilities/FastaManager.py -f get_stretch4 -coords tetra.contigs.10x.pilon_3.fa.all.gff_genesOnly.coord -fasta ../tetra.contigs.10x.pilon_3.fasta

module purge
module load BLAST/2.2.26-Linux_x86_64
formatdb -i tetra.contigs.10x.pilon_3.fa.all.gff_genesOnly.coord.fa -p F
blastall -p blastn -d tetra.contigs.10x.pilon_3.fa.all.gff_genesOnly.coord.fa -i tetra.contigs.10x.pilon_3.fa.all.gff_genesOnly.coord.fa -o blastnE00001 -e 0.00001 -a 8 -m 8
python ../../07_trin_10xpilonx3_2/augustus_retrain_strict/filter_genes_1.py -b blastnE00001 -f tetra.contigs.10x.pilon_3.fa.all.gff_genesOnly.coord.fa
python ~/GitHub/Utilities/FastaManager.py -f getseq2 -fasta tetra.contigs.10x.pilon_3.fa.all.gff_genesOnly.coord.fa -name tetra_genes_DupsRemoved.txt
```
*Note: Went from 29,962 genes for training to 7,949 genes *

B. Remove genes that aren't hits in GreenCut

formatdb -i tetra_genes_DupsRemoved.txt.fa -p F
```
formatdb -i ../proteins_CEG_GCcr_GCc64.fa -p T
blastall -p blastx -d ../proteins_CEG_GCcr_GCc64.fa -i tetra_genes_DupsRemoved.txt.fa -o blastx_greencut -e 0.01 -a 8 -m 8
python ../../07_trin_10xpilonx3_2/augustus_retrain_strict/filter_genes_2.py -b blastx_greencut -f1 tetra_genes_DupsRemoved.txt
python ../../07_trin_10xpilonx3_2/augustus_retrain_strict/filter_genes_3.py -gff tetra.contigs.10x.pilon_3.fa.all.gff -f tetra_genes_DupsRemoved_GreenCutHits.txt
```
*Note: Went from 7,949 genes for training to 983 genes*

C. Convert gff to genbank and split into training and testing (90/10)
```
perl ~/GitHub/Augustus/scripts/gff2gbSmallDNA.pl tetra.contigs.10x.pilon_3.fa.all.gff.filtered ~/02_Tetraselmis/06_MAKER/07_trin_10xpilonx3_2/tetra.contigs.10x.pilon_3.fasta 1000 tetra.contigs.10x.pilon_3.fa.all.gff.filtered.genebank
perl ~/GitHub/Augustus/scripts/randomSplit.pl tetra.contigs.10x.pilon_3.fa.all.gff.filtered.genebank 98
```
*Note 1000 is the max-size-gene-flanking (i.e. it will grab 1kb up and downstream of the gene to add to the genebank file)


#### Optimize hyper-parameters and then use to train/test augustus gene predictor
Prep configuration directory for tetraselmis
```
module load icc/2017.1.132-GCC-6.3.0-2.27  impi/2017.1.132
module load augustus
cd /mnt/home/azodichr/GitHub/Augustus/config/species/
mkdir tetraselmis_2sp
cp generic/* tetraselmis_2sp/.
cd tetraselmis_2sp/
for i in *; do mv $i $(echo $i | sed 's/generic/tetraselmis_2sp/'); done
sed -i 's/generic_/tetraselmis_2sp_/g' tetraselmis_2sp_parameters.cfg
```

Optimize & train using optimal parameters
```
optimize_augustus.pl --species=tetraselmis_2sp tetra.contigs.10x.pilon_3.fa.all.gff.filtered.genebank.train --metapars=/mnt/home/azodichr/GitHub/Augustus/config/species/tetraselmis_2sp/tetraselmis_2sp_metapars.cfg --AUGUSTUS_CONFIG_PATH=/mnt/home/azodichr/GitHub/Augustus/config/
```

**Training Results**

|Type|Sensitivity	|Specificity|
|---	|---	|---	|
|Nucleotide level|0.87|0.69|
|Transcript|0.104|0.082|

Apply to testing set
```
augustus --species=tetraselmis_2sp --AUGUSTUS_CONFIG_PATH=/mnt/home/azodichr/GitHub/Augustus/config/ tetra.contigs.10x.pilon_3.fa.all.gff.filtered.genebank.test > test_results.txt
```
**Testing Results**

|Type|Sensitivity|Specificity|
|---	|---	|---	|
|Nucleotide level|0.88 |0.71|
|Transcript|0.13|0.11|



## 3. Run MAKER with ab initio gene predictors
**Recycle the mapping of empicial evidence we have from the first MAKER round, so we don't have to perform all the BLASTs, etc. again**
```
awk '{ if ($2 == "est2genome") print $0 }' ../09_Pilon10x3_prot_3species/tetra.contigs.10x.pilon_3.fa.maker.output/tetra.contigs.10x.pilon_3.fa.all.gff > tetra10xpilonx3.maker.est2genome.gff
awk '{ if ($2 == "protein2genome") print $0 }' ../09_Pilon10x3_prot_3species/tetra.contigs.10x.pilon_3.fa.maker.output/tetra.contigs.10x.pilon_3.fa.all.gff > tetra10xpilonx3.maker.protein2genome.gff
awk '{ if ($2 ~ "repeat") print $0 }' ../09_Pilon10x3_prot_3species/tetra.contigs.10x.pilon_3.fa.maker.output/tetra.contigs.10x.pilon_3.fa.all.gff > tetra10xpilonx3.maker.repeats.gff
```
Parameters to adjust in the maker_opts.ctl file:
-est_gff=/mnt/home/azodichr/02_Tetraselmis/06_MAKER/11_Round2/tetra10xpilonx3.maker.est2genome.gff #aligned ESTs or mRNA-seq from an external GFF3 file
-protein_gff=/mnt/home/azodichr/02_Tetraselmis/06_MAKER/11_Round2/tetra10xpilonx3.maker.est2genome.gff  #aligned protein homology evidence from an external GFF3 file
-rm_gff=/mnt/home/azodichr/02_Tetraselmis/06_MAKER/11_Round2/tetra10xpilonx3.maker.repeats.gff #pre-identified repeat elements from an external GFF3 file
-snaphmm=/mnt/home/azodichr/02_Tetraselmis/06_MAKER/10_SNAP/tetra10xpilonx3.zff.length50_aed0.25.hmm #SNAP HMM file
-augustus_species=tetraselmis_2sp #Augustus gene prediction species model
-est2genome=0
-protein2genome=0


## 4. Iteratively run MAKER?



## 5. Assess quality using BUSCO 




# BUSCO Analysis



# Expression level time course

*Experiment Description:* Peter and Jake grew Tet in bioreactors and extracted mRNA at 0, 1, 2, 4, 6, 12, 24, and 48 hours. PE sequencing with two replicates each (different bioreactors). The goals are to look for cycling genes and use this data to aid in genome annotation once we get the PacBio data.

*Sequencing Details (from Kevin Carr):* Sequencing is complete for samples submitted to the RTSF Genomics Core, project NEO3896 (Tetraselmis 66_4). You submitted thirty (30) samples of total RNA for NGS library prep and sequencing.  Libraries were prepared using the Illumina TruSeq Stranded mRNA Library Preparation Kit. Completed libraries were QC'd using Qubit dsDNA HS, Caliper LabChipGX HS DNA and Kapa Biosystems Illumina Library Quantification qPCR assays. After quantitation, libraries were pooled in roughly equimolar amounts, 3 pools of 10 libraries each and each pool was loaded on one (1) lane of an Illumina HiSeq 2500 High Output flow cell (v2). Sequencing was carried out using HiSeq SBS reagents in a 2x125bp paired end format (PE125). Base calling was done by Illumina Real Time Analysis (RTA) v1.18.64 and output of RTA was demultiplexed and converted to FastQ format with Illumina Bcl2fastq v1.8.4. A summary of the output can be found on HPC at /mnt/home/azodichr/05_Tetraselmis/02_TimeCourse_Transcriptomics/20160527_SeqProduction_Kramer.xlsx


HPC LOCATION: /mnt/home/azodichr/05_Tetraselmis/02_TimeCourse_Transcriptomics/
Temporary storage: /mnt/scratch/azodichr/20160527_mRNASeq_PE/
Backup storage on calculon2: /home/azodi/RawData/20160527_mRNASeq_PE_Tetrasel


### Quality control and Trimming

**QC**
```
python ~shius/codes/qsub_hpc.py -f queue -u azodichr -c run_fastQC_1st.txt -w 230 -m 4 -n 200 -wd /mnt/home/azodichr/02_Tetraselmis/02_TimeCourse_Transcriptomics/01_DataProcessing/01_Fastq/ -mo fastqc
mv *.zip ../01_Fastq
mv *.html ../01_Fastq
```
*Note: GC Content graphs are bimodal. We presume the peak with the lower %GC content is made up of plastid sequence, and the higher GC content is the nuclear sequence. The decrease in the size of the lower %GC peak over time is consistent with Peter’s quality check results which show that over time, the plastids were degraded and less of that DNA was sequenced.*

**Trimming**
```
python ~shius/codes/qsub_hpc.py -f queue -u azodichr -c run_trimo_1.txt -w 30 -m 20 -n 200 -wd /mnt/home/azodichr/02_Tetraselmis/02_TimeCourse_Transcriptomics/01_DataProcessing/00_RawData/ -mo Trimmomatic
mv h* ../02_Trim
```

This will perform the following: [1] removes adapters (ILLUMINACLIP:TruSeq3-PE.fa:2:30:10), [2] removes leading low quality or N bases (below quality 3) (LEADING:3), [3] removes trailing low quality or N bases (below quality 3) (TRAILING:3), [4] scans the read with a 4-base wide sliding window, cutting when the average quality per base drops below 30 (SLIDINGWINDOW:4:30), [5] drops reads below the 36 bases long (MINLEN:36), and [6] cut 8 bases off the start of the read for the barcode (HEADCROP:8).

Example Command: 
```java -jar $TRIM/trimmomatic PE Tetra_h0_e1_GAGATTCC-GGCTCTGA_L006_R1_001.fastq.gz Tetra_h0_e1_GAGATTCC-GGCTCTGA_L006_R2_001.fastq.gz h0_e1_R1_PE h0_e1_R1_SE h0_e1_R2_PE h0_e1_R2_SE ILLUMINACLIP:$ADAPTOR:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:36 HEADCROP:8```

*Output:* 4 files for each sequence: R1_PE and R2_PE are for the paired output where both reads survived the processing, and R1_SE & R2_SE for the corresponding unpaired output where a read survived, but the partner read did not.


**QC after trimming**
```
python ~shius/codes/qsub_hpc.py -f queue -u azodichr -c run_fastQC_2nd.txt -w 60 -m 4 -n 200 -wd /mnt/home/azodichr/02_Tetraselmis/02_TimeCourse_Transcriptomics/01_DataProcessing/03_Fastq/ -mo fastqc
```

















######### BUSCO ##########

module load python3

!! If running in genome mode (--mode geno) need to set the paths to augustus
# export BUSCO_CONFIG_FILE="/mnt/home/azodichr/02_Tetraselmis/07_BUSCO/busco/config/config.ini"

export PATH="/mnt/home/azodichr/LocalPrograms/augustus.2.5.5/bin:$PATH"
export PATH="/mnt/home/azodichr/LocalPrograms/augustus.2.5.5/scripts:$PATH"
export AUGUSTUS_CONFIG_PATH="/mnt/home/azodichr/LocalPrograms/augustus.2.5.5/config/"


# Genome 
wd: /mnt/home/azodichr/02_Tetraselmis/07_BUSCO/00_Genome
python ../busco/scripts/run_BUSCO.py --i tetra.unitigs.fasta --out tetra_genome --lineage_path ../busco/lineages/eukaryota_odb9/ --mode geno
python ../busco/scripts/run_BUSCO.py --i tetra.unitigs.pilon.fasta --out tetra_pilon_genome --lineage_path ../busco/lineages/eukaryota_odb9/ --mode geno



# Protein sequences from MAKER annotation
wd: /mnt/home/azodichr/02_Tetraselmis/07_BUSCO/01_peptides
ln -s ../../06_MAKER/01_trinity_denovo/tetra.unitigs.maker.output/tetra.unitigs.all.maker.proteins.fasta
python ../busco/scripts/run_BUSCO.py --i tetra.unitigs.all.maker.proteins.fasta --out tetra_prot_TrinDeNovo --lineage_path ../busco/lineages/eukaryota_odb9/ --mode prot

ln -s ../../06_MAKER/02_trinity_guided/01_FromHISAT2/tetra.unitigs.maker.output/tetra.unitigs.all.maker.proteins.fasta tetra.unitigs.all.maker.TrinityHISAT.proteins.fasta
python ../busco/scripts/run_BUSCO.py --i tetra.unitigs.all.maker.TrinityHISAT.proteins.fasta --out tetra_prot_TrinGuid --lineage_path ../busco/lineages/eukaryota_odb9/ --mode prot

ln -s ../../06_MAKER/03_stringtie/01_FromHISAT2/tetra.unitigs.maker.output/tetra.unitigs.all.maker.proteins.fasta tetra.unitigs.all.maker.StringTieHISAT.proteins.fasta
python ../busco/scripts/run_BUSCO.py --i tetra.unitigs.all.maker.StringTieHISAT.proteins.fasta --out tetra_prot_StringTie --lineage_path ../busco/lineages/eukaryota_odb9/ --mode prot

ln -s ~/02_Tetraselmis/06_MAKER/04_trinity_pilon/tetra.unitigs.pilon.maker.output/tetra.unitigs.pilon.all.maker.proteins.fasta .
python ../busco/scripts/run_BUSCO.py --i tetra.unitigs.pilon.all.maker.proteins.fasta --out tetra_prot_TrinGuid_Pilon --lineage_path ../busco/lineages/eukaryota_odb9/ --mode prot

ln -s ../../06_MAKER/04_trinity_pilon2/tetra.unitigs.pilon.maker.output/tetra.unitigs.pilon.all.maker.proteins.fasta .
python ../busco/scripts/run_BUSCO.py --i tetra.unitigs.pilon.all.maker.proteins.fasta --out tetra_prot_TrinGuid_Pilon2 --lineage_path ../busco/lineages/eukaryota_odb9/ --mode prot


# Transcriptomes (testing to see if we're just not annotating them well)
wd: /mnt/home/azodichr/02_Tetraselmis/07_BUSCO/mkdir 02_transcriptomes

ln -s ~/02_Tetraselmis/02_TimeCourse_Transcriptomics/01_DataProcessing/06_TrinityGuided/03_Pilon/Tet_transcriptome_TrinGuid_Pilon_HISAT.fasta .
ln -s ~/02_Tetraselmis/02_TimeCourse_Transcriptomics/01_DataProcessing/07_StringTie/01_FromHISAT2/Tet_transcr_stringtie.fa .
ln -s ~/02_Tetraselmis/02_TimeCourse_Transcriptomics/01_DataProcessing/06_TrinityGuided/trinity_out_dir/Tet_transcr_guidedHISAT_trinity.fasta .
ln -s ~/02_Tetraselmis/02_TimeCourse_Transcriptomics/01_DataProcessing/06_TrinityGuided/03_Pilon/trinity_out_dir/Tet_transcriptome_TrinGuidPilon2_HISAT.fasta .

python ../busco/scripts/run_BUSCO.py --i Tet_transcr_guidedHISAT_trinity.fasta --out tetra_trans_TrinGuid --lineage_path ../busco/lineages/eukaryota_odb9/ --mode geno
python ../busco/scripts/run_BUSCO.py --i Tet_transcr_stringtie.fa --out tetra_trans_StringTie --lineage_path ../busco/lineages/eukaryota_odb9/ --mode geno
python ../busco/scripts/run_BUSCO.py --i Tet_transcriptome_TrinGuid_Pilon_HISAT.fasta --out tetra_trans_TrinGuid_Pilon --lineage_path ../busco/lineages/eukaryota_odb9/ --mode geno
python ../busco/scripts/run_BUSCO.py --i Tet_transcriptome_TrinGuidPilon2_HISAT.fasta --out tetra_trans_TrinGuid_Pilon2 --lineage_path ../busco/lineages/eukaryota_odb9/ --mode geno



# Plot BUSCO summaries:
wd: /mnt/home/azodichr/02_Tetraselmis/07_BUSCO/BUSCO_summaries 
cp ../00_genome/run_tetra_genome/short_summary_tetra_genome.txt .
cp ../00_Genome/run_tetra_pilon_genome/short_summary_tetra_pilon_genome.txt .
cp ../01_peptides/run_tetra_prot_TrinDeNovo/short_summary_tetra_prot_TrinDeNovo.txt .
cp ../01_peptides/run_tetra_prot_TrinGuid/short_summary_tetra_prot_TrinGuid.txt .
cp ../01_peptides/run_tetra_prot_StringTie/short_summary_tetra_prot_StringTie.txt .
cp ../01_peptides/run_tetra_prot_TrinGuid_Pilon2/short_summary_tetra_prot_TrinGuid_Pilon2.txt .
cp ../02_transcriptomes/run_tetra_trans_TrinGuid_Pilon/short_summary_tetra_trans_TrinGuid_Pilon.txt .
cp ../02_transcriptomes/run_tetra_trans_TrinGuid_Pilon2/short_summary_tetra_trans_TrinGuid_Pilon2.txt .
cp ../02_transcriptomes/run_tetra_trans_trinity/short_summary_tetra_trans_trinity.txt .
cp ../02_transcriptomes/run_tetra_trans_StringTie/short_summary_tetra_trans_StringTie.txt .

python busco/scripts/generate_plot.py --working_directory BUSCO_summaries/



# Run busco on unassembled reads to determine if the genes are there and just not assembling!!
wd: /mnt/home/azodichr/02_Tetraselmis/07_BUSCO/03_unassembled

ln -s /mnt/scratch/panchyni/TetraFASTA/tetra.unassembled.fasta .
python ../busco/scripts/run_BUSCO.py --i tetra.unassembled.fasta --out tetra_unassembled --lineage_path ../busco/lineages/eukaryota_odb9/ --mode geno































################ Transcriptome Assembly - Annotation - BUSCO - Quantification #######################






###################### MAKER Annotation - ab inito gene model predictions (i.e. Round #1) ###################### 

wd: /mnt/home/azodichr/02_Tetraselmis/06_MAKER/07_trin_10xpilonx3_2
ln -s /mnt/gs18/scratch/users/azodichr/10xgenom_2232_20180405/01_bam/Tet_10x_pilon_round3/tetra.contigs.10x.pilon_3.fa.fasta .
ln -s /mnt/gs18/scratch/users/azodichr/10xgenom_2232_20180405/02_Assembly/trinity_181228.Trinity.fasta .
ln -s /mnt/gs18/scratch/users/azodichr/10xgenom_2232_20180405/02_Assembly/trinity_181228.Trinity.fasta.stats .
ln -s ../annot_proteins.fa .
cp ../04_trinity_pilon/run_maker .

module purge
module load GCC/7.2.0-2.29 MPICH/3.2.1 
module load maker/2.31.9

maker -CTL
# Specify file names and set est2gene=1 in the "maker_opts.ctl" file

sbatch run_maker


# Check results:
cd /mnt/home/azodichr/02_Tetraselmis/06_MAKER/04_trinity_pilon/tetra.unitigs.pilon.maker.output
# Check datastore_index.log to make sure each contig finished:
grep "FAILED" tetra.contigs.10x.pilon_3.fa_master_datastore_index.log | wc -l
grep "SKIPPED" tetra.contigs.10x.pilon_3.fa_master_datastore_index.log | wc -l
grep "RETRY" tetra.contigs.10x.pilon_3.fa_master_datastore_index.log | wc -l

fasta_merge -d tetra.contigs.10x.pilon_3.fa_master_datastore_index.log
gff3_merge -d tetra.contigs.10x.pilon_3.fa_master_datastore_index.log

# Get number of genes
grep ">" tetra.contigs.10x.pilon_3.fa.all.maker.proteins.fasta | wc -l
## 42,219

python ~shius/codes/FastaManager.py -f get_sizes -fasta tetra.contigs.10x.pilon_3.fa.all.maker.proteins.fasta



#####################  Retraining ab initio gene models using AUGUSTUS   #######################3#


Following instructions from: http://augustus.gobics.de/binaries/retraining.html

module load icc/2017.1.132-GCC-6.3.0-2.27  impi/2017.1.132
module load augustus

wd: /mnt/home/azodichr/02_Tetraselmis/06_MAKER/07_trin_10xpilonx3_2/augustus_retrain
ln -s ../tetra.contigs.10x.pilon_3.fa.maker.output_EST_PROT_2Gene/tetra.contigs.10x.pilon_3.fa.all.gff .

1. Compile training and testing genes

A. Remove genes from ab initio MAKER round 1 that have >90% amino acid similarity (over 90% of the length of the gene) with other genes (cause redundant genes can cause overfitting).
grep -P "maker\tgene" tetra.contigs.10x.pilon_3.fa.all.gff > tetra.contigs.10x.pilon_3.fa.all.gff_genesOnly
python ~/GitHub/Utilities/FastaManager.py -f gff_to_coord -gff tetra.contigs.10x.pilon_3.fa.all.gff_genesOnly
python ~/GitHub/Utilities/FastaManager.py -f get_stretch4 -coords tetra.contigs.10x.pilon_3.fa.all.gff_genesOnly.coord -fasta ../tetra.contigs.10x.pilon_3.fasta

module purge
module load BLAST/2.2.26-Linux_x86_64
formatdb -i tetra.contigs.10x.pilon_3.fa.all.gff_genesOnly.coord.fa -p F
blastall -p blastn -d tetra.contigs.10x.pilon_3.fa.all.gff_genesOnly.coord.fa -i tetra.contigs.10x.pilon_3.fa.all.gff_genesOnly.coord.fa -o blastnE00001 -e 0.00001 -a 8 -m 8
python ../../07_trin_10xpilonx3_2/augustus_retrain_strict/filter_genes_1.py -b blastnE00001 -f tetra.contigs.10x.pilon_3.fa.all.gff_genesOnly.coord.fa
python ~/GitHub/Utilities/FastaManager.py -f getseq2 -fasta tetra.contigs.10x.pilon_3.fa.all.gff_genesOnly.coord.fa -name tetra_genes_DupsRemoved.txt

*RESULTS*
* 07_trin_10xpilonx3_2: went from 21,562 genes in the original gff to XX genes after removing the ones with >90% sequence identity on at least 90% of the gene length. *
* 08_Pilon10x3_remake_prot: went from 53,375 genes for training to 13345 *
* 09_Pilon10x3_prot_3species: went from 29,962 genes for training to 7949 genes *



B. Remove genes that aren't hits in GreenCut 
formatdb -i tetra_genes_DupsRemoved.txt.fa -p F

*07*
ln -s ../../
formatdb -i ../proteins_CEG_GCcr_GCc64.fa -p T
blastall -p blastx -d ../proteins_CEG_GCcr_GCc64.fa -i tetra_genes_DupsRemoved.txt.fa -o blastx_greencut -e 0.01 -a 8 -m 8

*08*
formatdb -i ../greencut_noPlastid_cd.fa -p T
blastall -p blastx -d ../greencut_noPlastid_cd.fa -i tetra_genes_DupsRemoved.txt.fa -o blastx_greencut -e 0.01 -a 8 -m 8


*09*
formatdb -i ../proteins_CEG_GCcr_GCc64.fa -p T
blastall -p blastx -d ../proteins_CEG_GCcr_GCc64.fa -i tetra_genes_DupsRemoved.txt.fa -o blastx_greencut -e 0.01 -a 8 -m 8


python ../../07_trin_10xpilonx3_2/augustus_retrain_strict/filter_genes_2.py -b blastx_greencut -f1 tetra_genes_DupsRemoved.txt


C. Remove genes and accompanying gene model info from gff file if not in tetra_genes_DupsRemoved_GreenCutHits.txt

python ../../07_trin_10xpilonx3_2/augustus_retrain_strict/filter_genes_3.py -gff tetra.contigs.10x.pilon_3.fa.all.gff -f tetra_genes_DupsRemoved_GreenCutHits.txt


*RESULTS*
* 07_trin_10xpilonx3_2: went from XX genes to 750 genes. *
* 08_Pilon10x3_remake_prot: went from 13,345 genes for training to 2275 genes*
* 09_Pilon10x3_prot_3species: went from 7,951 genes for training to 983 genes *



# Convert gff into genbank format for Augustus
perl ~/GitHub/Augustus/scripts/gff2gbSmallDNA.pl tetra.contigs.10x.pilon_3.fa.all.gff.filtered ~/02_Tetraselmis/06_MAKER/07_trin_10xpilonx3_2/tetra.contigs.10x.pilon_3.fasta 1000 tetra.contigs.10x.pilon_3.fa.all.gff.filtered.genebank


Note 1000 is the max-size-gene-flanking (i.e. it will grab 1kb up and downstream of the gene to add to the genebank file)

# Split genbank file into training and testing (use 10% of genes for testing to start)

perl ~/GitHub/Augustus/scripts/randomSplit.pl tetra.contigs.10x.pilon_3.fa.all.filtered.genbank 227  # 08
perl ~/GitHub/Augustus/scripts/randomSplit.pl tetra.contigs.10x.pilon_3.fa.all.gff.filtered.genebank 98  # 09


2. Create and optimize meta parameters file for tetraselmis, train model, and test model
module load icc/2017.1.132-GCC-6.3.0-2.27  impi/2017.1.132
module load augustus

# 08
mkdir tetraselmis
cp generic/* tetraselmis/.
cd tetraselmis/
for i in *; do mv $i $(echo $i | sed 's/generic/tetraselmis/'); done
sed -i 's/generic_/tetraselmis_/g' tetraselmis_parameters.cfg

optimize_augustus.pl --species=tetraselmis tetra.contigs.10x.pilon_3.fa.all.gff.filtered.genebank.train --metapars=/mnt/home/azodichr/GitHub/Augustus/config/species/tetraselmis/tetraselmis_metapars.cfg --AUGUSTUS_CONFIG_PATH=/mnt/home/azodichr/GitHub/Augustus/config/


# 09
mkdir tetraselmis_2sp
cp generic/* tetraselmis_2sp/.
cd tetraselmis_2sp/
for i in *; do mv $i $(echo $i | sed 's/generic/tetraselmis_2sp/'); done
sed -i 's/generic_/tetraselmis_2sp_/g' tetraselmis_2sp_parameters.cfg

optimize_augustus.pl --species=tetraselmis_2sp tetra.contigs.10x.pilon_3.fa.all.gff.filtered.genebank.train --metapars=/mnt/home/azodichr/GitHub/Augustus/config/species/tetraselmis_2sp/tetraselmis_2sp_metapars.cfg --AUGUSTUS_CONFIG_PATH=/mnt/home/azodichr/GitHub/Augustus/config/

* Note this optimization step also trains using the best parameters.*

Training results:
08: Transcript level: sensitivity = 0.075, specificity = 0.067. Nucleotide level: sensitivity = 0.86, specificity = 0.68
09: Transcript level: sensitivity = 0.104, specificity = 0.082. Nucleotide level: sensitivity = 0.87, specificity = 0.69


3. Test AUGUSTUS with selected parameters

# 08 
augustus --species=tetraselmis_2sp --AUGUSTUS_CONFIG_PATH=/mnt/home/azodichr/GitHub/Augustus/config/ tetra.contigs.10x.pilon_3.fa.all.gff.filtered.genbank.test > test_results.txt
---------------------------------------------\
                 | sensitivity | specificity |
---------------------------------------------|
nucleotide level |       0.828 |       0.692 |
---------------------------------------------/
----------------------------------------------------------------------------\
transcript | #pred | #anno |   TP |   FP |   FN | sensitivity | specificity |
----------------------------------------------------------------------------|
gene level |   267 |   227 |   23 |  244 |  204 |       0.101 |      0.0861 |
----------------------------------------------------------------------------/


# 09
augustus --species=tetraselmis_2sp --AUGUSTUS_CONFIG_PATH=/mnt/home/azodichr/GitHub/Augustus/config/ tetra.contigs.10x.pilon_3.fa.all.gff.filtered.genebank.test > test_results.txt

---------------------------------------------\
                 | sensitivity | specificity |
---------------------------------------------|
nucleotide level |       0.878 |       0.711 |
---------------------------------------------/
----------------------------------------------------------------------------\
transcript | #pred | #anno |   TP |   FP |   FN | sensitivity | specificity |
----------------------------------------------------------------------------|
gene level |   114 |    98 |   13 |  101 |   85 |       0.133 |       0.114 |
----------------------------------------------------------------------------/


4. Train SNAP
module purge
module load GCC/7.2.0-2.29 MPICH/3.2.1
module load maker/2.31.9

# export 'confident' gene models from MAKER and rename to something meaningful
maker2zff -x 0.25 -l 50 -d ../09_Pilon10x3_prot_3species/tetra.contigs.10x.pilon_3.fa.maker.output/tetra.contigs.10x.pilon_3.fa_master_datastore_index.log

rename 's/genome/tetra10xpilonx3.zff.length50_aed0.25/g' *
for i in *; do mv $i $(echo $i | sed 's/genome/tetra10xpilonx3.zff.length50_aed0.25/'); done


# gather some stats and validate
fathom tetra10xpilonx3.zff.length50_aed0.25.ann tetra10xpilonx3.zff.length50_aed0.25.dna -gene-stats > gene-stats.log
fathom tetra10xpilonx3.zff.length50_aed0.25.ann tetra10xpilonx3.zff.length50_aed0.25.dna -validate > validate.log

# collect the training sequences and annotations, plus 1000 surrounding bp for training
fathom tetra10xpilonx3.zff.length50_aed0.25.ann tetra10xpilonx3.zff.length50_aed0.25.dna -categorize 1000 > categorize.log
fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log

# create the training parameters
mkdir params
cd params
forge ../export.ann ../export.dna > ../forge.log
cd ..

# assembly the HMM
hmm-assembler.pl tetra10xpilonx3.zff.length50_aed0.25 params > tetra10xpilonx3.zff.length50_aed0.25.hmm




5. Re-run MAKER with augustus trained gene models

# Recycle the mapping of empicial evidence we have from the first MAKER round, so we don't have to perform all the BLASTs, etc. again
# transcript alignments
awk '{ if ($2 == "est2genome") print $0 }' ../09_Pilon10x3_prot_3species/tetra.contigs.10x.pilon_3.fa.maker.output/tetra.contigs.10x.pilon_3.fa.all.gff > tetra10xpilonx3.maker.est2genome.gff
# protein alignments
awk '{ if ($2 == "protein2genome") print $0 }' ../09_Pilon10x3_prot_3species/tetra.contigs.10x.pilon_3.fa.maker.output/tetra.contigs.10x.pilon_3.fa.all.gff > tetra10xpilonx3.maker.protein2genome.gff
# repeat alignments
awk '{ if ($2 ~ "repeat") print $0 }' ../09_Pilon10x3_prot_3species/tetra.contigs.10x.pilon_3.fa.maker.output/tetra.contigs.10x.pilon_3.fa.all.gff > tetra10xpilonx3.maker.repeats.gff





###################### MAKER Annotation - ab inito gene model predictions (i.e. Round #2) ###################### 

1. Copy MAKER opts from first round and add AUGUSTUS hmm file and tell it not to predict gene models only from ESTs or Proteins

maker_opts.ctl:
 augustus_species=tetraselmis
 est2genome=0
 protein2genome=0

sbatch run_maker





###### BUSCO Analysis ######

module load python3

!! If running in genome mode (--mode geno) need to set the paths to augustus

export PATH="/mnt/home/azodichr/LocalPrograms/augustus.2.5.5/bin:$PATH"
export PATH="/mnt/home/azodichr/LocalPrograms/augustus.2.5.5/scripts:$PATH"
export AUGUSTUS_CONFIG_PATH="/mnt/home/azodichr/LocalPrograms/augustus.2.5.5/config/"


# Genome 
ln -s /mnt/gs18/scratch/users/azodichr/10xgenom_2232_20180405/01_bam/Tet_10x_pilon_round3/tetra.contigs.10x.pilon_3.fa.fasta 00_Genome/.
export BUSCO_CONFIG_FILE="/mnt/home/azodichr/02_Tetraselmis/07_BUSCO/busco/config/config.ini"
python ../busco/scripts/run_BUSCO.py --i tetra.contigs.10x.pilon_3.fa.fasta --out tetra_pilon10xx3_genome --lineage_path ../busco/lineages/eukaryota_odb9/ --mode geno

# Transcriptomes (testing to see if we're just not annotating them well)
ln -s /mnt/gs18/scratch/users/azodichr/10xgenom_2232_20180405/02_Assembly/trinity_181228.Trinity.fasta 02_transcriptomes/.
python ../busco/scripts/run_BUSCO.py --i trinity_181228.Trinity.fasta --out tetra_trans_TrinGuid_Pilon10xx3 --lineage_path ../busco/lineages/eukaryota_odb9/ --mode tran

# Protein sequences from MAKER annotation
ln -s ~/02_Tetraselmis/06_MAKER/07_trin_10xpilonx3_2/tetra.contigs.10x.pilon_3.fa.maker.output/tetra.contigs.10x.pilon_3.fa.all.maker.proteins.fasta 01_peptides/.
python ../busco/scripts/run_BUSCO.py --i tetra.contigs.10x.pilon_3.fa.all.maker.proteins.fasta --out tetra_prot_TrinGuid_Pilon10xx3 --lineage_path ../busco/lineages/eukaryota_odb9/ --mode prot



# Plot BUSCO summaries:
wd: /mnt/home/azodichr/02_Tetraselmis/07_BUSCO/BUSCO_summaries 
cp ../00_Genome/run_tetra_pilon10x_genome/short_summary_tetra_pilon10x_genome.txt .
cp ../01_peptides/run_tetra_prot_TrinGuid_Pilon10x/short_summary_tetra_prot_TrinGuid_Pilon10x.txt .
cp ../02_transcriptomes/run_tetra_trans_TrinGuid_Pilon10x/short_summary_tetra_trans_TrinGuid_Pilon10x.txt .

python busco/scripts/generate_plot.py --working_directory BUSCO_summaries/




