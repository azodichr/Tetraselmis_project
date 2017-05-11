# Celera Resources
Help thread: https://www.biostars.org/p/115081/

Simple Control File: http://www.cbcb.umd.edu/software/PBcR/data/selfSampleData/pacbio.spec

RunCA Manual (PBcR uses the same options): http://wgs-assembler.sourceforge.net/wiki/index.php/RunCA


# Prep Celera on Calculon2
### Download and compile Celera into Calculon2
Log onto Calculon2 root: (root@calculon2.plantbiology.msu.edu).
Note: On Calculon2 you don’t have to change where you are installing like we did in HPC.

#### Canu
<pre><code>$ cd /share/apps/_tarArchives   #save program zip file into archive
# In a local terminal window:
$ sftp root@calculon2.plantbiology.msu.edu 
$ put Downloads/canu-1.4.Linux-amd64.tar.xz

# Back in Calculon2
$ mv root/canu-1.4.Linux-amd64.tar.xz share/apps/_tarArchives/
$ mkdir ../Canu
$ cp canu-1.4.Linux-amd64.tar.xz ../Canu/
$ tar xvfJ canu-1.4.Linux-amd64.tar.xz          # Decompress
$ mkdir 1.4 
$ mv canu-1.4 1.4/                     # Nick stores everything by name/version.
$ cd 1.4/canu-1.4/
    
# Test out canu to see if it's working
$ cd Linux-amd64/bin/
$ ./canu</code></pre>

#### Celera
<pre><code>$ cd apps/_tarArchives/               #save program zip file into archive
$ wget [url to program]
$ mkdir ../Celera
$ cp .bx2 ../Celera
$ mv ../Celera
$ tar xjf wgs-8.3rc2.tar.bz2
$ rm wgs-8.3rc2.tar.bz2
$ mkdir 8.3 
$ mv wgs-8.3 8.3/                     # Nick stores everything by name/version.
$ cd 8.3/kmer/
$ make install
$ cd ../src/
$ make                                # This creates the Linux-amd64 folder in the 8.3 dir, which contains the bin 

# Test out PBcR to see if it's working
$ cd ../Linux-amd64/bin/
$ ./PBcR</code></pre>


### Make Canu into a Calculon2 module
#### Prep Canu directory in Modules
<pre><code>$ cd /usr/share/Modules/modulefiles/bio
$ mkdir Canu
$ cd Canu
$ vim 1.4                           #Copy from another module file. Change the appfile, app, ver, and the path to the bin (the first prepend-path line - note any lines after that can be removed unless they are needed by the new module).
$ cp 1.4 canu
# Go back to modulefiles dir (cd ../../)
$ vim UpdateBioModules.sh           # Add "cp bio/Canu/canu ."
$ sh UpdateBioModules.sh            # ls and you should see canu in the list.</code></pre>

Pull Modules Here to Prep for Pushing them in the Comp Nodes
<pre><code>$ cd /share/modules
$ sh 1_PullCustomModules.sh</code></pre>

Load to the computing modules
<pre><code>$ ssh compute-0-1
$ cd /share/modules/
$ sh 2_PushCustomModules.sh
$ cd /usr/share/Modules/modulefiles
$ sh UpdateBioModules.sh               # This should make the canu file show up!

#Test
$ module load canu
$ canu                                 #or try $ ./canu
$ logout

# Repeate in compute-0-2:
$ ssh compute-0-2
$ ...
</code></pre>

### Make Celera into a Calculon2 module
#### Prep Celera directory in Modules
<pre><code>$ cd /usr/share/Modules/modulefiles/bio
$ mkdir Celera
$ cd Celera
$ vim 8.3                           #Copy format from another module file. Just be sure to change the appfile, app, ver!
$ cp 8.3 celera
# Go back to modulefiles dir.
$ vim UpdateBioModules.sh           # Add the new celera to the list
$ sh UpdateBioModules.sh            # ls and you should see celera in the list.</code></pre>

Pull Modules Here to Prep for Pushing them in the Comp Nodes
<pre><code>$ cd /share/modules
$ sh 1_PullCustomModules.sh</code></pre>

Load to the computing modules
<pre><code>$ ssh compute0-0
$ cd /share/modules/
$ sh 2_PushCustomModules.sh
$ cd ../../usr/share/Modules/modulefiles
$ sh UpdateBioModules.sh               # This should make the celera file show up!
#Test
$ module load celera
$ PBcR
$ logout
$ ssh compute0-1
#Repeat above steps in this compute node.</code></pre>


# Test run on E. coli dataset using Canu
Instructions from Canu: http://canu.readthedocs.io/en/latest/quick-start.html

<pre><code>$ curl -L -o p6.25x.fastq http://gembox.cbcb.umd.edu/mhap/raw/ecoli_p6_25x.filtered.fastq
$ canu  -p ecoli -d ecoli-auto  genomeSize=4.8m  -pacbio-raw p6.25x.fastq gnuplotTested=true useGrid=false</code></pre>

Note on useGrid: Canu wants to use the grid-engine (Sungrid something...) but ours isn't working so turn grid off (ok since we only have 2 compute nodes, Nick can just set it so that it uses all the resources on one node.)


# Run PBcR on Example E. coli dataset using Celera
The PacBio corrected Reads (PBcR) algorithm from Celera is designed for high-noise single molecule sequencing. It uses MinHash alignment process (MHAP) to increase speed and decrease memory usage.

ssh azodi@calculon2.plantbiology.msu.edu

### Get data onto Calculon2
In a local terminal window:
<pre><code>$ sftp azodi@calculon2.plantbiology.msu.edu 
$ sftp> put Downloads/ecoli_filtered.fastq.gz</code></pre>

Back to Calculon2
<pre><code>$ gunzip ecoli_filtered.fastq.gz
$ ssh compute-0-0
$ module load celera
$ cd 00_PacBio/00_Ecoli/</code></pre>

### Make celera.spec file. See examples
$ nohup PBcR -l ecoli -s celera.spec -pbCNS -fastq ecoli_filtered.fastq genomeSize=4650000 sgeName=ecoli "sge=-A ecoli” &

*Error that it is trying to submit a job!



# Assemble Tetraselmis using Canu 
Nick uploaded the data in the commons folder (/export/home/common/TetraPacBio/). Did a preliminary read quality check 81% of the reads had a quality score (PHRED) of 7 or better and 64% had a score of 10 or better (10% error rate). Total reads = 17741. 

Test by assembling just the first subreads fastq file:
<pre><code> $ cp /export/home/common/TetraPacBio/AllFASTQ/*.subreads.fastq .
$ nohup canu -p Tetra664 -d Tetra664_170112 genomeSize=680m errorRate=0.013 -pacbio-raw *subreads.fastq gnuplotTested=true useGrid=false & </code></pre>

The PLB planned generator test (1/16/17 at 9:00 am) shut down calculon2 and stopped the code. Restart on 1/19/17 at 1:30 pm.
<pre><code>$ nohup canu -p Tetra664 -d Tetra664_170112 genomeSize=680m errorRate=0.013 -pacbio-raw *subreads.fastq gnuplotTested=true useGrid=false & </code></pre>
Job ID: 11266

# Run Canu locally installed on Nicks's HPC
wd: /mnt/scratch/azodichr/Tet_assemb

ssh dev-intel16

module load Java/1.8.0_31

<pre><code> 
$ nohup /mnt/home/panchyni/bin/canu -correct -p Tetra664 -d Tetra664_170214 -pacbio-raw  /mnt/scratch/azodichr/Tet_assemb/00_Corrected/*subreads.fastq genomeSize=680m errorRate=0.013 gridOptionsJobName=TetAssembly_1 gridOptions="-l walltime=167:55:00,ncpus=28,mem=32GB" maxMemory=256 maxThreads=28 & 

$ nohup /mnt/home/panchyni/bin/canu -correct -p Tetra664 -d Tetra664_170214 -pacbio-raw  /mnt/scratch/azodichr/Tet_assemb/00_Corrected/*subreads.fastq genomeSize=680m errorRate=0.013 gridOptionsJobName=TetAssembly_1 gridOptions="-l walltime=3:55:00,mem=8GB -V" maxMemory=256 maxThreads=28 &

$ nohup /mnt/home/panchyni/bin/canu -correct -p Tetra664 -d Tetra664_170214 -pacbio-raw  /mnt/scratch/azodichr/Tet_assemb/00_Corrected/*subreads.fastq genomeSize=680m errorRate=0.013 gridOptionsJobName=TetAssembly_1 gridOptions="-V" maxMemory=256 maxThreads=28 &

$ /mnt/home/panchyni/bin/canu -correct -p Tetra664 -d Tetra664_170215 -pacbio-raw  /mnt/scratch/azodichr/Tet_assemb/00_Corrected/*subreads.fastq genomeSize=680m errorRate=0.013 gridOptionsJobName=TetAssembly_1 maxMemory=256 maxThreads=28 gridOptions="-l nodes=1:ppn=8,walltime=3:55:00,mem=4GB -V" gridOptionsCNS="-l nodes=1:ppn=8,walltime=48:00:00,mem=32GB -V" gridOptionsCOR="-l nodes=6:ppn=28,walltime=167:55:00,mem=32GB -V" gridOptionsOEA="-l nodes=6:ppn=28,walltime=167:55:00,mem=32GB -V" gridOptionsCORMHAP="-l nodes=6:ppn=28,walltime=167:55:00,mem=32GB -V" gridOptionsOBTMHAP="-l nodes=6:ppn=28,walltime=167:55:00,mem=32GB -V" gridOptionsOBTOVL="-l nodes=6:ppn=28,walltime=167:55:00,mem=32GB -V" gridOptionsUTGOVL="-l nodes=6:ppn=28,walltime=167:55:00,mem=32GB -V"

$ /mnt/home/panchyni/bin/canu -correct -p Tetra664 -d Tetra664_170221 -pacbio-raw  /mnt/scratch/azodichr/Tet_assemb/00_Corrected/*subreads.fastq genomeSize=680m errorRate=0.013 gridOptionsJobName=TetAssembly_2 maxMemory=256 maxThreads=28 gridOptions="-l nodes=1:ppn=8,walltime=3:59:00,mem=4GB -V"

$ /mnt/home/panchyni/bin/canu -correct -p Tetra664 -d Tetra664_170221_B -pacbio-raw  /mnt/scratch/azodichr/Tet_assemb/00_Corrected/*subreads.fastq genomeSize=680m errorRate=0.013 gridOptionsJobName=TetAssembly_2 maxMemory=256 maxThreads=28 gridOptions="-l nodes=1:ppn=8,walltime=3:59:00,mem=4GB -V" gridOptionsCORMHAP="-l mem=32GB"
</code></pre>

### Try again with adding more wall time and memory to the corrections mhap step. 
With the last runs there were a lot of walltime/memory errors. 
<pre><code>
$ /mnt/home/panchyni/bin/canu -correct -p Tetra664 -d Tetra664_170222 -pacbio-raw  /mnt/scratch/azodichr/Tet_assemb/00_Corrected/*subreads.fastq genomeSize=680m errorRate=0.013 gridOptionsJobName=170222 maxMemory=256 maxThreads=28 gridOptions="-l nodes=1:ppn=8,walltime=3:59:00,mem=4GB -V" gridOptionsCORMHAP="-l walltime=8:00:00,mem=32GB"
</code></pre>
1. Precompute step timed out on 1 job (39688266[170]). Reran it using the following code from scratch/azodichr/Tet_assemb/Tetra664_170222/canu-scripts/canu.01.out:
<pre><code>
      qsub \
        -l nodes=1:ppn=8,walltime=3:59:00,mem=4GB -V -l walltime=8:00:00,mem=32GB -l mem=18g -l nodes=1:ppn=9 \
        -d `pwd` -N "cormhap_Tetra664_170222" \
        -t 170 \
        -j oe -o /mnt/scratch/azodichr/Tet_assemb/Tetra664_170222/correction/1-overlapper/precompute.\$PBS_ARRAYID.out \
        /mnt/scratch/azodichr/Tet_assemb/Tetra664_170222/correction/1-overlapper/precompute.sh
</code></pre>

2. Correction MHAP had 11 job cancel (39713843[X] - X = 12, 23, 71, 194, 206, 207, 209, 212, 215, 218, 221. Reran it using the following code from scratch/azodichr/Tet_assemb/Tetra664_170222/canu-scripts/canu.02.out:
<pre><code>
qsub -l nodes=1:ppn=8,walltime=3:59:00,mem=4GB -V -l walltime=8:00:00,mem=32GB -l mem=18g -l nodes=1:ppn=9 -d `pwd` -N "cormhap_Tetra664_170222" -t X -j oe -o /mnt/scratch/azodichr/Tet_assemb/Tetra664_170222/correction/1-overlapper/mhap.\$PBS_ARRAYID.out /mnt/scratch/azodichr/Tet_assemb/Tetra664_170222/correction/1-overlapper/mhap.sh
</code></pre>

### Try again with adding even more wall time to the corrections mhap step. 
With the last runs there were a lot of walltime/memory errors. 
<pre><code>
/mnt/home/panchyni/bin/canu -correct -p Tetra664 -d Tetra664_170228 -pacbio-raw  /mnt/scratch/azodichr/Tet_assemb/00_Corrected/*subreads.fastq genomeSize=680m errorRate=0.013 gridOptionsJobName=170228 maxMemory=256 maxThreads=28 gridOptions="-l nodes=1:ppn=8,walltime=3:59:00,mem=4GB -V" gridOptionsCORMHAP="-l walltime=24:00:00,mem=32GB"
</code></pre>


### Nick IDed a subset of fastq files that might be causing the problem (m161028...)
Move "bad" fastq files to /mnt/scratch/azodichr/Tet_assemb/00_Corrected/00_Bad_m161028/
Have Shin-Han run: 
<pre><code>module load Java/1.8.0_31

/mnt/home/panchyni/bin/canu -correct -p Tetra664 -d Tetra664_170321 -pacbio-raw /mnt/scratch/azodichr/Tet_assemb/00_Corrected/*subreads.fastq genomeSize=680m errorRate=0.013 gridOptionsJobName=Tet_0321 maxMemory=256 maxThreads=28 gridOptions="-l nodes=1:ppn=8,walltime=3:59:00,mem=4GB -V" gridOptionsCORMHAP="-l walltime=12:00:00,mem=32GB"</code></pre>

Nick may have found more bad datasets (m161:022, 023, 024, 027, 029). 
Origional data size: 93 G
After removing m161028: 92 G
After removing the rest: 86 G

<pre><code>module load Java/1.8.0_31

/mnt/home/panchyni/bin/canu -correct -p Tetra664 -d Tetra664_170323 -pacbio-raw /mnt/scratch/azodichr/Tet_assemb/00_Corrected/*subreads.fastq genomeSize=680m errorRate=0.013 gridOptionsJobName=Tet_0323 maxMemory=256 maxThreads=28 gridOptions="-l nodes=1:ppn=8,walltime=3:59:00,mem=4GB -V" gridOptionsCORMHAP="-l walltime=18:00:00,mem=32GB"</code></pre>

4/6/17
<pre><code> /mnt/home/panchyni/bin/canu -correct -p Tetra664 -d Tetra664_170321 -pacbio-raw /mnt/scratch/azodichr/Tet_assemb/00_Corrected/*subreads.fastq genomeSize=680m errorRate=0.013 gridOptionsJobName=Tet_0406 maxMemory=256 maxThreads=28 gridOptions="-l nodes=1:ppn=8,walltime=3:59:00,mem=8GB -V" gridOptionsCORMHAP="-l walltime=12:00:00,mem=64GB"</code></pre>

4/7/17: From yesterday's run, it completed all of the mhap runs after 24 hours, but then the next "bucketize" step gave us segmentation faults for every run - try again with asking for the maximum amount of memory for each section.
<pre><code> /mnt/home/panchyni/bin/canu -correct -p Tetra664 -d Tetra664_170407 -pacbio-raw /mnt/scratch/azodichr/Tet_assemb/00_Corrected/*subreads.fastq genomeSize=680m errorRate=0.013 gridOptionsJobName=Tet_0406 maxMemory=256 maxThreads=28 gridOptions="-l nodes=1:ppn=8,walltime=3:59:00,mem=256GB -V" gridOptionsCORMHAP="-l walltime=12:00:00,mem=256GB"</code></pre>

4/7/17: Nick remembered he had that seg fault error before with the test set - it was a problem with canu that they corrected - nick downloaded the updated version, so I should run it with that before trying it with the huge memory requests
<pre><code> /mnt/home/panchyni/bin/canu0317 -correct -p Tetra664 -d Tetra664_170407 -pacbio-raw /mnt/scratch/azodichr/Tet_assemb/00_Corrected/*subreads.fastq genomeSize=680m errorRate=0.013 gridOptionsJobName=Tet_0406 maxMemory=256 maxThreads=28 gridOptions="-l nodes=1:ppn=8,walltime=3:59:00,mem=64GB -V" gridOptionsCORMHAP="-l walltime=12:00:00,mem=64GB"</code></pre>


4/14/17: Changing it to requesting 68 GB in the main gridOptions slowed it down A LOT! change back to 8BG.
<pre><code> /mnt/home/panchyni/bin/canu0317 -correct -p Tetra664 -d Tetra664_170414 -pacbio-raw /mnt/scratch/azodichr/Tet_assemb/00_Corrected/*subreads.fastq genomeSize=680m errorRate=0.013 gridOptionsJobName=Tet_0406 maxMemory=256 maxThreads=28 gridOptions="-l nodes=1:ppn=8,walltime=3:59:00,mem=8GB -V" gridOptionsCORMHAP="-l walltime=12:00:00,mem=64GB"</code></pre>


# Run assembly on the trimmed reads Nick completed on Calculon2
5/10/17
<pre><code>scp azodi@calculon2.plantbiology.msu.edu:/data/scratch/panchy/TetraselmisAssembly/Tetra664_04122017_TrimAssemble/Tetra664.TrimmedReads.fasta
module load Java/1.8.0_31
</code></pre>

First try:  ### Jobs started going over their wall time during the UTGOVL stage (first 70 jobs finished in time, then they started erroring out... set walltime limit to 20 hours.)
<pre><code>/mnt/home/panchyni/bin/canu0317 -assemble -p Tetra664 -d Tetra664_170511 genomeSize=680m correctedErrorRate=0.013 gridOptionsJobName=Tet_0510 -pacbio-corrected Tetra664.TrimmedReads.fasta merylMemory=8 gridOptions="-l nodes=1:ppn=8,walltime=3:59:00,mem=8GB -V" gridOptionsUTGOVL="-l walltime=20:00:00,mem=12GB"
</code></pre>

### Canu process notes:
Correction process (NDS:TSK:Req’dMemory:hrs)
1. Submits canu_X that runs for 10 minutes then “H” (1:8:4:4)
2. Submits meryl_X that runs for under 2 hours then finishes (1:9:18:4)
3. Submits cormhap_X (1:9:18:168)



# Down the road:
LastZ - compare genome across algae (multiZ)
Green Cut - known photosynthetic genes to use for training the maker annotation. Do these ones manually.
