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
<pre><code>
module load Java/1.8.0_31

/mnt/home/panchyni/bin/canu -correct -p Tetra664 -d Tetra664_170321 -pacbio-raw /mnt/scratch/azodichr/Tet_assemb/00_Corrected/*subreads.fastq genomeSize=680m errorRate=0.013 gridOptionsJobName=Tet_0321 maxMemory=256 maxThreads=28 gridOptions="-l nodes=1:ppn=8,walltime=3:59:00,mem=4GB -V" gridOptionsCORMHAP="-l walltime=12:00:00,mem=32GB"
</code></pre>

### Canu process notes:
Correction process (NDS:TSK:Req’dMemory:hrs)
1. Submits canu_X that runs for 10 minutes then “H” (1:8:4:4)
2. Submits meryl_X that runs for under 2 hours then finishes (1:9:18:4)
3. Submits cormhap_X (1:9:18:168)



# Down the road:
LastZ - compare genome across algae (multiZ)
Green Cut - known photosynthetic genes to use for training the maker annotation. Do these ones manually.
