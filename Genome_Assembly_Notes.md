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

