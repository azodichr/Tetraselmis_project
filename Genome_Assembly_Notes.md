# Celera Resources
Help thread: https://www.biostars.org/p/115081/

Simple Control File: http://www.cbcb.umd.edu/software/PBcR/data/selfSampleData/pacbio.spec

RunCA Manual (PBcR uses the same options): http://wgs-assembler.sourceforge.net/wiki/index.php/RunCA


# Prep Celera on Calculon2
### Download and compile Celera into Calculon2
Log onto Calculon2 root: (root@calculon2.plantbiology.msu.edu).
Note: On Calculon2 you don’t have to change where you are installing like we did in HPC.

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


### Make Celera into a Calculon2 module
Prep Celera directory in Modules
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



# Run PBcR on Example E. coli dataset
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
