# DNSS2
Deep learning architectures for protein secondary structure prediction (version 2)

Test Environment
--------------------------------------------------------------------------------------
Red Hat Enterprise Linux Server release 6.4 (Santiago)

Installation Steps
--------------------------------------------------------------------------------------

**(A) Download and Unzip DNSS2 source package**  

Create a working directory called 'DNSS2' where all scripts, programs and databases will reside:

Download the DNSS2 code:
```
cd ~/
git clone https://github.com/multicom-toolbox/DNSS2.git
cd DNSS2
```

**(B) Install tensorflow, Keras, and h5py and Update keras.json**  

(a) Create python virtual environment (if not installed)
```
virtualenv ~/python_virtualenv_DNSS2
source ~/python_virtualenv_DNSS2/bin/activate
pip install --upgrade pip
```

(b) Install Keras:
```
pip install keras==1.2.2
```

(c) Install theano, numpy, h5py and tensorflow: 
```
pip install numpy==1.12.1
pip install theano==0.9.0
pip install tensorflow==1.5
pip install h5py
```


(f) Add the entry [“image_dim_ordering": "tf”,] to your keras..json file at ~/.keras/keras.json. After the update, your keras.json should look like the one below:  
```
{
    "epsilon": 1e-07,
    "floatx": "float32",
    "image_dim_ordering":"tf",
    "image_data_format": "channels_last",
    "backend": "tensorflow"
}
```
**(C) Download programs**
```
cd programs
wget http://mirrors.vbi.vt.edu/mirrors/ftp.ncbi.nih.gov/blast/executables/blast%2B/2.2.25/ncbi-blast-2.2.25%2B-x64-linux.tar.gz
tar -zxf ncbi-blast-2.2.25%2B-x64-linux.tar.gz


wget https://github.com/soedinglab/hh-suite/releases/download/v3.0-beta.3/hhsuite-3.0-beta.3-Linux.tar.gz
tar -zxf hhsuite-3.0-beta.3-Linux.tar.gz

```

**(D) Download programs**

* DNSS2 requires non-redundent sequence database formated by blast(i.e., uniref90) and hhblits database(i.e., uniclust30_2017_10). If the two databases haven't been downloaded and formated, please try following steps:

```
cd database
mkdir uniref90
cd uniref90
wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
gzip -d uniref90.fasta.gz
../programs/ncbi-blast-2.2.25+/bin/makeblastdb -in  uniref90.fasta


mkdir uniclust30_2017_10
cd uniclust30_2017_10
wget http://wwwuser.gwdg.de/~compbiol/uniclust/2017_10/uniclust30_2017_10_hhsuite.tar.gz

```


**(E) Configuration**

* Update the database paths for the following variables in the script 'configure.pl'
```
$uniref90db = 'your_db_path/uniref90/uniref90.fasta';
$uniclust30db = 'your_db_path/uniclust30_2017_10/uniclust30_2017_10';
```

(b) run 'configure.pl'
```
perl configure.pl
```


**(F) Test**
* There are three ways to indicate the protein to predict:

(1) Predict from protein file:
```
   Usage:
   $ perl runDNSS2.pl -seq <file name>.fasta -file -out <output folder>

   Example:
   a)perl runDNSS2.pl -seq test/2SN3-A.fasta -file -out ./output/2SN3-A
```


(2) Directly give a sequence:
```
   Usage:
   $  perl runDNSS2.pl -seq <AA sequence> -name <Protein Name> -out   -out <output folder>

   Example:
   a). perl runDNSS2.pl -seq KEGYLVKKSDGAKYGXLKLGENEGCDTEDKAKNQGGSYGYXYAFACWDEGLPESTPTYPLPNKSA -name 2SN3-A  -out ./output/2SN3-A
```


(3) Predicting multiple proteins:

```
   Usage:
   $ perl runDNSS2.pl -indir <input directory> -out <output directory>

   Example:
   a) perl runDNSS2.pl -indir ./test/ -out ./output/
```



