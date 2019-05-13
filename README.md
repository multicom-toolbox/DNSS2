# DNSS2
Deep learning architectures for protein secondary structure prediction (version 2)

Test Environment
--------------------------------------------------------------------------------------
Red Hat Enterprise Linux Server release 6.4 (Santiago)

Installation Steps
--------------------------------------------------------------------------------------

**(A) Download and Unzip DNSS2 source package**  

Create a working directory called 'DNSS2' where all scripts, programs and databases will reside:
```
Download the DeepSF code:
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

**(E) Configuration**

```
perl configure.pl
```