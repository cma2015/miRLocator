## miRLocator: A python implementation and web server for predicting miRNAs from pre-miRNA sequences

<a href="http://mirlocator.omicstudio.cloud" target="_blank"><img src="https://img.shields.io/badge/Web_server-ready-red.svg"></a><a href="https://hub.docker.com/r/malab/mirlocator" target="_blank"><img src="https://img.shields.io/badge/Docker_image-ready-red.svg" target="_blank"></a><a href="https://hub.docker.com/r/malab/mirlocator" target="_blank"><img src="https://img.shields.io/docker/pulls/malab/mirlocator"></a><a href="https://github.com/cma2015/miRLocator" target="_blank"><img src="https://img.shields.io/badge/Source%20codes-support-blue"></a>

## Introduction

miRNAs are ~21-nucleotide noncoding RNAs, some of which are known to sit at the heart of regulating gene expression in plant growth, development and response to environmental biotic and abiotic stresses. The most commonly used method for discovering miRNAs is based on nest-generation sequencing (NGS) technologies. However, this type of experimental method requires the identification of expressed miRNAs and has a limited ability to detect miRNAs that exhibit low, linkage, stress, developmental and/or cell-specific expression. Therefore, computational tools are urgently required to locate the precise mature miRNAs from pre-miRNA sequences. Here, we present an ML-based system with random forest algorithm named miRLocator for the computational prediction of mature miRNAs within plant pre-miRNAs. We implemented miRLocator into a Docker image (<https://hub.docker.com/r/malab/mirlocator>) and web interface (<http://miRLocator.omicstudio.cloud>) to maximize its practicability. All source codes are availabel at: <https://github.com/cma2015/miRLocator>.

## miRLocator local server construction

#### System Requirement
* Ubuntu (>= 14.04)  

#### Dependencies  
* ViennaRNA V2.0:
  ```bash
  $ sudo apt-add-repository ppa:j-4/vienna-rna
  $ sudo apt-get update
  $ sudo apt-get install vienna-rna
  ```
* Python package: 
  ```bash
  $ pip install scikit-learn scikit-neuralnetwork numpy scipy flask flask-WTF
  ```

#### Building local miRLocator server
```bash
$ git clone https://github.com/cma2015/miRLocator.git
$ cd miRLocator
$ sudo python upload.py
```
Then miRLocator will be available at http://0.0.0.0:8080.

## miRLocator Docker image installation
### Docker installation and start
#### For Windows (Test on Windows 10 Enterprise version):
- Download [Docker](<https://download.docker.com/win/stable/Docker%20for%20Windows%20Installer.exe>) for windows
- Double click the EXE file to open it;
- Follow the wizard instruction and complete installation;
- Search docker, select ___Docker for Windows___ in the search results and clickit.
#### For Mac OS X (Test on macOS Sierra version 10.12.6 and macOS High Sierra version 10.13.3): ####
- Download [Docker](<https://download.docker.com/mac/stable/Docker.dmg>) for Mac os <br>
- Double click the DMG file to open it;
- Drag the docker into Applications and complete installation;
- Start docker from Launchpad by click it.
#### For Ubuntu (Test on Ubuntu 14.04 LTS and Ubuntu 16.04 LTS): ####
- Go to [Docker](<https://download.docker.com/linux/ubuntu/dists/>), choose your Ubuntuversion, browse to ___pool/stable___ and choose ___amd64, armhf, ppc64el or s390x.____ Download the ___DEB___ file for the Docker version you want to install;
- Install Docker, supposing that the DEB file is download into following path: `/home/docker-ce<version-XXX>~ubuntu_amd64.deb`
```bash
$ sudo dpkg -i /home/docker-ce<version-XXX>~ubuntu_amd64.deb      
$ sudo apt-get install -f
```
 ### Verify if Docker is installed correctly ### 
----------------------------------------
   Once Docker installation is completed, we can run `hello-world` image to verify if Docker is installed correctly. Open terminal in Mac OS X and Linux operating system and open CMD for Windows operating system, then type the following command:
```bash
$ docker run hello-world
```
   **<font color =red>Note</font>:** root permission is required for Linux operating system.
   **<font color =red>Note</font>:** considering that differences between different computers may exist, please refer to [official installation manual](https://docs.docker.com/install) if instructions above don’t work.

 ### miRLocator installation from Docker Hub ###
--------------------------------
  For Mac OS X and Linux operating systems, open the terminal, for Windows operating system, open CMD. Typing the following command:
```bash
# Pull miRLocator from Docker Hub
$ docker pull malab/mirlocator
```
### Quckly start Docker miRLocator ###
```bash
# Run miRLocator in Docker as it in local disk and mount a volume
$ docker run -it -v data_dir:/data malab/mirlocator
# enter into the file dictionary of miRLocator
$ cd miRLocator
# list all permissible command-line arguments
$ python miRLocator.py ‐h
```
#### Constructing miRNA prediction model
```bash
# Prediction model construction
$ python miRLocator.py ‐p training ‐i samples/trainingData.txt ‐o train_output ‐m prediction_model ‐k 0
# Prediction model construction with 5‐fold cross validation
$ python miRLocator.py ‐p training ‐i samples/trainingData.txt ‐o train_output ‐m prediction_model ‐k 5
```

#### Mature miRNA prediction
```bash
# Prediction mode
$ python miRLocator.py ‐p prediction ‐i samples/predictionData.txt ‐o predict_output ‐m samples/trained_prediction_model
# Evaluating the performance of the model
$ python miRLocator.py ‐p prediction ‐i samples/predictionData.txt ‐o predict_output ‐m samples/trained_prediction_model ‐a predictionData_Annotated.txt
```

## Citation

Zhang T, Ju L, Zhai J, Song Y, Song J, Ma C. (2019) [**miRLocator: A python implementation and web server for predicting miRNAs from pre-miRNA sequences.**](https://link.springer.com/protocol/10.1007/978-1-4939-9042-9_6) In: de Folter S. (eds) Plant MicroRNAs. Methods in Molecular Biology, vol 1932. Humana Press, New York, NY

Cui H, Zhai J, Ma C. (2015). [**miRLocator: machine learning-based prediction of mature microRNAs within plant pre-miRNA sequences.**](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0142753) PLoS One, 10(11), e0142753.


