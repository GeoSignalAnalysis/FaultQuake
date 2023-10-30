# FaultQuake

![FaultQuake](https://github.com/GeoSignalAnalysis/FaultQuake/blob/main/workflow_faultQuake15.png)




#  FaultQuake: An Open-Source Python Tool for Estimating Seismic Activity Rates in Faults
IPIML is primarily developed and tested on Debian-based Linux OS systems. Therefore, we suggest using IPIML on such environments for the best experience. While it's possible to use IPIML on Windows and macOS, there may be challenges during compiling and running the workflow due to potential compatibility issues.

We greatly value community contributions and are steadfastly committed to continuously addressing and resolving any bugs that arise in the repository. Should you encounter any issues, please don't hesitate to contact us.

We implement the IPIML workflow in six steps, using a FaultQuake conda environment:

## Installation
The installation guides for these environments are provided below:

# FaultQuake environment:
Create and activate a conda environment, IPIML for detecting the primary events:
If you want to process with CPU:
```bash
conda config --add channels conda-forge
conda config --set channel_priority strict
conda create -n IPIML python=3.9 tensorflow==2.11.1 keras==2.11.0 h5py obspy spyder pygmt matplotlib pyyaml pandas tqdm pyproj jupyter notebook basemap six numpy protobuf
conda activate IPIML
pip install keras-rectified-adam seisbench
```


## Usage 

Optimized Deep Learning (DL)-based workflows can improve the efficiency and accuracy of earthquake detection and location processes. IPIML is a six-step automated event detection, phase association, and earthquake location workflow, which integrates the state-of-the-art Pair-Input DL model and waveform Migration Location methods (IPIML). 

 

## To cite: 
Please cite the following paper in your documents if you use MALMI in your work. 

In case you utilize IPIML for processing your data, it would be appreciated if you cite the following paper(s):


BibTex:
```




```

## License 
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. For more details, see in the license file.

## Contributing
If you would like to contribute to the project or have any suggestions about the code, please feel free to create Pull Requests, raise issues and contact me. 
If you have any questions about the usage of this package or find bugs in the code, please also feel free to contact me.

## Contact information 
Copyright(C) 2023 Hamzeh Mohammadigheymasi 
Author: Hamzeh Mohammadigheymasi (hamzeh@ubi.pt), Peidong Shi (peidong.shi@sed.ethz.ch), and Xiao Zhuowei  (xiaozhuowei@mail.iggcas.ac.cn)



