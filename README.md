#  FaultQuake: An Open-Source Python Tool for Estimating Seismic Activity Rates in Faults


# OVCW workflow

![FaultQuake](https://github.com/GeoSignalAnalysis/FaultQuake/blob/main/drawing_conflation6.png)



# FaultQuake workflow



![FaultQuake](https://github.com/GeoSignalAnalysis/FaultQuake/blob/main/FaultQuake_workflow.png)


FaultQuake is primarily developed and tested on Debian-based Linux OS systems. Therefore, we suggest using FaultQuake in such environments for the best experience. While it's possible to use FaultQuake on Windows and macOS, there may be challenges during compiling and running the workflow due to potential compatibility issues.

We greatly value community contributions and are steadfastly committed to continuously addressing and resolving any bugs that arise in the repository. Should you encounter any issues, please don't hesitate to contact us.

We implement the FaultQuake workflow in six steps, using a FaultQuake conda environment:

## Installation
The installation guides for these environments are provided below:

# FaultQuake environment:
Create and activate a conda environment, FaultQuake for detecting the primary events:


```bash
conda create -n faultquake python=3.10
conda activate faultquake
pip install statsmodels
pip install numpy scipy matplotlib PyQt5

```


## How to run FaultQuake 
```bash
conda activate faultquake
python ./FaultQuake.py

```


## Usage 


 

## To cite: 




BibTex:
```
@inproceedings{ntavakolizadeh2023,
  title={FaultQuake: An Open-Source Python Tool for Estimating Seismic Activity Rates in Faults},
  author={Tavakolizadeh, Nasrin and Mohammadigheymasi, Hamzeh and Matias, Lu{\'\i}s and Silveira, Gra{\c{c}}a and Fernandes, Rui and Dolatabadi, Nima},
  booktitle={Computers and Geosciences},
  pages={xxx},
  year={xxx}
}


```

## License 
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. For more details, see in the license file.

## Contributing
If you would like to contribute to the project or have any suggestions about the code, please feel free to create Pull Requests, raise issues and contact me. 
If you have any questions about the usage of this package or find bugs in the code, please also feel free to contact me.

## Contact information 
Copyright(C) 2023 Nasrin Tavakolizadeh 
Author: Nasrin Tavakolizadeh (n.tavakolizadeh@ubi.pt), Hamzeh Mohammadigheymasi (hamzeh@ubi.pt), Francesco Visini (francesco.visini@ingv.it)



