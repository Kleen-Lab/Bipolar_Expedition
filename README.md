# Public code repository for the Journal of Engineering manuscript "Distance of bipolar re-referencing imparts nonlinear frequency-specific influences on intracranial recording signal measurements"

================================

## Package setup. 

### 1. Required software 
MATLAB (tested with Version 2023-2025).    
Chronux toolbox.   
R (tested with version 4.4).   
Python (version 3). 

### 2. Other dependencies 
None

### 3. Installation
To install the package and have a copy of the code to edit locally, navigate to where you would like to store the package code in your terminal. Clone the package.

```https://github.com/Kleen-Lab/bipolar_rereferencing.git```

You will also need the Chronux package. [Chronux](https://chronux.org/). We have packaged a version of Chronux under helpers 

Note that you will have to configure paths in the code to wherever your data is stored. The code is setup such that you have put the data in a specified structure in a ```/data``` folder inside of this code directory, and have a ```/results``` folder inside of this code directory as well. 

### 4. Download data and setup

Download data from : [10.17605/OSF.IO/CN8W4](https://osf.io/cn8w4/)

point your ```setup_bipolar.m``` to a ```data/``` folder. 

point your ```setup_bipolar.m``` to a ```results/``` folder.

Place the ```baseline-high-density-data``` data in a subfolder under the data directory. 

Place the ```imaging``` data in a subfolder under the data directory. 

## Files of interest 

Run the ```master_bipolar_2025.m``` file to generate all analyses and figures. Information to run all analysis and plots is included in the master_script and the functions it calls.

The supplementary figure 1 requires R to recreate. After running the linear analysis, save data into the ```tableoutput.xlsx``` file for plotting via ```supp_figure_1_stats.R```. 

The ```helpers``` folder has the associated helper functions that the scripts call to generate the analyses and figures. 

## License

BSD-3 License
