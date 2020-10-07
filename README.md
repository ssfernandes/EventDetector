# WINTENDED - Event Detector
This repository contains the source code for the manuscript: "Evolving Social Networks Analysis via Tensor Decompositions: From Global Event Detection 
Towards Local Pattern Discovery and Specification" by Sofia Fernandes, Hadi Fanaee-T and João Gama (2019)

## Requirements
Our code makes use of the [Tensor Toolbox for MATLAB](http://www.sandia.gov/~tgkolda/TensorToolbox/).
       
## Usage
For illustration purposes, we provide a synthetic time-evolving network (datasets/synthi.mat) and a  demo file (demo.m).

In the demo, we apply our method (```eventdetector```) with varying window lengths and number of factors and  then rank the instants abnormality according to the manuscript procedure using ```top_events function```.

For details on the usage of each function (namely, regarding input and output variables), please check the corresponding file.

