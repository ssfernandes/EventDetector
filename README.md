# EventDetector
Windowed Tensor Decomposition Event Detector for Time-evolving networks

SOURCE 
------
    Evolving Social Networks Analysis via Tensor Decompositions: From Global Event Detection 
    Towards Local Pattern Discovery and Specification
    Sofia Fernandes, Hadi Fanaee-T and João Gama
    International Conference on Discovery Science (2019)

DESCRIPTION
-----------
     MATLAB implementation of our Windowed Tensor Decomposition Event Detector for Time-evolving networks
       
FILES
-----
    For illustration purposes, we provide a sythetic time-evolving network (datasets/synthi.mat) and a 
    demo file (demo.m).
    In the demo, we apply our method (eventdetector) with varying window lengths and number of factors and 
    then rank the instants abnormality according to the manuscript procedure using top_events function.
 
DEPENDECIES
-----------
    Our code makes use of the Tensor Toolbox for MATLAB, which is avialable at
        http://www.sandia.gov/~tgkolda/TensorToolbox/.

CONTACT
-------
    To report any bug, comment or suggestion, contact me at sdsf@inesctec.pt.

