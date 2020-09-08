# Falcon
Using machine learning to build a fast detector simulator that maps parton jets to reco-jets.

## Requirements
 - [Root](https://root.cern/install/) (tested on version 6.22.2)
 - [FastJet](http://fastjet.fr/) (tested on version 3.3.4)
Using the makefile requires the following environmental variables to be set
(if you ran ```source bin/thisroot.sh``` in your root installation or
installed root in an active conda environment, ROOTSYS has already been set).
```
export ROOTSYS=path/to/your/root/installation 
export FASTJETSYS=path/to/your/fastJet/installation
```
