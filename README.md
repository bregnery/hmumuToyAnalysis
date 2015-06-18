hmumuToyAnalysis
================

H->mumu Ntuple example analyzer code

Data samples are available on the UF analysis cluster in

    /raid/raid8/jhugon/higgsSamples/stage1/

which physically resides on the melrose node.

Even better, you can access it on the /cms filesystem here:

    /cms/data/store/user/jhugon/hmumu/stage1/

The current analyzer.C includes:
*  Pileup reweighting
*  Higgs signal muon resolution smearing to match the MuSclefit data resolution
*  Muon ID, Iso, and selection cuts mirroring the full analysis

It does not include:  
*  Muon efficiency scale factors
*  HLT trigger matching or checks
*  MuScleFit muon corrections

Instructions
============

Log on to one of the UF Tier 2 servers.  
Optional: Create your own scratch directory:

    cd /scratch/osg/
    mkdir <your_directory>
    cd <your_directory>/

You will need a working copy of ROOT to analyze Ntuples.  The easiest way to get one is to create a CMSSW release, move into it and setup your environment.

    cmsrel CMSSW_7_4_2
    cd CMSSW_7_4_2
    cmsenv

Now you should be able to type the root command and get a working prompt, if so you can continue to setup the anlaysis package.

First, clone this git repository:

    git clone https://github.com/bregnery/hmumuToyAnalysis.git

Then recompile CMSSW:

    scram b
    
move into hmumuToyAnalysis.  analyzer.C is an script for analyzing our Ntuples.  To run it, and load the needed libraries, execute:

    root -b -q -x run.C

It spits out .root files with histograms inside.  If you would like to analyze a different set of input files, then you can change run.C

To make images from the histograms, just run:

    ./plot.py
    
Currently, plot.py stacks each MC histogram, after rescaling appropriately, and displays the data on top.

To display the histograms, type the command:

    display <filename>.png
    
Scale factors, monte carlo pathways, and the maximum number of events can be adjusted in the run file:

    emacs run.C
    
Make sure to run over all events when analyzing monte carlo data
    
The selection criteria can be adjusted in the analysis file:

    emacs analyzer.C
    
