This repository contains scripts used in the analysis of AmBe data in SNO+.

Scripts and directories of interest:
  - /utils/: contains the pre-processors used on ntuples prior to the data
    cleaning sacrifice analyses.  Use in the following order:
    1. Process all data with the subtupler program.  Will add the u.r and posr3
       variables to the newly made ntuple files.  The bash script with the
       subtupler is useful for running the subtupler code on an entire
       directory of ntuples at once.
    2. Process all data with find_prompt_delayed.py.  The script opens root files
       and creates prompt/delayed pairs based on nhit/time/fit criteria given as
       flags at runtime.
    /scripts/: Provides examples of how to use the AmBe data cleaning classes
    in lib.  
      * SacrificeEval.py:runs the sacrifice plot maker class on all data in the 
    directories specified inside of the SacrificeEval.py script.  You will have to open
    SacrificeEval.py in a text editor and point to the directories where your
    freshly made processed/production (data/MC) subtuples are.
      * NeutronSacrifice.py: Runs a script that calculates the 4.4 MeV gamma
    and neutron sacrifices using information from fits to the data's interevent
    time distribution.  
  - /lib/: Contains classes that can open SNO+ ntuples and estimate the data 
    cleaning sacrifice.  Given a config file (see /config/), different preliminary
    cuts can be applied prior to evaluating the data cleaning sacrifice.
  - /config/: Examples of configuration files read in by the classes defined in 
    files in the /lib/ directory.  A text file in this directory gives a rundown
    of how each toggle can be used when evaluating the DC sacrifice on data.
  - /uncertainties/: Old ript for calculating the uncertainty of the neutron
    capture's data cleaning sacrifice.  Deprecated at this point, but
    you may find use for it.  The data cleaning sacrifices and event
    rates in N16, AmBe, and background data must be put into a JSON format and
    fed into CalculateUncertainty.py (see values.json for example).  See Teal's
    DocDB slides on how to find each term, as well as where the equations come
    from.
