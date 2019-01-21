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
    SacrificeEval.py: runs the sacrifice plot maker class on all data in the 
    directories specified inside of SacrificeEval.py.  You will have to open
    SacrificeEval.py in a text editor and point to the directories where your
    freshly made processed/production (data/MC) subtuples are.
  - /plots/: functions that can prepare data output from find_prompt_delayed.py
    for plotting with matplotlib/seaborn.  Plots also perform fits and output
    deviation from flat for use in systematic uncertainty calculation.
  - /uncertainties/: Script for calculating the uncertainty of the neutron
    capture's data cleaning sacrifice.  Mostly deprecated at this point, but
    you may find use for it.  The data cleaning sacrifices and event
    rates in N16, AmBe, and background data must be put into a JSON format and
    fed into CalculateUncertainty.py (see values.json for example).  See Teal's
    DocDB slides on how to find each term, as well as where the equations come
    from.
