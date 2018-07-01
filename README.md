This repository contains scripts used in the analysis of AmBe data in SNO+.

Scripts and directories of interest:

  - find_prompt_delayed.py: File used to create prompt/delayed pairs based
    on nhit/time/fit criteria.  
  - /plots/: functions that can prepare data output from find_prompt_delayed.py
    for plotting with matplotlib/seaborn.  Plots also perform fits and output
    deviation from flat for use in systematic uncertainty calculation.
  - /uncertainties/: Script for calculating the uncertainty of the neutron
    capture's data cleaning sacrifice.  The data cleaning sacrifices and event
    rates in N16, AmBe, and background data must be put into a JSON format and
    fed into CalculateUncertainty.py (see values.json for example).  See Teal's
    DocDB slides on how to find each term, as well as where the equations come
    from.
