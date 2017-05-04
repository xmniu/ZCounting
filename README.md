# ZCounting
A package designed to do Z counting in DQM offline framework.

## Idea
Keep event info (number of primary vertices, lumi-section, etc) and track and muon collection on the fly(analyze), then fill histograms in the end of each run(endRun).

## Histograms
12 = 3 (3 efficiency categories) * 2 (2 eta region) * 2 (passing and failling) for efficiency calculation
1 for Z raw yields
1 for numbers of primary vertices, to be used in next step for MC simulation
 
