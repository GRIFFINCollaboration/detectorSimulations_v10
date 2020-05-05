
To simulate TRIFIC resolving power for an experiment follow these steps:
1. Set experimental parameters (beam Z/A/Energy and gas pressure) in trific.mac
2. run "./Griffinv10 trific.mac"
3. run "root LatestTrificSort.C"
4. mv TRIFIC.root TRIFIC_ZA.root (where ZA is replaced by the current ion)
5. repeat 1-4 for each ion of interest
6. run "hadd Sum.root TRIFIC_*.root"
Then view resultant spectra in Sum.root

Always worth running a few ions with vistrific.mac to ensure setup is as intended
In Griffinv10 gui can navigate to "Scene Tree" and turn off drawing of fTrificPipePhysical to see internals

If any changes to X/Y grid setup is made (after Oct 2019) then one must update DetectionSystemTrific.cc AND LatestTrificSort.C
