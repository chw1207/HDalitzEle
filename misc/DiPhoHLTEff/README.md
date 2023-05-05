# Diphoton30_22(18) HLT measurements
This trigger is the main trigger for merged category. To measure the scale factors of it, Please follow the following steps.

1. Produce the Tnp ntuples with the corresponding trigger filters and variables needed for merged electron ID
   - [Modified EgammaAnalysis-TnPTreeProducer](https://github.com/chw1207/EgammaAnalysis-TnPTreeProducer): please follow the instruction of this repository to produce the ntupls

2. Merge the ntuples and add the merged ID score to the ntuples.
   - Use `runMergeTnpNtuples.C` to execute the main script `mergeTnpNtuples.C`
   ```bash 
    root -l -b -q 'runMergeTnpNtuples.C("seed")'   # merge ntuples for seeded leg
    root -l -b -q 'runMergeTnpNtuples.C("unseed")' # merge ntuples for unseeded leg
   ```

3. Use Tag-and-Probe technic to measure the the efficiencies and scale factors
    - [Modified egm_tnp_analysis](https://github.com/chw1207/egm_tnp_analysis): please follow the instruction of this repository to perform the measurements.