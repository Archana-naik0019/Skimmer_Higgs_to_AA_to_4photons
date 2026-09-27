This skimmer processes input files and retains only the subset of NanoAOD branches (required for the H to AA to 4 photons analysis) specified in [`core/config_Data_updated.py`](https://github.com/Archana-naik0019/Skimmer_Higgs_to_AA_to_4photons/blob/main/core/config_Data_updated.py) for data and in [`core/config/updated.py`](https://github.com/Archana-naik0019/Skimmer_Higgs_to_AA_to_4photons/blob/main/core/config_updated.py) for MC samples. All skimmed output files will contain exclusively these configured branches.
Additionally, a second tier of skimming is performed via [`selection/event.py`](https://github.com/Archana-naik0019/Skimmer_Higgs_to_AA_to_4photons/blob/main/selection/event.py) by applying:
- Lumi-based filtering
- High-Level Trigger (HLT) filtering
- Basic photon selection cuts

#How to run

1. **Interactive Run**
  To run the skimmer interactively on a single ROOT file or a test sample, use the following command:
    ```bash
    python nano_reduce.py --input /path/to/input_file.root --output /path/to/output_file.root --is-data --lumimask <golden.json file name>

2. **Condor Submission**
  To submit jobs on HT Condor, use the script, use the following command:
   ```bash
   python3 condor_submitter.py --fileset <.json file with xrd links to NanoAODs>  --output-dir <outdir name> --jobflavour <>

Reference used to write the skimmer scripts: https://github.com/Bapi2X22/My_Analysis/tree/main/2024_efficiency_study/Backgrounds/Skimmer .
