# PRA2003 — Monitoring Bacterial Movement and Populations
**Leire Ugalde Salvá — i6384470**

## Overview and context 

This project analyzes bacterial tracking data to study the movement and proliferation of different bacterial strains under specific conditions, and how behavior differs between wild-type (WT) strains and their genetic variants (mutants). This matters because mutant strains (antibiotic resistance, capsule deficient, drug resistant...) can behave differently to their WT pair under the same conditions. These mutations can potentially affect bacterial growth or movement. The impact of bacterial mutations under these conditions are not only measured by survival, but also by the possible asymmetry between a WT/mutant and their momentum. For each event, the data is read and used to answer the three questions below.

## Research Questions

1. What are the average counts of each bacterial strain and their statistical uncertainties?
2. Is there any asymmetry between the normal (WT) and mutant strain of each pair?
3. Is there any asymmetry as a function of their momentum?

## Installation and usage
Needs R (tested with R 4.x), no additional packets are required.

To run:
1. Put all the data files (output-Set0.txt, output-Set1.txt ... output-Set10.txt) in the same folder as the scripts.
2. Open that folder in RStudio/VS Code, or make sure your working directory is set to it 
3. Depending on what you want to run:
   - single dataset, week 3 version: `source("PRA2003Week3Corrected.R")` will ask you to type a filename, or just press enter to use output-Set1.txt
   - all 10 datasets combined, week 4: `source("PRA2003Week4.R")` runs automatically, no input needed, just needs PRA2003Week3Corrected.R,  in the same folder
4. Output:
   - week 3 prints the results table straight to the console
   - week 4 prints the combined results and the WT/mutant comparison to the console, and also saves them as CSV files (results-PerDataset.csv, results-Combined.csv, results-PairwiseAsymmetry.csv) in the same folder, so they can be checked without re-running the script

## Files 
| File | Description |
|---|---|
| `PRA2003Week2.R` | Reads one event from the data file (`output-Set0.txt`) and calculates the momentum magnitude of each bacterium. |
| `PRA2003Week3.R` | First attempt at reading the average bacteria count per event and its statistical uncertainty (default `output-Set1.txt`, but allows user input). **Contained errors, not used in later weeks.** |
| `PRA2003Week3Corrected.R` | Corrected version: reads the average bacteria count per event and its statistical uncertainty (default `output-Set1.txt`). This is the version used in all later weeks. |
| `PRA2003Week4.R` | Runs the Week 3 analysis across all 10 datasets (`output-Set1.txt` to `output-Set10.txt`), combines the results into an overall average and uncertainty per strain, and tests whether there is asymmetry between each WT/mutant pair. |

Additional files will be added and documented in upcoming weeks.

### More detail on `PRA2003Week4.R`

- Runs the same strain-counting logic from `PRA2003Week3Corrected.R` on each of the 10 files individually.
- Combines the results by summing the raw counts and event totals across all 10 datasets before dividing, this pools the counting data.
- Computes an overall average count per event and uncertainty for each of the 12 strains.
- Determines whether each WT/mutant pair is symmetric. For each pair (e.g. *E. coli* WT (211) vs. *E. coli* mutant (−211)):
  - Calculates the difference between the two averages.
  - Calculates the combined uncertainty on that difference: `sqrt(uncertaintyA² + uncertaintyB²)`.
  - Calculates a z-score (`difference / combined uncertainty`) to determine how many sigmas apart the values are:
    - |z| ≥ 3 → statistically significant asymmetry
    - 2 ≤ |z| < 3 → weak evidence of asymmetry
    - |z| < 2 → no real asymmetry
- Output
    - `results-PerDataset.csv` per-strain averages for each of 10 datasets individually
    - `results-Combined.csv`pooled average and uncertainty per strain across 10 datasets
    - `results-PairwiseAsymmetry.csv` difference, combined uncertainty and z-score for each WT/mutant pair. Answers whether the pairs are symmetrical 


## Input data Format

- Data file `output-Set0.txt` is structured as:
  - Header line: `event_id  n_particles`
  - One line per particle: `px py pz code`
    - `px`, `py`, `pz` are momentum components.
    - `code` identifies the strain (normal or mutant).
- The remaining files (`output-Set1.txt`, `output-Set2.txt`, …) follow the same format.

## Results
Answering research questions above 

### 1. Average count per event and uncertainty (combined across all 10 datasets)

| Strain ID | Strain | Average count per event | Uncertainty |
|---|---|---|---|
| 211 | E. coli WT | 19.94951 | 2.08e-03 |
| -211 | E. coli mutant | 19.91721 | 2.08e-03 |
| 321 | Bacillus subtilis WT | 2.50915 | 7.37e-04 |
| -321 | Bacillus subtilis mutant | 2.50346 | 7.36e-04 |
| 2212 | Pseudomonas aeruginosa WT | 1.20803 | 5.11e-04 |
| -2212 | Pseudomonas aeruginosa antibiotic-resistant | 1.18416 | 5.06e-04 |
| 3122 | Streptococcus pneumoniae | 0.27660 | 2.45e-04 |
| -3122 | Capsule-deficient S. pneumoniae | 0.27170 | 2.43e-04 |
| 3312 | Mycobacterium tuberculosis | 0.03944 | 9.24e-05 |
| -3312 | Drug-resistant M. tuberculosis | 0.03900 | 9.19e-05 |
| 3334 | Salmonella enterica | 0.00119 | 1.60e-05 |
| -3334 | Salmonella mutant | 0.00115 |  1.58e-05 |

Results from `PRA2003Week4.R`

### 2. Asymmetry between WT and mutant strains

For each pair, the difference between the two averages was compared against the combined uncertainty on that difference, giving a z-score (number of sigmas apart the two values are). This is the basis for the asymmetry conclusion below. 

| Pair | Average A | Average B | Difference | Combined uncertainty | z-score | Conclusion |
|---|---|---|---|---|---|---|
| E. coli WT (211) vs. mutant (−211) | 19.94951 | 19.91721 | 0.03230 | 2.94e-03 | 10.99 | Asymmetry is significant |
| Bacillus subtilis WT (321) vs. mutant (−321) | 2.50915 | 2.50346 | 0.00569 | 1.04e-03 | 5.46 | Asymmetry is significant |
| Pseudomonas aeruginosa WT (2212) vs. resistant (−2212) | 1.20803 | 1.18416 | 0.02387 | 7.20e-04 | 33.17 | Asymmetry is significant |
| S. pneumoniae (3122) vs. capsule-deficient (−3122) | 0.27660 | 0.27170 | 0.00490 | 3.45e-04 | 14.23 | Asymmetry is significant |
| M. tuberculosis (3312) vs. drug-resistant (−3312) | 0.03944 | 0.03900 | 0.00044 | 1.30e-04 | 3.38 | Asymmetry is significant |
| Salmonella enterica (3334) vs. mutant (−3334) | 0.00119 | 0.00115 | 0.00004 | 2.25e-05 | 1.58 | Asymmetry is significant |



### 3. Asymmetry as a function of momentum
To be completed in the following week. 









