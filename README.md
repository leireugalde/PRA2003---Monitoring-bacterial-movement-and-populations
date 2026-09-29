# PRA2003 — Monitoring Bacterial Movement and Populations
**Leire Ugalde Salvá — i6384470**

## Overview and context 

This project analyzes bacterial tracking data to study the movement and proliferation of different bacterial strains under specific conditions, and how behavior differs between wild-type (WT) strains and their genetic variants (mutants). This matters because mutant strains don't always behave the same as their WT counterpart under the same conditions, resistance or other mutations can come at a cost to growth or movement, or not affect it at all. Looking at survival alone doesn't capture that, so this project also checks for asymmetry between WT/mutant pairs and their momentum.

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
| `PRA2003Week4.R` | Runs the Week 3 analysis across all 10 datasets (`output-Set1.txt` to `output-Set10.txt`) and computes the mean and standard deviation of each strain's count across the 10 sub-samples. |

Additional files will be added and documented in upcoming weeks.

**Method: sub-sampling**

The full sample (5M events) is split into 10 sub-samples of 500K events each (`output-Set1.txt` to `output-Set10.txt`). Each quantity is worked out separately within every subsample. The final result is just the mean of the 10 subsample results, which comes out equal to the result for the full sample. Its statistical uncertainty is taken as the standard deviation across those 10 subsample results.

- Runs the same strain-counting logic from `PRA2003Week3Corrected.R` on each of the 10 files individually.
- Computes, per strain, the mean count across the 10 subsamples and its SD.
- Nothing beyond this is calculated by the script: the WT/mutant asymmetry comparison (see Question 2 below) was worked out manually in Excel, using the mean and SD values the script produces as the starting point.

**Output**
- `results-PerDataset.csv` - per strain averages for each of the 10 datasets individually
- `results-Combined.csv` - mean and standard deviation per strain across the 10 subsamples


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
| 211 | E. coli WT | 19.94951 | 3.27e-02 |
| -211 | E. coli mutant | 19.91721 | 3.19e-02 |
| 321 | Bacillus subtilis WT | 2.50915 | 4.77e-03 |
| -321 | Bacillus subtilis mutant | 2.50346 | 5.50e-03 |
| 2212 | Pseudomonas aeruginosa WT | 1.20803 | 1.90e-03 |
| -2212 | Pseudomonas aeruginosa antibiotic-resistant | 1.18416 | 2.41e-03 |
| 3122 | Streptococcus pneumoniae | 0.27660 | 1.07e-03 |
| -3122 | Capsule-deficient S. pneumoniae | 0.27170 | 9.85e-04 |
| 3312 | Mycobacterium tuberculosis | 0.03944 | 2.84e-04 |
| -3312 | Drug-resistant M. tuberculosis | 0.03900 | 4.02e-04 |
| 3334 | Salmonella enterica | 0.00119 | 4.17e-05 |
| -3334 | Salmonella mutant | 0.00115 |  5.08e-05 |

Results from `PRA2003Week4.R`

### 2. Asymmetry between WT and mutant strains

To check for asymmetry, A had to be worked out separately for each of the 10 sub-samples first. The following formula: 

**A = (N_WT − N_mutant) / (N_WT + N_mutant)**

A's uncertainty is just the standard deviation of its 10 values across the sub-samples. Working out A per subsample instead of on the pooled totals matters because it keeps WT and mutant counts paired to the same events the whole time. They're measured together, not as two separate things. A pair is called asymmetric if |A| is at least 3 times bigger than that uncertainty (the 3σ cutoff), and symmetric otherwise.

All of the values were calculated manually in Excel by applying the formula individually. 

| Pair | Difference (WT − mutant) | Asymmetry A (%) | Significance | Result |
|---|---|---|---|---|
| E. coli | 0.032 ± 0.005 | 0.08 ± 0.01 | 7.2σ | Asymmetric |
| B. subtilis | 0.006 ± 0.003 | 0.11 ± 0.07 | 1.7σ | Symmetric |
| P. aeruginosa | 0.024 ± 0.002 | 1.0 ± 0.1 | 10σ | Asymmetric |
| S. pneumoniae | 0.0049 ± 0.0006 | 0.9 ± 0.1 | 8.5σ | Asymmetric |
| M. tuberculosis | 0.0004 ± 0.0005 | 0.6 ± 0.6 | 0.9σ | Symmetric |
| Salmonella | 0.00004 ± 0.00007 | 1.6 ± 2.9 | 0.5σ | Symmetric |

**Conclusion:** Three pairs come out asymmetric at the 3σ threshold: E. coli, P. aeruginosa and S. pneumoniae, all with WT more abundant than the mutant. The other three: B. subtilis, M. tuberculosis and Salmonella, come out symmetric.


### 3. Asymmetry as a function of momentum
To be completed in the following week. 









