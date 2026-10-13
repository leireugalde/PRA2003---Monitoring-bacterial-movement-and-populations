# PRA2003 — Monitoring Bacterial Movement and Populations
**Leire Ugalde Salvá — i6384470**

## Overview and context 

This project analyzes bacterial tracking data to study the movement and proliferation of different bacterial strains under specific conditions, and how behavior differs between wild-type (WT) strains and their genetic variants (mutants). This matters because mutant strains don't always behave the same as their WT counterpart under the same conditions, resistance or other mutations can come at a cost to growth or movement, or not affect it at all. Looking at survival alone doesn't capture that, so this project also checks for asymmetry between WT/mutant pairs and their momentum. This can show whether a mutation can have an impact on the amount of bacteria and their movement (momentum), and whether there is a dependence. 

## Research Questions

1. What are the average counts of each bacterial strain and their statistical uncertainties?
2. Is there any asymmetry between the normal (WT) and mutant strain of each pair?
3. Is there any asymmetry as a function of their momentum?

## Installation and usage
Needs R (tested with R 4.x), no additional packets are required.
The momentum analysis (Week 5/6) also requires Python 3 with numpy, pandas and matplotlib installed. 

To run:
1. Put all the data files (output-Set0.txt, output-Set1.txt ... output-Set10.txt) in the same folder as the scripts.
2. Open that folder in RStudio/VS Code, or make sure your working directory is set to it 
3. Depending on what you want to run:
 - single dataset, week 3 version: `source("PRA2003Week3Corrected.R")` will ask you to type a filename, otherwise, press enter to use output-Set1.txt
   - all 10 datasets combined, week 4: `source("PRA2003Week4.R")` runs automatically, no input needed, just needs PRA2003Week3Corrected.R,  in the same folder
   - momentum analysis (week 5/6) following two steps: 
     1. `source("deliverable6F_momentum.R")` reads all 10 datasets and makes `momenta-Set1.csv` ... `momenta-Set10.csv`.
        - This takes a few minutes as the files are heavy, but it only has to be run once
     2. `python deliverable6F_momentum.py` asks for the WT ID of a pair: 211, 321, 2212, 3122, 3312 or 3334.
        - Run it once per pair, it only takes a few seconds
4. Output:
   - week 3 prints the results table straight to the console
   - week 4 prints the combined results and the WT/mutant comparison to the console, and then saves them as CSV files (`results-PerDataset.csv´, `results-Combined.csv´, `results-PairwiseAsymmetry.csv´) in the same folder, so they can be checked without re-running the script
   - the momentum R script prints how many bacteria of each ID it found per dataset (as a check, these should match the week 4 counts) and saves the `momenta-SetN.csv` files
   - the momentum Python script saves, for the chosen pair: `asymmetry-<ID>.png` (asymmetry vs momentum), `spectrum-<ID>.png` (dN/dP of WT and mutant vs momentum) and `asymmetry-<ID>.csv` (all numbers per momentum bin). It also prints the momentum percentiles, the counts per sub-sample, the table per bin, the overall asymmetry and the line fit to the terminal.

## Files 
| File | Description |
|---|---|
| `PRA2003Week2.R` | Reads one event from the data file (`output-Set0.txt`) and calculates the momentum magnitude of each bacterium. |
| `PRA2003Week3.R` | First attempt at reading the average bacteria count per event and its statistical uncertainty (default `output-Set1.txt`, but allows user input). **Contained errors, not used in later weeks.** |
| `PRA2003Week3Corrected.R` | Corrected version: reads the average bacteria count per event and its statistical deviation as uncertainty (default `output-Set1.txt`). This is the version used in all later weeks. |
| `PRA2003Week4.R` | Runs the Week 3 analysis across all 10 datasets (`output-Set1.txt` to `output-Set10.txt`) and computes the mean and standard deviation of each strain's count across the 10 sub-samples. |
| `deliverable6F_momentum` | Goes through 10 datasets, keeps only the 12 strain IDs (6 WT/mutant pairs), calculates the momentum of every single bacteria and saves ID + momentum in `momenta-Set1.csv` ... `momenta-Set10.csv`. |
| `deliverable6_momentum.py` | Reads the momenta of one WT/mutant pair, puts them in momentum bins, calculates the asymmetry and dN/dP per bin for every sub-sample, fits a straight line to A vs momentum and makes the plots. |
| `asymmetry-<ID>.png` / `.csv`, `spectrum-<ID>.png` | Output of `deliverable6_momentum.py` for each of the 6 pairs. |


**Method: sub-sampling**

The full sample (5M events) is split into 10 sub-samples of around 500K events each (`output-Set1.txt` to `output-Set10.txt`). Each quantity is worked out separately within every subsample. The final result is just the mean of the 10 subsample results, which comes out equal to the result for the full sample. Its statistical uncertainty is taken as the standard deviation across those 10 subsample results.

- Runs the same strain-counting logic from `PRA2003Week3Corrected.R` on each of the 10 files individually.
- Computes, per strain, the mean count across the 10 subsamples and its SD.
- Beyond this the asymmetry for bacterial count was worked out manually in Excel, using the values produced by the code as the starting point. 

**Method: asymmetry (relative formula)**
To compare a WT strain with its mutant, the relative asymmetry is used:

**A = (N_WT − N_mutant) / (N_WT + N_mutant)**

The relative formula is used instead of the plain difference N_WT − N_mutant because:
- The pairs have very different numbers of bacteria (ex. E. coli vs Salmonella), which can impact the results by rather showing which pair is more common, rather than the actual asymmetry, therefore skewing the results.
- The relative formula is used therefore used to normalize the variance caused by individual size, strength, or total scale
- A has no unit and always lies between −1 and +1: A = 0 means equal numbers, A > 0 means more WT, A < 0 means more mutant.

*A is always worked out separately for each of the 10 sub-samples first and then averaged, so that WT and mutant counts stay paired to the same events.*

**Output**
- `results-PerDataset.csv` - per strain averages for each of the 10 datasets individually
- `results-Combined.csv` - mean and standard deviation per strain across the 10 subsamples

**Method: asymmetry as a function of momentum (week 5/6)**
 
**Momentum** 
For every bacterium the total momentum is the length of its momentum vector, described as p = √(px² + py² + pz²), in units of 10⁻²⁰ kg·m/s.This focuses on how much momentum the bacterium has. This is done in R (`deliverable6_momentum.R`) because reading the raw files requires longer. The files for each pair can be done much faster on Python, as these use significantly smaller CSV files.

**Bins**
The momenta of a pair is organized into different momentum bins. This creates a histogram per sub-sample. Each pair had the same bin number and widths, as the asymmetry (A) compares these. Different pairs will have different bins, as they have different amounts of data. The bin edges, as in what determines the intervals,  was chosen based on the percentiles of each pair, which is visible in teh script. Most bacteria have lower momnetums, so there are more and more narrow bins in those regions. As the momentum increases the bins get wider. This allows each bin to have enough bacteria to have a useful A value. The first bin edge is 0, and everything in the last edge goes into the tailbin. This is so no bacteria gets left out, and is observed as the grey shaded area. 

*This shows the bins that were used for each pair* 

| Pair | Bin edges (10⁻²⁰ kg·m/s) | Number of bins (incl. overflow) |
|---|---|---|
| E. coli (211) | 0, 1, 1.5, 2, 2.5, 3, 4, 6, 8, 12, 16, 24, 40 | 13 |
| B. subtilis (321) | 0, 1, 2, 3, 5, 8, 16 | 7 |
| P. aeruginosa (2212) | 0, 1, 1.5, 2, 2.5, 3, 4, 6, 8, 12, 16 | 11 |
| S. pneumoniae (3122) | 0, 1, 1.5, 2, 2.5, 3, 4, 6, 8, 12, 16 | 11 |
| M. tuberculosis (3312) | 0, 1, 1.5, 2, 2.5, 3, 4, 6, 8, 12, 16 | 11 |
| Salmonella (3334) | 0, 1.5, 3, 6 | 4 |
 
E. coli has two extra bins as it has a lot more data. Salmonella on the other hand has significantly less, so it is only given 4 bins. 
 
**Asymmetry per bin.** 
Same idea as Question 2 (check in answering of Question 2 why the relative formula is used). (see there for why the relative formula is used, mainly: the number of bacteria changes a lot with momentum, and A is not affected by that), but now per momentum bin: A = (N_WT − N_mutant) / (N_WT + N_mutant) is calculated in every bin for every one of the 10 sub-samples, and then averaged over the sub-samples. The average is weighted by the number of bacteria the sub-sample has in that bin, because a sub-sample with only a few bacteria in a bin gives a very noisy A and shouldn't count as much as one with thousands. The error bars are the (weighted) standard deviation of A across the 10 sub-samples, the same type of uncertainty as in Questions 1 and 2.

**Asymmetry per bin** 
The asymmetry A (see *Method: asymmetry*) is calculated in every momentum bin, separately for each of the 10 subsamples, and then averaged over the subsamples. This gives one A value per bin, which is what is plotted against momentum. The average is then weighted by how many bacteria each subsample has in that specific bin. If a subsample only has a few bacteria in a bin, it will give a noisy A, so it should count less than a subsample that has thousands. The same relative formula is used as in answering Question 2. 

**dN/dP** 
The script also makes `spectrum-<ID>.png`, which shows how the WT and mutant are spread over the momenta: dN/dP = (number in the bin) / (bin width) / 500,000 events. Dividing by the width is very necessary to correct as the bins are not equally wide, otherwise wide bins would look big just because they are wide. A doesn't need this because the width is in both the top and bottom of the fraction and cancels. The overflow bin has no upper edge (no width), so it is not in this plot.
 
**Answering whether asymmetry is momentum dependent?**
To test this a straight line A = a·P + b is fitted through the points, excluding the tailbin. If the slope a is consistent with zero, A doesn't depend on momentum; if it is more than 3σ away from zero, it does. For only the line fit  the uncertainty of each point is the standard error of the mean, SD/√10, instead of the SD. The fit asks whether the average A changes between bins, and the uncertainty on an average of 10 sub-samples is SD/√10. The error bars in the plots are still the SD, to show how much the sub-samples spread. The straight line is just a simple test for "does A change with momentum", it is not a model of the exact shape.

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

The asymmetry A was calculated for each pair using the relative formula explained in *Method: asymmetry*, separately for each of the 10 sub-samples.A's uncertainty is just the standard deviation of its 10 values across the sub-samples. Working out A per subsample instead of on the pooled totals matters because it keeps WT and mutant counts paired to the same events the whole time. They're measured together, not as two separate things. A pair is called asymmetric if |A| is at least 3 times bigger than that uncertainty (the 3σ cutoff), and symmetric otherwise. 

All of the values were calculated manually in Excel by applying the formula individually. 

| Pair | Difference (WT − mutant) | Asymmetry A (%) | Significance | Result |
|---|---|---|---|---|
| E. coli | 0.032 ± 0.005 | 0.08 ± 0.01 | 7.2σ | Asymmetric |
| B. subtilis | 0.006 ± 0.003 | 0.11 ± 0.07 | 1.7σ | Symmetric |
| P. aeruginosa | 0.024 ± 0.002 | 1.0 ± 0.1 | 10σ | Asymmetric |
| S. pneumoniae | 0.0049 ± 0.0006 | 0.9 ± 0.1 | 8.5σ | Asymmetric |
| M. tuberculosis | 0.0004 ± 0.0005 | 0.6 ± 0.6 | 0.9σ | Symmetric |
| Salmonella | 0.00004 ± 0.00007 | 1.6 ± 2.9 | 0.5σ | Symmetric |

**Conclusion:** 
Three pairs come out asymmetric at the 3σ threshold: E. coli, P. aeruginosa and S. pneumoniae, all with WT more abundant than the mutant. The other three: B. subtilis, M. tuberculosis and Salmonella, come out symmetric.


### 3. Asymmetry as a function of momentum
Results from `deliverable6_momentum.py` (method explained above).

**How to read the graphs**
- green points: average A of the 10 sub-samples in that momentum bin.
   - Vertical bar = SD across the sub-samples
   - Horizontal bar = the width of the bin
- pink band: overall A of the pair (all momenta together) ± SD (same value as in Question 2)
- red dashed line: A = 0, equal numbers of WT and mutant --> not asymmetric
- blue line: straight-line fit, the legend gives the slope and how many σ it is away from zero
- grey area: the tailbin bin (all bacteria above the last bin edge), not used in the fit

**Summary**
| Pair | Overall A (± SD) | Slope a (per 10⁻²⁰ kg·m/s) | Slope significance | Does A depend on momentum? |
|---|---|---|---|---|
| E. coli | 0.00081 ± 0.00011 | (9.2 ± 1.4) × 10⁻⁵ | 6.4σ | Yes, A increases with P |
| B. subtilis | 0.00114 ± 0.00066 | (3.0 ± 0.6) × 10⁻⁴ | 5.2σ | Yes, A increases with P |
| P. aeruginosa | 0.00998 ± 0.00100 | (3.7 ± 0.8) × 10⁻⁴ | 4.5σ | Yes, A increases with P |
| S. pneumoniae | 0.00894 ± 0.00106 | (5.1 ± 1.6) × 10⁻⁴ | 3.2σ | Borderline, just above 3σ |
| M. tuberculosis | 0.00562 ± 0.00623 | (5.0 ± 4.4) × 10⁻⁴ | 1.1σ | No |
| Salmonella | 0.015 ± 0.029 | (−1.3 ± 0.8) × 10⁻² | 1.8σ | Not enough data |

#### E. coli (211 / −211)
PICTURE 

All 13 points are above zero, so WT is more common than the mutant at every momentum, never the other way around. The asymmetry is small (below 1% everywhere). In the lowest bins (P < 3) A is only about 0.0003–0.0008 and the error bars of single bins reach zero, but from P ≈ 3 upward most points are clearly above zero, and the highest bins have the largest values (about 0.0025 between 16 and 40, and about 0.008 above 40). So the asymmetry found in Question 2 (7.2σ) is present across the whole momentum range, and it is largest for the fastest bacteria. The overall A (0.00081) is close to the low-momentum values because most E. coli have a low momentum, so those bins dominate the total.




