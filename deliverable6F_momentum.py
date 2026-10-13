# PRA2003 - asymmetry and dN/dP vs momentum for one WT/mutant pair

# input:  momenta-Set1..10.csv (from R, columns ID and p)
# output: asymmetry-ID.png  -> A against momentum + a straight-line fit
#         spectrum-ID.png   -> dN/dP of WT and mutant against momentum
#         asymmetry-ID.csv  -> the numbers behind both plots
#
# rough idea of what happens:
#   1. read momenta of the WT and mutant of one pair, separately for each of the 10 sets
#   2. put them in momentum bins (histogram per set)
#   3. per bin: calculate A for every set, average 10 sets, error = standard deviation (SD) of the 10 sets
#   4. fit line through A vs P -> is the slope different from 0? (= does A depend on momentum)
#   5. averaging for dN/dP, and make the plots

import numpy as np               
import pandas as pd                
import matplotlib.pyplot as plt    

n_sets = 10                  
events_per_set = 500_000     


#   error bars= SD  = standard deviation of the 10 subsamples -> how much one subsample scatters
use_sem = False
error_div = np.sqrt(n_sets) if use_sem else 1   # what the SD gets divided by
err_name = "SEM" if use_sem else "SD"           # for labels and printouts

#  WT has the positive ID and the mutant the same number but negative
# abs() so it still works if someone types -2212 by accident (protection)
wt = abs(int(input("WT ID (e.g. 2212): ")))
mut = -wt


# bins
# bin edges per pair. 
# WT and mutant have the same bins to compare asymmetry
# different pairs have different bins due to different amount of data
# narrow bins at low p, that's where most bacteria are
# wider bins in tail, less bacterial data 
# everything above the last  goes into one extra overflow bin
# number of bins chosen based on percentiles
bins_by_id = {
    211:  [0, 1, 1.5, 2, 2.5, 3, 4, 6, 8, 12, 16, 24, 40],  # E. coli
    321:  [0, 1, 2, 3, 5, 8, 16],                          # B. subtilis
    2212: [0, 1, 1.5, 2, 2.5, 3, 4, 6, 8, 12, 16],         # P. aeruginosa 
    3122: [0, 1, 1.5, 2, 2.5, 3, 4, 6, 8, 12, 16],         # S. pneumoniae
    3312: [0, 1, 1.5, 2, 2.5, 3, 4, 6, 8, 12, 16],         # M. tuberculosis 
    3334: [0, 1.5, 3, 6],                                  # Salmonella: less bacteria so less bins 
}
bins = np.array(bins_by_id[wt], dtype=float)
print("bins:", bins.tolist(), "+ overflow above", bins[-1])


# reading data
# p_wt[0] = all WT momenta in Set1, p_wt[1] = Set2, etc. (same for the mutant)
p_wt, p_mut = [], []
for i in range(1, n_sets + 1):
    df = pd.read_csv(f"momenta-Set{i}.csv")
    p_wt.append(df.loc[df["ID"] == wt, "p"].to_numpy())    # focus only the rows with  WT ID
    p_mut.append(df.loc[df["ID"] == mut, "p"].to_numpy())  # focus only the rows with the mutant ID

# all momenta of pair in one array
all_p = np.concatenate(p_wt + p_mut)
# percentiles: 90th = the momentum that 90% of the bacteria are below
# decide where bins should be and where tail starts
print("percentiles 1,5,10,25,50,75,90,95,99:", np.percentile(all_p, [1, 5, 10, 25, 50, 75, 90, 95, 99]).round(2))
print("max p:", all_p.max())

# how many WT and mutant per set + the median momentum of  set
# Set 10 has the same number of bacteria but a higher momentum (we found that earlier). 
# makes the spread and error bars a bit bigger
n_wt = np.array([len(x) for x in p_wt])
n_mut = np.array([len(x) for x in p_mut])
print("\nset   N WT   N mut   median p")
for i in range(n_sets):
    print(i + 1, n_wt[i], n_mut[i], round(np.median(np.concatenate([p_wt[i], p_mut[i]])), 3))


def weighted_mean_std(values, weights):
    # weighted mean + weighted std over the 10 sets (rows), for every bin (columns) at once
    #  weight = number of bacteria set has in specific bin
    total = np.maximum(weights.sum(axis=0), 1)        # sum of the weights (at least 1)
    n = np.maximum((weights > 0).sum(axis=0), 2)      # how many sets actually have bacteria in the bin
    mean = (weights * values).sum(axis=0) / total     
    # weighted variance. mean from same data so divide by n-1 
    var = (weights * (values - mean) ** 2).sum(axis=0) / total * n / (n - 1)
    return mean, np.sqrt(var)


# histograms
# making overflow bin (P > last edge)
edges = np.append(bins, all_p.max() + 1)
# np.histogram counts how many momenta fall in each bin
# table: 10 rows (sets) x number of bins (columns), last column = overflow
counts_wt = np.array([np.histogram(x, bins=edges)[0] for x in p_wt])
counts_mut = np.array([np.histogram(x, bins=edges)[0] for x in p_mut])


# asymmetry per bin 
# A = (N_WT - N_mut) / (N_WT + N_mut)
#   A > 0 -> more WT than mutant,  A = 0 -> equal,  A < 0 -> more mutant
# calculated for every set and every bin at the same time (numpy)
N = counts_wt + counts_mut
A_sets = (counts_wt - counts_mut) / np.maximum(N, 1)   # max(N,1)
A, A_std = weighted_mean_std(A_sets, N)                # average of the 10 sets, per bin
# A_std =  error bar
# error_div is 1 (we use the SD)
A_err = A_std / error_div

# overall asymmetry: all momenta together
# reshape(-1, 1) makes it a column so it fits in the same function as above
A_tot_sets = ((n_wt - n_mut) / (n_wt + n_mut)).reshape(-1, 1)
A_tot, A_tot_std = weighted_mean_std(A_tot_sets, (n_wt + n_mut).reshape(-1, 1))
A_tot, A_tot_err = A_tot[0], A_tot_std[0] / error_div   # [0] to get a normal number out of the array


# does A depending on momentum? use a straight line A = a*P + b through points
#   slope a = 0 -> no momentum dependence
#   slope a clearly not 0 -> A changes with momentum
#   the overflow bin is left out: it has no upper edge, so it has no real "middle" momentum
#   for the fit chi^2 use the SEM (SD / sqrt(10))
#   check if a flat line (A = overall A everywhere) fits, with chi^2
#   chi^2 = sum of ((point - flat line) / error)^2
#   if  flat line fits, every point is 1 error away -> chi^2 ≈ number of points
#   chi^2 a lot bigger than the number of points -> flat doesn't fit -> A depends on P
width = np.diff(bins)                       # width of every normal bin
centers = (bins[:-1] + bins[1:]) / 2        # middle of every normal bin (x-value of the points)
A_fit = A[:-1]                              # ignore overflow bin
err_fit = A_std[:-1] / np.sqrt(n_sets)      # SEM of every bin, only for the fit
ok = err_fit > 0                            # skip empty bins
slope = None
if ok.sum() >= 3:                           # need at least 3 points 
    (slope, intercept), cov = np.polyfit(centers[ok], A_fit[ok], 1, w=1 / err_fit[ok], cov="unscaled")
    slope_err = np.sqrt(cov[0, 0])          # error on the slope = sqrt of the variance on the diagonal
    print(f"\nline fit: slope = {slope:.2e} +- {slope_err:.2e} per 10^-20 kg m/s"
          f"  -> {abs(slope) / slope_err:.1f} sigma from 0")
chi2_flat = np.sum(((A_fit[ok] - A_tot) / err_fit[ok]) ** 2)
print(f"flat line (A = {A_tot:.5f} everywhere): chi^2 = {chi2_flat:.1f} for {ok.sum()} points"
      f"  (should be around {ok.sum()} if A doesn't depend on P)")


# dN/dP = number of bacteria in a bin / width of  bin, per event
# divide by width because the bins aren't equally wide
# convert from a frequency to a density 
dNdP_wt = counts_wt[:, :-1] / width / events_per_set
dNdP_mut = counts_mut[:, :-1] / width / events_per_set
# normal mean (every set counts equally)
wt_mean, wt_err = dNdP_wt.mean(axis=0), dNdP_wt.std(axis=0, ddof=1) / error_div
mut_mean, mut_err = dNdP_mut.mean(axis=0), dNdP_mut.std(axis=0, ddof=1) / error_div


# table
# one row per bin. A_over_err = how many error bars A is away from 0
# rule of thumb we use: > 3 means a real asymmetry, < 3 could just be chance
labels = [f"{bins[b]:g}-{bins[b + 1]:g}" for b in range(len(bins) - 1)] + [f">{bins[-1]:g}"]
table = pd.DataFrame({
    "P_range": labels,
    "dNdP_WT": np.append(wt_mean, np.nan), "dNdP_WT_err": np.append(wt_err, np.nan),
    "dNdP_mutant": np.append(mut_mean, np.nan), "dNdP_mutant_err": np.append(mut_err, np.nan),
    "A": A, "A_err": A_err, "A_over_err": A / A_err,
})
print()
print(table.round(5).to_string(index=False))
table.to_csv(f"asymmetry-{wt}.csv", index=False)

#ANSWERING ALL RQ
# research question 1: average number of each strain per subsample, with its error
print(f"\naverage count per set: WT {n_wt.mean():.1f} +- {n_wt.std(ddof=1) / error_div:.1f}, "
      f"mutant {n_mut.mean():.1f} +- {n_mut.std(ddof=1) / error_div:.1f}  ({err_name})")
# research question 2: the overall asymmetry (if bigger than 3sigma significant)
print(f"overall A = {A_tot:.5f} +- {A_tot_err:.5f} ({err_name})"
      f"  -> {A_tot / (A_tot_std[0] / np.sqrt(n_sets)):.1f} sigma from 0 (using SEM)")

names = {211: "E. coli", 321: "B. subtilis", 2212: "P. aeruginosa",
         3122: "S. pneumoniae", 3312: "M. tuberculosis", 3334: "Salmonella"}
unit = "$10^{-20}$ kg m/s"   # write -20 as a superscript


# plot 1: A vs momentum (RQ 3)
# overflow bin has no end, so just as big as previous bin (in grey)

last_w = bins[-1] - bins[-2]
x = np.append(centers, bins[-1] + last_w / 2)    # x position of every point (incl. overflow)
xerr = np.append(width, last_w) / 2              # horizontal bar = how wide the bin is (half each side)

fig, ax = plt.subplots(figsize=(9, 6))
ax.axvspan(bins[-1], bins[-1] + last_w, color="0.93")   # grey area for the overflow bin
ax.text(bins[-1] + last_w / 2, 0.97, f"P > {bins[-1]:g}", transform=ax.get_xaxis_transform(),
        ha="center", va="top", color="0.4")             # x in momentum, y as fraction of the plot height
# band = overall A ± its error. 5 decimals
ax.axhspan(A_tot - A_tot_err, A_tot + A_tot_err, color="orchid", alpha=0.3,
           label=f"overall A = {A_tot:.5f} ± {A_tot_err:.5f}")
ax.axhline(0, color="red", ls="--", label="A = 0 (no asymmetry)")   
#  fitted line, only over the normal bins 
if slope is not None:
    p_line = np.array([bins[0], bins[-1]])
    ax.plot(p_line, slope * p_line + intercept, color="tab:blue", lw=1.5,
            label=f"line fit: slope = {slope:.1e} ± {slope_err:.1e} ({abs(slope) / slope_err:.1f}σ)")
# the actual points: vertical bar = SD of the 10 sets, horizontal bar = bin width
ax.errorbar(x, A, yerr=A_err, xerr=xerr, fmt="o", capsize=4, color="tab:green",
            label=f"mean A of 10 subsamples ± {err_name}")
ax.set_xticks(bins[bins % 1 == 0])       # only whole numbers on the axis to avoid overlap
lo, hi = ax.get_ylim()
ax.set_ylim(lo, hi + 0.4 * (hi - lo))    # extra space on top so no legend covering
ax.set_xlabel(f"momentum P ({unit})")
ax.set_ylabel("A = (N$_{WT}$ − N$_{mut}$) / (N$_{WT}$ + N$_{mut}$)")
ax.set_title(f"{names[wt]}: WT ({wt}) vs mutant ({mut})")
ax.legend(loc="upper left")
fig.tight_layout()
fig.savefig(f"asymmetry-{wt}.png", dpi=200)


# plot 2: dN/dP of WT and mutant 
# shows how the bacteria are spread over the momenta
# log scale to avoid data skewness (a lot more for smaller values)

fig2, bx = plt.subplots(figsize=(9, 6))
bx.errorbar(centers, wt_mean, yerr=wt_err, xerr=width / 2, fmt="o", capsize=4,
            color="orchid", label=f"WT ({wt})")
# mutant = hollow and a bit bigger, because WT and mutant are almost on top of each other
bx.errorbar(centers, mut_mean, yerr=mut_err, xerr=width / 2, fmt="o", ms=10, mfc="none", capsize=4,
            color="tab:green", label=f"mutant ({mut})")
bx.set_yscale("log")
bx.set_xticks(bins[bins % 1 == 0])
bx.set_xlabel(f"momentum P ({unit})")
bx.set_ylabel(f"dN/dP (bacteria per event per {unit})")
bx.set_title(f"{names[wt]}: dN/dP, mean of 10 subsamples ± {err_name}")
bx.legend()
fig2.tight_layout()
fig2.savefig(f"spectrum-{wt}.png", dpi=200)

plt.show()
