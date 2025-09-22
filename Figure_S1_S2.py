import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

##########
# This script plots violin plots of the VAF data for Figures S1 and S2.
# VAF data files are included in the data/ folder as well as in the Supplementary Files in the publication.
##########

# List the Group IDs
groups = ["A1", "A2", "A3", "A4", "A5", "B1", "B2", "B3", "B4", "B5"]
# Path to the VAF data. Assumes that group ID is in the {}.
data_path = "data/VAF/{}_VAF.csv"
# File name template in which to save the figures. Assumes that group ID is in the {}.
save_fig_path = "figures/{}_vaf.svg" 

# Loop through the groups and plot
for group in groups:
    df = pd.read_csv(data_path.format(group))
    # If there are any columns that don't have data, set them to zero for the plotting.
    for c in df.columns:
        if df[c].isnull().all():
            df[c] = 0
    # Plot the violin plot
    plt.figure()
    sns.violinplot(df, inner="quart", bw_adjust=0.5, fill=True, linecolor="black", linewidth=.5)
    plt.tight_layout()
    # Save and show
    if save_fig_path:
        plt.savefig(save_fig_path.format(group))
    plt.show()