import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import matplotlib.ticker as ticker

#data shift
yshift = 1.34e6

# -------------------------------
# Load .dat file (auto-detect columns)
# -------------------------------
def load_dat_file(filename):
    try:
        data = np.loadtxt(filename, comments=['#', '%'])
        ncols = data.shape[1]

        x = data[:, 0]
        y = data[:, 1]

        xerr = None
        yerr = None

        if ncols == 3:
            yerr = data[:, 2]

        elif ncols >= 4:
            xerr = data[:, 2]
            yerr = data[:, 3]

        return x, y, xerr, yerr

    except Exception as e:
        print(f"Error reading {filename}: {e}")
        return None, None, None, None

# -------------------------------
# Main
# -------------------------------
def main():
    parser = argparse.ArgumentParser(description="Plot .dat files")
    files = ["Parent_tf1.dat","Daughter_tf1.dat","BetaN_Daughter_tf1.dat","expbg_tf1.dat","Graph_gr.dat","fit_tf1.dat"] 
    labels = [r"$^{44}\mathrm{S}$", r"$^{44}\mathrm{Cl}$", r"$^{43}\mathrm{Cl}$", "Background", "Exp", "Total fit"] 
    color = ["blue", "orange", "green", "brown", "red"]
    parser.add_argument('--logx', action='store_true')
    parser.add_argument('--logy', action='store_true')
    parser.add_argument('--output', default=None)

    args = parser.parse_args()

    plt.figure()
    i = 0
    for file in files:
        x, y, xerr, yerr = load_dat_file(file)
        if x is None:
            continue

        label = labels[i]
        # -----------------------
        # Decide plot type
        # -----------------------
        if(file == "Graph_gr.dat" or file == "fit_tf1.dat"):
            if yerr is not None:
                 plt.errorbar(x, y, yerr = yerr, fmt='.',linestyle='none',alpha= 1, markersize=7, zorder=1, capsize=4, label=label, color = color[i])
            else:
                plt.plot(x, y, label=label)  
        else:
            if(label == "Background"):
                plt.plot(x, y, label= label, linestyle = 'dotted', color = color[i])
            else:
                plt.plot(x, y + yshift, label= label, color = color[i])
        i += 1

    # -------------------------------
    # Styling
    # -------------------------------
    # labels_map = {
    # "Graph_gr.dat": r"$^{44}\mathrm{S}$ data",
    # "run1.dat": "Background",
    # "run2.dat": "Gate 1"
    # }
    plt.xlabel("Time (ms)")
    plt.ylabel("Decay/ms")
    plt.title(r"$^{44}\mathrm{S}$ data")
    plt.legend(loc="lower right")
    
    ax = plt.gca()

    ax.minorticks_on()

    ax.tick_params(which='both', direction='in', top=False, right=False)
    ax.tick_params(which='major', length=6)
    ax.tick_params(which='minor', length=3)
    ax.text(
        0.98, 0.95,
        r"$^{44}\mathrm{S}$",
        transform=ax.transAxes,
        ha="right",
        va="top",
	fontsize=10,
        #bbox=dict(boxstyle="round", alpha=0.6)
     )

    ax.text(
        0.98, 0.85,
        r"$T_{1/2}$ = 120.8(15)ms",
        transform=ax.transAxes,
        ha="right",
        va="top",
	fontsize=8,
       #bbox=dict(boxstyle="round", alpha=0.6)
    )

    plt.yscale('log', nonpositive='clip')

    plt.xlim(0, 4000)
    plt.ylim(1.339e6, 1.42e6)

    if args.output:
        plt.savefig(args.output, dpi=600)
        print(f"Saved plot to {args.output}")
    else:
        plt.show()


# -------------------------------
if __name__ == "__main__":
    main()
