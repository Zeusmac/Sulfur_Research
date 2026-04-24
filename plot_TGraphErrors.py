import ROOT
import numpy as np
import matplotlib.pyplot as plt
import argparse


# -------------------------------
# Extract TGraphErrors from ROOT
# -------------------------------
def get_graph(root_file, name):
    obj = root_file.Get(name)

    if not obj:
        raise ValueError(f"Object {name} not found")

    if not obj.InheritsFrom("TGraphErrors"):
        raise TypeError(f"{name} is not a TGraphErrors")

    n = obj.GetN()

    x = np.array([obj.GetPointX(i) for i in range(n)])
    y = np.array([obj.GetPointY(i) for i in range(n)])
    ex = np.array([obj.GetErrorX(i) for i in range(n)])
    ey = np.array([obj.GetErrorY(i) for i in range(n)])

    return x, y, ex, ey


# -------------------------------
# Extract TF1 objects
# -------------------------------
def get_tf1s(root_file):
    tf1s = []

    for key in root_file.GetListOfKeys():
        obj = key.ReadObj()
        if obj.InheritsFrom("TF1"):
            tf1s.append(obj)

    return tf1s


# -------------------------------
# Compute T1/2 from exponential
# assumes lambda is parameter index 1
# -------------------------------
def compute_half_life(tf1):
    try:
        lam = tf1.GetParameter(1)
        lam_err = tf1.GetParError(1)

        if lam <= 0:
            return None

        t12 = np.log(2) / lam
        t12_err = (np.log(2) / lam**2) * lam_err

        return t12, t12_err

    except:
        return None


# -------------------------------
# Main
# -------------------------------
def main():
    parser = argparse.ArgumentParser(description="PyROOT plot with TF1 overlays")

    parser.add_argument("file", help="ROOT file")
    parser.add_argument("graph", help="TGraphErrors name")

    parser.add_argument("--xlabel", default="X Axis")
    parser.add_argument("--ylabel", default="Y Axis")
    parser.add_argument("--title", default="Fit with TF1")
    parser.add_argument("--logy", action="store_true")

    args = parser.parse_args()

    root_file = ROOT.TFile.Open(args.file)

    # -------------------------------
    # Graph
    # -------------------------------
    x, y, ex, ey = get_graph(root_file, args.graph)

    plt.figure()

    plt.errorbar(
        x, y,
        xerr=ex,
        yerr=ey,
        fmt='o',
        label=args.graph,
        capsize=2
    )

    # -------------------------------
    # TF1 overlays
    # -------------------------------
    tf1s = get_tf1s(root_file)

    text_lines = []

    x_dense = np.linspace(np.min(x), np.max(x), 1000)

    for tf1 in tf1s:
        # Draw function
        y_fit = np.array([tf1.Eval(xx) for xx in x_dense])
        plt.plot(x_dense, y_fit, label=tf1.GetName())

        # χ2 / NDF
        chi2 = tf1.GetChisquare()
        ndf = tf1.GetNDF()

        text_lines.append(f"{tf1.GetName()}: χ²/NDF = {chi2:.2f}/{ndf}")

        # Half-life (if exponential-like)
        t12_result = compute_half_life(tf1)
        if t12_result:
            t12, t12_err = t12_result
            text_lines.append(f"{tf1.GetName()}: T1/2 = {t12:.3f} ± {t12_err:.3f}")

    # -------------------------------
    # Styling
    # -------------------------------
    plt.xlabel(args.xlabel)
    plt.ylabel(args.ylabel)
    plt.title(args.title)
    plt.legend()

    if args.logy:
        plt.yscale("log")

    plt.grid(True)

    # -------------------------------
    # Add annotation box
    # -------------------------------
    text = "\n".join(text_lines)

    plt.gca().text(
        0.65, 0.95,
        text,
        transform=plt.gca().transAxes,
        fontsize=10,
        verticalalignment='top',
        bbox=dict(boxstyle="round", alpha=0.8)
    )

    plt.tight_layout()
    plt.show()


# -------------------------------
if __name__ == "__main__":
    main()
