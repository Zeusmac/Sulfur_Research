import uproot
import numpy as np
import matplotlib.pyplot as plt
import argparse

def get_histogram(file, hist_name):
    """
    Extract TH1 histogram data from ROOT file
    """
    hist = file[hist_name]

    values = hist.values()
    edges = hist.axis().edges()

    # Convert to bin centers
    x = 0.5 * (edges[:-1] + edges[1:])
    y = values

    return x, y, edges


def get_tf1_functions(file):
    """
    Extract all TF1 objects from ROOT file
    """
    tf1_list = []

    for key in file.keys():
        obj = file[key]

        # uproot identifies TF1 as "TF1"
        if obj.classname == "TF1":
            tf1_list.append(obj)

    return tf1_list


def evaluate_tf1(tf1, x):
    """
    Evaluate TF1 manually using its formula and parameters
    """
    formula = tf1.member("fFormula")
    params = tf1.member("fParams")

    # Replace ROOT-style parameters [0], [1], ... with actual values
    expr = formula
    for i, p in enumerate(params):
        expr = expr.replace(f"[{i}]", str(p))

    # Evaluate safely
    y = eval(expr, {"x": x, "np": np, "exp": np.exp, "sin": np.sin, "cos": np.cos})

    return y


def main():
    parser = argparse.ArgumentParser(description="Plot ROOT histogram with TF1 overlays")

    parser.add_argument("file", help="ROOT file")
    parser.add_argument("hist", help="Histogram name")

    parser.add_argument("--xlabel", default="ms, 1ms per bin")
    parser.add_argument("--ylabel", default="Counts")
    parser.add_argument("--title", default="44S")
    parser.add_argument("--logy", action="store_true")

    args = parser.parse_args()

    root_file = uproot.open(args.file)

    # --- Histogram ---
    x, y, edges = get_histogram(root_file, args.hist)

    plt.figure()

    plt.step(x, y, where='mid', label=args.hist)

    # --- TF1 functions ---
    tf1s = get_tf1_functions(root_file)

    x_dense = np.linspace(edges[0], edges[-1], 1000)

    for tf1 in tf1s:
        try:
            y_fit = evaluate_tf1(tf1, x_dense)
            name = tf1.name
            plt.plot(x_dense, y_fit, label=name)
        except Exception as e:
            print(f"Skipping {tf1.name}: {e}")

    # --- Labels ---
    plt.xlabel(args.xlabel)
    plt.ylabel(args.ylabel)
    plt.title(args.title)
    plt.legend()

    if args.logy:
        plt.yscale("log")

    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    main()
