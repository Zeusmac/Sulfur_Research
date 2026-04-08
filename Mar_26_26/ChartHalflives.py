import os
import re
import pandas as pd
import matplotlib.pyplot as plt
folder_path = "BestFits"

results = []

for filename in os.listdir(folder_path):
    if not filename.endswith(".txt"):
        continue

    filepath = os.path.join(folder_path, filename)

    with open(filepath, "r", encoding="utf-8") as f:
        text = f.read()

    section_match = re.search(
    	r"-+\s*Half-lives\s*-+(.*?)(?:\n\s*-{5,}|\Z)",
    	text,
   	re.DOTALL
	)
    if not section_match:
        print(f"❌ Section not found in {filename}")
        continue

    section = section_match.group(1)
    match = re.search(
        r"lambda_parent\s+([0-9.eE+-]+)\s+([-0-9.eE+]+)\s+\+\-([0-9.eE+-]+)",
        section
    )

    if match:
        value = float(match.group(1))
        err_minus = float(match.group(2))
        err_plus = float(match.group(3))

        results.append((filename, value, err_minus, err_plus))
    else:
        print(f"❌ lambda_p not found in {filename}")

    # Convert to DataFrame
    df = pd.DataFrame(results, columns=["File Name", "lambda_p"])

    # Convert to DataFrame
    df = pd.DataFrame(results, columns=["File Name", "lambda_p"])

    # Print table
    print(df)

    # Save to CSV
    df.to_csv("lambda_p_summary.csv", index=False)

    # Plot
    plt.figure()
    plt.bar(df["File Name"], df["lambda_p"])
    plt.xticks(rotation=90)
    plt.xlabel("File Name")
    plt.ylabel("lambda_p (ms)")
    plt.title("lambda_p from files")
    plt.tight_layout()
    plt.show()
