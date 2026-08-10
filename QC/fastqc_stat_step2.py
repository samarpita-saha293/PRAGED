import os
import csv
import pandas as pd

# Current directory containing all FastQC HTML files
html_dir = "."

# Output file in the same directory
output_file = "fastqc_summary_paired.tsv"

samples = {}

# Find all FastQC HTML files in the current directory
html_files = [
    f for f in os.listdir(html_dir)
    if f.endswith(".html")
]

for html_file in sorted(html_files):

    html_path = os.path.join(html_dir, html_file)

    try:
        # Read all tables from the FastQC HTML
        tables = pd.read_html(html_path)

    except Exception as e:
        print(f"Could not read {html_file}: {e}")
        continue

    # Find the Basic Statistics table
    basic_stats = None

    for table in tables:
        # FastQC Basic Statistics table has "Measure" and "Value"
        if "Measure" in table.columns and "Value" in table.columns:
            basic_stats = table
            break

    if basic_stats is None:
        print(f"Basic Statistics table not found: {html_file}")
        continue

    # Convert the table into a dictionary
    stats = dict(
        zip(
            basic_stats["Measure"].astype(str).str.strip(),
            basic_stats["Value"].astype(str).str.strip()
        )
    )

    # Get the filename reported by FastQC
    filename = stats.get("Filename", "")

    if not filename:
        print(f"Filename not found in: {html_file}")
        continue

    # Determine R1/R2
    if "_R1.fastq.gz" in filename:
        sample = filename.replace("_R1.fastq.gz", "")

        samples.setdefault(sample, {})

        samples[sample]["R1_Total_Sequences"] = stats.get(
            "Total Sequences", ""
        )

        samples[sample]["R1_Sequence_Length"] = stats.get(
            "Sequence length", ""
        )

    elif "_R2.fastq.gz" in filename:
        sample = filename.replace("_R2.fastq.gz", "")

        samples.setdefault(sample, {})

        samples[sample]["R2_Total_Sequences"] = stats.get(
            "Total Sequences", ""
        )

        samples[sample]["R2_Sequence_Length"] = stats.get(
            "Sequence length", ""
        )

    else:
        print(f"Skipping non-R1/R2 file: {filename}")


# Write paired summary
with open(output_file, "w", newline="") as outfile:

    writer = csv.writer(outfile, delimiter="\t")

    writer.writerow([
        "Sample",
        "R1_Total_Sequences",
        "R2_Total_Sequences",
        "R1_Sequence_Length",
        "R2_Sequence_Length"
    ])

    for sample in sorted(samples):

        writer.writerow([
            sample,
            samples[sample].get("R1_Total_Sequences", ""),
            samples[sample].get("R2_Total_Sequences", ""),
            samples[sample].get("R1_Sequence_Length", ""),
            samples[sample].get("R2_Sequence_Length", "")
        ])

print(f"\nDone!")
print(f"Output written to: {output_file}")
