import os
import csv

# ---------------------------------------------------------
# Ask whether the data is Genome or Exome
# ---------------------------------------------------------

while True:
    choice = input(
        "Is this Genome or Exome data?\n"
        "A. Genome\n"
        "B. Exome\n"
        "Enter A or B: "
    ).strip().upper()

    if choice in ["A", "B"]:
        break

    print("Invalid choice. Please enter A or B.")

if choice == "A":
    target_region_size = 3_000_000_000
    data_type = "Genome"
else:
    target_region_size = 37_453_133
    data_type = "Exome"


# ---------------------------------------------------------
# Input files
# ---------------------------------------------------------

adapter_file = "adapter_summary.tsv"
fastqc_file = "fastqc_summary_paired.tsv"
sample_list_file = "sample_list_paired.tsv"

output_file = "final_summary.tsv"


# ---------------------------------------------------------
# Check files
# ---------------------------------------------------------

for file in [adapter_file, fastqc_file, sample_list_file]:

    if not os.path.isfile(file):
        print(f"Error: {file} not found in the current directory.")
        exit()


# ---------------------------------------------------------
# Function to read TSV into dictionary using Sample
# ---------------------------------------------------------

def read_tsv(filename):

    data = {}

    with open(filename, "r", newline="") as f:

        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:

            sample = row["Sample"]

            data[sample] = row

    return data


# ---------------------------------------------------------
# Read all three files
# ---------------------------------------------------------

adapter_data = read_tsv(adapter_file)
fastqc_data = read_tsv(fastqc_file)
sample_list_data = read_tsv(sample_list_file)


# ---------------------------------------------------------
# Merge all samples
# ---------------------------------------------------------

all_samples = set()

all_samples.update(adapter_data.keys())
all_samples.update(fastqc_data.keys())
all_samples.update(sample_list_data.keys())


# ---------------------------------------------------------
# Define column order
# ---------------------------------------------------------

adapter_columns = [
    "Adapter_Beginning(x,y)",
    "Adapter_End(x,y)"
]

fastqc_columns = [
    "R1_Total_Sequences",
    "R2_Total_Sequences",
    "R1_Sequence_Length",
    "R2_Sequence_Length"
]

sample_list_columns = [
    "R1_File",
    "R1_Size",
    "R2_File",
    "R2_Size"
]


# ---------------------------------------------------------
# Write final TSV
# ---------------------------------------------------------

with open(output_file, "w", newline="") as out:

    writer = csv.writer(out, delimiter="\t")

    writer.writerow([
        "Sample",
        *adapter_columns,
        *fastqc_columns,
        *sample_list_columns,
        "Coverage"
    ])

    for sample in sorted(all_samples):

        adapter = adapter_data.get(sample, {})
        fastqc = fastqc_data.get(sample, {})
        sample_list = sample_list_data.get(sample, {})

        # ---------------------------------------------
        # Fetch FastQC values for coverage calculation
        # ---------------------------------------------

        try:

            r1_sequences = float(
                fastqc.get("R1_Total_Sequences", 0)
            )

            r2_sequences = float(
                fastqc.get("R2_Total_Sequences", 0)
            )

            r1_length = float(
                fastqc.get("R1_Sequence_Length", 0)
            )

            r2_length = float(
                fastqc.get("R2_Sequence_Length", 0)
            )

            # Total bases sequenced
            total_bases = (
                r1_sequences * r1_length
                +
                r2_sequences * r2_length
            )

            # Coverage
            coverage = total_bases / target_region_size

            coverage = f"{coverage:.2f}"

        except (ValueError, TypeError, ZeroDivisionError):

            coverage = ""

        # ---------------------------------------------
        # Write row
        # ---------------------------------------------

        writer.writerow([
            sample,

            adapter.get("Adapter_Beginning(x,y)", ""),
            adapter.get("Adapter_End(x,y)", ""),

            fastqc.get("R1_Total_Sequences", ""),
            fastqc.get("R2_Total_Sequences", ""),
            fastqc.get("R1_Sequence_Length", ""),
            fastqc.get("R2_Sequence_Length", ""),

            sample_list.get("R1_File", ""),
            sample_list.get("R1_Size", ""),
            sample_list.get("R2_File", ""),
            sample_list.get("R2_Size", ""),

            coverage
        ])


print("\nData type:", data_type)
print("Target region size:", target_region_size)
print("Written:", output_file)
