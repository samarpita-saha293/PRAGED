import csv

input_file = "sample_list.tsv"
output_file = "sample_list_paired.tsv"

samples = {}

with open(input_file) as infile:
    for line in infile:
        line = line.strip()

        if not line:
            continue

        fields = line.split()

        # Ignore any line that doesn't look like ls -lh output
        if len(fields) < 9:
            print("Skipping:", line)
            continue

        filename = fields[-1]
        size = fields[4]

        if filename.endswith("_R1.fastq.gz"):
            sample = filename[:-len("_R1.fastq.gz")]
            samples.setdefault(sample, {})
            samples[sample]["R1_File"] = filename
            samples[sample]["R1_Size"] = size

        elif filename.endswith("_R2.fastq.gz"):
            sample = filename[:-len("_R2.fastq.gz")]
            samples.setdefault(sample, {})
            samples[sample]["R2_File"] = filename
            samples[sample]["R2_Size"] = size

with open(output_file, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["Sample","R1_File","R1_Size","R2_File","R2_Size"])

    for sample in sorted(samples):
        writer.writerow([
            sample,
            samples[sample].get("R1_File",""),
            samples[sample].get("R1_Size",""),
            samples[sample].get("R2_File",""),
            samples[sample].get("R2_Size","")
        ])

print("Done!")
