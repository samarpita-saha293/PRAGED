import os
import zipfile
import csv
import re

output_file = "adapter_summary.tsv"

with open(output_file, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow([
        "Sample",
        "Adapter_Beginning(x,y)",
        "Adapter_End(x,y)"
    ])

    for zip_file in os.listdir("."):
        if not zip_file.endswith("_fastqc.zip"):
            continue

        with zipfile.ZipFile(zip_file) as z:

            data_file = [f for f in z.namelist()
                         if f.endswith("fastqc_data.txt")][0]

            content = z.read(data_file).decode()

        sample = os.path.basename(zip_file).replace("_fastqc.zip", "")

        in_adapter = False
        positions = []
        values = []

        for line in content.splitlines():

            if line.startswith(">>Adapter Content"):
                in_adapter = True
                continue

            if in_adapter and line.startswith(">>END_MODULE"):
                break

            if in_adapter:

                if line.startswith("#"):
                    continue

                cols = line.split("\t")

                if len(cols) < 2:
                    continue

                try:
                    x = cols[0]

                    y = max(float(v) for v in cols[1:])

                    positions.append(x)
                    values.append(y)

                except:
                    pass

        first = None

        for x, y in zip(positions, values):
            if y > 0:
                first = (x, y)
                break

        if first:

            last = (positions[-1], values[-1])

            writer.writerow([
                sample,
                f"({first[0]},{first[1]})",
                f"({last[0]},{last[1]})"
            ])

print("Written:", output_file)
