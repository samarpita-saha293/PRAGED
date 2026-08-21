## Raw QC stats sheet 

Step_1: Navigate to the directory containing all R1 and R2 FASTQ.gz files received from the sequencing centre for the batch. Use the below mentioned command to extract file information in a .tsv files.

```
ls -lh > sample_list.tsv
```
Now, move sample_list.tsv into QC directory. 

Step_2: [file_metrics_step1.py](file_metrics_step1.py) - Extracst R1 and R2 files sizes from 
- Input: sample_list.tsv
- Output: sample_list_paired.tsv

Step_3: [fastqc_stat_step2.py](fastqc_stat_step2.py) - Fetches and returns values of Total Sequences and Sequence Lengths From Basic stats table for each sample. Make sure that the fastqc.html files of all the required samples are present in the source directory.
- Input: All html from source directory
- Output: fastqc_paired_summary.tsv

Step_4: [adapter_content_step3.py](adapter_content_step3.py) - Extracts adapter-content metrics from FASTQC html reports and generates a tsv summary containing the first and last detected coordinates for each sample
- Input: All html from source directory
- Output: adapter_summary.tsv

Step_5: [merge_info_step4.py](merge_info_step4.py) - Merges sample metadata from previously generated outputs and calculates per-sample sequencing coverage based on the selected genome or exome target region size.
- Input: sample_list_paired.tsv, fastqc_paired_summary.tsv, adapter_summary.tsv
- Output: final_summary.tsv


Additionally, you may use this command to fetch Total Sequences:
```
grep -oP '(?<=Total Sequences</td><td>)[0-9]+' fastqc_report.html
```
