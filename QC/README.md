## Raw QC stats sheet 

Step_1: Navigate to the directory containing all R1 and R2 FASTQ.gz files received from the sequencing centre for the batch. Use the below mentioned command to extract file information in a .tsv files.

```
ls -lh > sample_list.tsv
```

Step_2: [file_metrics_step1.py](file_metrics_step1.py) - Run the script to extract R1 and R2 files sizes from sample_list.tsv into sample_list_paired.tsv.

Step_3: [fetch_file_info_step1.py](fetch_file_info_step1.py) - Fetches and returns values of Total Sequences and Sequence Lengths From Basic stats table for each sample. Make sure that the fastqc.html files of all the required samples are present in the source directory.
