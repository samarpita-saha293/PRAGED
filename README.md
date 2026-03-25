# PRAGED
__A complete guide for analysing PRAGED data on a routine basis__

Please use the latest versions of these scripts on the HPC. All the scripts are stored in the mentioned below directory:
```
/scratch/diag_ashwin/NGS_Databases_Tools/scripts
```

## WGS and WES data Preprocessing

__Please compile in order:__
1. [MissionProjectFileRename.sh](MissionProjectFileRename.sh) - The renaming script enforces filename conventions according to PRAGED standards.
2. [MissionProjectFileTransfer.sh](MissionProjectFileTransfer.sh) - This script moves files into Mission project directories and separates PRAGED data from unrelated files.
3. [exomescript.sh](exomescript.sh)/[genomescript.sh](genomescript.sh) - Utility functions for evaluation

__In place of exome_script you may use the following script as per your requirement:__

- [genomescript_withadapterseq.sh](genomescript_withadapterseq.sh) -
- [genomescript_withadapterseq.sh](genomescript_withadapterseq.sh) -
- [MissionExomeBatchRun.sh](MissionExomeBatchRun.sh) -
- [MissionExomeBatchRunwithadapterseq.sh](MissionExomeBatchRunwithadapterseq.sh) -
- [MissionExomeBatchRun_withouttrimgalore.sh](MissionExomeBatchRun_withouttrimgalore.sh) -

## Data upload to AWS

4. [train.py](train.py) - Trains, fine-tunes the model and saves it

## Data upload to Franklin


5. [eval.py](eval.py) - Loads the trained model and reports metrics
6. [eval_met.py](eval_met.py) - Prints out the model metrics
7. [eval_miscl.py](eval_miscl.py) - Generates misclassified images


