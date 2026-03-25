# PRAGED
Step by step instruction for analysing WGS and WES data on a routine basis are stored

Please make sure to use the latest versions of these scripts on HPC. The pipeline uses a fine-tuned ResNet-50 model on the [PneumoniaMNIST dataset](https://www.kaggle.com/datasets/rijulshr/pneumoniamnist/data) and is designed to work efficiently even on CPU setups.

### WGS and WES data Preprocessing
We use the PneumoniaMNIST dataset in .npz format with NumPy arrays for training, validation, and testing. To Obtain the .npz File go to the  [PneumoniaMNIST dataset](https://www.kaggle.com/datasets/rijulshr/pneumoniamnist/data) and place it in the root directory of this project.

__Please compile in order :__
1. [MissionProjectFileRename.sh](MissionProjectFileRename.sh) - The renaming script enforces filename conventions according to PRAGED standards.
2. [MissionProjectFileTransfer.sh](MissionProjectFileTransfer.sh) - This script moves files into Mission project directories and separates PRAGED data from unrelated files.
3. [utils.py](utils.py) - Utility functions for evaluation
4. [train.py](train.py) - Trains, fine-tunes the model and saves it
5. [eval.py](eval.py) - Loads the trained model and reports metrics
6. [eval_met.py](eval_met.py) - Prints out the model metrics
7. [eval_miscl.py](eval_miscl.py) - Generates misclassified images
8. [pr_curve.py](pr_curve.py) - Generates the Precision-Recall curve (Highlights performance under class imbalance by focusing on positive class predictions)
9. [eval_roc.py](eval_roc.py) - Generates the ROC curve

Lastly, __saved_model__ stores the trained model checkpoint.

## Code Description


First, load train/val/test images and labels from a NumPy zipped file.

```
data = np.load("pneumoniamnist.npz")
```
