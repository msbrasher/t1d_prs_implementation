# PRS Ancestry Calibration
### Overview
This is an R-function to regress the effects of genetic principal components out of the mean and variance of calculated polygenic risk scores, based on methods from the eMERGE network.

---

### Input & Output
**Input**

The function requires a data frame in R with at least a numeric column of PRS values and any number of columns of numeric genetic principal component values. The data frame can contain other columns, and columns can be in any order.

For example:

<img width="378" height="71" alt="Screenshot 2025-07-28 at 11 58 34 AM" src="https://github.com/user-attachments/assets/5523966e-24d2-4509-aa48-d58e50865df1" />

**Output**

The function will output a data frame that is the same as the input data frame, but with two additional columns appended. 

*calib1* = contains the PRS values with PCs regressed out of the mean

*calib2* = contains the PRS values with PCs regressed out of both the mean and variance



---

### Example
```
source("ancestry_calib_function.R")

data <- read.table("input_file.txt", header = T)

data_calib <- ancestry_calib(data, pc_prefix = "pc", num_pcs = 10, score_col = "SCORE")
```

---

### Citation

To cite this work please cite our recent publication: 

Brasher, M. S. et al. Enabling reproducible type 1 diabetes polygenic risk scoring for clinical and translational applications. Preprint at https://doi.org/10.1101/2025.07.15.25331523 (2025).

As well as the original work proposing this method:

Khera, A. V. et al. Whole-Genome Sequencing to Characterize Monogenic and Polygenic Contributions in Patients Hospitalized With Early-Onset Myocardial Infarction. Circulation 139, 1593–1602 (2019).

Khan, A. et al. Genome-wide polygenic score to predict chronic kidney disease across ancestries. Nat. Med. 28, 1412–1420 (2022).
