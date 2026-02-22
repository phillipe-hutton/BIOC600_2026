# Data Analysis in R and RStudio

Once you have processed your files, you will need to analyze them to be able to interpret your results. Python and R are the two most popular coding languages used by bioinformaticians and data scientists in general. Python is typically used for general-purpose programming and for building machine learning models, whereas R is used for statistical analysis and data visualization. In this workshop, we will be using R and the Bioconductor repository to analyze our processed genomics data.

## The RStudio IDE

To begin, open RStudio. The IDE looks like this:

---

![RStudioIDE](https://github.com/user-attachments/assets/d06742f0-587d-4a70-86b9-232d391acf67)

---

- The top-left panel is where you will view R scripts and data frames (data frames are tables of data that we will manipulate in R).
- The top-right panel shows variables and data frames as well as previous commands
- The bottom-right panel shows files and directories, packages, plots, and documentation for commands and packages
- The bottom-left panel is where the R console and the terminal are located (this is where the scripts are run)

## Data Wrangling with the tidyverse package

To begin, you will first install and load the tidyverse package:

```
#R

install.packages("tidyverse") # install tidyverse
library(tidyverse) # load the package
```

Here is a cheet sheet 
