# MinLinMo
Fast and efficient *n<<p* variable selection for linear models

## Compiling and installing MinLinMo
MinLinMo was written in C++ version 14. It has been developed for both Intel and ARM
processors. For Intel processors, both multi-threading and AVX2 are employed for maximum
performance. For ARM processors, NEON is used instead of AVX2. As AVX2 supports 256 bit registers (NEON currently only supports 128 bit), parts
of MinLinMo will likely run faster on Intel processors. During the model building phase, which is
also the most time consuming phase, multi-threading did not improve performance.

**More detailed instructions for compiling MinLinMo on a specific platform can be found in the corresponding directory above.**

To compile MinLinMo for Ubuntu-type Linux on **Intel CPUs with AVX2** make sure that the **GNU Scientific Library** (GSL) is
installed. This can be done by typing the following in the command line:

 ```sudo apt install libgsl-dev```

MinLinMo can now be compiled with:

 ```g++ MinLinMo.cpp -march=core-avx2 -lgsl -lblas -lpthread -O3 -o MinLinMo -std=c++14```

 For **ARM Linux** (the GSL library can be installed with the same command as above on Ubuntu-type Linux):

 ```g++ MinLinMo.cpp -lpthread -std=c++14 -O3 -lgsl -lblas -lpthread -march=armv8-a -o MinLinMo```

For **Apple silicon OS X**, use **Homebrew** to install GSL:

 ```brew install gsl```

MinLinMo can now be compiled using:

 ```g++ MinLinMo.cpp -lpthread -std=c++14 -O3 -lgsl -lblas -I /opt/homebrew/Cellar/gsl/2.7.1/include/ -L /opt/homebrew/Cellar/gsl/2.7.1/lib/ -march=armv8-a -o MinLinMo```

For **Intel OS X** follow the instructions for Intel-based Linux above. This will hopefully work, but the SIMD library could be different (has not been tested):

 ```g++ MinLinMo.cpp -lpthread -std=c++14 -O3 -lgsl -lblas -I /opt/homebrew/Cellar/gsl/2.7.1/include/ -L /opt/homebrew/Cellar/gsl/2.7.1/lib/ -march= core-avx2```

For **Windows (and Intel-based Macs with Boot camp)**, the easiest way to run MinLinMo is to download the _.exe_ and _.dll_ files from the Windows directory and make sure the files are in the same directory.

To compile MinLinMo from scratch, it will likely be easiest to install GSL through Visual Studio 2017 or later. This can be done
as follows:

- Create a new C++ project, right click the project name under Solution Explorer and choose
“Manage Nuget Packages”
- Click the Browser tab and search with the keywords “Microsoft.gsl”, which will filter out the Microsoft GSL version. It can now be installed by clicking the install button

MinLinMo can now be compiled with the Microsoft C++ compiler from Visual Studio, see Windows section for full details.

## Testing MinLinMo ##

MinLinMo default settings are: predictors must correlate 0.1 with the outcome, the linear model must improve with at least 1% (R<sup>2 </sup>) variance explaned and an added predictor must correlate 0.1 with current model residuals to be considered in the linear model. 

To test MinLinMo on a dataset, first download the included _mtcars_ dataset (from the R package _data(mtcars)_). That is, the files _mpg.txt_ (outcome) and _mtcars.csv_ (predictors). For Linux and OS X, run MinLinMo from a path to the directory it was installed in or in the current directory:

 ```
~$ ./MinLinMo mpg.txt mtcars.csv
=============================
Simple variable selector v1.0
=============================
Number of arguments:3
Outcome vector filename entered:mpg.txt
Matrix filename entered:mtcars.csv
Correlation entered:0.1
Delta Rsq entered:0.01
Residual correlation entered:0.1
Loaded outcome:
Number of rows:32
Loaded matrix data:
Number of rows:32	Number of columns:10
Calculating correlation matrix
Time taken (sec/millisec):0.00213961
Number of predictors correlating with outcome:10
Computing model...
39.6863 intercept
-3.19097 wt
-1.50779 cyl
Final model R2:0.818519
Number of predictors:2
Time taken (sec/millisec):0.000248568
```

 For Windows, the _.dll_ files must be in the same directory as the _MinLinMo.exe_ file. MinLinMo can then be run with a path to the directory it was installed, or in the current by typing:

 ```MinLinMo.exe mpg.txt mtcars.csv```

 To adjust only the parameter for Pearson correlation between predictors and outcome type:

 ```./MinLinMo mpg.txt mtcars.csv 0.1```

 To also adjust R<sup>2</sup> to 0.1% variance explaned increase:

 ```./MinLinMo mpg.txt mtcars.csv 0.2 0.001```

 The residual correlation can be adjusted as well:

 ```
 ./MinLinMo mpg.txt mtcars.csv 0.2 0.001 0.5
=============================
Simple variable selector v1.0
=============================
Number of arguments:6
Outcome vector filename entered:mpg.txt
Matrix filename entered:mtcars.csv
Correlation entered:0.2
Delta Rsq entered:0.001
Residual correlation entered:0.5
Loaded outcome:
Number of rows:32
Loaded matrix data:
Number of rows:32	Number of columns:10
Calculating correlation matrix
Time taken (sec/millisec):0.00214258
Number of predictors correlating with outcome:10
Computing model...
19.7462 intercept
-5.04798 wt
0.929198 qsec

Final model R2:0.814445

Number of predictors:2
Time taken (sec/millisec):0.000223395
```

## Preparing datasets for MinLinMo ##

To analyse datasets, MinLinMo will require a file with the outcome which is just a list of numbers with a title on top. In addition, MinLinMo will require a comma-separated (_.csv_-type file) prediction matrix with column names and without rownames.

Here is an example of a random dataset consisting of 50,000 columns and 1,000 rows created with **Python** and the **Pandas** library:
```
import numpy as np
import pandas as pd
numcols = 50000
numrows = 1000
matrix = np.random.normal(0, 1, (numrows, numcols))
headers = ["V" + str(i) for i in range(0, numcols)]
df = pd.DataFrame( matrix, columns=headers)
df.to_csv("random_matrix.csv", index=False)

# Outcome vector
vector = np.random.normal(0,1,numrows)
header = ["Outcome"]
outc = pd.DataFrame(vector, columns=header)
outc.to_csv("Outcome.csv", index=False)
```
Or in **R**, using the **data.table** library for fast storage:

```
library(data.table)
numcols <- 50000
numrows <- 1000
matr <- matrix( rnorm( numcols * numrows ), ncol=numcols, nrow=numrows)
headers <- paste("pred",1:numcols, sep="")
colnames(matr) <- headers
df <- as.data.frame( matr )
fwrite(matr, file="random_matrix.csv", quote=F)

vect <- rnorm( numrows )
outc <- as.data.frame( vect )
names(outc) <- "Outcome"
fwrite(outc, file="Outcome.csv", quote=F)
```

Testing the dataset with MinLinMo:

```
$ ./MinLinMo Outcome.csv random_matrix.csv 
=============================
Simple variable selector v1.0
=============================
Number of arguments:3
Outcome vector filename entered:Outcome.csv
Matrix filename entered:random_matrix.csv
Correlation entered:0.1
Delta Rsq entered:0.01
Residual correlation entered:0.1
Loaded outcome:
Number of rows:1000
Loaded matrix data:
Number of rows:1000	Number of columns:50000
Calculating correlation matrix
Time taken (sec/millisec):0.036771
Number of predictors correlating with outcome:74
Computing model...
-0.0281351 intercept
-0.0966533 V49059
0.099808 V23010
0.123745 V19734
-0.0958096 V28082
-0.120931 V47759
0.122009 V42874
0.115737 V18680
0.129641 V16339
0.100093 V494
0.128163 V16052
-0.099705 V23107
0.105727 V15277
-0.0890509 V23105
0.106072 V25399
-0.106206 V21016
0.109707 V19123
0.122541 V26187
0.104949 V17961

Final model R2:0.218852

Number of predictors:18
Time taken (sec/millisec):0.152509
```

## Prediction with a MinLinMo trained model ##

MinLinMo's selected variables can either be employed in a separate analysis or they can be used, together with the estimated coefficients, to predict the outcome on a test dataset. To do the latter, make a _.txt_ file with the coefficients and their corresponding variables by copying them from the screen and pasting them into a text editor. Then add columns titles, such as _est_ and _id_ (space separated), as seen below:

```
est id
-0.0281351 intercept
-0.0966533 V49059
0.099808 V23010
0.123745 V19734
-0.0958096 V28082
-0.120931 V47759
0.122009 V42874

```
Save the estimates and call the text file _ests.txt_. Prediction is likely easier in R, so start by loading the MinLinMo estimates:

```
ests <- read.table("ests.txt", header=T, sep=" ")
```
Make sure that the column names are strings and not factors:

```
ests$id <- as.character(ests$id)
```

Load your test dataset, containing the predictor matrix (e.g. _mldat_), into R and make sure that the variable names of your previously loaded estimates (_ests$id_) are the same as those in the predictor matrix (except _intercept_, of course).

The outcome can now be predicted as follows:

```
pred <- cbind( 1, as.matrix( mldat[,ests$id[-1]] ) ) %*% ests$est
```
The variable _pred_ should now contain a vector of the predicted values that can be compared to the given values if present.
