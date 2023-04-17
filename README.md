# lablaster: Laser Ablation BLASt Through Endpoint in R

<!-- badges: start -->
<!-- badges: end -->

## Description

Detect when the laser ablates through the target sample in a laser ablation mass spectrometry time resolved analysis.

This function imports a data frame containing a single time resolved laser ablation mass spectrometry analysis of a foraminifera (or other carbonate shell).
It assumes that the first row contains the signal of the target sample and that background correction has already been applied.

Column names referencing the time stamps and target signal are specified as function arguments.
Then the signal column is evaluated by smoothing the values using a moving average, then scaled between 0-1 and the rate of change over a number of observations.

The function identifies the maximum rate of signal change, locates the corresponding time stamp and then subtracts the time it took for the laser to blast through the target.
The result (endTime) is the time stamp of the last observation while the laser is still focussed the desired target.
Subsequently the function removes the rows of data that occur between the final laser shot of high signal counts and the maximum rate of change resulting with a data frame containing only the desired geochemical target. Visualisation tools for manual validation are also included.

## Installation

The easiest method to install lablaster is from CRAN

```r
install.packages("lablaster")
```
Or install the development version from GitHub

``` r
#install.packages("devtools")
#library(devtools)
install_github("alexsb1/lablaster")
```

## Arguments

**detectDf** A data frame containing a single time resolved analysis, with a column referencing time and another with the corresponding measured counts data.

**dt** An integer that controls the number of observations (rows) are used in calculating a rolling lagged difference in 44Ca signal. Using a lower value for a faster blast through of chamber wall can improve end point sensitivity. Default = 10.

**smoothing** Controls the length of the moving average filter to apply over the dt period. Default = 5.

**timeCol** The column title in the data frame identifying the time stamp of the time resolved analysis. Default = "Time".

**signalCol** The column title in the data frame identifying the numerical data used to identify the laser ablation blast through endpoint. This is often <sup>44</sup>Ca but could be any column of numerical data that you want to detect the endpoint. Default = "Ca44".

**profile** Logical. A visualisation of the endpoint detection mechanism as a ggplot2 object. To make this profile plot, set the argument to TRUE, otherwise set it to FALSE. Specifying to not make a plot can save a substantial amount of time. Default = "FALSE".

**timeUnits** The units that the time resolved analysis is measured in. This is the units of the timeCol. This argument is a string and is only necessary if the argument `profile = "TRUE"`. Default = "seconds".

## Output
The function returns a data frame containing the columns:

`dfReturn$df` contains a data frame with only the observations between the first data frame row and the endTime.

`dfReturn$startTime` contains the earliest time in your TRA as a numerical value.

`dfReturn$endTime` contains the last time step while the laser is still focussed the desired target in your TRA before the as a numerical value.

`dfReturn$profile` contains a visualisation of your TRA identifying where the laser ablated through the carbonate shell as a ggplot object. This is only available if a profile was generated using profile = "TRUE".


## Examples

* Example with every argument explicitly specified
``` r
library(lablaster)
# Example data: An antepenultimate chamber of Menardella exilis foraminifera 72, identified hereon as "foram72shot3".
endPoint(detectDf = foram72shot3, dt = 10, smoothing = 5, timeCol = "Time", signalCol = "Ca44", profile = "TRUE",  timeUnits = "seconds")
```

* The simplest usage, using all the default values. \
`endPoint(data1)`

* Specifying custom dt and smoothing arguments. \
`endPoint(data1, dt = 15, smoothing = 10)`

* Specifying a different dataframe column to use for detection by passing the "Mg24" dataframe heading as substitute for Ca44 and using default dt and smoothing. \
`endPoint(data1, signalCol = "Mg24")`

* Specifying to make a profile plot. This can take a substantial amount of time. \
_Note that if you use a time unit other than seconds you must specify these here._ \
`endPoint(data1, profile = "TRUE", timeUnits = "seconds")`


## Details
This function was designed for detecting when a laser had ablated through a foraminifera test and keeping only the outputted rows of a time resolved analysis that were relevant by returning a start time and end time of desired data.

Major and trace elements in the foraminifera test were analysed using a New Wave UP193 laser ablation system (ArF source, 30 µm spot diameter, 22% power and 5 Hz pulse rate) coupled to an Agilent 8900 Triple Quad (MS-QQQ) in single quadrupole mode using a He and Ar gas mixture (900 ml/min) at the University of Southampton. Each laser spot pattern was sequenced with a 15 s warmup, 50 s laser pulse and 15 s washout. The default values within the function are based on outputs from this setup.

## Reference

You should cite this function as

`Searle-Barnes et al., "Laser Ablation BLASt Through Endpoint in R", 2023, https://github.com/alexsb1/lablaster`

## Acknowledgements
Thomas Ezard, J Andy Milton, Chris Standish and Gavin Foster.

School of Ocean and Earth Science, National Oceanography Centre Southampton, European Way, Southampton, SO14 3ZH, UK

This work was supported by Natural Environment Research Council award NE/P019269/1.

## Licence
This is released under GNU General Public License v3.0
