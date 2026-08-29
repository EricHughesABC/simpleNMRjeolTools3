#A simple R script to invert spectrum point-by-point
#
# The rhdf5 library can be installed by running the following from the R prompt:
# > install.packages("BiocManager")
# > BiocManager::install("rhdf5")
# 
#How to use:
#In Jason, add "External command" to the processing list
#Set for external command parameters:
#Set "cmd" to "R" or to the full path to Rscript.exe (e.g. C:\Program Files\R\R-4.2.1\bin\Rscript.exe)
#Set "arguments" to the path of the script (e.g. "C:\Users\<username>\Desktop\invert.R")
#Set "Data file" to "Spectrum as JJH5"
#Press "Apply"

library("rhdf5")

dataset = "/JasonDocument/DataPoints/"

args = commandArgs(trailingOnly=TRUE)
filename = args[1]

dsetre = h5read(filename, paste(dataset, "0", sep=""))
dsetim = h5read(filename, paste(dataset, "1", sep=""))

print(paste("dataset size:", length(dsetre)))
h5write(-dsetre, filename, paste(dataset, "0", sep=""))
h5write(-dsetim, filename, paste(dataset, "1", sep=""))

print("dataset changed")
