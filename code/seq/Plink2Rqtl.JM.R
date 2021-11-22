#!/bin/R
args = commandArgs(trailingOnly=TRUE)

library(qtl)
source("code/seq/PLINK2RQTL.f2.R")

dir <- args[1]
name <- args[2]
tag <- args[3]

setwd(dir)

ped <- paste(dir, name, ".ped", sep = "")
map <- paste(dir, name, ".map", sep = "")
PLINKtoCSVR(ped = ped, map = map, out = paste(dir, name, tag , sep = ""))

