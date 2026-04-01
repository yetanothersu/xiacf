# R/globals.R

# This file is for declaring global variables to avoid "no visible binding for global variable" notes during R CMD check.
utils::globalVariables(c(
    "Lag",
    "Xi_Threshold_95",
    "ACF",
    "CCF",
    "Xi",
    "i",
    "y_vals",
    "Direction"
))
