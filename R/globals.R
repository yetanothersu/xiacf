# R/globals.R

# 検査ツール（R CMD check）を黙らせるためのグローバル変数宣言
utils::globalVariables(c(
    "Lag",
    "Xi_Threshold_95",
    "ACF",
    "Xi",
    "i"
))
