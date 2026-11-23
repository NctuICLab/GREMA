# Export the SOS data of the Bioconductor package SPEM (data/sos.rda) to a CSV file.
# Usage: Rscript export_sos.R <path to sos.rda> <output csv>
# (the rda holds an ExpressionSet; the S4 object is read through its attributes, so Biobase is not needed)
args <- commandArgs(trailingOnly = TRUE)
e <- new.env(); load(args[1], e)
m <- get("exprs", attr(get("sos", e), "assayData"))
write.csv(m, args[2])
