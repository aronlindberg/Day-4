# Load necessary libraries & datasets
library(TraMineR)
data(biofam)

# Parse and label the data
biofam$chort <- cut(biofam$birthyr, c(1900,1930,1940,1950,1960), labels=c("1900-1929", "1930-1939", "1940-1949", "1950-1959"), right=FALSE)

bf.states <- c("Parent", "Left", "Married", "Left/Married", "Child", "Left/Child", "Left/Married/Child", "Divorced")

bf.shortlab <- c("P", "L", "M", "LM", "C", "LC", "LMC", "D")

# Create a sequence object
biofam.seq <- seqdef(biofam[,10:25], states=bf.shortlab, labels=bf.states, weights=biofam$wp00tbgs)

# Computing longitudinal characteristics
bf.lgth <- seqlength(biofam.seq)
bf.tran <- seqtransn(biofam.seq)
bf.sseq <- seqsubsn(biofam.seq)
bf.entr <- seqient(biofam.seq)
bf.turb <- seqST(biofam.seq)
bf.ici <- seqici(biofam.seq)

# Merge all the characteristics into a single table
