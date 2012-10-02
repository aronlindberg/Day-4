# Load necessary libraries & datasets
library(TraMineR)
data(biofam)

# Parse and label the data
biofam$cohort <- cut(biofam$birthyr, c(1900,1930,1940,1950,1960), labels=c("1900-1929", "1930-1939", "1940-1949", "1950-1959"), right=FALSE)

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
bf.long <- data.frame(bf.lgth, bf.tran, bf.sseq, bf.entr, bf.turb, bf.ici)

# Renaming the columns
names(bf.long) <- c("lgth", "tran", "sseq", "entr", "turb", "ici")

# Create a histogram
par(mfrow = c(2, 3))
col=c("red", "green3", "blue", "cyan", "magenta", "yellow", "grey")
for (i in 2:ncol(bf.long)){
  hist(bf.long[,i], col=col[i], main=names(bf.long)[i], xlab="")
}

# Extract distinct subsequences and durations
bf.dss <-seqdss(biofam.seq)
tail(bf.dss)
bf.dur <- seqdur(biofam.seq)
tail(bf.dur)

# Mean and variance of time in distinct successive states
meant <- apply(bf.dur, 1, mean, na.rm=TRUE)
summary(meant)

vart <- apply(bf.dur, 1, var, na.rm=TRUE)
summary(vart)

# Create scatterplots-matrices for entropy, turbulence, and complexity
plot(bf.long[,c("entr", "turb","ici")])

# Boxplots of complexity by birth cohorts
boxplot(bf.long$ici ~ biofam$cohort, col="magenta")

# Regress complexity on birth cohort, sex, and language
lm.ici <- lm(bf.long$ici ~ cohort + sex + plingu02, data=biofam)
summary(lm.ici)