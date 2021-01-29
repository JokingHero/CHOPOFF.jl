rm(list = ls(all.names = TRUE))
gc(reset = TRUE)

library(data.table)
library(Biostrings)
library(ggplot2)

# CRISPRofftargetHunter
hunter <- fread("~/.julia/dev/CRISPRofftargetHunter/test/sample_data/guides_not_norm.csv")

# crispritz
cz <- fread("soft/crispritz/guides_cr.output.targets.txt")
cz$guide <- substr(gsub("-", "", cz$crRNA), start = 1, stop = 20)
cz <- cz[Total <= 3,
         .(czD0 = sum(Total == 0),
           czD1 = sum(Total == 1),
           czD2 = sum(Total == 2),
           czD3 = sum(Total == 3)), by = "guide"]
cz <- cz[order(czD0, czD1, czD2, czD3)]

ma <- function(x, n = 10){filter(x, rep(1 / n, n), sides = 2)}
ranks <- match(hunter$guide, cz$guide)

ggplot(data = data.frame(x = 1:length(ranks), y = ranks), aes(x = x, y = y)) +
  geom_line()
