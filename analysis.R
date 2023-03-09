library(ggplot2)

generations <- 10

ancestry_mut <- read.csv(file = paste0(getwd(), "/output/mut_1.csv"))
ancestry_true <- read.csv(file = paste0(getwd(), "/output/sticklebacks_ancestryproportions_1.csv"))


corresult <- cor.test(ancestry_mut$P1AncestryProportion , ancestry_true$P1AncestryProportion)
plot(ancestry_mut$P1AncestryProportion ~ ancestry_true$P1AncestryProportion, main = paste("P1 proportion ", "Correlation r = " , round(corresult$estimate, 5), " in Gen ", generations))


corresult <- cor.test(ancestry_mut$P2AncestryProportion , ancestry_true$P2AncestryProportion)
plot(ancestry_mut$P2AncestryProportion ~ ancestry_true$P2AncestryProportion, main = paste("P2 proportion ", "Correlation r = " , round(corresult$estimate, 5), " in Gen ", generations))


corresult <- cor.test(ancestry_mut$P3AncestryProportion , ancestry_true$P3AncestryProportion)
plot(ancestry_mut$P3AncestryProportion ~ ancestry_true$P3AncestryProportion, main = paste("P3 proportion ", "Correlation r = " , round(corresult$estimate, 5), " in Gen ", generations))


corresult <- cor.test(ancestry_mut$P4AncestryProportion , ancestry_true$P4AncestryProportion)
plot(ancestry_mut$P4AncestryProportion ~ ancestry_true$P4AncestryProportion, main = paste("P4 proportion ", "Correlation r = " , round(corresult$estimate, 5), " in Gen ", generations))

Freq_dropout <-  rowSums(1*as.matrix((ancestry_true > 0) & (ancestry_mut == 0)))
sum(Freq_dropout>0) / length(Freq_dropout)

Num_categories <- 2^generations
fractions <- 1/Num_categories
Actual_ancestry_bin <- round(ancestry_true*Num_categories)/Num_categories
Inferred_ancestry_bin <- round(ancestry_mut*Num_categories)/Num_categories
#table(Inferred_ancestry_bin$P1AncestryProportion, Actual_ancestry_bin$P1AncestryProportion)
Nwrong_bins <- rowSums(1*as.matrix(Actual_ancestry_bin != Inferred_ancestry_bin))
sum(Nwrong_bins > 0)
miscat <- sum(Nwrong_bins > 0)
miscat / 1000