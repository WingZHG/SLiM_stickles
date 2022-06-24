library(ggplot2)

generations <- 5

ancestry_mut <- read.csv(file = paste0(getwd(), "/output/mut.csv"))
ancestry_true <- read.csv(file = paste0(getwd(), "/output/sticklebacks_ancestryproportions.csv"))


corresult <- cor.test(ancestry_mut$P1AncestryProportion , ancestry_true$P1AncestryProportion)
plot(ancestry_mut$P1AncestryProportion ~ ancestry_true$P1AncestryProportion, main = paste("P1 proportion ", "Correlation r = " , round(corresult$estimate, 5), " in Gen ", generations))


corresult <- cor.test(ancestry_mut$P2AncestryProportion , ancestry_true$P2AncestryProportion)
plot(ancestry_mut$P2AncestryProportion ~ ancestry_true$P2AncestryProportion, main = paste("P2 proportion ", "Correlation r = " , round(corresult$estimate, 5), " in Gen ", generations))


corresult <- cor.test(ancestry_mut$P3AncestryProportion , ancestry_true$P3AncestryProportion)
plot(ancestry_mut$P3AncestryProportion ~ ancestry_true$P3AncestryProportion, main = paste("P3 proportion ", "Correlation r = " , round(corresult$estimate, 5), " in Gen ", generations))


corresult <- cor.test(ancestry_mut$P4AncestryProportion , ancestry_true$P4AncestryProportion)
plot(ancestry_mut$P4AncestryProportion ~ ancestry_true$P4AncestryProportion, main = paste("P4 proportion ", "Correlation r = " , round(corresult$estimate, 5), " in Gen ", generations))
