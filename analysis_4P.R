
file_prefix <- "sticklebacks_ancestryproportions"
file_suffix <- ".csv"

# create an empty list to store the data frames
ancestry_true_list <- list()

# loop through each file and read it into R
for (i in 1:iteration) {
  filename <- paste0(file_prefix, "_", i, file_suffix)
  filepath <- paste0(getwd(), "/output/", filename)
  ancestry_true_list[[i]] <- read.csv(filepath)
}

# combine the data frames into one
ancestry_true <- do.call(rbind, ancestry_true_list)
ancestry_true <- ancestry_true[, -1]
colnames(ancestry_true) <- c("FG", "TL", "WB", "WT")
col_means <- colMeans(ancestry_true)

col_means
