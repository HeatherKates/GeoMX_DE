# Load necessary libraries
library(ggplot2)
library(viridis)
library(dplyr)
library(gridExtra)

INS_corrected_ <- data.frame(cbind(colnames(spe@assays@data$corrected_counts),spe@assays@data$corrected_counts["INS",]))
colnames(INS_corrected_counts) <- c("AOI","INS_corrected_counts")
INS_corrected_counts <- data.frame(cbind(colnames(corrected_counts_subset),corrected_counts_subset["INS",]))
colnames(INS_corrected_counts) <- c("AOI","INS_corrected_counts")
AOI_target <- data.frame(cbind(spe@colData$sample_id,as.character(spe@colData$AOI_target)))
INS_AOI_target <- merge(INS_corrected_counts,AOI_target,by.x="AOI",by.y="X1")
# Plot
INS_AOI_target$INS_corrected_counts <- as.numeric(INS_AOI_target$INS_corrected_counts)
INS_AOI_target %>%
  ggplot( aes(x=X2, y=INS_corrected_counts, fill=X2)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("A boxplot with jitter") +
  xlab("")
