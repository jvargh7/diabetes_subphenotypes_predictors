# Load required library
library(ggplot2)

# Create data
subtype_data <- data.frame(
  subtype = c("MARD", "MOD", "SIDD", "SIRD"),
  count = c(455, 838, 28, 241),
  percent = c(29.1, 53.6, 1.8, 15.4)
)

# Define colors
cluster_colors <- c("MOD" = "#F8BDA4", "SIRD" = "#A1C3AC", "SIDD" = "#ACD9EA", "MARD" = "#D0ACC9")

# Create label with subtype and percent
subtype_data$label <- paste0(subtype_data$subtype, "\n", subtype_data$percent, "%")

p <- ggplot(subtype_data, aes(x = "", y = count, fill = subtype)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  theme_void() +
  scale_fill_manual(values = cluster_colors) +
  guides(fill = "none")

# 保存高清图像（300 dpi，适合出版）
ggsave(
  filename = "subtype_pie_chart.png",  # 可改为 .tiff 或 .pdf 等格式
  plot = p,
  width = 6,       # 宽度（英寸）
  height = 6,      # 高度（英寸）
  dpi = 300        # 分辨率
)

