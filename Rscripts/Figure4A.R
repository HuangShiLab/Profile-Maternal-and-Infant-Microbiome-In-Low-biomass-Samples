#######################################
# Monther-infant Study
#
# Figure 4A
#
# Author: JIANG Yuesong
#######################################

library(ggplot2)
library(dplyr)
# library(coin)

# data import
setwd("~/Downloads/data")

data <- read.table("alpha-diversity.tsv", header = TRUE, sep = "\t")

# # 执行 Wilcoxon 检验
wilcox_result <- wilcox_test(shannon_entropy ~ time, data = data) %>%
  add_significance()  %>% add_xy_position(x = "time")

# 绘制箱线图，设置填充轮廓线的颜色和粗细
boxplot <- ggplot(data, aes(x = time, y = shannon_entropy)) +
  geom_boxplot(aes(fill = time),
               color = "black",  # 设置箱线图的轮廓颜色为黑色
               fill = c("#87CEEB","#7CCD7C")) + 
  geom_beeswarm(alpha = 0.6, color = "black") +
  labs(x = "", y = "Shannon Diversity") +  # 设置标签
  stat_pvalue_manual(wilcox_result, label = "p.signif", tip.length = 0.02) +  
  theme_minimal() +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.3, color = "black"))

# 自动换行
boxplot <- boxplot + scale_x_discrete(labels = function(x) sapply(strwrap(x, width = 15, simplify = FALSE), paste, collapse = "\n"))

# 打印箱线图
print(boxplot)


# 保存图像为PDF格式，指定尺寸
ggsave("./plots/4A.png", boxplot, width = 2.5, height = 3)

