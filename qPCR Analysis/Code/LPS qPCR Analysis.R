############################################################
### LPS.PBS E3.E4 QPCR ANALYSIS
### Ayush Noori
############################################################

# load packages
library(data.table)
library(stringr)
library(tidyverse)

# read and write to Excel files
library(openxlsx)

# statistical analysis
library(psych)

# plot data
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

# caption text and display significance
# devtools::install_github("wilkelab/ggtext")
library(ggtext)
library(ggsignif)

# heatmap
library(pheatmap)


###################################
### SETUP
###################################

# set working directory
homedir = "" # change to your directory
lps_dir = file.path(homedir, "Results", "LPS Plots")
setwd(homedir)

lps = read.xlsx("Data/qPCR Data.xlsx", sheet = 1)
lps_RQ = read.xlsx("Data/qPCR Data.xlsx", sheet = 3)
if(file.exists("LPS ANOVA Results.txt")) file.remove("LPS ANOVA Results.txt")


###################################
### FILTER DATA
###################################

# filter LPS
rownames(lps) = lps$X1; lps = lps[, -1]
colnames(lps)[1:3] = c("Genotype", "Treatment", "Sex")

# append RQ data
colnames(lps_RQ) = paste0(colnames(lps_RQ), ".RQ")
lps = cbind(lps, lps_RQ[, 5:19])

# convert to factor
lps$Genotype = factor(lps$Genotype, levels = c(3, 4))
lps$Treatment = factor(lps$Treatment, levels = c(2, 1))
lps$Sex = factor(lps$Sex, levels = c(1, 2))

# remap factors
lps$Genotype = plyr::revalue(lps$Genotype, replace = c(`3` = "APOE3", `4` = "APOE4"))
lps$Treatment = plyr::revalue(lps$Treatment, replace = c(`2` = "PBS", `1` = "LPS"))
lps$Sex = plyr::revalue(lps$Sex, replace = c(`1` = "Male", `2` = "Female"))

# get gene indices
idx = c(4:7, 9:12, 14:17)
genes = colnames(lps)[idx]
idx_RQ = idx+15
genes_RQ = colnames(lps)[idx_RQ]

# initialize plot list
plots = vector(mode = "list", length = length(genes))


for (i in 1:length(idx)) {
  
  ###################################
  ### STATISTICS
  ###################################
  
  # select gene for this iteration
  gene = genes[i]
  gene_RQ = genes_RQ[i]
  cat(paste0("Analysis #", i, ": ", gene, "\n"))
  
  # perform two-way ANOVA
  # use the * rather than + to add interaction term
  res_aov = aov(get(gene) ~ Genotype * Treatment, data = lps) # excluding Sex
  stat = summary(res_aov)[[1]]
  
  # print to output
  sink(file = "LPS ANOVA Results.txt", append = TRUE)
  if(i != 1) cat("\n\n")
  cat(paste0(gene, " ANOVA:\n"))
  print(as.matrix(stat))
  sink()
  
  # calculate summary statistics
  dat = describeBy(lps[[gene_RQ]], list(lps$Genotype, lps$Treatment), mat=TRUE, digits=2)
  names(dat)[names(dat) == "group1"] = "Genotype"
  names(dat)[names(dat) == "group2"] = "Treatment"
  dat$Genotype = factor(dat$Genotype, levels = c("APOE3", "APOE4"))
  dat$Treatment = factor(dat$Treatment, levels = c("PBS", "LPS"))
  
  # make sure that error bars don't become negative
  for (r in 1:nrow(dat)) {
    if(dat$mean[r] - dat$sd[r] < 0) dat$sd[r] = dat$mean[r]
  }
  
  # define error bar limits
  limits = aes(ymax = mean + sd, ymin = mean - sd)
  
  
  ###################################
  ### PLOT DATA
  ###################################
  
  stat_plot = data.frame(stat$`Pr(>F)`[1:3])
  colnames(stat_plot) = c("PValue")
  stat_plot[["Label"]] = c("Genotype", "Treatment", "Interaction")
  stat_plot[["Significance"]] = stat_plot$PValue < 0.05
  
  stat_plot$PValue = round(stat_plot$PValue, 4)
  
  # plot caption
  caption_labs = stat_plot$PValue
  
  # conditional formating for caption label
  for (k in 1:length(caption_labs)){
    if(caption_labs[k] == 0) {
      caption_labs[k]  = "<span style = 'color:#C33C54;'>*p* < 0.0001</span>"
    } else if (caption_labs[k] < 0.05) {
      caption_labs[k]  = paste0("<span style = 'color:#C33C54;'>", "*p* = ", caption_labs[k], "</span>")
    } else if (caption_labs[k] < 0.1) {
      caption_labs[k]  = paste0("<span style = 'color:#D4AFB9;'>", "*p* = ", caption_labs[k], "</span>")
    } else if(caption_labs[k] >= 0.1) {
      caption_labs[k]  = "<span style = 'color:#41414E;'>N.S.</span>"
    }
  }
  
  # summary caption label
  plot_caption = paste0("<span style = 'color:#41414E;'>**Genotype:** </span>", caption_labs[1], "<br>",
                        "<span style = 'color:#41414E;'>**Treatment:** </span>", caption_labs[2], "<br>",
                        "<span style = 'color:#41414E;'>**Interaction:** </span>", caption_labs[3], "<br>")
 
  # make barplot
  p = ggplot(dat, aes(x = Genotype, y = mean, color = Genotype, group = Treatment)) +
    geom_bar(stat='identity', fill = "white", size = 1, width = 0.5) +
    scale_color_manual(values = c("#785D37", "#62BBA5")) +
    geom_errorbar(limits, color = "#41414E", width=0.25) +
    geom_jitter(data = lps, aes(x = Genotype, y = .data[[gene_RQ]], color = Genotype, shape = Sex),
                size = 2, width=0.15) +
    facet_wrap(~ Treatment) +
    ggtitle(gene) +
    ylab("Expression Fold-Change") +
    xlab("") +
    scale_y_continuous(expand = expansion(mult = c(-0.01, .1))) +
    labs(caption = plot_caption) +
    theme_linedraw() +
    theme(plot.title = element_text(size = 20, hjust = 0.5, color = "#41414E", face = "bold.italic"),
          plot.caption = element_markdown(hjust = 0.5, size = 16),
          axis.text.x.bottom = element_text(size = 12, face = "bold.italic", color = "#41414E"),
          axis.title.y = element_text(size = 14, face = "bold", color = "#41414E"),
          strip.background = element_rect(fill = "#5D5D6F"),
          strip.text = element_text(size=14, color = "white", face = "bold"),
          panel.background = element_rect(fill = "#FDF5ED"),
          axis.text.x = element_text(size=12, color = "#41414E", face = "italic"),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank()
    )
  
  # save plot
  ggsave(file.path(lps_dir, paste0(gene, " Plot.pdf")), p, width = 6, height = 6)
  
  # add plot to total list
  plots[[i]] = p
  
}

p_legend = ggplot(lps, aes(x = Genotype, y = gene_RQ, color = Genotype, shape = Sex)) +
  scale_color_manual(values = c("#785D37", "#62BBA5")) +
  geom_jitter(size = 2, width=0.1) + 
  theme(legend.title = element_text(size = 18, face = "bold", color = "#41414E"),
        legend.text = element_text(size = 14, color = "#41414E"),
        legend.position = "top")

my_legend = get_legend(p_legend)

# combine and arrange all plots using NULL values as placeholders to arrange plots
all_plots = ggarrange(NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                      plots[[7]], NULL,
                      plots[[8]] + rremove("ylab"),  NULL,
                      plots[[5]] + rremove("ylab"),  NULL,
                      plots[[6]] + rremove("ylab"),
                      plots[[2]],  NULL,
                      plots[[1]] + rremove("ylab"),  NULL,
                      plots[[4]] + rremove("ylab"),  NULL,
                      plots[[3]] + rremove("ylab"),
                      plots[[9]],  NULL,
                      plots[[12]] + rremove("ylab"),  NULL,
                      plots[[10]] + rremove("ylab"),  NULL,
                      plots[[11]] + rremove("ylab"),
                      ncol = 7,
                      nrow = 4,
                      widths = c(1, 0.1, 1, 0.1, 1, 0.1, 1),
                      heights = c(0.05, 1, 1, 1),
                      legend.grob = my_legend,
                      legend = "bottom")

ggsave("LPS-PBS e3-e4 qPCR Figure.pdf", all_plots, width = 16, height = 18)
# ggsave("LPS-PBS e3-e4 qPCR Figure.png", all_plots, width = 16, height = 18, units = "in", dpi = 900)



###################################
### CORRELATION HEATMAPS
###################################

cor_dat = lps[, 4:18]
print(sapply(cor_dat, function(x) agostino.test(as.numeric(x), alternative = "two.sided")$p.value))

lps_cor = cor(cor_dat, method = "spearman")
lps_cor = lps_cor[-c(5, 10, 15), -c(5, 10, 15)]

newnames = lapply(rownames(lps_cor),  function(x) bquote(italic(.(x))))

hm = pheatmap(lps_cor,
              color = colorRampPalette(c("#313695", "#8FC3DC", "#E3E0DD", "#F88D52", "#A50026"))(100),
              breaks = seq(from = -1, to = 1, length.out = 101),
              cellwidth = 25, cellheight = 25,
              cluster_cols = TRUE, cluster_rows = TRUE,
              border_color = NA,
              show_colnames = TRUE,
              show_rownames = TRUE,
              silent = TRUE,
              labels_row = as.expression(newnames), labels_col = as.expression(newnames))

ggsave("LPS-PBS e3-e4 qPCR Correlation Heatmap.pdf", hm, width = 8, height = 7.5)
# ggsave("LPS-PBS e3-e4 qPCR Correlation Heatmap.png", hm, width = 8, height = 7.5, units = "in", dpi = 900)



###################################
### CORRELATION MATRIX
###################################

wb = createWorkbook()

hs = createStyle(fontColour = "#FAF3DD", fgFill = "#5171A5", fontName = "Arial Black",
                 halign = "left", valign = "center", textDecoration = "Bold",
                 border = "Bottom", borderStyle = "thick", fontSize = 10)

addWorksheet(wb, sheetName = "LPS.PBS-E2.3.4")
writeDataTable(wb, "LPS.PBS-E2.3.4", x = data.frame(lps_cor), tableStyle = "TableStyleMedium15", headerStyle = hs, rowNames = TRUE)
addStyle(wb, "LPS.PBS-E2.3.4", createStyle(fontColour = "#363635", fgFill = "#FAF3DD", fontName = "Arial", fontSize = 10), rows = 2:16, cols = 2:16, gridExpand = T)
addStyle(wb, "LPS.PBS-E2.3.4", hs, rows = 2:16, cols = 1, gridExpand = T)

saveWorkbook(wb, "LPS-PBS e3-e4 Correlation Matrix.xlsx", overwrite = TRUE)