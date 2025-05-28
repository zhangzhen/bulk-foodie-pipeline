library(data.table)    # 1.17.0
library(dplyr)    # 1.1.4
library(ggplot2)    # 3.4.4

datasets <- list.dirs(path = '../02-merge-outputs/ratio_TSS', recursive = F)

pdf('conv_ratio_HK_TSS.pdf', width = 7, height = 7)

for (f in datasets) {
    prefix <- unlist(strsplit(f, '/'))[4]
    fname <- paste(f, paste0(prefix, '_TSS_HK.tsv'), sep = '/')

    TSS.HK <- fread(fname, col.names = c('Conv', 'Unconv', 'Rel_pos'))

    TSS.HK.sum <- TSS.HK %>% group_by(Rel_pos) %>% summarise(Conv = sum(Conv), Unconv = sum(Unconv), Ratio = Conv / (Conv + Unconv), Count = n()) %>%
        filter(Count >= 100)

    ratio.max <- quantile((TSS.HK.sum %>% filter(Rel_pos > -150, Rel_pos <= 50))$Ratio, probs = 0.95)
    ratio.min <- quantile((TSS.HK.sum %>% filter(Rel_pos > 50, Rel_pos <= 250))$Ratio, probs = 0.05)

    TSS.HK.sum.plot <- TSS.HK.sum %>% ggplot(aes(x = Rel_pos, y = Ratio)) +
        geom_line(linewidth = 1, color = '#6495ED') +
        coord_cartesian(xlim = c(-500, 500), ylim = c(0, 1)) +
        geom_hline(yintercept = ratio.max, linewidth = 1, color = '#FF7F50') +
        geom_hline(yintercept = ratio.min, linewidth = 1, color = '#FF7F50') +
        labs(x = 'Relative Position to TSS (bp)', y = 'Converted Ratio',
             title = paste0(prefix, ', Max - Min = ', format(ratio.max - ratio.min, digits = 4))) +
        theme_classic() +
        theme(axis.text = element_text(size = 15, family = 'sans', color = 'black', face = 'bold'),
              axis.title = element_text(size = 15, family = 'sans', color = 'black', face = 'bold'),
              plot.title = element_text(size = 15, family = 'sans', color = 'black', face = 'bold', hjust = 0.5))

    print(TSS.HK.sum.plot)
}

dev.off()

