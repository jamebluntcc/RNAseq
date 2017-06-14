library(ggplot2)
library(reshape2)
library(scales)
library(argparser)

options(bitmapType='cairo')

theme_onmath <- function(base_size = 14){
  theme_bw()+
    theme(panel.background = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = 'black'),
          axis.text = element_text(color = 'black',face = 'bold'),
          axis.title = element_text(face = 'bold',size = base_size),
          axis.title.x = element_text(vjust = -0.2),
          axis.title.y = element_text(angle=90,vjust =2),
          plot.title = element_text(face = "bold",size = rel(1.2), hjust = 0.5),
          legend.key = element_blank(),
          legend.title = element_text(face = 'italic'),
          strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
          strip.text = element_text(face="bold")
    )
}

scale_color_onmath <- function(...){
  library(scales)
  discrete_scale("colour","onmath",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}
scale_fill_onmath <- function(...){
  library(scales)
  discrete_scale("colour","onmath",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}

p <- argparser::arg_parser('this is a R script to plot quantification saturation')
p <- add_argument(p,'file_path',help = 'saturation file path',type = "character")
p <- add_argument(p,'out_prefix',help = 'putout saturation plot prefix',type = "character")
argv <- parse_args(p)


saturation_df = read.delim(argv$file_path)

saturation_plot_df = melt(saturation_df)

saturation_plot_df$rpkm <- factor(saturation_plot_df$rpkm, levels = c('0-1', '1-5', '5-10', '10-50', '50-100', '>100'))
pdf(paste(argv$out_prefix,'saturation.plot.pdf', sep = '.'), width = 8, height = 6)
ggplot(saturation_plot_df, aes(x = variable, y = value, group = rpkm, color = rpkm)) + geom_point() + geom_line()+
     xlab('Mapped reads') + ylab('Fraction of transcripts within \n15% of real value') + theme_onmath()+scale_color_onmath() + scale_y_continuous(labels = scales::percent) +
     scale_x_discrete(breaks = paste0('X',seq(from=5,to=95,by=5)),
                      labels = scales::percent(seq(from=0.05,to=0.95,by=0.05)))
dev.off()

png(paste(argv$out_prefix,'saturation.plot.png', sep = '.'), width = 8, height = 6, units = 'in', res = 300)
ggplot(saturation_plot_df, aes(x = variable, y = value, group = rpkm, color = rpkm)) + geom_point() + geom_line()+
  xlab('Mapped reads') + ylab('Fraction of transcripts within \n15% of real value') + theme_onmath()+scale_color_onmath() + scale_y_continuous(labels = scales::percent) +
  scale_x_discrete(breaks = paste0('X',seq(from=5,to=95,by=5)),
                   labels = scales::percent(seq(from=0.05,to=0.95,by=0.05)))
dev.off()

