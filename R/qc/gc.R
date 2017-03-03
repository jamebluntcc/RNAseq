#2017-03-01
options(warn = -1)
options(stringsAsFactors = F)
suppressMessages(require('tidyverse',quietly = T))
suppressMessages(require('reshape2',quietly = T))
suppressMessages(require('scales',quietly = T))
suppressMessages(require('argparser',quietly = T))
#----theme set----
theme_Publication <- function(base_size=14, base_family="helvetica") {
  suppressMessages(require('grid',quietly = T))
  suppressMessages(require('ggthemes',quietly = T))
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}
#----plot----
p <- arg_parser('test argparser')
p <- add_argument(p,'--gcfile',help = 'gc stats',default = NULL)
p <- add_argument(p,'--prefix',help = 'path to save plot',default = './')
argv <- parse_args(parser = p)

data <- read.delim(argv$gcfile,header = T)
data[,2:dim(data)[2]] <- data[,2:dim(data)[2]] / 100
data.m <- melt(data,id='X.Base')
gg <- ggplot(data.m,aes(x = X.Base,y = value,colour = variable))+
  geom_line()+geom_vline(xintercept = 150,linetype = 2)+
  scale_x_continuous(breaks = seq(from = 0,to = 300,by = 25),labels = seq(from = 0,to = 300,by = 25))+
  scale_y_continuous(breaks = seq(0,0.7,by = 0.1),labels = percent(seq(0,0.7,by = 0.1)))+
  xlab("Postion")+ylab("Percent(%)")+guides(color = guide_legend(title = ""))+theme_Publication()+
  scale_colour_Publication()
ggsave(paste(argv$prefix,'gc.png',sep = "."),plot = gg,width = 8,height = 6,dpi = 300,type = 'cairo')
ggsave(paste(argv$prefix,'gc.pdf',sep = "."),plot = gg,width = 8,height = 6,device = cairo_pdf)

