#####
##### This is just for plotting LFDs by date, or month and year
##### Using ggplot (in tidyverse)

library(fishmethods)
library(tidyverse)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

stbr=read.csv("stbrasLF-trunc.csv", header = T) #stbrasLF is Meob trunc is 2004-2009

view(stbr)

stbr %>% ggplot(aes(Length_CM, Frequency))+geom_col()+facet_grid(DATE~.) #Too big?
stbr %>% ggplot(aes(Length_CM, Frequency))+geom_col()+facet_grid(Month~.)

stbr %>% ggplot(aes(Length_CM, Frequency)) + geom_col(color = "black", fill = "light blue") + facet_grid(Year~.) + #other option fill=navy
  scale_x_continuous(name="Fork length (cm)")+
  
 theme(
  panel.background = element_blank(),
  axis.title.y=element_text(color="black", size = 12, face = "bold"),
  axis.text.y=element_text(color="black", size = 10),
  axis.title.x = element_text(size=12, color = "black", face = "bold"),
  axis.text.x = element_text(angle =90 ,
                             hjust = 1,
                             vjust = 0.5,
                             size = 10, 
                             colour = "black"),
  axis.line = element_line(size = 0.5, colour = "black", linetype=1),
  axis.ticks = element_line(size = 0.5, color="black") , 
  axis.ticks.length = unit(0.20, "cm"),
  panel.grid = element_line(color = "grey",
                            size = 0.2,
                            linetype = 1)
 )
