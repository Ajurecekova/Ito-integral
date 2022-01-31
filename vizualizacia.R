library(tidyverse)
datatable <- read.table("Trajectories.dsv",sep='|',header=TRUE)

ggplot(data =datatable,aes(x=aes_x,y=value,group=aes_y))+
  geom_line(color="red")


datatable <- read.table("Ito.dsv",sep='|',header=TRUE)

ggplot(data =datatable,aes(x=aes_x,y=value,group=aes_y))+
  geom_line(color = '#E69F00')
