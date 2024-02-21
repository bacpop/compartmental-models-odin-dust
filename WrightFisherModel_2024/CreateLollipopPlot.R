library(ggplot2)
library(gridExtra)
library(grid)

lollipop_cluster_freqs <- function(plot_title = "Generic Plot Title",data1, data2, data3, model1, model2, model3){
  lollipop_data_1 <- data.frame(
    x=1:length(data1),
    model_1=model1,
    data_1=as.numeric(data1)
  )
  lollipop_data_2 <- data.frame(
    x=1:length(data2),
    model_2=model2,
    data_2=as.numeric(data2)
  )
  lollipop_data_3 <- data.frame(
    x=1:length(data3),
    model_3=model3,
    data_3=as.numeric(data3)
  )
  
  # Change baseline
  lollipop_plot_1 <- ggplot(lollipop_data_1) +
    geom_segment( aes(x=x, xend=x, y=model_1, yend=data_1), color="grey") +
    geom_point( aes(x=x, y=model_1, color="Model non-VT"), size=3 ) +
    geom_point( aes(x=x, y=data_1, color="Data non_VT"), size=3 ) +
    scale_color_manual(values = c("#E69F00", "#56B4E9"),
                       guide  = guide_legend(), 
                       name   = "Group") +
    coord_flip()+
    #theme_ipsum() +
    theme(legend.position = "none") +
    ggtitle("2001") +
    ylab("Frequency") +
    xlab("Clusters") +
    theme(axis.title  = element_text(size = 20), axis.text = element_text(size = 20), plot.title = element_text(size = 25,hjust = 0.5))  +
    ylim(0, max(max(lollipop_data_1$model_1),max(lollipop_data_1$data_1),max(lollipop_data_2$model_2),max(lollipop_data_2$data_2),max(lollipop_data_3$model_3),max(lollipop_data_3$data_3)))
  lollipop_plot_2 <- ggplot(lollipop_data_2) +
    geom_segment( aes(x=x, xend=x, y=model_2, yend=data_2), color="grey") +
    geom_point( aes(x=x, y=model_2, color="Model"), size=3 ) +
    geom_point( aes(x=x, y=data_2, color="Data"), size=3 ) +
    scale_color_manual(values = c("#E69F00", "#56B4E9"),
                       guide  = guide_legend(), 
                       name   = "Group") +
    coord_flip()+
    #theme_ipsum() +
    theme(legend.position = "none") +
    ggtitle("2004") +
    ylab("Frequency") +
    xlab("Clusters") +
    theme(axis.title  = element_text(size = 20), axis.text = element_text(size = 20), plot.title = element_text(size = 25,hjust = 0.5))  +
    ylim(0, max(max(lollipop_data_1$model_1),max(lollipop_data_1$data_1),max(lollipop_data_2$model_2),max(lollipop_data_2$data_2),max(lollipop_data_3$model_3),max(lollipop_data_3$data_3)))
  lollipop_plot_3 <- ggplot(lollipop_data_3) +
    geom_segment( aes(x=x, xend=x, y=model_3, yend=data_3), color="grey") +
    geom_point( aes(x=x, y=model_3, color="Model"), size=3 ) +
    geom_point( aes(x=x, y=data_3, color="Data"), size=3 ) +
    scale_color_manual(values = c("#E69F00", "#56B4E9"),
                       guide  = guide_legend(), 
                       name   = "Legend") +
    coord_flip()+
    #theme_ipsum() +
    theme(legend.position = c(.8,.8),legend.text = element_text(size = 20),legend.title = element_text(size = 20)) +
    ggtitle("2007") +
    ylab("Frequency") +
    xlab("Clusters") +
    theme(axis.title  = element_text(size = 20), axis.text = element_text(size = 20), plot.title = element_text(size = 25,hjust = 0.5))  +
    ylim(0, max(max(lollipop_data_1$model_1),max(lollipop_data_1$data_1),max(lollipop_data_2$model_2),max(lollipop_data_2$data_2),max(lollipop_data_3$model_3),max(lollipop_data_3$data_3)))
  
  grid.arrange(lollipop_plot_1 + scale_y_continuous(limits = c(NA,0.17)) + theme(plot.margin = unit(c(.5,0.5,1,0.5), "cm"),axis.text.y = element_blank()) ,lollipop_plot_2 + scale_y_continuous(limits = c(NA,0.17)) +  theme(plot.margin = unit(c(.5,0.5,1,0.5), "cm"),axis.text.y = element_blank()), lollipop_plot_3 + scale_y_continuous(limits = c(NA,0.17)) + theme(plot.margin = unit(c(.5,0.5,1,0.5), "cm"),axis.text.y = element_blank()), ncol = 3, nrow=1, top = textGrob(plot_title,gp=gpar(fontsize=20,font=3)))
}