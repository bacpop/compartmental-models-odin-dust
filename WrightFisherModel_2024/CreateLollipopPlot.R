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
  
  grid.arrange(lollipop_plot_1 + scale_y_continuous(limits = c(NA,0.2)) + theme(plot.margin = unit(c(.5,0.5,1,0.5), "cm"),axis.text.y = element_blank()) ,lollipop_plot_2 + scale_y_continuous(limits = c(NA,0.2)) +  theme(plot.margin = unit(c(.5,0.5,1,0.5), "cm"),axis.text.y = element_blank()), lollipop_plot_3 + scale_y_continuous(limits = c(NA,0.2)) + theme(plot.margin = unit(c(.5,0.5,1,0.5), "cm"),axis.text.y = element_blank()), ncol = 3, nrow=1, top = textGrob(plot_title,gp=gpar(fontsize=20,font=3)))
}

### Creating a lollipop plot with three points per line

lollipop_cluster_freqs_3points <- function(year = "year unknown", plot_title = "Generic Plot Title",data1, model_name_1 ="Model 1", model1, model_name_2 ="Model 2", model2){
  lollipop_data_1 <- data.frame(
    x=1:length(data1),
    model_1=model1,
    data_1=as.numeric(data1),
    model_2 = model2
  )
  # Change baseline
  lollipop_plot_1 <- ggplot(lollipop_data_1) +
    geom_segment( aes(x=x, xend=x, y=model_1, yend=data_1), color="grey") +
    geom_segment( aes(x=x, xend=x, y=model_2, yend=data_1), color="grey") +           
    geom_point( aes(x=x, y=model_1, color=model_name_1), size=3, shape = 1, stroke = 2, alpha = 0.7) +
    geom_point( aes(x=x, y=data_1, color="Data"), size=3, shape = 1, stroke = 2, alpha = 0.7) +
    geom_point(aes(x=x, y=model_2, color=model_name_2), size=3, shape = 1, stroke = 2, alpha = 0.7) +
    #geom_point( aes(x=x, y=model_1, color=model_name_1), size=5, alpha = 0.7) +
    #geom_point( aes(x=x, y=data_1, color="Data"), size=5, alpha = 0.7) +
    #geom_point(aes(x=x, y=model_2, color=model_name_2), size=5, alpha = 0.7) +
    scale_color_manual(values = c("#E69F00","#56B4E9","#CC79A7"),
                       guide  = guide_legend(), 
                       name   = "Group") +
    coord_flip()+
    #theme_ipsum() +
    theme(legend.position = c(.8,.8),legend.text = element_text(size = 20),legend.title = element_text(size = 20)) +
    ggtitle(year) +
    ylab("Frequency") +
    xlab("Clusters") +
    theme(axis.title  = element_text(size = 20), axis.text = element_text(size = 20), plot.title = element_text(size = 25,hjust = 0.5))  +
    ylim(0, max(max(lollipop_data_1$model_1)))
  grid.arrange(lollipop_plot_1 + scale_y_continuous(limits = c(NA,0.2)) + theme(plot.margin = unit(c(.5,0.5,1,0.5), "cm"),axis.text.y = element_blank()), ncol = 1, nrow=1, top = textGrob(plot_title,gp=gpar(fontsize=20,font=3)))
}

lollipop_cluster_freqs_2x2points <- function(year = "year unknown", plot_title = "Generic Plot Title",data1, model_name_1 ="Model 1", model1, model_name_2 ="Model 2", model2){
  lollipop_data_1 <- data.frame(
    x=1:(length(data1)),
    model_1=model1,
    data_1=as.numeric(data1),
    data_2=as.numeric(data1),
    model_2 =model2
  )
  # Change baseline
  lollipop_plot_1 <- ggplot(lollipop_data_1) +
    geom_segment( aes(x=x, xend=x, y=model_1, yend=data_1), color="grey") +
    geom_segment( aes(x=x-0.3, xend=x-0.3, y=model_2, yend=data_2), color="grey") +
    geom_point( aes(x=x, y=data_1, color="Data"), size=5) +
    geom_point( aes(x=x-0.3, y=data_2, color="Data"), size=5) +
    geom_point( aes(x=x, y=model_1, color=model_name_1), size=5) +
    geom_point( aes(x=x-0.3, y=model_2, color=model_name_2), size=5) +
    #geom_point( aes(x=x, y=model_1, color=model_name_1), size=5, alpha = 0.7) +
    #geom_point( aes(x=x, y=data_1, color="Data"), size=5, alpha = 0.7) +
    #geom_point(aes(x=x, y=model_2, color=model_name_2), size=5, alpha = 0.7) +
    scale_color_manual(values = c("#E69F00","#56B4E9","#CC79A7"),
                       guide  = guide_legend(), 
                       name   = "Group") +
    coord_flip()+
    #theme_ipsum() +
    theme(legend.position = c(.8,.8),legend.text = element_text(size = 20),legend.title = element_text(size = 20)) +
    ggtitle(year) +
    ylab("Frequency") +
    xlab("Clusters") +
    theme(axis.title  = element_text(size = 20), axis.text = element_text(size = 20), plot.title = element_text(size = 25,hjust = 0.5))  +
    ylim(0, max(max(lollipop_data_1$model_1)))
  grid.arrange(lollipop_plot_1 + scale_y_continuous(limits = c(NA,0.2)) + theme(plot.margin = unit(c(.5,0.5,1,0.5), "cm"),axis.text.y = element_blank()), ncol = 1, nrow=1, top = textGrob(plot_title,gp=gpar(fontsize=20,font=3)))
}


# Creating a joint Lollipop plot for VTs and NVTs
lollipop_cluster_freqs_VTandNVT <- function(year = "year unknown", plot_title = "Generic Plot Title",data1, model_name_1 ="Model 1", model1, VT_vec){
  lollipop_data_1 <- data.frame(
    x=1:length(data1),
    model_1=model1,
    data_1=as.numeric(data1),
    Vaccine_Types <- sapply(VT_vec, function(x) if(x==1){"VT"}else{"NVT"})
  )
  # Change baseline
  lollipop_plot_1 <- ggplot(lollipop_data_1) +
    geom_segment( aes(x=x, xend=x, y=model_1, yend=data_1), color="grey") +
    geom_point( aes(x=x, y=model_1, color=model_name_1, shape = Vaccine_Types), size=3, stroke = 2) +
    geom_point( aes(x=x, y=data_1, color="Data",shape = Vaccine_Types), size=3,  stroke = 2) +
    #geom_point( aes(x=x, y=model_1, color=model_name_1), size=5, alpha = 0.7) +
    #geom_point( aes(x=x, y=data_1, color="Data"), size=5, alpha = 0.7) +
    #geom_point(aes(x=x, y=model_2, color=model_name_2), size=5, alpha = 0.7) +
    scale_color_manual(values = c("#E69F00","#56B4E9","#CC79A7"),
                       guide  = guide_legend(), 
                       name   = "Group") +
    scale_shape_manual(values=c(19, 1))+
    coord_flip()+
    #theme_ipsum() +
    theme(legend.position = c(.8,.8),legend.text = element_text(size = 20),legend.title = element_text(size = 20)) +
    ggtitle(year) +
    ylab("Frequency") +
    xlab("Clusters") +
    theme(axis.title  = element_text(size = 20), axis.text = element_text(size = 15), plot.title = element_text(size = 25,hjust = 0.5))  +
    ylim(0, max(max(lollipop_data_1$model_1)))
  grid.arrange(lollipop_plot_1 + scale_y_continuous(limits = c(NA,0.2)) + theme(plot.margin = unit(c(.5,0.5,1,0.5), "cm"),axis.text.y = element_blank()), ncol = 1, nrow=1, top = textGrob(plot_title,gp=gpar(fontsize=20,font=3)))
}

# Creating a joint Lollipop plot for VTs and NVTs
lollipop_cluster_freqs_VTandNVT_labelSero <- function(year = "year unknown", plot_title = "Generic Plot Title",data1, model_name_1 ="Model 1", model1, VT_vec, SeroLabel){
  lollipop_data_1 <- data.frame(
    x=1:length(data1),
    model_1=model1,
    data_1=as.numeric(data1),
    Vaccine_Types <- sapply(VT_vec, function(x) if(x==1){"VT"}else{"NVT"})
  )
  # Change baseline
  lollipop_plot_1 <- ggplot(lollipop_data_1) +
    geom_segment( aes(x=x, xend=x, y=model_1, yend=data_1), color="grey") +
    geom_vline(xintercept = 10, color = "white", linewidth = 1.5) +
    geom_vline(xintercept = 20, color = "white", linewidth = 1.5) +
    geom_vline(xintercept = 30, color = "white", linewidth = 1.5) +
    geom_vline(xintercept = 40, color = "white", linewidth = 1.5) +
    geom_vline(xintercept = 50, color = "white", linewidth = 1.5) +
    geom_vline(xintercept = 60, color = "white", linewidth = 1.5) +
    geom_point( aes(x=x, y=model_1, color=model_name_1, shape = Vaccine_Types), size=3, stroke = 2) +
    geom_point( aes(x=x, y=data_1, color="Data",shape = Vaccine_Types), size=3,  stroke = 2) +
    #geom_point( aes(x=x, y=model_1, color=model_name_1), size=5, alpha = 0.7) +
    #geom_point( aes(x=x, y=data_1, color="Data"), size=5, alpha = 0.7) +
    #geom_point(aes(x=x, y=model_2, color=model_name_2), size=5, alpha = 0.7) +
    scale_color_manual(values = c("#E69F00","#56B4E9","#CC79A7"),
                       guide  = guide_legend(), 
                       name   = "Group") +
    scale_shape_manual(values=c(19, 1))+
    coord_flip()+
    #theme_ipsum() +
    theme(legend.position = c(.8,.8),legend.text = element_text(size = 20),legend.title = element_text(size = 20)) +
    ggtitle(year) +
    ylab("Frequency") +
    xlab("Clusters") +
    theme(axis.title  = element_text(size = 20), axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 20), plot.title = element_text(size = 25,hjust = 0.5))  +
    ylim(0, max(max(lollipop_data_1$model_1)))
  grid.arrange(lollipop_plot_1 + scale_y_continuous(limits = c(NA,0.2))  + scale_x_continuous(breaks = 1:length(data1), labels = 1:length(data1), sec.axis = dup_axis(name = "Serotypes", labels = SeroLabel))+ theme(plot.margin = unit(c(.5,0.5,1,0.5), "cm")), ncol = 1, nrow=1, top = textGrob(plot_title,gp=gpar(fontsize=20,font=3)))
}


# Creating a joint Lollipop plot for VTs and NVTs and two models
lollipop_cluster_freqs_2x2_VTandNVT <- function(year = "year unknown", plot_title = "Generic Plot Title",data1, model_name_1 ="Model 1", model1, model_name_2 ="Model 2", model2, VT_vec){
  lollipop_data_1 <- data.frame(
    x=1:length(data1),
    model_1=model1,
    data_1=as.numeric(data1),
    Vaccine_Types <- sapply(VT_vec, function(x) if(x==1){"VT"}else{"NVT"}),
    data_2=as.numeric(data1),
    model_2 =model2
  )
  # Change baseline
  lollipop_plot_1 <- ggplot(lollipop_data_1) +
    geom_segment( aes(x=x, xend=x, y=model_1, yend=data_1), color="grey") +
    geom_segment( aes(x=x-0.3, xend=x-0.3, y=model_2, yend=data_2), color="grey") +
    geom_point( aes(x=x, y=data_1, color="Data",shape = Vaccine_Types), size=3, stroke = 2) +
    geom_point( aes(x=x-0.3, y=data_2, color="Data",shape = Vaccine_Types), size=3, stroke = 2) +
    geom_point( aes(x=x, y=model_1, color=model_name_1,shape = Vaccine_Types), size=3, stroke = 2) +
    geom_point( aes(x=x-0.3, y=model_2, color=model_name_2,shape = Vaccine_Types), size=3, stroke = 2) +
    #geom_segment( aes(x=x, xend=x, y=model_1, yend=data_1), color="grey") +
    #geom_point( aes(x=x, y=model_1, color=model_name_1, shape = Vaccine_Types), size=3, stroke = 2) +
    #geom_point( aes(x=x, y=data_1, color="Data",shape = Vaccine_Types), size=3,  stroke = 2) +
    #geom_point( aes(x=x, y=model_1, color=model_name_1), size=5, alpha = 0.7) +
    #geom_point( aes(x=x, y=data_1, color="Data"), size=5, alpha = 0.7) +
    #geom_point(aes(x=x, y=model_2, color=model_name_2), size=5, alpha = 0.7) +
    scale_color_manual(values = c("#E69F00","#56B4E9","#CC79A7"),
                       guide  = guide_legend(), 
                       name   = "Group") +
    scale_shape_manual(values=c(19, 1))+
    coord_flip()+
    #theme_ipsum() +
    theme(legend.position = c(.8,.8),legend.text = element_text(size = 20),legend.title = element_text(size = 20)) +
    ggtitle(year) +
    ylab("Frequency") +
    xlab("Clusters") +
    theme(axis.title  = element_text(size = 20), axis.text = element_text(size = 20), plot.title = element_text(size = 25,hjust = 0.5))  +
    ylim(0, max(max(lollipop_data_1$model_1)))
  grid.arrange(lollipop_plot_1 + scale_y_continuous(limits = c(NA,0.2)) + theme(plot.margin = unit(c(.5,0.5,1,0.5), "cm"),axis.text.y = element_blank()), ncol = 1, nrow=1, top = textGrob(plot_title,gp=gpar(fontsize=20,font=3)))
}

# Creating a joint Lollipop plot for VTs and NVTs and two models
lollipop_cluster_freqs_2x2_VTandNVT_labelSero <- function(year = "year unknown", plot_title = "Generic Plot Title",data_name = "Data", data1, model_name_1 ="Model 1", model1, model_name_2 ="Model 2", model2, VT_vec, SeroLabel){
  lollipop_data_1 <- data.frame(
    x=1:length(data1),
    model_1=model1,
    data_1=as.numeric(data1),
    Vaccine_Types <- sapply(VT_vec, function(x) if(x==1){"VT"}else{"NVT"}),
    data_2=as.numeric(data1),
    model_2 =model2
  )
  # Change baseline
  lollipop_plot_1 <- ggplot(lollipop_data_1) +
    geom_segment( aes(x=x, xend=x, y=model_1, yend=data_1), color="grey") +
    geom_segment( aes(x=x-0.3, xend=x-0.3, y=model_2, yend=data_2), color="grey") +
    geom_point( aes(x=x, y=data_1, color=data_name,shape = Vaccine_Types), size=3, stroke = 2) +
    geom_point( aes(x=x-0.3, y=data_2, color=data_name,shape = Vaccine_Types), size=3, stroke = 2) +
    geom_point( aes(x=x, y=model_1, color=model_name_1,shape = Vaccine_Types), size=3, stroke = 2) +
    geom_point( aes(x=x-0.3, y=model_2, color=model_name_2,shape = Vaccine_Types), size=3, stroke = 2) +
    #geom_segment( aes(x=x, xend=x, y=model_1, yend=data_1), color="grey") +
    #geom_point( aes(x=x, y=model_1, color=model_name_1, shape = Vaccine_Types), size=3, stroke = 2) +
    #geom_point( aes(x=x, y=data_1, color="Data",shape = Vaccine_Types), size=3,  stroke = 2) +
    #geom_point( aes(x=x, y=model_1, color=model_name_1), size=5, alpha = 0.7) +
    #geom_point( aes(x=x, y=data_1, color="Data"), size=5, alpha = 0.7) +
    #geom_point(aes(x=x, y=model_2, color=model_name_2), size=5, alpha = 0.7) +
    scale_color_manual(values = c("#E69F00","#56B4E9","#CC79A7"),
                       guide  = guide_legend(), 
                       name   = "Group") +
    scale_shape_manual(values=c(19, 1))+
    coord_flip()+
    #theme_ipsum() +
    theme(legend.position = c(.8,.8),legend.text = element_text(size = 20),legend.title = element_text(size = 20)) +
    ggtitle(year) +
    ylab("Frequency") +
    xlab("Clusters") +
    theme(axis.title  = element_text(size = 20), axis.text = element_text(size = 20), plot.title = element_text(size = 25,hjust = 0.5))  +
    ylim(0, max(max(lollipop_data_1$model_1)))
  grid.arrange(lollipop_plot_1 + scale_y_continuous(limits = c(NA,0.22)) + scale_x_continuous(breaks = 1:length(data1), labels = 1:length(data1), sec.axis = dup_axis(name = "Serotypes", labels = SeroLabel))+ theme(plot.margin = unit(c(.5,0.5,1,0.5), "cm")), ncol = 1, nrow=1, top = textGrob(plot_title,gp=gpar(fontsize=20,font=3)))
}
