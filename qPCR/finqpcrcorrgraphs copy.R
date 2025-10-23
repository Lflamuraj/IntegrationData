#clears everything on R
rm(list=ls())

library(ggplot2)
library(ggpubr)


data <- read.csv("All_data_qpcr_NEWSMAD7.csv")

#allom correstion
Man_corrected <- lm(Man ~ standard_length, data = data)
summary(Man_corrected)

Hyd_corrected <- lm(Hyd ~ standard_length, data = data)
summary(Hyd_corrected)

CB1_corrected <- lm(CB1 ~ standard_length, data = data)
summary(CB1_corrected)

CB5_corrected <- lm(CB5 ~ standard_length, data = data)
summary(CB5_corrected)

Clth_corrected <- lm(Clth ~ standard_length, data = data)
summary(Clth_corrected)

test<- lm(Hyd_corrected ~ Man_corrected)

#Man_residual <- residuals(Man_corrected)
#Hyd_residual <- residuals(Man_corrected)
#CB1_residual <- residuals(CB1_corrected)
#CB5_residual <- residuals(CB5_corrected)
#Clth_residual <- residuals(Clth_corrected)

#ManHyd_model <- lm(Hyd_residual ~ Man_residual)
#summary(ManHyd_model)
#plot(ManHyd_model)



#smad7 mand on x-axis
ggplot(data, aes(x= Man, y= Hyd))+ 
  geom_point(aes(shape = factor(Genus)), size = 4) +
  scale_shape_manual(values=c(15, 16, 17, 18, 6)) +
  geom_smooth(method= lm, se=FALSE, color= "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(color = "black",
                                                                        fill = NA,
                                                                        size = 1))

Manhydsmad7 <- lm(Hyd ~ Man, data = data)
summary(Manhydsmad7)

ggplot(data, aes(x= Man, y= CB1))+ 
  geom_point(aes(shape = factor(Genus)), size = 4) +
  scale_shape_manual(values=c(15, 16, 17, 18, 6)) +
  geom_smooth(method= lm, se=FALSE, color= "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(color = "black",
                                                                        fill = NA,
                                                                        size = 1))

Mancb1smad7 <- lm(CB1 ~ Man, data = data)
summary(Mancb1smad7)

ggplot(data, aes(x= Man, y= CB5))+ 
  geom_point(aes(shape = factor(Genus)), size = 4) +
  scale_shape_manual(values=c(15, 16, 17, 18, 6)) +
  geom_smooth(method= lm, se=FALSE, color= "black") +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), panel.border = element_rect(color = "black",
                                                                               fill = NA,
                                                                               size = 1))

Mancb5smad7 <- lm(CB5 ~ Man, data = data)
summary(Mancb5smad7)

ggplot(data, aes(x= Man, y= Clth))+ 
  geom_point(aes(shape = factor(Genus)), size = 4) +
  scale_shape_manual(values=c(15, 16, 17, 18, 6)) +
  geom_smooth(method= lm, se=FALSE, color= "black") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(color = "black",
                                                                        fill = NA,
                                                                        size = 1))
  
  Manclthsmad7 <- lm(Clth ~ Man, data = data)
  summary( Manclthsmad7)

#smad7 clth on x-axis
ggplot(data, aes(x= Clth, y= Man))+ 
  geom_point(aes(shape = factor(Genus)), size = 4) +
  scale_shape_manual(values=c(15, 16, 17, 18, 6)) +
  geom_smooth(method= lm, se=FALSE, color= "black") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(color = "black",
                                                                        fill = NA,
                                                                        size = 1))

Clthmandsmad7 <- lm(Man ~ Clth, data = data)
summary(Clthmandsmad7)

ggplot(data, aes(x= Clth, y= Hyd))+ 
  geom_point(aes(shape = factor(Genus)), size = 4) +
  scale_shape_manual(values=c(15, 16, 17, 18, 6)) +
  geom_smooth(method= lm, se=FALSE, color= "black") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(color = "black",
                                                                        fill = NA,
                                                                        size = 1))

Clthhydsmad7 <- lm(Hyd ~ Clth, data = data)
summary(Clthhydsmad7)


ggplot(data, aes(x= Clth, y= CB1))+ 
  geom_point(aes(shape = factor(Genus)), size =4) +
  scale_shape_manual(values=c(15, 16, 17, 18, 6)) +
  geom_smooth(method= lm, se=FALSE, color= "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(color = "black",
                                                                        fill = NA,
                                                                        size = 1))

ClthCB1smad7 <- lm(CB1 ~ Clth, data = data)
summary(ClthCB1smad7)


ggplot(data, aes(x= Clth, y= CB5))+ 
  geom_point(aes(shape = factor(Genus)), size = 4) +
  scale_shape_manual(values=c(15, 16, 17, 18, 6)) +
  geom_smooth(method= lm, se=FALSE, color= "black")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(color = "black",
                                                                        fill = NA,
                                                                        size = 1))

ClthCB5smad7 <- lm(CB5 ~ Clth, data = data)
summary(ClthCB5smad7)
------

sp7data <- read.csv("sp7_data_qpcr.csv")  
  
#sp7 mand on x-axis
  ggplot(sp7data, aes(x= Man, y= Hyd))+ 
    geom_point(aes(shape = factor(Genus)), size = 4) +
    scale_shape_manual(values=c(15, 16, 17, 18, 6)) +
  geom_smooth(method= lm, se=FALSE, color= "black") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), panel.border = element_rect(color = "black",
                                                                          fill = NA,
                                                                          size = 1))

Manhydsp7 <- lm(Hyd ~ Man, data = sp7data)
summary(Manhydsp7)

ggplot(sp7data, aes(x= Man, y= CB1))+ 
  geom_point(aes(shape = factor(Genus)), size = 4) +
  scale_shape_manual(values=c(15, 16, 17, 18, 6)) +
  geom_smooth(method= lm, se=FALSE, color= "black") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(color = "black",
                                                                        fill = NA,
                                                                        size = 1))

Mancb1sp7 <- lm(CB1 ~ Man, data = sp7data)
summary(Mancb1sp7)

ggplot(sp7data, aes(x= Man, y= CB5))+ 
  geom_point(aes(shape = factor(Genus)), size = 4) +
  scale_shape_manual(values=c(15, 16, 17, 18, 6)) +
  geom_smooth(method= lm, se=FALSE, color= "black") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(color = "black",
                                                                        fill = NA,
                                                                        size = 1))

Mancb5sp7 <- lm(CB5 ~ Man, data = sp7data)
summary(Mancb5sp7)

ggplot(sp7data, aes(x= Man, y= Clth))+ 
  geom_point(aes(shape = factor(Genus)), size = 4) +
  scale_shape_manual(values=c(15, 16, 17, 18, 6)) +
  geom_smooth(method= lm, se=FALSE, color= "black") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(color = "black",
                                                                        fill = NA,
                                                                        size = 1))


Manclthsp7 <- lm(Clth ~ Man, data = sp7data)
summary( Manclthsp7)

#sp7 clth on x-axis
ggplot(sp7data, aes(x= Clth, y= Man))+ 
  geom_point(aes(shape = factor(Genus)), size = 4) +
  scale_shape_manual(values=c(15, 16, 17, 18, 6)) +
  geom_smooth(method= lm, se=FALSE, color= "black") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(color = "black",
                                                                        fill = NA,
                                                                        size = 1))


Clthmandsp7 <- lm(Man ~ Clth, data = sp7data)
summary(Clthmandsp7)

ggplot(sp7data, aes(x= Clth, y= Hyd))+ 
  geom_point(aes(shape = factor(Genus)), size = 4) +
  scale_shape_manual(values=c(15, 16, 17, 18, 6)) +
  geom_smooth(method= lm, se=FALSE, color= "black") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(color = "black",
                                                                        fill = NA,
                                                                        size = 1))

Clthhydsp7 <- lm(Hyd ~ Clth, data = sp7data)
summary(Clthhydsp7)

ggplot(sp7data, aes(x= Clth, y= CB1))+ 
  geom_point(aes(shape = factor(Genus)), size = 4) +
  scale_shape_manual(values=c(15, 16, 17, 18, 6)) +
  geom_smooth(method= lm, se=FALSE, color= "black") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(color = "black",
                                                                        fill = NA,
                                                                        size = 1))


ClthCB1sp7 <- lm(CB1 ~ Clth, data = sp7data)
summary(ClthCB1sp7)

ggplot(sp7data, aes(x= Clth, y= CB5))+ 
  geom_point(aes(shape = factor(Genus)), size = 4) +
  scale_shape_manual(values=c(15, 16, 17, 18, 6)) +
  geom_smooth(method= lm, se=FALSE, color= "black") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(color = "black",
                                                                        fill = NA,
                                                                        size = 1))

ClthCB5sp7 <- lm(CB5 ~ Clth, data = sp7data)
summary(ClthCB5sp7)
------
#check of smad and sp7
  
checkdata <- read.csv("crosscheck_sp7_and_smad7.csv")

#mand
ggplot(checkdata, aes(x= Man_Sp7, y= Man_smad7))+ 
  geom_point(aes(color = factor(Spp)), size = 2.5) +
  geom_smooth(method= lm, se=FALSE, color= "black")

Man_gene <- lm(Man_smad7 ~ Man_Sp7, data = checkdata)
summary(Man_gene)

#hyd
ggplot(checkdata, aes(x= Hyd_Sp7, y= Hyd_smad7))+ 
  geom_point(aes(color = factor(Spp)), size = 2.5) +
  geom_smooth(method= lm, se=FALSE, color= "black")

Hyd_gene <- lm(Hyd_smad7 ~ Hyd_Sp7, data = checkdata)
summary(Hyd_gene)

#cb1
ggplot(checkdata, aes(x= CB1_Sp7, y= CB1_smad7))+ 
  geom_point(aes(color = factor(Spp)), size = 2.5) +
  geom_smooth(method= lm, se=FALSE, color= "black")

CB1_gene <- lm(CB1_smad7 ~ CB1_Sp7, data = checkdata)
summary(CB1_gene)

#cb5
ggplot(checkdata, aes(x= CB5_Sp7, y= CB5_smad7))+ 
  geom_point(aes(color = factor(Spp)), size = 2.5) +
  geom_smooth(method= lm, se=FALSE, color= "black")

CB5_gene <- lm(CB5_smad7 ~ CB5_Sp7, data = checkdata)
summary(CB5_gene)

#clth
ggplot(checkdata, aes(x= Clth_Sp7, y= Clth_smad7))+ 
  geom_point(aes(color = factor(Spp)), size = 2.5) +
  geom_smooth(method= lm, se=FALSE, color= "black")

Clth_gene <- lm(Clth_smad7 ~ Clth_Sp7, data = checkdata)
summary(Clth_gene)


------
#animals removed 
  
data2 <-read.csv("All_data_qpcr_NEWSMAD7_OLRemoved.csv")

ggplot(data2, aes(x= Man, y= Hyd))+ 
  geom_point(aes(color = factor(Spp)), size = 2.5) +
  geom_smooth(method= lm, se=FALSE, color= "black")

OLRManhydsmad7 <- lm(Hyd ~ Man, data = data2)
summary(OLRManhydsmad7)

ggplot(data2, aes(x= Man, y= CB1))+ 
  geom_point(aes(color = factor(Spp)), size = 2.5) +
  geom_smooth(method= lm, se=FALSE, color= "black")

OLRMancb1smad7 <- lm(CB1 ~ Man, data = data2)
summary(OLRMancb1smad7)

ggplot(data2, aes(x= Man, y= CB5))+ 
  geom_point(aes(color = factor(Spp)), size = 2.5) +
  geom_smooth(method= lm, se=FALSE, color= "black")

OLRMancb5smad7 <- lm(CB5 ~ Man, data = data2)
summary(OLRMancb5smad7)

ggplot(data2, aes(x= Man, y= Clth))+ 
  geom_point(aes(color = factor(Spp)), size = 2.5) +
  geom_smooth(method= lm, se=FALSE, color= "black")
       
OLRManclthsmad7 <- lm(Clth ~ Man, data = data2)
summary(OLRManclthsmad7)
       
#smad7 clth on x-axis
ggplot(data2, aes(x= Clth, y= Man))+ 
  geom_point(aes(color = factor(Spp)), size = 2.5) +
  geom_smooth(method= lm, se=FALSE, color= "black")
       
OLRClthmandsmad7 <- lm(Man ~ Clth, data = data2)
 summary(OLRClthmandsmad7)
       
ggplot(data2, aes(x= Clth, y= Hyd))+ 
  geom_point(aes(color = factor(Spp)), size = 2.5) +
  geom_smooth(method= lm, se=FALSE, color= "black")
       
OLRClthhydsmad7 <- lm(Hyd ~ Clth, data = data2)
summary(OLRClthhydsmad7)
       
ggplot(data2, aes(x= Clth, y= CB1))+ 
  geom_point(aes(color = factor(Spp)), size = 2.5) +
  geom_smooth(method= lm, se=FALSE, color= "black")
       
OLRClthCB1smad7 <- lm(CB1 ~ Clth, data = data2)
summary(OLRClthCB1smad7)
       
ggplot(data2, aes(x= Clth, y= CB5))+ 
  geom_point(aes(color = factor(Spp)), size = 2.5) +
  geom_smooth(method= lm, se=FALSE, color= "black")
       
OLRClthCB5smad7 <- lm(CB5 ~ Clth, data = data2)
summary(OLRClthCB5smad7)
------
  ###heat map###
  data_heat <- data.frame(
    dims = c("Mandible", "Cleithrum"),
    Mandible = c(0, 0.6928),
    Hyoid = c(0.5498, 0.673),
    CB1 = c(0.6007, 0.5443),
    CB5 = c(0.5125, 0.2817),
    Cleithrum = c(0.6928, 0)
  )



library(reshape2)

# Reshape the data for heatmap
data_melted <- melt(data_heat, id.vars = "dims", variable.name = "variables", value.name = "value")


library(ggplot2)
  
  # Create the heatmap
  ggplot(data_melted, aes(x = dims, y = variables, fill = value)) +
  geom_tile() +
    scale_fill_gradient2(low = "cornflowerblue", mid = "khaki1", high = "palevioletred2", midpoint = 0.5, 
                       limits = c(0.0, 1.0), oob = scales::squish) +
  labs(title = "Heatmap of Data", fill = "R-squared") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45,         # Angle of text (rotate text)
                               hjust = 1,          # Horizontal justification
                               size = 12,          # Font size of the x-axis labels
                               color = "black",     # Color of the x-axis labels
                               family = "Arial"),  # Keep the angle for x-axis text
    axis.text.y = element_text(angle = 45,         # Angle of text (rotate text)
                               hjust = 1,          # Horizontal justification
                               size = 12,          # Font size of the x-axis labels
                               color = "black",     # Color of the x-axis labels
                               family = "Arial"),  # Keep the angle for y-axis text
    axis.title.x = element_blank(),  # Remove x-axis label
    axis.title.y = element_blank(),  # Remove y-axis label
    axis.ticks = element_blank()     # Remove axis ticks
  ) 

------
  ###heat map###
    data_heat2 <- data.frame(
      dims = c("Mandible", "Cleithrum"),
      Mandible = c(0, 0.0333),
      Hyoid = c(0.03078, 0.0528),
      CB1 = c(0.06079, 0.03298),
      CB5 = c(0.09883, 0.5946),
      Cleithrum = c(0.03333, 0)
    )
  
  
  
  library(reshape2)
  
  # Reshape the data for heatmap
  data_melted2 <- melt(data_heat2, id.vars = "dims", variable.name = "variables", value.name = "value")
  
  
  library(ggplot2)
  
  # Create the heatmap
  ggplot(data_melted2, aes(x = dims, y = variables, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "cornflowerblue", mid = "khaki1", high = "palevioletred2", midpoint = 0.5,
                         limits = c(0.0, 1.0), oob = scales::squish) +
    labs(title = "Heatmap of Data", fill = "R-squared") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45,         # Angle of text (rotate text)
                                 hjust = 1,          # Horizontal justification
                                 size = 12,          # Font size of the x-axis labels
                                 color = "black",     # Color of the x-axis labels
                                 family = "Arial"),  # Keep the angle for x-axis text
      axis.text.y = element_text(angle = 45,         # Angle of text (rotate text)
                                 hjust = 1,          # Horizontal justification
                                 size = 12,          # Font size of the x-axis labels
                                 color = "black",     # Color of the x-axis labels
                                 family = "Arial"),  # Keep the angle for y-axis text
      axis.title.x = element_blank(),  # Remove x-axis label
      axis.title.y = element_blank(),  # Remove y-axis label
      axis.ticks = element_blank()     # Remove axis ticks
    ) 
  
