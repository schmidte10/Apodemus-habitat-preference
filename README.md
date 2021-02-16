# Apodemus-habitat-preference
Description of statistical analysis run for 'investigating habitat preferences of Apodemus spp. in a traditionally managed agricultural landscape using geospatial analysis'. Code for only one species is shown, however in total three different species were examined within this species. 
## Author: Elliott Schmidt 
## Corresponding author:  
## Date: 15/02/2021 

## Set path
```{set working directory, warning=F, message=F}
setwd("C://PATH")
```
## load the necessary packages
```{r install packages, warning=F, message=F}
library(tidyr)
library(plyr)
library(dplyr)
library(readr)
library(vcd)
library(readxl)
library(janitor)
library(reshape2)
library(vegan)
library(ggrepel)
library(ggplot2)
library(GGally)
library(data.table)
library(cowplot)
library(gmodels)
library(car)
library(sn)
library(emmeans)
```
## Data importing, cleaning, and organizing
```{r data_organization, warning=F, message=F}
SmallMammalR_input <- read_excel("C://PATH/FOLDER/SmallMammalR_input.xlsx", 
                                 col_types = c("text", "text", "numeric", 
                                               "text", "numeric", "date", "numeric", 
                                               "text", "text", "text", "text", "text", 
                                               "text", "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "numeric", 
                                               "numeric", "numeric", "text", "skip"))

temp1 <- SmallMammalR_input[SmallMammalR_input$n_r=='N',] %>%             #extract only new captures (1 observation per individual)
  clean_names()%>%  #clean column names (super useful function)
  subset(species=='AA'|species=='AF'|species=='AS') %>%                   #extract Apodemus species
  unite("siteid2", village, year, habitat, date, sep="", remove = F) %>%  #making unique code for each site (step 1 of 2)
  unite("siteid", village, year, habitat, sep="", remove = F)             #making unique code for each site (step 2 of 2)

temp2 <- temp1 %>% reshape2::dcast(siteid~species, length)                #reformatting data relevant for species counts 
```

We need to find out how many traps were function for trapping nights in each location. For example traps not set, disturbed or stolen (all included as trap disturbances) were not included in trap nights because animals were not able to be captured in them

```
temp3 <-   temp1 %>% dplyr::select(siteid, siteid2, habitat,number_of_traps) %>%  #extract relevant columns
  distinct(siteid2, .keep_all = T)%>%                                             #removing duplicates
  group_by(siteid)%>%                                                             #grouping by the variable (siteid) - necessary for subsequent lines
  mutate(total_traps=sum(number_of_traps))%>%                                     #calculate the total number of traps for each traping period (3-6 consecutive days)
  distinct(siteid, .keep_all=T)%>%                                                #removing duplicates
  ungroup()%>%                                                                    #ungroup data
  full_join(temp2, by='siteid')                                                   #join with data frame temp2 - species captured and the number of traps used
  
temp4 <- SmallMammalR_input[SmallMammalR_input$n_r=="trap malfunction",] %>%      #extract samples that were recorded as trap malfunctions
  clean_names()%>%                                                                #we have seen this before - look at 'temp1' code
  unite("siteid2", village, year, habitat, date, sep="", remove = F)%>%           #code is very similar to 'temp1' refer back to previous code if confused
  unite("siteid", village, year, habitat, sep="", remove = F)%>% 
  dplyr::select(siteid, siteid2, habitat,number_of_traps,n_r)  

temp5 <- ddply(temp4,.(siteid2,n_r),nrow) %>%   #calculate the number of trap disturbances 
  full_join(temp4, by='siteid2')%>%             #joining data frame to 'temp4' so trap distrubances can be substraced from traps used
  distinct(siteid2, .keep_all=T)%>%             #removing duplicates
  group_by(siteid)%>%
  mutate(trap_mal=sum(V1))%>%                   #creating new column (just renaming V1 column that contains the number of trap dist.)
  dplyr::select(siteid,trap_mal)%>%             #extract relevant columns
  distinct(siteid, .keep_all=T)                 #remove duplicates

temp6 <- temp3 %>% full_join(temp5, by='siteid')%>%                 #joining 'temp3' with 'temp5' by 'siteid'
  replace_na(list(trap_mal=0)) %>%                                  #some 'NAs' are in the dataframe -sites where no animals were captured - NAs removed
  mutate(funct_trap=total_traps-trap_mal)%>%                        #creating functional trap column (number of working traps)
  mutate_each(funs(./funct_trap*100), matches("AS|AA|AF"))%>%       #transforming abundances to captures per hundred trap nights
  drop_na() %>%  mutate(habitat = case_when((habitat=="RE") ~ "WE", #some NA's are back because they were not removed from 'temp3'
                             (habitat=="LNV2")~ "LNV",              #correcting some data entry errors 'RE' should be 'WE' - LNV2 and MNV2 sites change to 
                             (habitat=="MNV2") ~ "MNV",             #LNV and MNV (somtime there were 2 LNV or MNV samples at a location)
                             TRUE ~ as.character(.$habitat)))%>% 
  subset(!habitat=="SCRUB")                                         #Scrub was a habitat used in earlier years but then abandoned importing land use data
```
Importing geospatial data 
```
temp7 <- read_excel("C://PATH/FOLDER/SmallMammal_coordinates_data (Salvat automat).xlsx", 
                                                           col_types = c("skip", "text", "numeric", 
                                                                         "numeric", "numeric", "text", "numeric", 
                                                                         "numeric", "numeric", "numeric"))%>% 
  setnames(old=c("Habitat type","distance to forest (<20 meters)","distance of forest (>20 meters)",   
                 "distance of nearest road/track(m)",
                 "distance to village center(km)"), new=c('habitat','tl_edge','for_edge', 
                                                          'dist_road','dist_village'))%>%
  clean_names()%>%
  unite("siteid", village, year, habitat, sep="", remove = F)%>% 
  dplyr::select(siteid, village, longitude, latitude, tl_edge, for_edge, dist_road, dist_village)%>% 
  mutate(dist_village=dist_village*1000)%>% 
  drop_na(longitude)
  
temp8 <- read_excel("C://PATH/FOLDER/SmallMammalBufferDataR.xlsx")%>% 
  setnames(old=c("Habitat type"), new=c('habitat'))%>%  #renaming columns
  clean_names()%>%                                                        #cleaning names
   unite("siteid", village, year, habitat, sep="", remove = F) %>%        #creating new unique variable
  full_join(temp7, by='siteid')%>%                                        #joining dataframes
  mutate(habitat = case_when((habitat=="RE") ~ "WE",                      #correcting some data entry errors
                             TRUE ~ as.character(.$habitat)))%>% 
  subset(!habitat=="Scrub")%>%                                            #removing 'scrub' habitat - see above
  full_join(temp6, by='siteid') %>%                                       #joining land-use dataframe to small mammal abundance information
  setnames(old=c("habitat.x","village.x"), new=c('habitat','village'))%>% #renaming columns
  drop_na(village)%>%                                                     #removing samples where NA's for villages - sites where no animals were captured
  drop_na(AA)%>%                                                          #removing samples where NA's for species - land use data available but no species data
  mutate(across(77:79, round, 0)) %>%                                     #species abundances were rounded so that a quassipoisson distribution could be used
  mutate(forest= con_for+mix_for+broad_leaved_for)%>%                     #all forest types combined into one variable
  mutate(forest_400= con_for_400+mix_for_400+broad_leaved_for_400)%>%     
  mutate(forest_1000= con_for_1000+mix_for_1000+broad_leaved_for_1000)
```
Reshaping data for difference analysis
```
temp9 <- temp8 %>% reshape2::melt(id.vars=c("siteid","habitat","year","village"), measure.vars=c("AF","AS","AA"))   #reshaping data for difference analysis
```
### Correlation plots
```{r correlation_plots, warning=F, message=F}
ggpairs(temp8[,c(9,11,13,15,17:22,26,82)], upper = list(continuous = wrap("cor", method = "spearman"))) 
```
### CCA
```{r CCA example,warning=F, message=F}
names(temp8[,c(9,11,15,17:22,26,82)])
names(temp8[,c(77:79)])
temp10 <- temp8 %>% rowwise() %>% filter(sum(c_across(77:79))!=0)          #dataframe to be used in cca ordination (remove rows with no captures)
set.seed(123)
ord <- cca(temp10[,c(77:79)] ~ elev_mean+twi_mean+lai_mean+disc_urb_fab+   #cca ordination
             non_irr_arable+perm_crops+complex_cult_patterns+
             agri_natural_veg+forest+trans_wood_shrub, data = temp10)
vif.cca(ord)                                                               #double checking for collinearity - 'pastures' was initially included but removed due to collinearity
summary(ord)[1]
summary(eigenvals(ord)) 
fit <- envfit(ord,  temp10[,c(9,11,15,17:19,21,22,26,82)], perm=999, display = "lc")
scores(fit, 'vectors')

cca_1 <- data.frame(summary(ord)[[1]])
cca_2 <- scores(fit, "vectors") %>% rbind(scores(fit, "factor")) %>% as.data.frame()

colors<- c("grey19","grey19","grey19")
sp <- c("A.agr","A.flav","A.sly")
var <- c("ELEV","LAI","TWI","DUF","NIA","PC","CCP","ANV","TWS","FOR")

cca.biplot_200 <- ggplot(cca_1, aes(x=CCA1, y=CCA2))+theme_classic()+theme(legend.position = "none",
  axis.text = element_text(size=12),axis.title.x = element_text(size=12),
  axis.title.y = element_text(size=12))+
  geom_hline(yintercept = 0, linetype="solid",color="grey81")+ 
  geom_vline(xintercept = 0, linetype="solid", color="grey81")+
  geom_segment(data=cca_2, aes(x=0, y=0, xend=CCA1, yend=CCA2), size=1, color="forestgreen",alpha =0.45)+
  geom_text_repel(data = cca_2,aes(x=CCA1,y=CCA2),label=var,hjust=0.5,vjust=-1)+
  geom_text_repel(aes(label=sp),fontface="italic",max.iter = 3000,hjust=0,vjust=0)+
  geom_point(data=cca_2, aes(x=CCA1, y=CCA2),shape=2,color="black",size=3)+ 
  geom_point(aes(color="black"),size=4,shape=16,fill="black",stroke=1.5)+
  scale_color_manual("sp",values = colors)+
  xlim(-1.5,1.5)+xlab("CCA1 [22.2%]")+ylab("CCA2 [4.8%]");cca.biplot_200 
```
### General linerized models 
```{r glm, warning=F, message=F} 
dat.qp.glm <- glm(AF ~ elev_mean+lai_mean+twi_mean+disc_urb_fab+perm_crops+
                    complex_cult_patterns+agri_natural_veg+trans_wood_shrub+forest, 
                  family="quasipoisson", data = temp8) 
par(mfrow=c(2,2))
plot(dat.qp.glm)
vif(dat.qp.glm)
fit.val <- data.frame(obs=temp8$AF, fitted=fitted(dat.qp.glm))
ggplot(fit.val, aes(y=fitted, x=obs))+geom_point()
summary(dat.qp.glm)
ci(dat.qp.glm)
```


### Interpopulation differences
```{r interpopulation_differences, warning=F, message=F} 
woodland_edge <- glm(value ~ variable, 
                  family="quasipoisson", data = temp9[(temp9$habitat=="WE"),])
par(mfrow=c(2,2))
plot(woodland_edge)
summary(woodland_edge)
ci(woodland_edge)
marginal = emmeans(woodland_edge, ~ variable)
pairs(marginal)
```


### Body condition analysis 
```{r body condition analysis, warning=F, message=F}
bc <- temp1 %>% 
  subset(age =="A") %>% 
  subset(reprod_cond=="TL"|reprod_cond=="TS"|reprod_cond=="TA"|
           reprod_cond=="IMP/NC"| reprod_cond=="NC"|reprod_cond=="IMP") %>% 
  drop_na(tl,fl,hb,weight)%>% 
  subset(weight>10) %>% #removing outliers (adults should not weight less than 20g 
  subset(hb>40) %>%  #head and body length should be great than 40mm - mistake in the data. outlier removed.
  subset(fl>15.5) %>%
  mutate(habitat = case_when((habitat=="RE") ~ "WE", 
                                  (habitat=="LNV2")~ "LNV", 
                                  (habitat=="MNV2") ~ "MNV", 
                                  TRUE ~ as.character(.$habitat)))%>% 
  subset(!habitat=="SCRUB") 
af <- bc %>% filter(species =="AF")
as <- bc %>% filter(species =="AS")
aa <- bc %>% filter(species =="AA")

## AF #####
pca_af <- princomp(~ tl+fl+hb, af, cor=T)
loadings(pca_af)
summary(pca_af)
screeplot(pca_af)
eigenvals(pca_af)
body_size_af <- (pca_af$scores[,1])
head(body_size_af)


bodycond_corr_af <- lm(af$weight~body_size_af)
summary(bodycond_corr_af)
plot(body_size_af ~ af$weight)
af$bodycondition <- bodycond_corr_af$residuals
shapiro.test(af$bodycondition)
glm <- selm(bodycondition~habitat, data=af)
summary(glm)
slot(glm, "opt.method")$convergence #sucessful if equals zero

sn.mple(y=af$bodycondition,opt.method = "nlminb")$cp
cp2dp(sn.mple(y=af$bodycondition,opt.method = "nlminb")$cp, family = "SN")

par(mfrow=c(1,2))
x <- rsn(n=100, xi=-7.83, omega = 10.616, alpha = 2.30)
hist(x)
hist(af$bodycondition)


ggplot(af, aes(y=bodycondition, x=habitat))+geom_boxplot()
```

## PLOTS 


### Plot 1 (see CCA section above)
```{r plot1}
plot_grid(cca.biplot_200, 
          nrow=1, 
          rel_widths = c(1))
```



### Plot2 
```{r plot2}
par(mfrow=c(2,2))
dat.qp.glm <- glm(AF ~ forest, 
                  family="quasipoisson", data = temp8) 
par(mar=c(4,5,4,4))
plot(AF ~ forest, data= temp8, type='n', ann=F, axes=F, ylim=c(0,40), xlim=c(0,50))
points(AF ~ forest, data= temp8,pch=16)
xs <- seq(0,50,l=1000)
ys <- predict(dat.qp.glm, newdata=data.frame(forest=xs), type = "response",se=T)
points(ys$fit ~ xs, col = "red", type = "l")
lines(ys$fit - 1 * ys$se.fit ~ xs, col = "black", type = "l", lty = 2)
lines(ys$fit + 1 * ys$se.fit ~ xs, col = "black", type = "l", lty = 2)
axis(1)
mtext("% forest cover", 1, cex = 1.5, line = 3)
axis(2, las = 2)
mtext("A. flavicollis/100 trapnights", 2, cex = 1.5, line = 3)
box(bty = "l")
text(40,32, "t-value =2.28", cex=1)
text(40,30, "p-value =0.027", cex=1)
text(0,40, "A", cex=2)

dat.qp.glm <- glm(AS ~ forest, 
                  family="quasipoisson", data = temp8) 
plot(AS ~ forest, data= temp8, type='n', ann=F, axes=F, ylim=c(0,40), xlim=c(0,50))
points(AS ~ forest, data= temp8,pch=16) 
xs <- seq(0,50,l=1000)
ys <- predict(dat.qp.glm, newdata=data.frame(forest=xs), type = "response",se=T)
points(ys$fit ~ xs, col = "grey", type = "l")
lines(ys$fit - 1 * ys$se.fit ~ xs, col = "black", type = "l", lty = 2)
lines(ys$fit + 1 * ys$se.fit ~ xs, col = "black", type = "l", lty = 2)
axis(1)
mtext("% forest cover", 1, cex = 1.5, line = 3)
axis(2, las = 2)
mtext("A. sylvaticus/100 trapnights", 2, cex = 1.5, line = 3)
box(bty = "l")
text(40,32, "t-value =1.92", cex=1)
text(40,30, "p-value =0.06", cex=1)
text(0,40, "B", cex=2)

dat.qp.glm <- glm(AA ~ twi_mean, 
                  family="quasipoisson", data = temp8) 
par(mar = c(4,5,4,4))
plot(AA ~ twi_mean, data= temp8, type='n', ann=F, axes=F, ylim=c(0,40), xlim=c(5,10.5))
points(AA ~ twi_mean, data= temp8,pch=16) 
xs <- seq(0,50,l=1000)
ys <- predict(dat.qp.glm, newdata=data.frame(twi_mean=xs), type = "response",se=T)
points(ys$fit ~ xs, col = "red", type = "l")
lines(ys$fit - 1 * ys$se.fit ~ xs, col = "black", type = "l", lty = 2)
lines(ys$fit + 1 * ys$se.fit ~ xs, col = "black", type = "l", lty = 2)
axis(1)
mtext("Topographical wetness index (TWI)", 1, cex = 1.5, line = 3)
axis(2, las = 2)
mtext("A. agrarius/100 trapnights", 2, cex = 1.5, line = 3)
box(bty = "l")
text(9.5,32, "t-value =2.69", cex=1)
text(9.5,30, "p-value =0.0097", cex=1)
text(5,40, "C", cex=2)

dat.qp.glm <- glm(AF ~ twi_mean, 
                  family="quasipoisson", data = temp8) 
par(mar = c(4,5,4,4))
plot(AF ~ twi_mean, data= temp8, type='n', ann=F, axes=F, ylim=c(0,40), xlim=c(5,10.5))
points(AF ~ twi_mean, data= temp8,pch=16) 
xs <- seq(0,50,l=1000)
ys <- predict(dat.qp.glm, newdata=data.frame(twi_mean=xs), type = "response",se=T)
points(ys$fit ~ xs, col = "grey", type = "l")
lines(ys$fit - 1 * ys$se.fit ~ xs, col = "black", type = "l", lty = 2)
lines(ys$fit + 1 * ys$se.fit ~ xs, col = "black", type = "l", lty = 2)
axis(1)
mtext("Topographical wetness index (TWI)", 1, cex = 1.5, line = 3)
axis(2, las = 2)
mtext("A. flavicollis/100 trapnights", 2, cex = 1.5, line = 3)
box(bty = "l")
text(9.5,32, "t-value =-2.00", cex=1)
text(9.5,30, "p-value =0.051", cex=1)
text(5,40, "D", cex=2)
```



### Plot3 
```{r plot3}
temp9$habitat<-  factor(temp9$habitat, levels=c("HNV", "MNV", "LNV", "WE"))

fig <- ggplot(temp9, aes(variable, value,color=habitat,fill=habitat))+ 
  geom_boxplot()+theme_classic()+ylim(0,45)+
  ylab("captures/100 trap night")+xlab("Species")+ 
  scale_x_discrete(labels=c("AF" = "A. flavicollis", "AS" = "A. sylvaticus", 
                            "AA" = "A. agrarius"))+
  scale_fill_manual(name= "Habitat", values = c("grey98", "grey77","grey40","grey18"), 
                    labels=c("High \n nature value","Medium \n nature value","Low \n nature value","Woodland edge"))+
  scale_color_manual(name = "Habitat", values = c("black", "black", "black","black"),guide=F)+
theme( 
  axis.title.x = element_text(size = 18), 
  axis.title.y = element_text(size = 18), 
  axis.text.x = element_text(size = 15), 
  axis.text.y = element_text(size = 15), 
  legend.position = "bottom",legend.title = element_blank(),legend.text = element_text(size=14))+
  annotate("text", x=0.72,y=21,label="A", size=6)+ 
  annotate("text", x=0.91,y=18,label="B", size=6)+
  annotate("text", x=1.10,y=14,label="B", size=6)+
  annotate("text", x=1.28,y=41,label="A", size=6)+ 
  
  annotate("text", x=1.72,y=15,label="B", size=6)+ 
  annotate("text", x=1.91,y=16,label="AB", size=6)+
  annotate("text", x=2.10,y=15,label="B", size=6)+
  annotate("text", x=2.28,y=30,label="A", size=6)+ 
  
  annotate("text", x=2.72,y=14,label="AB", size=6)+ 
  annotate("text", x=2.91,y=13,label="AB", size=6)+
  annotate("text", x=3.10,y=30,label="B", size=6)+
  annotate("text", x=3.28,y=15,label="A", size=6);fig
```



### Plot4
```{r plot4}
par(mfrow=c(1,2))
dat.qp.glm <- glm(AF ~ tl_edge, 
                  family="quasipoisson", data = temp8) 
par(mar=c(4,5,4,4))
plot(AF ~ tl_edge, data= temp8, type='n', ann=F, axes=F)
points(AF ~ tl_edge, data= temp8,pch=16)
xs <- seq(0,500,l=1000)
ys <- predict(dat.qp.glm, newdata=data.frame(tl_edge=xs), type = "response",se=T)
points(ys$fit ~ xs, col = "red", type = "l")
lines(ys$fit - 1 * ys$se.fit ~ xs, col = "black", type = "l", lty = 2)
lines(ys$fit + 1 * ys$se.fit ~ xs, col = "black", type = "l", lty = 2)
axis(1)
mtext("Distance from forest edge", 1, cex = 1.5, line = 3)
axis(2, las = 2)
mtext("A. flavicollis/100 trapnights", 2, cex = 1.5, line = 3)
box(bty = "l")
text(325,32, "t-value =-2.02", cex=1)
text(325,30, "p-value =0.049", cex=1)

dat.qp.glm <- glm(AF ~ dist_road, 
                  family="quasipoisson", data = temp8) 
par(mar = c(4,5,4,4))
plot(AF ~ dist_road, data= temp8, type='n', ann=F, axes=F)
points(AF ~ dist_road, data= temp8,pch=16) 
xs <- seq(0,550,l=1000)
ys <- predict(dat.qp.glm, newdata=data.frame(dist_road=xs), type = "response",se=T)
points(ys$fit ~ xs, col = "red", type = "l")
lines(ys$fit - 1 * ys$se.fit ~ xs, col = "black", type = "l", lty = 2)
lines(ys$fit + 1 * ys$se.fit ~ xs, col = "black", type = "l", lty = 2)
axis(1)
mtext("Distance from road", 1, cex = 1.5, line = 3)
axis(2, las = 2)
mtext("", 2, cex = 1.5, line = 3)
box(bty = "l")
text(450,32, "t-value =-2.48", cex=1)
text(450,30, "p-value =0.016", cex=1)
```
