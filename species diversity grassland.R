library(data.table)
library(ggplot2)
library(dplyr)
library(terra)
library(cowplot)
library(effects)
library(jtools)
library(piecewiseSEM)
library(lmerTest)
library(grfxtools)
library(glmm.hp)
library(ggsci)

#read sample Data##
rawData <- fread('C:\\Users\\Administrator\\Desktop\\awData.csv') 
df_scale <- fread('C:\\Users\\Administrator\\Desktop\\df_scale.csv')
df_scale2 <- fread('C:\\Users\\Administrator\\Desktop\\df_scale2.csv') 
#remote data
grassData <- rast('E:\\grassAllData.tif')
names(grassData) <- c('pr','srad','tmmn','tmmx','vpd','vs',
                      'spNum','agb','bgb','soc','TN','soilBulk')
phenoData <- rast('E:\\grassPhenoData.tif')
names(phenoData) <- c('NDVIsos','EVIsos','NIRvsos','NDVIeos','EVIeos',
                      'NIRveos','SIFsos','SIFeos')
`add<-`(grassData,phenoData)
crs <- '+proj=longlat +datum=WGS84'
globalRaster <- rast(vals=1:518400,nrows=360, ncols=1440,xmin=-180, xmax=180,ymin=0, ymax=90,crs=crs)
grassData2 <- resample(grassData,globalRaster)
######Fig.1b###########
#eos~div model
model_eos_div <- lmer(NIRveos~spNum+(1|grassType),data=rawData)
rsquared(model_eos_div)$Conditional %>% round(2)
out <- summary(model_eos_div)
trend_unscaled <- out$coefficients["spNum","Estimate"]/ sd(rawData$spNum)
error_unscaled <- out$coefficients["spNum","Std. Error"]/ sd(rawData$spNum)
trend_unscaled
error_unscaled

parEOSDiv <- partialize(model_eos_div,"spNum")
out_model_eos_div <- allEffects(model_eos_div,partial.residuals = TRUE)
# windowsFonts(Times_New_Roman=windowsFont("Times New Roman"))
ggplot_eos_div <- function(x){
  df <- tibble(upper = x$`spNum`$upper[,1],
               lower = x$`spNum`$lower[,1],
               off   = x$`spNum`$fit[,1],
               spNum = x$`spNum`$x[,1])
  gg <- ggplot() +
    geom_hex(data = parEOSDiv, aes(spNum,NIRveos)) +
    scale_fill_gradientn("",colours = alpha(colorRampPalette( c("#2878b5","#ff8884","#f3d266"))(5),0.7),
                         trans = "log", limits=c(1,6),breaks=c(1,2,4,6)) +
    geom_ribbon(data = df, aes(x = spNum, ymin = lower, ymax = upper), alpha = 0.3, fill = "#2878b5") +
    geom_line(data = df, aes(spNum, off), col = "#2878b5",linewidth=1) +
    theme_bw() + 
    theme(
      panel.grid=element_blank(),
      legend.position = c(0.9,0.8),
      legend.key.size = unit(0.3,"cm"),
      legend.background = element_rect('transparent'),
      legend.text = element_text(size=5),
      axis.text = element_text(size=8),
      axis.title = element_text(size=8),
    )+
    xlab("Diversity")+
    ylab("EOS (DOY)")
  return(gg)
}

lab1 <- "paste(\" EOS ~ Diversity \",italic(p), \" < 0.001*\")"

fig1b<-ggplot_eos_div(out_model_eos_div)+
  annotate("text", x = 14, y = 310, label = lab1,parse = TRUE,
           size = 2)
#ggsave("C:\\Users\\Administrator\\Desktop\\草地多样性\\fig1b.pdf", fig1b,width=8, height=7, units="cm")

######Fig.1c###########
# library(RColorBrewer)
library(paletteer)
polarPlot <- function(polarData,x,y,cor,statData){
  Polar <- ggpolar(pole = "N", max.lat = 80, min.lat = 25,
                   max.lon = 180, min.lon = -180,
                   lat.ax.vals = c(30, 60),
                   longitude.spacing = 60,
                   ax.labs.size = 1,
                   lat.ax.labs.pos = -150,
                   plt.lat.labels=FALSE,
                   plt.lon.labels=FALSE,
                   size.axes = 0.1,
                   size.outer = 0.1,
                   size=0.2,
                   clip = "off",
                   land.fill.colour = "#F2F2F2",#"transparent",#透明色#ffffcc #E8DDDC
                   land.outline.colour = "gray1")+
    theme(
      axis.text.y = element_blank(),
      axis.ticks=element_blank(),
    )+
    geom_tile(data = fcor2,aes_string(x=x,y=y,fill=cor))+
    # geom_point(data=polarData, aes_string(x = x, y = y, color = cor), size = 1)+
    scale_fill_gradientn(colors=paletteer_c("ggthemes::Orange-Blue-White Diverging", 4),limits=c(-1,1))+
    cowplot::theme_cowplot()+
    theme(text = element_text(size = 2),
          # plot.title = element_text(hjust = 0.5, vjust = -60, size = 5),
          legend.position = c(1,0.55),legend.direction="vertical",
          legend.key.size = unit(0.2,"cm"),
          legend.text=element_text(size=4),legend.title=element_text(size=6),
          plot.margin = margin(t = 0.1,r = 0.3, b = 0.1, l = 0.1, unit = 'inch')
          )
  
  Hist <- ggplot(data = statData)+
    geom_bar(mapping = aes(x =type, y = val, fill = type),
             stat = 'identity', width = 0.5,position = 'dodge2',show.legend = F)+
    scale_y_continuous(breaks = seq(0,60,30), labels = paste0(seq(0,60,30),'%'))+
    #scale_fill_viridis_d(direction = -1,option = "H")+
    scale_fill_manual(values=c("#9E3D22","#2B5C8A"))+
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 5, color = 'black'),
          axis.line = element_line(linetype = 1, size =0.5),
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"),
          panel.grid = element_blank()
          # axis.ticks.length.x = unit(-0.15, "cm")
    )
  
  Plot <- ggdraw()+
    draw_plot(Polar, 0,0,1,1)+
    draw_plot(Hist, 0.1, 0.05,0.45,0.3)
}

eos_diversity_tif <- grassData2[[c('NIRveos','spNum')]] %>% terra::scale()
fcor <- focalReg(eos_diversity_tif,9,"ols",na.rm=T) %>% fortify() %>% na.omit()

names(fcor)[4] <- 'slope'

stat <- rbind(data.table(type='positive',val=sum(fcor$slope>0)*100/nrow(fcor)),
              data.table(type='negtative',val=sum(fcor$slope<=0)*100/nrow(fcor)))

divEOSplot <- polarPlot(fcor,'x','y','slope',stat)
# divEOSplot


######Fig2a-d#############
##fig2a####
#eos~div+agb model
model_eos_agb_div <- lmer(NIRveos~agb+spNum+(1|grassType),data=rawData)
out <- summary(model_eos_agb_div)
rsquared(model_eos_agb_div)$Conditional %>% round(2)

trend_unscaled <- out$coefficients["agb","Estimate"]/ sd(rawData$agb)
error_unscaled <- out$coefficients["agb","Std. Error"]/ sd(rawData$agb)
trend_unscaled
error_unscaled

parDivAgbEOS <- partialize(model_eos_agb_div,"agb")
out_model_eos_agb_div <- allEffects(model_eos_agb_div,partial.residuals = TRUE)
ggplot_eos_agb_div <- function(x){
  df <- tibble(upper = x$`agb`$upper[,1],
               lower = x$`agb`$lower[,1],
               off   = x$`agb`$fit[,1],
               agb = x$`agb`$x[,1])
  gg <- ggplot() +
    geom_hex(data = parDivAgbEOS, aes(agb,NIRveos),show.legend = F) +
    scale_fill_gradientn("",colours = alpha(colorRampPalette( c("#2878b5","#ff8884","#f3d266"))(5),0.7),
                         trans = "log", limits=c(1,6),breaks=c(1,2,4,6)) +
    geom_ribbon(data = df, aes(x = agb, ymin = lower, ymax = upper), alpha = 0.3, fill = "#2878b5") +
    geom_line(data = df, aes(agb, off), col = "#2878b5",linewidth=1) +
    theme_bw() + 
    theme(
      panel.grid=element_blank(),
      legend.position = c(0.9,0.8),
      legend.key.size = unit(0.3,"cm"),
      legend.background = element_rect('transparent'),
      legend.text = element_text(size=5),
      axis.text = element_text(size=8),
      axis.title = element_text(size=8),
    )+
    xlab(expression(AGB~("g m"^ -1)))+
    ylab("EOS (DOY)")
  return(gg)
}
lab1 <- "paste(\" EOS ~ AGB + Diversity \",italic(p), \" = 0.701\")"


fig32<-ggplot_eos_agb_div(out_model_eos_agb_div)+
  annotate("text", x = 400, y = 310, label = lab1,parse = TRUE,
           size = 2)


##fig2b####
model_eos_bgb_div <- lmer(NIRveos~bgb+spNum+(1|grassType),data=rawData)
out <- summary(model_eos_bgb_div)
# rsquared(model_eos_bgb_div)$Conditional %>% round(2)
trend_unscaled <- out$coefficients["bgb","Estimate"]/ sd(rawData$bgb)
error_unscaled <- out$coefficients["bgb","Std. Error"]/ sd(rawData$bgb)
trend_unscaled
error_unscaled

parEOSbgbDiv <- partialize(model_eos_bgb_div,"bgb")
out_model_eos_bgb_div <- allEffects(model_eos_bgb_div,partial.residuals = TRUE)
ggplot_eos_bgb_div <- function(x){
  df <- tibble(upper = x$`bgb`$upper[,1],
               lower = x$`bgb`$lower[,1],
               off   = x$`bgb`$fit[,1],
               bgb = x$`bgb`$x[,1])
  gg <- ggplot() +
    geom_hex(data = parEOSbgbDiv, aes(bgb,NIRveos),show.legend = F) +
    scale_fill_gradientn("",colours = alpha(colorRampPalette( c("#2878b5","#ff8884","#f3d266"))(5),0.7),
                         trans = "log", limits=c(1,6),breaks=c(1,2,4,6)) +
    geom_ribbon(data = df, aes(x = bgb, ymin = lower, ymax = upper), alpha = 0.3, fill = "#2878b5") +
    geom_line(data = df, aes(bgb, off), col = "#2878b5",linewidth=1) +
    theme_bw() + 
    theme(
      panel.grid=element_blank(),
      legend.position = c(0.9,0.8),
      legend.key.size = unit(0.3,"cm"),
      legend.background = element_rect('transparent'),
      legend.text = element_text(size=5),
      axis.text = element_text(size=8),
      axis.title = element_text(size=8),
    )+
    xlab(expression(BGB~("g m"^ -1)))+
    ylab("EOS (DOY)")
  return(gg)
}
lab1 <- "paste(\" EOS ~ BGB + Diversity \",italic(p), \" < 0.001*\")"
fig35<-ggplot_eos_bgb_div(out_model_eos_bgb_div)+
  annotate("text", x = 3500, y = 310, label = lab1,parse = TRUE,
           size = 2)


#####fig2c####
eos_agb <- grassData2[[c('NIRveos','agb','spNum')]] %>% terra::scale()
fcor <- focalReg(eos_agb,9,na.rm=T) %>% fortify(.) %>% na.omit(.)

names(fcor)[4] <- 'slope'

agbstat <- rbind(data.table(type='positive',val=sum(fcor$slope>=0)*100/nrow(fcor)),
                 data.table(type='negtative',val=sum(fcor$slope<0)*100/nrow(fcor)))



agbEOSplot <- polarPlot(fcor,'x','y','slope',agbstat)


##fig2d####
eos_bgb <- grassData[[c('NIRveos','bgb','spNum')]] %>% terra::scale()
fcor <- focalReg(eos_bgb,9,na.rm=T) %>% fortify(.) %>% na.omit(.)
names(fcor)[4] <- c('slope')

bgbstat <- rbind(data.table(type='positive',val=sum(fcor$slope>0)*100/nrow(fcor)),
                 data.table(type='negtative',val=sum(fcor$slope<=0)*100/nrow(fcor)))

bgbPlot <- polarPlot(fcor,'x','y','slope',bgbstat)
bgbPlot



####Fig.2e##########
df_forest <- df_scale[,c('NIRveos','agb','bgb' ,'spNum','soilMoist','soc','n10','pr_growing','soilWeight5','vpd_growing','srad_growing','tmean_growing','grassType')]

colnames(df_forest) <- c('EOS','AGB','BGB','Diversity','SM','SOC','Soil_N','Precipitation','BD','VPD','SRAD','Temperature','grassType')
lme_SMF<-lmer(EOS ~ AGB + BGB +Diversity+SM+SOC+Soil_N+Precipitation+BD+VPD+SRAD+Temperature+(1|grassType),data=df_forest)
lmeResult <- summary(lme_SMF)
result <- lmeResult$coefficients %>% as.data.table()
result$index <- rownames(lmeResult$coefficients)
result <- subset(result,result$index!='(Intercept)')

df_forest$biomass <- df_forest$AGB+df_forest$BGB
df_forest$climate <- df_forest$Precipitation+df_forest$SRAD+df_forest$VPD+df_forest$Temperature
df_forest$soil <- df_forest$Soil_N+df_forest$SOC+df_forest$SM+df_forest$BD

glmLme<-lmer(EOS~soil+climate+biomass+Diversity+(1|grassType),data=df_forest)
lmerCoef <- summary(glmLme)$coef
allSum <- sum(abs(lmerCoef[2:5,1]))
soil <- abs(lmerCoef[2,1])/allSum*100
climate <- abs(lmerCoef[3,1])/allSum*100
biomass <- abs(lmerCoef[4,1])/allSum*100
diversity <- abs(lmerCoef[5,1])/allSum*100
# glmmmResults <- glmm.hp(glmLme)


result[which(result$index %in% c("AGB","BGB")),'type'] <- 'biomass'
result[which(result$index %in% c( "Precipitation","SRAD","VPD",'Temperature')),'type'] <- 'climate'
result[which(result$index %in% c( "Soil_N","BD",'SOC','SM')),'type'] <- 'soil'
result[which(result$index %in% c( "Diversity")),'type'] <- 'Diversity'

result <- result[order(result$type, result$ Estimate), ]
result$index <- factor(result$index, levels = c('Diversity','AGB','BGB','SM','SOC','Soil_N','BD','Precipitation','VPD','SRAD','Temperature'))
result$type <- factor(result$type, levels = c('Diversity','biomass','soil','climate'))

# windowsFonts(A=windowsFont("Times New Roman"),
#              B=windowsFont("Times New Roman"))
p1 <- ggplot(result, aes(x = index, y = Estimate, color = type)) +
  geom_point(size = 6) + 
  geom_errorbar(aes(ymin = Estimate - `Std. Error`, ymax = Estimate + `Std. Error`), width = 0.2, size =1,show.legend = F) +  
  scale_color_manual(values = c("#ff8884","#f3d266","#91D1C2FF","#2878b5")) +
  coord_flip() +  
  geom_hline(yintercept = 0, linetype = 2,size=1) + 
  labs(x = '', y = 'estimate', color = '') +  
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.line.x = element_line(), axis.ticks.y = element_blank(), 
        legend.key = element_blank()) +
  scale_x_discrete(breaks = result$index, labels = result$index, position = 'top') +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.5, 0.5))+
  theme(text=element_text(size=20))

p1
lfly <- data.frame(fre =c(climate,soil,biomass,diversity) , diet = c(''), beh = c('Diversity','biomass','soil','climate')) 
lfly$beh <- factor(lfly$beh,levels = lfly$beh)
p2 <- ggplot(lfly, aes(diet, fre, fill = beh))+
  geom_bar(stat = 'identity', position = 'fill', width = 0.3)+
  geom_col(width = 0.9) +
  labs(x = bquote(R^2==0.37), y = 'relative estimates (%)')+
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.ticks.x = element_blank(), axis.line.y = element_line(color = 'gray30'))+
  scale_fill_manual(values = c("#2878b5","#91D1C2FF","#f3d266","#ff8884"))+
  scale_x_discrete(position = 'top') +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))+
  theme(text=element_text(size=20))

p2
library(patchwork)
p2 <- p2 + theme(legend.position = 'none')
p2 + p1 + plot_layout(ncol = 2, widths = c(1, 2))


####Fig.2f##########
library(tidyverse)
library(dplyr)

###sem###
dataClimate <-  df_scale[,c('pr_growing','srad_growing','vpd_growing','tmean_growing')]
dataClimate %>%
  prcomp()->df_pca
result_pca <-data.frame(df_pca$x) 
summary(df_pca)
df_pca$rotation
climate_pca<-result_pca[,c(1:2)]

dataSoil <-  df_scale[,c('n','soilWeight5','soilMoist','soc')]
dataSoil %>%
  prcomp()->df_pca
result_pca <-data.frame(df_pca$x) 
summary(df_pca)
df_pca$rotation
soil_pca<-result_pca[,c(1:2)]

df_scale$pc1_clim <- climate_pca$PC1
df_scale$pc2_clim <- climate_pca$PC2
df_scale$pc1_soil <- soil_pca$PC1
df_scale$pc2_soil <- soil_pca$PC2

df_scale2 <- fread('C:\\Users\\Administrator\\Desktop\\草地多样性\\df_scale2.csv')

psem_lm2 <- psem(
  lm(agb~spNum+pc1_clim+pc2_clim+pc1_soil+pc2_soil,data = df_scale2),
  lm(bgb~spNum+pc1_clim+pc2_clim+pc1_soil+pc2_soil,data = df_scale2),
  # lm(soc~spNum+pc1_clim+pc2_clim+pc1_soil+(1|grassType),data = df_scale2),
  lm(NIRveos~agb+bgb+spNum,data = df_scale2),
  data=df_scale2
)
summary(psem_lm2)


#####supplement########
#######Fig1##########
library(ggplot2)
library(data.table)
library(terra)
library(tidyterra)
library(ggspatial)
library(cowplot)
library(stringr)
library(dplyr)
crs <- '+proj=longlat +datum=WGS84'
chinaProvinceBorder <- vect('E:\\doctor thesis\\data\\border\\bou2_4p.shp',crs=crs)
border1 <- vect('E:\\doctor thesis\\data\\border\\SouthSea\\九段线.shp',crs=crs)
dem <- rast('E:\\doctor thesis\\data\\enviroment\\elev.tif')
chinaDem <- crop(dem,chinaProvinceBorder,mask=T)
point <- fread('E:\\草地多样性\\allData.csv')


chinaMap <- ggplot()+
  geom_spatvector(data=border1)+
  geom_spatraster(data=chinaDem)+
  scale_fill_whitebox_c(palette = "deep",direction=-1,
                        guide = guide_legend(title = 'dem'))+
  geom_spatvector(data=chinaProvinceBorder,color='white',fill=NA,linewidth = 0.5,show.legend = F,na.rm=T)+
  geom_point(data = point,aes(lon,lat,color=grassType))+
  scale_color_discrete(name='grassType')+
  coord_sf(ylim = c(15,55),crs=crs)+ #限制地图分布范围
  annotation_scale(location = "bl")+ #比例尺
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering)+ #指北针
  #修改图例
  # scale_fill_discrete(name = "图例", labels = c('四川省'))+
  theme_bw()+
  theme(text = element_text(family = "Times_New_Roman",face='bold'),
        axis.text = element_text(family = 'Times_New_Roman',size = 14,face = 'bold'),
        axis.title.x = element_text(family = 'Times_New_Roman',size = 16,face = 'bold'),
        axis.title.y = element_text(family = 'Times_New_Roman',size = 16,face = 'bold'),
        #修改刻度线内
        axis.ticks.length=unit(0.2, "cm"),
        #加宽图边框
        #panel.border = element_rect(size=1),
        plot.background = element_rect(color = "white"),
        # axis.line = element_line(size = .8),
        axis.ticks = element_line(linewidth = .8),
        #添加网格线
        panel.grid = element_line(colour = "grey80",linewidth=.2),
        # 图例格式
        legend.title.align = 0.5, #居中
        legend.position = 'right')#设置图例位置
chinaMap
nine_map <- ggplot() +geom_spatvector(data=chinaProvinceBorder,linewidth = 1.5,show.legend = F,na.rm=T)+
  geom_spatvector(data=border1)+
  coord_sf(ylim = c(0,25),xlim = c(105,125),crs=crs)+
  theme(
    #aspect.ratio = 1.25, #调节长宽比
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(fill=NA,color="grey10",linetype=1,linewidth =1.),
    plot.margin=unit(c(0,0,0,0),"mm"))
nine_map
gg_inset_map <-  ggdraw() +
  draw_plot(chinaMap) +
  draw_plot(nine_map, x = 0.65, y = 0.16, width = 0.1, height = 0.25)
gg_inset_map


######figS3##
library(lme4)
library(modelr)
library(broom)
library(broom.mixed)
library(lmerTest)
library(tidyr)
library(ggplot2)
library(ggpubr)
# Model with varying slope and intercept??

m1 <- lmerTest::lmer(NIRveos ~ spNum + (1 + spNum | grassType), data = rawData)
m1
broom.mixed::tidy(m1, effects = "fixed")
broom.mixed::tidy(m1, effects = "ran_coefs")

p1<-rawData %>%
  add_predictions(m1) %>%
  ggplot(aes(
    x = spNum, y = NIRveos, group = grassType,
    colour = grassType
  )) +
  geom_point(size=2,alpha=0.4,show.legend = F) +
  theme_bw() + 
  geom_line(aes(x=spNum,y=pred),linewidth=1,show.legend = F)+
  stat_cor(method = "pearson",label.x.npc = 0.65,show.legend = F)+
  theme(axis.text = element_text(colour = "black"),
        axis.line = element_line(size=0.5, colour = 'black'),
        panel.grid=element_blank(),panel.background=element_blank())+
  labs(x = "Diversity", y = "EOS") +
  scale_colour_discrete("grassType")
p1
########
m2 <- lmerTest::lmer(agb ~ spNum + (1 + spNum | grassType), data = rawData)
m2

broom.mixed::tidy(m2, effects = "fixed")
broom.mixed::tidy(m2, effects = "ran_coefs")

p2<-rawData %>%
  add_predictions(m2) %>%
  ggplot(aes(
    x =  spNum, y =agb, group = grassType,
    colour = grassType
  )) +
  geom_point(size=2,alpha=0.4,show.legend = F) +
  theme_bw() +
  geom_line(aes(x=spNum,y=pred),linewidth=1,show.legend = F)+
  stat_cor(method = "pearson",label.x.npc = 0.65,show.legend = F)+
  theme(axis.text = element_text(colour = "black"),
        axis.line = element_line(size=0.5, colour = 'black'),
        panel.grid=element_blank(),panel.background=element_blank())+
  theme(legend.position = "none")+
  labs(x = "Diversity", y = "AGB") +
  scale_colour_discrete("grassType")
p2
####
m3 <- lmerTest::lmer( bgb ~ spNum + (1 + spNum | grassType), data = rawData)
m3
broom.mixed::tidy(m3, effects = "fixed")
broom.mixed::tidy(m3, effects = "ran_coefs")

p3<-rawData %>%
  add_predictions(m3) %>%
  ggplot(aes(
    x = spNum, y = bgb, group = grassType,
    colour = grassType
  )) +
  geom_point(size=2,alpha=0.4,show.legend = F) +
  theme_bw()+ 
  geom_line(aes(x=spNum,y=pred),linewidth=1,show.legend = F)+
  stat_cor(method = "pearson",label.x.npc = 0.65,show.legend = F)+
  theme(axis.text = element_text(colour = "black"),
        axis.line = element_line(size=0.5, colour = 'black'),
        panel.grid=element_blank(),panel.background=element_blank())+
  labs(x = "Diversity", y = "BGB") +
  scale_colour_discrete("grassType")
p3

m4 <- lmerTest::lmer( NIRveos ~ agb + (1 + agb | grassType),data = rawData)
m4
broom.mixed::tidy(m4, effects = "fixed")
broom.mixed::tidy(m4, effects = "ran_coefs")

p4<-rawData %>%
  add_predictions(m4) %>%
  ggplot(aes(
    x =  agb, y =NIRveos, group = grassType,
    colour = grassType
  )) +
  geom_point(size=2,alpha=0.4) +
  theme_bw() + 
  geom_line(aes(x=agb,y=pred),linewidth=1)+
  stat_cor(method = "pearson",label.x.npc = 0.65)+
  theme(axis.text = element_text(colour = "black"),
        axis.line = element_line(size=0.5, colour = 'black'),
        panel.grid=element_blank(),panel.background=element_blank())+
  theme(legend.position = "none")+
  labs(x = "AGB", y = "EOS") +
  scale_colour_discrete("grassType")
p4

m5 <- lmer( NIRveos ~bgb  + (1 + bgb | grassType), data = rawData)
m5
broom.mixed::tidy(m5, effects = "fixed")
broom.mixed::tidy(m5, effects = "ran_coefs")

p5<-rawData %>%
  add_predictions(m5) %>%
  ggplot(aes(
    x = bgb, y = NIRveos, group = grassType,
    colour = grassType
  )) +
  geom_point(size=2,alpha=0.4) +
  geom_line(aes(x=bgb,y=pred),linewidth=1)+
  stat_cor(method = "pearson",label.x.npc = 0.65)+
  theme_bw() + 
  theme(legend.position = c(1.5,0.5))+
  theme(axis.text = element_text(colour = "black"),
        axis.line = element_line(size=0.5, colour = 'black'),
        panel.grid=element_blank(),panel.background=element_blank())+
  labs(x = "BGB", y = "EOS") +
  scale_colour_discrete("grassType")
p5

p <- cowplot::plot_grid(p1, p2, p3, p4, p5, rel_widths = c(1,1),nrow = 3)
p
ggsave("C:\\Users\\Administrator\\Desktop\\草地多样性\\figS2.pdf",p, width=25, height=25, units="cm")

####Figs4 bootstrap####多重比较###data没更改####
library(lme4)
library(ggplot2)
library(agricolae)
library(ggsignif)
figsData <- fread('E:\\草地多样性\\data_pca.csv')
names(figsData)[32:35] <- c('climateP1','climateP2','soilP1','soilP2')
figsData2 <- dplyr::select(figsData,c(site,grassType,lat,lon,spNum,NIRveos,
                                      agb,bgb,climateP1,climateP2,soilP1))
figsDataScale <- apply(figsData2[,5:11],2,scale) %>% cbind(.,figsData2[,1:4]) %>% .[-c(99,39,53),]

######multiCompare##########
multiCompare <- function(data,y,x,title){
  lmerModel <- lmer(get(y) ~ get(x) + (1+ get(x)|grassType), data, REML = FALSE)
  lmerBoot <- bootMer(lmerModel,FUN=function(.){
    slope <- broom.mixed::tidy(.,effects = "ran_coefs")
    slope$estimate[7:12]
  },nsim=500,seed=1234)
  #####整理数据#######
  lmerPlot <- data.table(lmerBoot$t) %>%
    `colnames<-`(.,c('低地盐化草甸类','高寒草甸',
                     '荒漠草原','山地草甸','温性草甸草原','温性草原')) %>%
    melt.data.table(.,id=1:6,measure=1:6) %>%
    dplyr::select(.,-c(1:6)) %>% `colnames<-`(.,c('grassType','value'))
  ####多重比较######
  mod <-  aov(value ~ grassType, data=lmerPlot)
  re <-  LSD.test(mod,"grassType",alpha = 0.05)$groups
  
  ggplot(lmerPlot,aes(grassType,value,color=grassType))+ geom_boxplot()+
    annotate('text',x = 1:6, y = rep(1,6),label=re$groups)+ggtitle(title)
}

EOSdivComapre <- multiCompare(figsDataScale,'NIRveos','spNum','eos-div')
EOSdivComapre

agbdivComapre <- multiCompare(figsDataScale,'agb','spNum','agb-div')
agbdivComapre

bgbdivComapre <- multiCompare(figsDataScale,'bgb','spNum','bgb-div')
bgbdivComapre

agbEOSComapre <- multiCompare(figsDataScale,'agb','NIRveos','agb-eos')
agbEOSComapre

bgbEOSComapre <- multiCompare(figsDataScale,'bgb','NIRveos','bgb-eos')
bgbEOSComapre

#######figS5###########
partialGrass <- df7[,c('site','grassType','lat','lon','spNum','NIRveos','agb',
                       'bgb','pc1_clim','pc2_clim','pc1_soil')]
scaleData <- scale(partialGrass[,c(5:11)]) %>% as.data.table()

climatePartial <- ppcor::pcor.test(scaleData[,'NIRveos'],scaleData[,'spNum'],
                                   scaleData[,c('pc1_clim','pc2_clim')])
soilPartial <- ppcor::pcor.test(scaleData[,'NIRveos'],scaleData[,'spNum'],
                                scaleData[,'pc1_soil'])
climateSoilPartial <- ppcor::pcor.test(scaleData[,'NIRveos'],scaleData[,'spNum'],
                                       scaleData[,c('pc1_clim','pc2_clim','pc1_soil')])
partialData <- data.table(
  name=c('climatePartial','soilPartial','climateSoilPartial'),
  pcor=c(climatePartial$estimate,soilPartial$estimate,climateSoilPartial$estimate)
)

pbar<-ggplot(data=partialData,mapping=aes(x=reorder(name,pcor),y=pcor))+
  geom_bar(stat="identity",position="dodge",fill="#79CDCD")+
  theme_bw() + theme(panel.grid=element_blank())+theme(legend.position = "none")+
  xlab("Envrionment factors")+
  ylab("Pcor")
pbar




#####fig.s7############
library(ggpubr)
rawData <- fread('C:\\Users\\Administrator\\Desktop\\草地多样性\\rawData.csv') 
p2 <- ggplot(rawData, aes( x = spNum, y = agb )) +
  geom_point( size=2,alpha=0.7, color="#8FB3CF")+
  geom_smooth(method = 'lm', formula = y ~ x, se = T,color="#8FB3CF",linewidth=0.5) +
  stat_cor(method = "pearson")+
  theme_bw() + theme(panel.grid=element_blank())+theme(legend.position = "none")+
  xlab("Diversity")+
  ylab("AGB")
p2
ggsave("C:\\Users\\Administrator\\Desktop\\草地多样性\\Fig.s6agb.pdf",p2, width=8, height=7, units="cm")

p2 <- ggplot(rawData, aes( x = spNum, y = bgb )) +
  geom_point( size=2,alpha=0.7, color="#8FB3CF")+
  geom_smooth(method = 'lm', formula = y ~ x, se = T,color="#8FB3CF",linewidth=0.5) +
  stat_cor(method = "pearson")+
  theme_bw() + theme(panel.grid=element_blank())+theme(legend.position = "none")+
  xlab("Diversity")+
  ylab("BGB")
p2
ggsave("C:\\Users\\Administrator\\Desktop\\草地多样性\\Fig.s6bgb.pdf",p2, width=8, height=7, units="cm")


########tableS1##############
df2 <- fread('E:\\草地多样性\\allData.csv') %>%
  mutate(.,tmean=(tmmn+tmmx)*0.1/2,tmean_growing=(tmmn_growing+tmmx_growing)*0.1/2)
df3 <- subset(df2,agb<800) %>% dplyr::select(.,c(site,grassType,lat,lon,spNum,
                                                 soilWeight5,agb,bgb,carbon,n,p,soc,
                                                 NIRveos,pr_growing,srad_growing,
                                                 vpd_growing,tmean_growing))
grassDataClean <- df3[-which(site%in%c('ALT-116','ALT-041','ALT-057')),]
library(pastecs)
stat.desc(grassDataClean, basic=TRUE, desc=TRUE, norm=FALSE, p=0.95)
grassDataClean$grassType[which(grassDataClean$grassType%in%c('温性草原化荒漠类','温性荒漠','温性荒漠草原'))] <- '荒漠草原'
grassDataClean$grassType <- dplyr::case_when(
  grassDataClean$grassType == '低地盐化草甸类'~'Lowland meadows',
  grassDataClean$grassType == '高寒草甸'~'Mountain meadows',
  grassDataClean$grassType == '荒漠草原'~'Desert steppe',
  grassDataClean$grassType == '山地草甸'~'Mountain meadows',
  grassDataClean$grassType == '温性草甸草原'~'Temperate meadow steppe',
  grassDataClean$grassType == '温性草原'~'Temperate steppe'
)
fwrite(subset(grassDataClean,select=-c(carbon,p)),'E:\\草地多样性\\allData.csv')
allSummary <- data.table()
for (type in unique(grassDataClean$grassType)) {
  df <- subset(grassDataClean,grassType==type)
  sta <- stat.desc(df, basic=TRUE, desc=TRUE, norm=FALSE, p=0.95) %>%
    dplyr::mutate(.,type=type,name=rownames(.))
  allSummary <- rbind(allSummary,sta)
}





