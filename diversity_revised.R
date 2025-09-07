
# packages
library(lmerTest); 
library(emmeans); 
library(piecewiseSEM); 
library(broom.mixed); 
library(dplyr)

df <- fread('./df_final.csv')

# Extract standardized coefficients (multiply by sd of predictor / sd of response)
library(parameters)

#Stratified analysis: split by soil moisture median
df <- df %>% mutate(pr_bin = ifelse(pr_growing < median(pr_growing, na.rm=TRUE),'dry','wet'))

m_dry <- lmer(AGB_BGB_0_10 ~ spNum + pr_growing+soilMoist + (1|grassType), data=filter(df, pr_bin=='dry'))
m_wet <- lmer(AGB_BGB_0_10 ~ spNum + pr_growing+soilMoist + (1|grassType), data=filter(df, pr_bin=='wet'))
# summary(m_dry); 
# summary(m_bgb_10)

sum_dry <- summary(m_dry)
sum_wet <- summary(m_wet)


outDf <- data.frame(
  term=c('Intercept','Rich','Precip','SM','AIC','BIC','logLik','N','R2 margin','R2 condition'),
  Modeldry=c(
    paste0(round(sum_dry$coefficients[,1],2),' (',round(sum_dry$coefficients[,2],2),')**'),
    round(AIC(m_dry),2),round(BIC(m_dry),2),round(logLik(m_dry),2),
    nobs(m_dry),round(performance::r2(m_dry)$R2_marginal,2),
    round(performance::r2(m_dry)$R2_conditional,2)
  ),
  ModeldryP=c(
    ifelse(sum_dry$coefficients[,5]<0.01,'<0.01',
           ifelse(sum_dry$coefficients[,5]>0.05,round(sum_dry$coefficients[,5],2),'<0.05')),
    rep('-',6)
  ),
  Modelwet=c(
    paste0(round(sum_wet$coefficients[,1],2),' (',round(sum_wet$coefficients[,2],2),')**'),
    round(AIC(m_wet),2),round(BIC(m_wet),2),round(logLik(m_wet),2),
    nobs(m_wet),round(performance::r2(m_wet)$R2_marginal,2),
    round(performance::r2(m_wet)$R2_conditional,2)
  ),
  ModelwetP=c(
    ifelse(sum_wet$coefficients[,5]<0.01,'<0.01',
           ifelse(sum_wet$coefficients[,5]>0.05,round(sum_wet$coefficients[,5],2),'<0.05')),
    rep('-',6)
    
  )
)


plotDf <- data.frame(
  type=c(rep('dry',3),rep('wet',3)),
  Parameter=rep(c('Richness','Precipitation','Soil moisture'),2),
  slope=c(sum_dry$coefficients[2:4,1],sum_wet$coefficients[2:4,1])
  # sd=c(sum_dry$coefficients[2:4,2],sum_wet$coefficients[2:4,2])
)

ggplot(plotDf, aes(x = Parameter, y = slope, fill = type)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  # geom_errorbar(
  #   aes(ymin = slope - sd, ymax = slope + sd),
  #   position = position_dodge(width = 0.7),  # 使用相同的dodge宽度
  #   width = 0.2,
  #   size = 0.5  # 可选：调整误差线粗细
  # ) +
  scale_fill_manual(values = c("dry" = "#FDA35C", "wet" = "#76C1DF")) +
  labs(x = "", y = "Slope", fill = "") +
  theme_bw()


# Required libraries
library(lmerTest); 
library(emmeans); 
library(ggplot2); 
library(dplyr)
# library(parameters)    # for standardized coefficients
library(knitr)
library(gridExtra)

# full model
# 1) Standardized coefficients barplot
m1 <- lmer(NIRveos ~ spNum + pr_growing + soilMoist + (1|grassType), data=df)
sum_m <- summary(m1)

outDf <- data.frame(
  term=c('Intercept','Rich','Precip','SM','AIC','BIC','logLik','N','R2 margin','R2 condition'),
  Model=c(
    paste0(round(sum_m$coefficients[,1],2),' (',round(sum_m$coefficients[,2],2),')**'),
    round(AIC(m1),2),round(BIC(m1),2),round(logLik(m1),2),
    nobs(m1),round(performance::r2(m1)$R2_marginal,2),
    round(performance::r2(m1)$R2_conditional,2)
  ),
  ModelP=c(
    ifelse(sum_m$coefficients[,5]<0.01,'<0.01',
           ifelse(sum_m$coefficients[,5]>0.05,round(sum_m$coefficients[,5],2),'<0.05')),
    rep('-',6)
  )
)
# stdDF$Label <- c("Richness","Precipitation","Soil moisture")


# ggsave("Fig_std_coefs.png", p1, width=6, height=3.5, dpi=300)


# 2) Interaction plot: predicted NIRveos vs Richness at low/med/high SoilMoist
m2 <- lmer(NIRveos ~ spNum * soilMoist + pr_growing  + (1|grassType), data=df)
m3 <- lmer(NIRveos ~ spNum * pr_growing + soilMoist  + (1|grassType), data=df)


sum_m2 <- summary(m2)
sum_m3 <- summary(m3)

out_m2 <- data.frame(
  term=c(rownames(sum_m2$coefficients),c('AIC','BIC','logLik','N','R2 margin','R2 condition')),
  Model=c(
    paste0(round(sum_m2$coefficients[,1],2),' (',round(sum_m2$coefficients[,2],2),')**'),
    round(AIC(m2),2),round(BIC(m2),2),round(logLik(m2),2),
    nobs(m2),round(performance::r2(m2)$R2_marginal,2),
    round(performance::r2(m2)$R2_conditional,2)
  ),
  ModelP=c(
    ifelse(sum_m2$coefficients[,5]<0.01,'<0.01',
           ifelse(sum_m2$coefficients[,5]>0.05,round(sum_m2$coefficients[,5],2),'<0.05')),
    rep('-',6)
  )
)

out_m3 <- data.frame(
  term=c(rownames(sum_m3$coefficients),c('AIC','BIC','logLik','N','R2 margin','R2 condition')),
  Model=c(
    paste0(round(sum_m3$coefficients[,1],2),' (',round(sum_m3$coefficients[,2],2),')**'),
    round(AIC(m3),2),round(BIC(m3),2),round(logLik(m3),2),
    nobs(m3),round(performance::r2(m3)$R2_marginal,2),
    round(performance::r2(m3)$R2_conditional,2)
  ),
  ModelP=c(
    ifelse(sum_m3$coefficients[,5]<0.01,'<0.01',
           ifelse(sum_m3$coefficients[,5]>0.05,round(sum_m3$coefficients[,5],2),'<0.05')),
    rep('-',6)
  )
)


plotDf <- data.frame(
  type=c(rep('Full model',3),rep('Rich:SM',3),rep('Rich:Precip',3)),
  Parameter=rep(c('Richness','Precipitation','Soil moisture'),3),
  slope=c(sum_m$coefficients[2:4,1],sum_m2$coefficients[c(2,4,3),1],sum_m3$coefficients[2:4,1])
  # sd=c(sum_dry$coefficients[2:4,2],sum_wet$coefficients[2:4,2])
)

ggplot(plotDf, aes(x = Parameter, y = slope, fill = type)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("Full model" = "#76C1DF", "Rich:Precip" = "#FDA35C","Rich:SM" = "#9C0824")) +
  labs(x = "", y = "Slope", fill = "") +
  theme_bw()



# 3) Stratified slopes: dry vs bgb_10 (use your split)
# df <- df %>% mutate(sm_bin = ifelse(soilMoist < median(soilMoist, na.rm=TRUE), "dry","bgb_10"))
m_dry <- lm(NIRveos ~ spNum + pr_growing+ soilMoist , data=filter(df, pr_bin=="dry"))
m_wet <- lm(NIRveos ~ spNum + pr_growing +soilMoist, data=filter(df, pr_bin=="wet"))



sum_dry <- summary(m_dry)
sum_wet <- summary(m_wet)


outDf <- data.frame(
  term=c('Intercept','Rich','Precip','SM','AIC','BIC','logLik','N','R2 adjust'),
  Modeldry=c(
    paste0(round(sum_dry$coefficients[,1],2),' (',round(sum_dry$coefficients[,2],2),')**'),
    round(AIC(m_dry),2),round(BIC(m_dry),2),round(logLik(m_dry),2),
    nobs(m_dry),round(performance::r2(m_dry)$R2_adjusted,2)
  ),
  ModeldryP=c(
    ifelse(sum_dry$coefficients[,4]<0.01,'<0.01',
           ifelse(sum_dry$coefficients[,4]>0.05,round(sum_dry$coefficients[,4],2),'<0.05')),
    rep('-',5)
  ),
  Modelwet=c(
    paste0(round(sum_wet$coefficients[,1],2),' (',round(sum_wet$coefficients[,2],2),')**'),
    round(AIC(m_wet),2),round(BIC(m_wet),2),round(logLik(m_wet),2),
    nobs(m_wet),round(performance::r2(m_dry)$R2_adjusted,2)
  ),
  ModelwetP=c(
    ifelse(sum_wet$coefficients[,4]<0.01,'<0.01',
           ifelse(sum_wet$coefficients[,4]>0.05,round(sum_wet$coefficients[,4],2),'<0.05')),
    rep('-',5)
    
  )
)


plotDf <- data.frame(
  type=c(rep('dry',3),rep('wet',3)),
  Parameter=rep(c('Richness','Precipitation','Soil moisture'),2),
  slope=c(sum_dry$coefficients[2:4,1],sum_wet$coefficients[2:4,1])
  # sd=c(sum_dry$coefficients[2:4,2],sum_wet$coefficients[2:4,2])
)

ggplot(plotDf, aes(x = Parameter, y = slope, fill = type)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  # geom_errorbar(
  #   aes(ymin = slope - sd, ymax = slope + sd),
  #   position = position_dodge(width = 0.7),  # 使用相同的dodge宽度
  #   width = 0.2,
  #   size = 0.5  # 可选：调整误差线粗细
  # ) +
  scale_fill_manual(values = c("dry" = "#FDA35C", "wet" = "#76C1DF")) +
  labs(x = "", y = "Slope", fill = "") +
  theme_bw()



df_stack <- df %>%
  select(site, BGB_0_10, BGB_10_20, BGB_20_30) %>% na.omit() %>% 
  tidyr::pivot_longer(cols = starts_with("BGB_"), names_to = "depth", values_to = "root_g_m2") %>%
  mutate(depth = case_when(
    depth == "BGB_0_10" ~ "0-10 cm",
    depth == "BGB_10_20" ~ "10-20 cm",
    depth == "BGB_20_30" ~ "20-30 cm",
    TRUE ~ depth
  ))

# order sites by total BGB descending
site_order <- df %>% arrange(desc(BGB_0_30)) %>% pull(site)
# site_order2 <- str_split_i(site_order,'-',2) 

# df_stack$site2 <- str_split_i(df_stack$site,'-',2)
df_stack$site <- factor(df_stack$site, levels = site_order)

ggplot(df_stack, aes(x = site, y = root_g_m2, fill = depth)) +
  geom_col(position = "stack") +
  theme_bw() +
  scale_fill_manual(values = c("#F7C45E",'#F3892D','#9E3A26'))+
  scale_x_discrete(labels = df_stack$site_label)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = 'bottom') +
  labs(x = "Site (ordered by total BGB)", y = expression(Root~biomass~(g~m^{-2})), fill = "Depth")

df$prop_top10 <- df$BGB_0_10/df$BGB_0_30

ggplot(df, aes(y = prop_top10)) +
  geom_boxplot() +
  labs(x = "", y = "Proportion") +
  theme_bw()


# df <- df %>% mutate(eos_group = case_when(
#   !is.na(NIRveos) & is.numeric(NIRveos) ~ ntile(NIRveos, 3) %>% as.character(),
#   TRUE ~ as.character(NIRveos)
# ))





#compare ratio
compareDf <- df[df$BGB_0_30>0,]
formula <- y ~ x
ggplot(compareDf, aes(x = AGB_BGB_0_30, y = AGB_BGB_0_10)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  ggpmisc::stat_poly_eq(
    aes(label = paste(after_stat(eq.label), 
                      after_stat(rr.label), 
                      after_stat(p.value.label), 
                      sep = "~~~")),
    formula = formula,
    parse = TRUE,
    label.x = "right", 
    label.y = "bottom" 
  ) +
  labs(x = "BGB:AGB (0-30 cm)", y = "BGB:AGB (0-10 cm)") +
  theme_bw()


longDf <- compareDf %>%
  tidyr::pivot_longer(
    cols = c(AGB_BGB_0_10, AGB_BGB_0_30),
    names_to = "depth",
    values_to = "AGB_BGB"
  ) %>%
  mutate(
    depth = factor(depth, levels = c("AGB_BGB_0_10", "AGB_BGB_0_30"),
                   labels = c("0_10", "0_30"))
  )


mod_full <- lmer(AGB_BGB ~ spNum * depth + (1 | grassType), data = longDf, REML = FALSE)

mod_noInt <- lmer(AGB_BGB ~ spNum + depth + (1 | grassType), data = longDf, REML = FALSE)


sum_full <- summary(mod_full)
sum_noInt <- summary(mod_noInt)

out_full <- data.frame(
  term=c('Intercept','Rich','depth0_30','spNum:depth0_30','AIC','BIC','logLik','N','R2 marginal','R2 condition'),
  Modelfull=c(
    paste0(round(sum_full$coefficients[,1],2),' (',round(sum_full$coefficients[,2],2),')**'),
    round(AIC(mod_full),2),round(BIC(mod_full),2),round(logLik(mod_full),2),
    nobs(mod_full),round(performance::r2(mod_full)$R2_marginal,2),
    round(performance::r2(mod_full)$R2_conditional,2)
  ),
  ModelfullP=c(
    ifelse(sum_full$coefficients[,5]<0.01,'<0.01',
           ifelse(sum_full$coefficients[,5]>0.05,round(sum_full$coefficients[,5],2),'<0.05')),
    rep('-',6)
  )
)

out_noInt <- data.frame(
  term=c('Intercept','Rich','depth0_30','AIC','BIC','logLik','N','R2 marginal','R2 condition'),
  ModelnoInt=c(
    paste0(round(sum_noInt$coefficients[,1],2),' (',round(sum_noInt$coefficients[,2],2),')**'),
    round(AIC(mod_noInt),2),round(BIC(mod_noInt),2),round(logLik(mod_noInt),2),
    nobs(mod_noInt),round(performance::r2(mod_noInt)$R2_marginal,2),
    round(performance::r2(mod_noInt)$R2_conditional,2)
  ),
  ModelnoIntP=c(
    ifelse(sum_noInt$coefficients[,5]<0.01,'<0.01',
           ifelse(sum_noInt$coefficients[,5]>0.05,round(sum_noInt$coefficients[,5],2),'<0.05')),
    rep('-',6)
  )
)




# log-transform ratios to stabilize variance; remove zeros/negatives
mod_df <- df %>% filter(is.finite(AGB_BGB_0_30) & AGB_BGB_0_30 > 0 & is.finite(AGB_BGB_0_10) & AGB_BGB_0_10 > 0)
# mod_df2 <- subset(mod_df,mod_df$AGB_BGB_0_30<100)
# rawData$AGB_BGB_0_30 <- rawData$bgb/rawData$agb

m_bgb_30 <- lmer(AGB_BGB_0_30 ~ spNum+(1|grassType), data = mod_df)
m_bgb_10 <- lm(AGB_BGB_0_10 ~ spNum, data = mod_df)
# m_bgb_10 <- lm(log(AGB_BGB_0_10) ~ spNum, data = mod_df)



sum_bgb_30 <- summary(m_bgb_30)
sum_bgb_10 <- summary(m_bgb_10)


outDf <- data.frame(
  term=c('Intercept','Rich','Precip','SM','AIC','BIC','logLik','N','R2 margin','R2 condition'),
  Modelbgb_30=c(
    paste0(round(sum_bgb_30$coefficients[,1],2),' (',round(sum_bgb_30$coefficients[,2],2),')**'),
    round(AIC(m_bgb_30),2),round(BIC(m_bgb_30),2),round(logLik(m_bgb_30),2),
    nobs(m_bgb_30),round(performance::r2(m_bgb_30)$R2_marginal,2),
    round(performance::r2(m_bgb_30)$R2_conditional,2)
  ),
  Modelbgb_30P=c(
    ifelse(sum_bgb_30$coefficients[,5]<0.01,'<0.01',
           ifelse(sum_bgb_30$coefficients[,5]>0.05,round(sum_bgb_30$coefficients[,5],2),'<0.05')),
    rep('-',6)
  ),
  Modelbgb_10=c(
    paste0(round(sum_bgb_10$coefficients[,1],2),' (',round(sum_bgb_10$coefficients[,2],2),')**'),
    round(AIC(m_bgb_10),2),round(BIC(m_bgb_10),2),round(logLik(m_bgb_10),2),
    nobs(m_bgb_10),round(performance::r2(m_bgb_10)$R2_marginal,2),
    round(performance::r2(m_bgb_10)$R2_conditional,2)
  ),
  Modelbgb_10P=c(
    ifelse(sum_bgb_10$coefficients[,5]<0.01,'<0.01',
           ifelse(sum_bgb_10$coefficients[,5]>0.05,round(sum_bgb_10$coefficients[,5],2),'<0.05')),
    rep('-',6)
    
  )
)


plotDf <- data.frame(
  type=c(rep('bgb_30',3),rep('bgb_10',3)),
  Parameter=rep(c('Richness','Precipitation','Soil moisture'),2),
  slope=c(sum_bgb_30$coefficients[2:4,1],sum_bgb_10$coefficients[2:4,1])
  # sd=c(sum_bgb_30$coefficients[2:4,2],sum_bgb_10$coefficients[2:4,2])
)

ggplot(plotDf, aes(x = Parameter, y = slope, fill = type)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("bgb_30" = "#FDA35C", "bgb_10" = "#76C1DF")) +
  labs(x = "", y = "Slope", fill = "") +
  theme_bw()





library(effects)
library(jtools)
library(grfxtools)
library(cowplot)

#####fig.2a########
df <- fread('./df_final.csv')

model_eos_div <- lmer(NIRveos~spNum+(1|grassType),data=df)
rsquared(model_eos_div)$Conditional %>% round(2)
out <- summary(model_eos_div)
trend_unscaled <- out$coefficients["spNum","Estimate"]/ sd(df$spNum)
error_unscaled <- out$coefficients["spNum","Std. Error"]/ sd(df$spNum)
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
      legend.text = element_text(size=8),
      axis.text = element_text(size=8),
      axis.title = element_text(size=12),
    )+
    xlab("Diversity")+
    ylab("EOS (DOY)")
  return(gg)
}

lab1 <- "paste(\" EOS ~ Diversity \",italic(p), \" < 0.01**\")"

ggplot_eos_div(out_model_eos_div)+
  annotate("text", x = 12, y = 330, label = lab1,parse = TRUE,
           size = 4)

######fig 2c######
model_eos_div <- lmer(AGB_BGB_0_10~spNum+(1|grassType),data=df)
rsquared(model_eos_div)$Conditional %>% round(2)
out <- summary(model_eos_div)
trend_unscaled <- out$coefficients["spNum","Estimate"]/ sd(df$spNum)
error_unscaled <- out$coefficients["spNum","Std. Error"]/ sd(df$spNum)
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
    geom_hex(data = parEOSDiv, aes(spNum,AGB_BGB_0_10)) +
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
      legend.text = element_text(size=8),
      axis.text = element_text(size=8),
      axis.title = element_text(size=12),
    )+
    xlab("Diversity")+
    ylab("Ratio of BGB:AGB")
  return(gg)
}

lab1 <- "paste(\" Ratio of BGB:AGB ~ Diversity \",italic(p), \" < 0.01**\")"

ggplot_eos_div(out_model_eos_div)+
  annotate("text", x = 10, y = 140, label = lab1,parse = TRUE,
           size = 4)

#######fig2e############
model_eos_div <- lmer(NIRveos~AGB_BGB_0_10+spNum+(1|grassType),data=df)
rsquared(model_eos_div)$Conditional %>% round(2)
out <- summary(model_eos_div)
trend_unscaled <- out$coefficients["AGB_BGB_0_10","Estimate"]/ sd(df$spNum)
error_unscaled <- out$coefficients["AGB_BGB_0_10","Std. Error"]/ sd(df$spNum)
trend_unscaled
error_unscaled

parEOSDiv <- partialize(model_eos_div,"AGB_BGB_0_10")
out_model_eos_div <- allEffects(model_eos_div,partial.residuals = TRUE)
# windowsFonts(Times_New_Roman=windowsFont("Times New Roman"))
ggplot_eos_div <- function(x){
  df <- tibble(upper = x$`AGB_BGB_0_10`$upper[,1],
               lower = x$`AGB_BGB_0_10`$lower[,1],
               off   = x$`AGB_BGB_0_10`$fit[,1],
               AGB_BGB_0_10 = x$`AGB_BGB_0_10`$x[,1])
  gg <- ggplot() +
    geom_hex(data = parEOSDiv, aes(AGB_BGB_0_10,NIRveos)) +
    scale_fill_gradientn("",colours = alpha(colorRampPalette( c("#2878b5","#ff8884","#f3d266"))(5),0.7),
                         trans = "log", limits=c(1,6),breaks=c(1,2,4,6)) +
    geom_ribbon(data = df, aes(x = AGB_BGB_0_10, ymin = lower, ymax = upper), alpha = 0.3, fill = "#2878b5") +
    geom_line(data = df, aes(AGB_BGB_0_10, off), col = "#2878b5",linewidth=1) +
    theme_bw() + 
    theme(
      panel.grid=element_blank(),
      legend.position = c(0.9,0.8),
      legend.key.size = unit(0.3,"cm"),
      legend.background = element_rect('transparent'),
      legend.text = element_text(size=8),
      axis.text = element_text(size=8),
      axis.title = element_text(size=12),
    )+
    xlab("Ratio of BGB:AGB")+
    ylab("EOS")
  return(gg)
}

lab1 <- "paste(\" EOS  ~ Ratio of BGB:AGB \",italic(p), \" < 0.01**\")"

ggplot_eos_div(out_model_eos_div)+
  annotate("text", x = 100, y = 340, label = lab1,parse = TRUE,
           size = 4)



######fig.2b########


grassCleanData <- rast('./grassClean.tif')
polarPlot <- function(polarData,x,y,cor){
  ggpolar(pole = "N", max.lat = 80, min.lat = 25,
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
    geom_tile(data = polarData,aes_string(x=x,y=y,fill=cor))+
    # geom_point(data=polarData, aes_string(x = x, y = y, color = cor), size = 1)+
    scale_fill_gradientn(colors=paletteer_c("ggthemes::Orange-Blue-White Diverging", 4),limits=c(-1,1))+
    theme(
      # axis.text.y = element_blank(),
      # axis.ticks=element_blank(),
      # plot.title = element_text(hjust = 0.5, vjust = -60, size = 5),
      legend.position = 'right',
      legend.direction="vertical",
      legend.key.size = unit(0.5,"cm"),
      legend.text=element_text(size=8),
      legend.title=element_text(size=10),
      plot.margin = margin(t = 0.1,r = 0.3, b = 0.1, l = 0.1, unit = 'inch')
    )
  
}

eos_diversity_tif <- grassCleanData[[c('NIRveos','spNum')]] %>% terra::scale()
fcor <- focalReg(eos_diversity_tif,9,"ols",na.rm=T) %>% fortify() %>% na.omit()

names(fcor)[4] <- 'slope'

stat <- rbind(data.table(type='positive',val=sum(fcor$slope>0)*100/nrow(fcor)),
              data.table(type='negtative',val=sum(fcor$slope<=0)*100/nrow(fcor)))

polarPlot(fcor,'x','y','slope')

ggplot(data = stat)+
  geom_bar(mapping = aes(x =type, y = val, fill = type),
           stat = 'identity', width = 0.5,position = 'dodge2',show.legend = F)+
  scale_y_continuous(breaks = seq(0,60,30), labels = paste0(seq(0,60,30),'%'))+
  #scale_fill_viridis_d(direction = -1,option = "H")+
  scale_x_discrete(labels = c('Negative','Positve'))+
  scale_fill_manual(values=c("#9E3D22","#2B5C8A"))+
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 10, color = 'black'),
        axis.line = element_line(linetype = 1, size =0.5),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        panel.grid = element_blank()
        # axis.ticks.length.x = unit(-0.15, "cm")
  )

######fig.2d########
eos_bgb_ratio <- grassCleanData$bgb/grassCleanData$agb
eos_diversity_tif <- c(eos_bgb_ratio,grassCleanData$spNum) %>% terra::scale()
fcor <- focalReg(eos_diversity_tif,9,"ols",na.rm=T) %>% fortify() %>% na.omit()

names(fcor)[4] <- 'slope'

stat <- rbind(data.table(type='positive',val=sum(fcor$slope>0)*100/nrow(fcor)),
              data.table(type='negtative',val=sum(fcor$slope<=0)*100/nrow(fcor)))

polarPlot(fcor,'x','y','slope')

ggplot(data = stat)+
  geom_bar(mapping = aes(x =type, y = val, fill = type),
           stat = 'identity', width = 0.5,position = 'dodge2',show.legend = F)+
  scale_y_continuous(breaks = seq(0,60,20), labels = paste0(seq(0,60,20),'%'))+
  #scale_fill_viridis_d(direction = -1,option = "H")+
  scale_x_discrete(labels = c('Negative','Positve'))+
  scale_fill_manual(values=c("#9E3D22","#2B5C8A"))+
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 10, color = 'black'),
        axis.line = element_line(linetype = 1, size =0.5),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        panel.grid = element_blank()
        # axis.ticks.length.x = unit(-0.15, "cm")
  )
######fig.2f########
eos_bgb_ratio <- grassCleanData$bgb/grassCleanData$agb
eos_diversity_tif <- c(grassCleanData[['NDVIeos']],eos_bgb_ratio) %>% terra::scale()
fcor <- focalReg(eos_diversity_tif,9,"ols",na.rm=T) %>% fortify() %>% na.omit()

names(fcor)[4] <- 'slope'

stat <- rbind(data.table(type='positive',val=sum(fcor$slope>0)*100/nrow(fcor)),
              data.table(type='negtative',val=sum(fcor$slope<=0)*100/nrow(fcor)))

polarPlot(fcor,'x','y','slope')

ggplot(data = stat)+
  geom_bar(mapping = aes(x =type, y = val, fill = type),
           stat = 'identity', width = 0.5,position = 'dodge2',show.legend = F)+
  scale_y_continuous(breaks = seq(0,60,30), labels = paste0(seq(0,60,30),'%'))+
  #scale_fill_viridis_d(direction = -1,option = "H")+
  scale_x_discrete(labels = c('Negative','Positve'))+
  scale_fill_manual(values=c("#9E3D22","#2B5C8A"))+
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 10, color = 'black'),
        axis.line = element_line(linetype = 1, size =0.5),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        panel.grid = element_blank()
        # axis.ticks.length.x = unit(-0.15, "cm")
  )
#########fig3a################
df <- fread('./df_final.csv')
df_forest <- df[,c('NIRveos','AGB_BGB_0_10','spNum','soilMoist','soc','n','pr_growing','soilBD','vpd_growing','srad_growing','tmean_growing')] %>% 
  apply(.,2,scales::rescale)  %>% as.data.table() %>% cbind(.,grassType=df$grassType)

colnames(df_forest) <- c('EOS','Ratio','Diversity','SM','SOC','Soil_N','Precipitation','BD','VPD','SRAD','Temperature','grassType')
lme_SMF<-lm(EOS ~ Ratio +Diversity+SM+SOC+Soil_N+Precipitation+BD+VPD+SRAD+Temperature,data=df_forest)
lmeResult <- summary(lme_SMF)
result <- lmeResult$coefficients %>% as.data.table()
result$index <- rownames(lmeResult$coefficients)
result <- subset(result,result$index!='(Intercept)')

df_forest$biomass <- df_forest$Ratio
df_forest$climate <- df_forest$Precipitation+df_forest$SRAD+df_forest$VPD+df_forest$Temperature
df_forest$soil <- df_forest$Soil_N+df_forest$SOC+df_forest$SM+df_forest$BD

glmLme<-lm(EOS~soil+climate+biomass+Diversity,data=df_forest)
lmerCoef <- summary(glmLme)$coef
allSum <- sum(abs(lmerCoef[2:5,1]))
soil <- abs(lmerCoef[2,1])/allSum*100
climate <- abs(lmerCoef[3,1])/allSum*100
biomass <- abs(lmerCoef[4,1])/allSum*100
diversity <- abs(lmerCoef[5,1])/allSum*100
# glmmmResults <- glmm.hp(glmLme)


result[which(result$index %in% c("Ratio")),'type'] <- 'Biomass'
result[which(result$index %in% c( "Precipitation","SRAD","VPD",'Temperature')),'type'] <- 'Climate'
result[which(result$index %in% c( "Soil_N","BD",'SOC','SM')),'type'] <- 'Soil'
result[which(result$index %in% c( "Diversity")),'type'] <- 'Diversity'

result <- result[order(result$type, result$ Estimate), ]
result$index <- factor(result$index, levels = c('Diversity','Ratio','SM','SOC','Soil_N','BD','Precipitation','VPD','SRAD','Temperature'))
result$type <- factor(result$type, levels = c('Diversity','Biomass','Soil','Climate'))

# windowsFonts(A=windowsFont("Times New Roman"),
#              B=windowsFont("Times New Roman"))
p1 <- ggplot(result, aes(x = index, y = Estimate, color = type)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Estimate - `Std. Error`, ymax = Estimate + `Std. Error`), width = 0.2, size =1,show.legend = F) +
  scale_color_manual(values = c("#ff8884","#f3d266","#91D1C2FF","#2878b5")) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = 2,size=1) +
  labs(x = '', y = 'Estimate', color = '') +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line.x = element_line(), axis.ticks.y = element_blank(),
        legend.key = element_blank()) +
  scale_x_discrete(breaks = result$index, labels = result$index, position = 'top') +
  scale_y_continuous(limits = c(-0.7, 0.7))+
  theme(text=element_text(size=13))

p1
lfly <- data.frame(
  fre =c(climate,soil,biomass,diversity) , 
  diet = c(''),
  beh = c('Climate','Soil','Biomass','Diversity'))
lfly$beh <- factor(lfly$beh,levels = lfly$beh)
p2 <- ggplot(lfly, aes(diet, fre, fill = beh))+
  geom_bar(stat = 'identity', position = 'fill', width = 0.3)+
  geom_col(width = 0.9) +
  labs(x = bquote(R^2==0.37), y = 'Relative estimates (%)')+
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.ticks.x = element_blank(), axis.line.y = element_line(color = 'gray30'))+
  scale_fill_manual(values = c("#2878b5","#91D1C2FF","#f3d266","#ff8884"))+
  scale_x_discrete(position = 'top') +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0))+
  theme(text=element_text(size=13))

p2
library(patchwork)
p2 <- p2 + theme(legend.position = 'none')
p2 + p1 + plot_layout(ncol = 2, widths = c(1, 2))


####Fig.3b##########
library(tidyverse)
library(dplyr)

df_forest <- fread('./df_forest.csv')

result <- data.frame()
varNames1 <- lapply(c('Ratio','SM','SOC','Soil_N','BD'),FUN = function(x){
  return(c('Diversity',x))
}) 

varNames2 <- lapply(c('Diversity','Ratio','SM','SOC','Soil_N','BD'),FUN = function(x){
  return(c('EOS',x))
}) 

varNames3 <- lapply(c('SM','SOC','Soil_N','BD'),FUN = function(x){
  return(c('Ratio',x))
}) 

varNames4 <- c(varNames1,varNames2,varNames3) 

climateNames <- c('Precipitation','VPD','Temperature','SRAD')
pcorResult <- data.table()
for (i in 1:length(varNames4)) {
  corName1 <- df_forest[[varNames4[[i]][1]]]
  corName2 <- df_forest[[varNames4[[i]][2]]]
  corName3 <- df_forest[,..climateNames]
  
  pcorDf <- ppcor::pcor.test(corName1,corName2,corName3)
  pcroDf2 <- data.table(
    corNames1=varNames4[[i]][1],
    corNames2=varNames4[[i]][2],
    pcorDf
  )
  pcorResult <- rbind(pcorResult,pcroDf2)
  print(i)
}

pcorResult$pvalue <- ifelse(pcorResult$p.value<0.01,'<0.01',
                            ifelse(pcorResult$p.value>0.05,round(pcorResult$p.value,2),'<0.05'))


# psem_lm <- psem(
psem_lm2 <- psem(
  # lm(EOS~Diversity,data=df_forest),
  lm(EOS~Diversity+BD+Ratio+SM+SOC+Soil_N,data=df_forest),
  lm(Ratio~Diversity+BD+SOC+SM+Soil_N,data=df_forest),
  lm(BD~Diversity,data=df_forest),
  lm(SOC~Diversity,data=df_forest),
  lm(Soil_N~Diversity,data=df_forest),
  lm(SM~Diversity,data=df_forest)
)
sum_mod <- summary(psem_lm2)
plot(psem_lm2)



#########sup###############
#####fig.s1###########
library(ggpmisc)
df <- fread('./df_final.csv')

plotDf <- df[,{sum_mod <- summary(lm(NIRveos~spNum));
.(slope=sum_mod$coef[2,1],
  p =ifelse(sum_mod$coef[2,4]<0.01,'p<0.01',ifelse(sum_mod$coef[2,4]>0.01,sum_mod$coef[2,4],'p<0.05')))
},.(grassType)]


mod_full <- lm(NIRveos ~ spNum * grassType , data = df)

mod_noInt <- lm(NIRveos ~ spNum + grassType, data = df)


sum_full <- summary(mod_full)
sum_noInt <- summary(mod_noInt)

out_full <- data.frame(
  term=c(rownames(sum_full$coefficients),c('AIC','BIC','logLik','N','R2 adjust')),
  Modelfull=c(
    paste0(round(sum_full$coefficients[,1],2),' (',round(sum_full$coefficients[,2],2),')**'),
    round(AIC(mod_full),2),round(BIC(mod_full),2),round(logLik(mod_full),2),
    nobs(mod_full),round(performance::r2(mod_full)$R2_adjusted,2)
  ),
  ModelfullP=c(
    ifelse(sum_full$coefficients[,4]<0.01,'<0.01',
           ifelse(sum_full$coefficients[,4]>0.05,round(sum_full$coefficients[,4],2),'<0.05')),
    rep('-',5)
  )
)

out_noInt <- data.frame(
  term=c(rownames(sum_noInt$coefficients),c('AIC','BIC','logLik','N','R2 adjust')),
  ModelnoInt=c(
    paste0(round(sum_noInt$coefficients[,1],2),' (',round(sum_noInt$coefficients[,2],2),')**'),
    round(AIC(mod_noInt),2),round(BIC(mod_noInt),2),round(logLik(mod_noInt),2),
    nobs(mod_noInt),round(performance::r2(mod_full)$R2_adjusted,2)
  ),
  ModelnoIntP=c(
    ifelse(sum_noInt$coefficients[,4]<0.01,'<0.01',
           ifelse(sum_noInt$coefficients[,4]>0.05,round(sum_noInt$coefficients[,4],2),'<0.05')),
    rep('-',5)
  )
)

#######fig.s2richness 和eneness############
df_div <- fread('./df_div.csv')
# df_div <- subset(df_div,df_div$evenness_mean>0)
df_div2 <- df_div[,-c('AGB_BGB_0_10','NIRveos','grassType')]
colnames(df_div2) <- c('richness','shannon_index','simpson_index','evenness')

library(GGally)  

ggpairs(
  df_div2,
  upper = list(continuous = wrap("cor", size = 4)),
  lower = list(continuous = wrap("smooth", alpha = 0.4, se = FALSE)),
  diag  = list(continuous = "densityDiag")
)


## richness-only
m_rich <- lm(AGB_BGB_0_10 ~ spNum, data = df_div)
# evenness-only
m_even <- lm(AGB_BGB_0_10 ~ shannon_mean, data = df_div)


anova(m_rich, m_even)


sum_rich <- summary(m_rich)
sum_even <- summary(m_even)


outDf <- data.frame(
  term=c('Intercept','Rich','AIC','BIC','logLik','N','R2 adjust'),
  Modelrich=c(
    paste0(round(sum_rich$coefficients[,1],2),' (',round(sum_rich$coefficients[,2],2),')**'),
    round(AIC(m_rich),2),round(BIC(m_rich),2),round(logLik(m_rich),2),
    nobs(m_rich),round(performance::r2(m_rich)$R2_adjusted,2)
  ),
  ModelrichP=c(
    ifelse(sum_rich$coefficients[,4]<0.01,'<0.01',
           ifelse(sum_rich$coefficients[,4]>0.05,
                  round(sum_rich$coefficients[,4],2),'<0.05')),
    rep('-',5)
  ),
  Modeleven=c(
    paste0(round(sum_even$coefficients[,1],2),' (',round(sum_even$coefficients[,2],2),')**'),
    round(AIC(m_even),2),round(BIC(m_even),2),round(logLik(m_even),2),
    nobs(m_even),round(performance::r2(m_even)$R2_adjusted,2)
  ),
  ModelevenP=c(
    ifelse(sum_even$coefficients[,4]<0.01,'<0.01',
           ifelse(sum_even$coefficients[,4]>0.05,
                  round(sum_even$coefficients[,4],2),'<0.05')),
    rep('-',5)
    
  )
)




## richness-only
m_rich <- lm(NIRveos ~ spNum, data = df_div)
# evenness-only
m_even <- lm(NIRveos ~ shannon_mean, data = df_div)


sum_rich <- summary(m_rich)
sum_even <- summary(m_even)


outDf <- data.frame(
  term=c('Intercept','Rich','AIC','BIC','logLik','N','R2 adjust'),
  Modelrich=c(
    paste0(round(sum_rich$coefficients[,1],2),' (',round(sum_rich$coefficients[,2],2),')**'),
    round(AIC(m_rich),2),round(BIC(m_rich),2),round(logLik(m_rich),2),
    nobs(m_rich),round(performance::r2(m_rich)$R2_adjusted,2)
  ),
  ModelrichP=c(
    ifelse(sum_rich$coefficients[,4]<0.01,'<0.01',
           ifelse(sum_rich$coefficients[,4]>0.05,
                  round(sum_rich$coefficients[,4],2),'<0.05')),
    rep('-',5)
  ),
  Modeleven=c(
    paste0(round(sum_even$coefficients[,1],2),' (',round(sum_even$coefficients[,2],2),')**'),
    round(AIC(m_even),2),round(BIC(m_even),2),round(logLik(m_even),2),
    nobs(m_even),round(performance::r2(m_even)$R2_adjusted,2)
  ),
  ModelevenP=c(
    ifelse(sum_even$coefficients[,4]<0.01,'<0.01',
           ifelse(sum_even$coefficients[,4]>0.05,
                  round(sum_even$coefficients[,4],2),'<0.05')),
    rep('-',5)
    
  )
)


#######fig.s6############
df <- fread('./df_forest.csv')

partialGrass <- df[,c('Diversity','EOS','Ratio',
                      'pc1_clim','pc2_clim','pc1_soil','pc2_soil')]
scaleData <- scale(partialGrass) %>% as.data.table()

#eos-diversity
climatePartial <- ppcor::pcor.test(scaleData[,'EOS'],scaleData[,'Diversity'],
                                   scaleData[,c('pc1_clim','pc2_clim')])
soilPartial <- ppcor::pcor.test(scaleData[,'EOS'],scaleData[,'Diversity'],
                                scaleData[,c('pc1_soil','pc2_soil')])
climateSoilPartial <- ppcor::pcor.test(scaleData[,'EOS'],scaleData[,'Diversity'],
                                       scaleData[,c('pc1_clim','pc2_clim','pc1_soil','pc2_soil')])
partialData <- data.table(
  name=c('Model 1','Model 2','Model 3'),
  pcor=c(climatePartial$estimate,soilPartial$estimate,climateSoilPartial$estimate)
)

ggplot(data=partialData,mapping=aes(x=reorder(name,pcor),y=pcor))+
  geom_bar(width=0.5,stat="identity",position="dodge",fill="#79CDCD")+
  theme_bw() + 
  theme(panel.grid=element_blank())+theme(legend.position = "none")+
  xlab("")+
  ylab("Cofficient")

#eos-Ratio
climatePartial <- ppcor::pcor.test(scaleData[,'EOS'],scaleData[,'Ratio'],
                                   scaleData[,c('pc1_clim','pc2_clim')])
soilPartial <- ppcor::pcor.test(scaleData[,'EOS'],scaleData[,'Ratio'],
                                scaleData[,c('pc1_soil','pc2_soil')])
climateSoilPartial <- ppcor::pcor.test(scaleData[,'EOS'],scaleData[,'Ratio'],
                                       scaleData[,c('pc1_clim','pc2_clim','pc1_soil','pc2_soil')])
partialData <- data.table(
  name=c('Model 1','Model 2','Model 3'),
  pcor=c(climatePartial$estimate,soilPartial$estimate,climateSoilPartial$estimate)
)
partialData$name <- factor(partialData$name,levels = c('Model 1','Model 2','Model 3'))
ggplot(data=partialData,mapping=aes(x=name,y=pcor))+
  geom_bar(width=0.5,stat="identity",position="dodge",fill="#79CDCD")+
  theme_bw() + 
  theme(panel.grid=element_blank())+theme(legend.position = "none")+
  xlab("")+
  ylab("Cofficient")

