library(lme4)
library(car)
library(emmeans)
library(Rmisc)
library(tidyverse)

df<-read.csv('../../output/bursts_c.csv')
df$new_subject<-as.character(df$subject)
subjs_12m<-unique(df$subject[df$age=='12m'])
subj_idx<-100
for(subj in subjs_12m) {
  df$new_subject[df$subject==subj & df$age=='12m']<-paste0('sub-',subj_idx)
  subj_idx<-subj_idx+1
}
subjs_adult<-unique(df$subject[df$age=='adult'])
subj_idx<-200
for(subj in subjs_adult) {
  df$new_subject[df$subject==subj & df$age=='adult']<-paste0('sub-',subj_idx)
  subj_idx<-subj_idx+1
}
df$subject<-as.factor(df$new_subject)
df$age<-as.factor(df$age)
df$condition<-as.factor(df$condition)
df$epoch<-as.factor(df$epoch)
c_df<-df
c_df$cluster<-'C'

df<-read.csv('../../output/bursts_p.csv')
df$new_subject<-as.character(df$subject)
subjs_12m<-unique(df$subject[df$age=='12m'])
subj_idx<-100
for(subj in subjs_12m) {
  df$new_subject[df$subject==subj & df$age=='12m']<-paste0('sub-',subj_idx)
  subj_idx<-subj_idx+1
}
subjs_adult<-unique(df$subject[df$age=='adult'])
subj_idx<-200
for(subj in subjs_adult) {
  df$new_subject[df$subject==subj & df$age=='adult']<-paste0('sub-',subj_idx)
  subj_idx<-subj_idx+1
}
df$subject<-as.factor(df$new_subject)
df$age<-as.factor(df$age)
df$condition<-as.factor(df$condition)
df$epoch<-as.factor(df$epoch)
p_df<-df
p_df$cluster<-'P'

df<-rbind(c_df, p_df)
df$cluster<-as.factor(df$cluster)

type3 <- list(age = contr.sum, cluster = contr.sum)
model <- lmer(peak_amp_base ~ age*cluster+(1|subject),
              data = df,
              contrasts=type3,
              control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
results<-Anova(model, type = 3)
print(results)
emmeans(model, pairwise ~ age|cluster, infer=TRUE)
emmeans(model, pairwise ~ cluster|age, infer=TRUE)
for(cluster in c('C','P')) {
  for(age in c('9m','12m','adult')) {
    m=mean(df$peak_amp_base[df$age==age & df$cluster==cluster])
    sd=sd(df$peak_amp_base[df$age==age & df$cluster==cluster])
    print(paste0(age,' ',cluster,':  M=',m,', SD=',sd))
  }
}


type3 <- list(age = contr.sum, cluster = contr.sum)
model <- lmer(fwhm_freq ~ age*cluster+(1|subject),
              data = df,
              contrasts=type3,
              control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
results<-Anova(model, type = 3)
print(results)
emmeans(model, pairwise ~ age|cluster, infer=TRUE)
emmeans(model, pairwise ~ cluster|age, infer=TRUE)
for(cluster in c('C','P')) {
  for(age in c('9m','12m','adult')) {
    m=mean(df$fwhm_freq[df$age==age & df$cluster==cluster])
    sd=sd(df$fwhm_freq[df$age==age & df$cluster==cluster])
    print(paste0(age,' ',cluster,':  M=',m,', SD=',sd))
  }
}




df$dur_cycles=df$fwhm_time/(1000/df$peak_freq)
type3 <- list(age = contr.sum, cluster = contr.sum)
model <- lmer(dur_cycles ~ age*cluster+(1|subject),
              data = df,
              contrasts=type3,
              control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
results<-Anova(model, type = 3)
print(results)
emmeans(model, pairwise ~ age|cluster, infer=TRUE)
emmeans(model, pairwise ~ cluster|age, infer=TRUE)
for(cluster in c('C','P')) {
  for(age in c('9m','12m','adult')) {
    m=mean(df$dur_cycles[df$age==age & df$cluster==cluster])
    sd=sd(df$dur_cycles[df$age==age & df$cluster==cluster])
    print(paste0(age,' ',cluster,':  M=',m,', SD=',sd))
  }
}
