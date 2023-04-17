library(lme4)
library(car)
library(emmeans)
library(Rmisc)
library(tidyverse)
library(ggplot2)

sink('../../output/burst_stats_out.txt')

df<-read.csv('../../output/bursts.csv')
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
df$epoch<-as.factor(df$epoch)
df$peak_amp_base<-df$peak_amp_base*1e6

type3 <- list(age = contr.sum)
model <- lmer(peak_amp_base ~ age+(1|subject),
              data = df,
              contrasts=type3,
              control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
print(model)
results<-Anova(model, type = 2)
print(results)
print(emmeans(model, pairwise ~ age, infer=TRUE))
print(paste0('9m: M=',mean(df$peak_amp_base[df$age=='9m']),' SD=',sd(df$peak_amp_base[df$age=='9m'])))
print(paste0('12m: M=',mean(df$peak_amp_base[df$age=='12m']),' SD=',sd(df$peak_amp_base[df$age=='12m'])))
print(paste0('adult: M=',mean(df$peak_amp_base[df$age=='adult']),' SD=',sd(df$peak_amp_base[df$age=='adult'])))

print('')
print('')
print('')

type3 <- list(age = contr.sum)
model <- lmer(fwhm_freq ~ age+(1|subject),
              data = df,
              contrasts=type3,
              control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
print(model)
results<-Anova(model, type = 2)
print(results)
print(emmeans(model, pairwise ~ age, infer=TRUE))
print(paste0('9m: M=',mean(df$fwhm_freq[df$age=='9m']),' SD=',sd(df$fwhm_freq[df$age=='9m'])))
print(paste0('12m: M=',mean(df$fwhm_freq[df$age=='12m']),' SD=',sd(df$fwhm_freq[df$age=='12m'])))
print(paste0('adult: M=',mean(df$fwhm_freq[df$age=='adult']),' SD=',sd(df$fwhm_freq[df$age=='adult'])))

print('')
print('')
print('')


df$dur_cycles=df$fwhm_time/(1000/df$peak_freq)
type3 <- list(age = contr.sum)
model <- lmer(dur_cycles ~ age+(1|subject),
              data = df,
              contrasts=type3,
              control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
print(model)
results<-Anova(model, type = 2)
print(results)
print(emmeans(model, pairwise ~ age, infer=TRUE))
print(paste0('9m: M=',mean(df$dur_cycles[df$age=='9m']),' SD=',sd(df$dur_cycles[df$age=='9m'])))
print(paste0('12m: M=',mean(df$dur_cycles[df$age=='12m']),' SD=',sd(df$dur_cycles[df$age=='12m'])))
print(paste0('adult: M=',mean(df$dur_cycles[df$age=='adult']),' SD=',sd(df$dur_cycles[df$age=='adult'])))

print('')
print('')
print('')


exe_df<-df[df$cluster=='ipsi' | df$cluster=='contra',]

type3 <- list(age = contr.sum, cluster=contr.sum)
model <- lmer(peak_amp_base ~ age*cluster+(1|subject),
              data = exe_df,
              contrasts=type3,
              control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
print(model)
results<-Anova(model, type = 3)
print(results)
print(emmeans(model, pairwise ~ age|cluster, infer=TRUE))
print(emmeans(model, pairwise ~ cluster|age, infer=TRUE))

print(ks.test(exe_df$peak_amp_base[exe_df$age=='9m' & exe_df['cluster']=='ipsi'],exe_df$peak_amp_base[exe_df$age=='9m' & exe_df['cluster']=='contra']))
print(ks.test(exe_df$peak_amp_base[exe_df$age=='12m' & exe_df['cluster']=='ipsi'],exe_df$peak_amp_base[exe_df$age=='12m' & exe_df['cluster']=='contra']))
print(ks.test(exe_df$peak_amp_base[exe_df$age=='adult' & exe_df['cluster']=='ipsi'],exe_df$peak_amp_base[exe_df$age=='adult' & exe_df['cluster']=='contra']))

print(paste0('9m ipsi: M=',mean(exe_df$peak_amp_base[exe_df$age=='9m' & exe_df['cluster']=='ipsi']),' SD=',sd(exe_df$peak_amp_base[exe_df$age=='9m' & exe_df['cluster']=='ipsi'])))
print(paste0('9m contra: M=',mean(exe_df$peak_amp_base[exe_df$age=='9m' & exe_df['cluster']=='contra']),' SD=',sd(exe_df$peak_amp_base[exe_df$age=='9m' & exe_df['cluster']=='contra'])))
print(paste0('12m ipsi: M=',mean(exe_df$peak_amp_base[exe_df$age=='12m' & exe_df['cluster']=='ipsi']),' SD=',sd(exe_df$peak_amp_base[exe_df$age=='12m' & exe_df['cluster']=='ipsi'])))
print(paste0('12m contra: M=',mean(exe_df$peak_amp_base[exe_df$age=='12m' & exe_df['cluster']=='contra']),' SD=',sd(exe_df$peak_amp_base[exe_df$age=='12m' & exe_df['cluster']=='contra'])))
print(paste0('adult ipsi: M=',mean(exe_df$peak_amp_base[exe_df$age=='adult' & exe_df['cluster']=='ipsi']),' SD=',sd(exe_df$peak_amp_base[exe_df$age=='adult' & exe_df['cluster']=='ipsi'])))
print(paste0('adult contra: M=',mean(exe_df$peak_amp_base[exe_df$age=='adult' & exe_df['cluster']=='contra']),' SD=',sd(exe_df$peak_amp_base[exe_df$age=='adult' & exe_df['cluster']=='contra'])))

print('')
print('')
print('')


type3 <- list(age = contr.sum, cluster=contr.sum)
model <- lmer(fwhm_freq ~ age*cluster+(1|subject),
              data = exe_df,
              contrasts=type3,
              control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
print(model)
results<-Anova(model, type = 3)
print(results)
print(emmeans(model, pairwise ~ age|cluster, infer=TRUE))
print(emmeans(model, pairwise ~ cluster|age, infer=TRUE))

print(ks.test(exe_df$fwhm_freq[exe_df$age=='9m' & exe_df['cluster']=='ipsi'],exe_df$fwhm_freq[exe_df$age=='9m' & exe_df['cluster']=='contra']))
print(ks.test(exe_df$fwhm_freq[exe_df$age=='12m' & exe_df['cluster']=='ipsi'],exe_df$fwhm_freq[exe_df$age=='12m' & exe_df['cluster']=='contra']))
print(ks.test(exe_df$fwhm_freq[exe_df$age=='adult' & exe_df['cluster']=='ipsi'],exe_df$fwhm_freq[exe_df$age=='adult' & exe_df['cluster']=='contra']))

print(paste0('9m ipsi: M=',mean(exe_df$fwhm_freq[exe_df$age=='9m' & exe_df['cluster']=='ipsi']),' SD=',sd(exe_df$fwhm_freq[exe_df$age=='9m' & exe_df['cluster']=='ipsi'])))
print(paste0('9m contra: M=',mean(exe_df$fwhm_freq[exe_df$age=='9m' & exe_df['cluster']=='contra']),' SD=',sd(exe_df$fwhm_freq[exe_df$age=='9m' & exe_df['cluster']=='contra'])))
print(paste0('12m ipsi: M=',mean(exe_df$fwhm_freq[exe_df$age=='12m' & exe_df['cluster']=='ipsi']),' SD=',sd(exe_df$fwhm_freq[exe_df$age=='12m' & exe_df['cluster']=='ipsi'])))
print(paste0('12m contra: M=',mean(exe_df$fwhm_freq[exe_df$age=='12m' & exe_df['cluster']=='contra']),' SD=',sd(exe_df$fwhm_freq[exe_df$age=='12m' & exe_df['cluster']=='contra'])))
print(paste0('adult ipsi: M=',mean(exe_df$fwhm_freq[exe_df$age=='adult' & exe_df['cluster']=='ipsi']),' SD=',sd(exe_df$fwhm_freq[exe_df$age=='adult' & exe_df['cluster']=='ipsi'])))
print(paste0('adult contra: M=',mean(exe_df$fwhm_freq[exe_df$age=='adult' & exe_df['cluster']=='contra']),' SD=',sd(exe_df$fwhm_freq[exe_df$age=='adult' & exe_df['cluster']=='contra'])))

print('')
print('')
print('')


type3 <- list(age = contr.sum, cluster=contr.sum)
model <- lmer(dur_cycles ~ age*cluster+(1|subject),
              data = exe_df,
              contrasts=type3,
              control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
print(model)
results<-Anova(model, type = 3)
print(results)
print(emmeans(model, pairwise ~ age|cluster, infer=TRUE))
print(emmeans(model, pairwise ~ cluster|age, infer=TRUE))

print(ks.test(exe_df$dur_cycles[exe_df$age=='9m' & exe_df['cluster']=='ipsi'],exe_df$dur_cycles[exe_df$age=='9m' & exe_df['cluster']=='contra']))
print(ks.test(exe_df$dur_cycles[exe_df$age=='12m' & exe_df['cluster']=='ipsi'],exe_df$dur_cycles[exe_df$age=='12m' & exe_df['cluster']=='contra']))
print(ks.test(exe_df$dur_cycles[exe_df$age=='adult' & exe_df['cluster']=='ipsi'],exe_df$dur_cycles[exe_df$age=='adult' & exe_df['cluster']=='contra']))


print(paste0('9m ipsi: M=',mean(exe_df$dur_cycles[exe_df$age=='9m' & exe_df['cluster']=='ipsi']),' SD=',sd(exe_df$dur_cycles[exe_df$age=='9m' & exe_df['cluster']=='ipsi'])))
print(paste0('9m contra: M=',mean(exe_df$dur_cycles[exe_df$age=='9m' & exe_df['cluster']=='contra']),' SD=',sd(exe_df$dur_cycles[exe_df$age=='9m' & exe_df['cluster']=='contra'])))
print(paste0('12m ipsi: M=',mean(exe_df$dur_cycles[exe_df$age=='12m' & exe_df['cluster']=='ipsi']),' SD=',sd(exe_df$dur_cycles[exe_df$age=='12m' & exe_df['cluster']=='ipsi'])))
print(paste0('12m contra: M=',mean(exe_df$dur_cycles[exe_df$age=='12m' & exe_df['cluster']=='contra']),' SD=',sd(exe_df$dur_cycles[exe_df$age=='12m' & exe_df['cluster']=='contra'])))
print(paste0('adult ipsi: M=',mean(exe_df$dur_cycles[exe_df$age=='adult' & exe_df['cluster']=='ipsi']),' SD=',sd(exe_df$dur_cycles[exe_df$age=='adult' & exe_df['cluster']=='ipsi'])))
print(paste0('adult contra: M=',mean(exe_df$dur_cycles[exe_df$age=='adult' & exe_df['cluster']=='contra']),' SD=',sd(exe_df$dur_cycles[exe_df$age=='adult' & exe_df['cluster']=='contra'])))

print('')
print('')
print('')


type3 <- list(age = contr.sum, cluster=contr.sum)
model <- lmer(peak_freq ~ age*cluster+(1|age/subject),
              data = exe_df,
              contrasts=type3,
              control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
print(model)
results<-Anova(model, type = 3)
print(results)
print(emmeans(model, pairwise ~ age|cluster, infer=TRUE))
print(emmeans(model, pairwise ~ cluster|age, infer=TRUE))
print(paste0('9m ipsi: M=',mean(exe_df$peak_freq[exe_df$age=='9m' & exe_df['cluster']=='ipsi']),' SD=',sd(exe_df$peak_freq[exe_df$age=='9m' & exe_df['cluster']=='ipsi'])))
print(paste0('9m contra: M=',mean(exe_df$peak_freq[exe_df$age=='9m' & exe_df['cluster']=='contra']),' SD=',sd(exe_df$peak_freq[exe_df$age=='9m' & exe_df['cluster']=='contra'])))
print(paste0('12m ipsi: M=',mean(exe_df$peak_freq[exe_df$age=='12m' & exe_df['cluster']=='ipsi']),' SD=',sd(exe_df$peak_freq[exe_df$age=='12m' & exe_df['cluster']=='ipsi'])))
print(paste0('12m contra: M=',mean(exe_df$peak_freq[exe_df$age=='12m' & exe_df['cluster']=='contra']),' SD=',sd(exe_df$peak_freq[exe_df$age=='12m' & exe_df['cluster']=='contra'])))
print(paste0('adult ipsi: M=',mean(exe_df$peak_freq[exe_df$age=='adult' & exe_df['cluster']=='ipsi']),' SD=',sd(exe_df$peak_freq[exe_df$age=='adult' & exe_df['cluster']=='ipsi'])))
print(paste0('adult contra: M=',mean(exe_df$peak_freq[exe_df$age=='adult' & exe_df['cluster']=='contra']),' SD=',sd(exe_df$peak_freq[exe_df$age=='adult' & exe_df['cluster']=='contra'])))

# type3 <- list(age = contr.sum)
# model <- lmer(dur_cycles ~ peak_amp_base*age+(1|subject),
#               data = df,
#               contrasts=type3,
#               control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
# results<-Anova(model, type = 2)
# print(results)
# emtrends(model, pairwise ~ age, infer=TRUE, var='peak_amp_base')
# 
# type3 <- list(age = contr.sum)
# model <- lmer(dur_cycles ~ fwhm_freq*age+(1|subject),
#               data = df,
#               contrasts=type3,
#               control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
# results<-Anova(model, type = 2)
# print(results)
# emtrends(model, pairwise ~ age, infer=TRUE, var='fwhm_freq')
# 
# type3 <- list(age = contr.sum)
# model <- lmer(dur_cycles ~ peak_freq*age+(1|subject),
#               data = df,
#               contrasts=type3,
#               control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
# results<-Anova(model, type = 2)
# print(results)
# emtrends(model, pairwise ~ age, infer=TRUE, var='peak_freq')
# 
# type3 <- list(age = contr.sum)
# model <- lmer(peak_amp_base ~ fwhm_freq*age+(1|subject),
#               data = df,
#               contrasts=type3,
#               control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
# results<-Anova(model, type = 2)
# print(results)
# emtrends(model, pairwise ~ age, infer=TRUE, var='fwhm_freq')
# 
# type3 <- list(age = contr.sum)
# model <- lmer(peak_amp_base ~ peak_freq*age+(1|subject),
#               data = df,
#               contrasts=type3,
#               control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
# results<-Anova(model, type = 2)
# print(results)
# emtrends(model, pairwise ~ age, infer=TRUE, var='peak_freq')
# 
# type3 <- list(age = contr.sum)
# model <- lmer(fwhm_freq ~ peak_freq*age+(1|subject),
#               data = df,
#               contrasts=type3,
#               control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
# results<-Anova(model, type = 2)
# print(results)
# emtrends(model, pairwise ~ age, infer=TRUE, var='peak_freq')

sink()