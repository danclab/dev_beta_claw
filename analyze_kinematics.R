library(lme4)
library(car)
library(lsmeans)

sink('../../output/kinematics_stats_out.txt')

df_9m_exe<-read.csv('../../data/9m/behavior_exe.csv')
df_9m_exe$bi<-FALSE
df_9m_exe$Fthand[df_9m_exe$Fthand=='R ']<-'R'
df_9m_exe<-df_9m_exe[df_9m_exe$FTGhand=='L' | df_9m_exe$FTGhand=='R',]
df_9m_exe$Age<-'9m'
df_9m_exe$Subject<-paste0(df_9m_exe$Subject,'-9m')

df_12m_exe<-read.csv('../../data/12m/behavior_exe.csv')
df_12m_exe$bi<-FALSE
df_12m_exe$Fthand[df_12m_exe$Fthand=='R ']<-'R'
df_12m_exe$Fthand[df_12m_exe$Fthand=='L ']<-'L'
df_12m_exe$bi[!(df_12m_exe$FTGhand=='L' | df_12m_exe$Fthand=='R')]<-TRUE
df_12m_exe<-df_12m_exe[df_12m_exe$FTGhand!='',]
df_12m_exe$Age<-'12m'
df_12m_exe$Subject<-paste0(df_12m_exe$Subject,'-12m')

df_adult_exe<-read.csv('../../data/adult/behavior_exe.csv')
df_adult_exe$bi<-FALSE
df_adult_exe$bi[!(df_adult_exe$FTGhand=='L' | df_adult_exe$Fthand=='R')]<-TRUE
df_adult_exe<-df_adult_exe[df_adult_exe$FTGhand!='',]
df_adult_exe$Age<-'adult'
df_adult_exe$Subject<-paste0(df_adult_exe$Subject,'-adult')

df_exe<-rbind(df_9m_exe,df_12m_exe,df_adult_exe)
type3 <- list(Age = contr.sum)
model <- glmer(bi ~ Age+(1|Subject),
              data = df_exe,
              contrasts=type3,
              control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)),
              family = binomial)
print(model)
results<-Anova(model, type = 2)
print(results)
print(paste0('9m: M=',mean(df_exe$bi[df_exe$Age=='9m'],na.rm=TRUE),' SD=',sd(df_exe$bi[df_exe$Age=='9m'],na.rm = TRUE)))
print(paste0('12m: M=',mean(df_exe$bi[df_exe$Age=='12m'],na.rm = TRUE),' SD=',sd(df_exe$bi[df_exe$Age=='12m'],na.rm = TRUE)))
print(paste0('adult: M=',mean(df_exe$bi[df_exe$Age=='adult'],na.rm = TRUE),' SD=',sd(df_exe$bi[df_exe$Age=='adult'],na.rm = TRUE)))


print('')
print('')
print('')

df_exe<-df_exe[df_exe$Fthand=='L' | df_exe$FTGhand=='R',]
df_exe$right<-TRUE
df_exe$right[df_exe$FTGhand=='L']<-FALSE
type3 <- list(Age = contr.sum)
model <- glmer(right ~ Age+(1|Subject),
               data = df_exe,
               contrasts=type3,
               control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)),
               family = binomial)
print(model)
results<-Anova(model, type = 2)
print(results)
print(emmeans(model, pairwise ~ Age, infer=TRUE))
print(paste0('9m: M=',mean(df_exe$right[df_exe$Age=='9m'],na.rm=TRUE),' SD=',sd(df_exe$right[df_exe$Age=='9m'],na.rm = TRUE)))
print(paste0('12m: M=',mean(df_exe$right[df_exe$Age=='12m'],na.rm = TRUE),' SD=',sd(df_exe$right[df_exe$Age=='12m'],na.rm = TRUE)))
print(paste0('adult: M=',mean(df_exe$right[df_exe$Age=='adult'],na.rm = TRUE),' SD=',sd(df_exe$right[df_exe$Age=='adult'],na.rm = TRUE)))


print('')
print('')
print('')


df_9m<-read.csv('../../data/9m/derivatives/NEARICA_behav/processed_kinematics.csv')
df_9m$Age<-'9m'
df_9m$Subject<-paste0(df_9m$Subject,'-9m')
df_9m_exe=df_9m[df_9m$Condition=='exe',]
df_9m_obs=df_9m[df_9m$Condition=='obs',]

df_12m<-read.csv('../../data/12m/derivatives/NEARICA_behav/processed_kinematics.csv')
df_12m$Age<-'12m'
df_12m$Subject<-paste0(df_12m$Subject,'-12m')
df_12m_exe=df_12m[df_12m$Condition=='exe',]
df_12m_obs=df_12m[df_12m$Condition=='obs',]

df_adult<-read.csv('../../data/adult/derivatives/NEARICA_behav/processed_kinematics.csv')
df_adult$Age<-'adult'
df_adult$Subject<-paste0(df_adult$Subject,'-adult')
df_adult_exe=df_adult[df_adult$Condition=='exe',]
df_adult_obs=df_adult[df_adult$Condition=='obs',]

df_exe<-rbind(df_9m_exe,df_12m_exe,df_adult_exe)
type3 <- list(Age = contr.sum)
model <- lmer(ReachDur ~ Age+(1|Subject),
              data = df_exe,
              contrasts=type3,
              control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
print(model)
results<-Anova(model, type = 2)
print(results)
print(emmeans(model, pairwise ~ Age, infer=TRUE))
print(paste0('9m: M=',mean(df_exe$ReachDur[df_exe$Age=='9m'],na.rm=TRUE),' SD=',sd(df_exe$ReachDur[df_exe$Age=='9m'],na.rm = TRUE)))
print(paste0('12m: M=',mean(df_exe$ReachDur[df_exe$Age=='12m'],na.rm = TRUE),' SD=',sd(df_exe$ReachDur[df_exe$Age=='12m'],na.rm = TRUE)))
print(paste0('adult: M=',mean(df_exe$ReachDur[df_exe$Age=='adult'],na.rm = TRUE),' SD=',sd(df_exe$ReachDur[df_exe$Age=='adult'],na.rm = TRUE)))


print('')
print('')
print('')


type3 <- list(Age = contr.sum)
model <- lmer(GraspDur ~ Age+(1|Subject),
              data = df_exe,
              contrasts=type3,
              control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
print(model)
results<-Anova(model, type = 2)
print(results)
print(emmeans(model, pairwise ~ Age, infer=TRUE))
print(paste0('9m: M=',mean(df_exe$GraspDur[df_exe$Age=='9m'],na.rm=TRUE),' SD=',sd(df_exe$GraspDur[df_exe$Age=='9m'],na.rm = TRUE)))
print(paste0('12m: M=',mean(df_exe$GraspDur[df_exe$Age=='12m'],na.rm = TRUE),' SD=',sd(df_exe$GraspDur[df_exe$Age=='12m'],na.rm = TRUE)))
print(paste0('adult: M=',mean(df_exe$GraspDur[df_exe$Age=='adult'],na.rm = TRUE),' SD=',sd(df_exe$GraspDur[df_exe$Age=='adult'],na.rm = TRUE)))


print('')
print('')
print('')


type3 <- list(Age = contr.sum)
model <- lmer(EndDur ~ Age+(1|Subject),
              data = df_exe,
              contrasts=type3,
              control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
print(model)
results<-Anova(model, type = 2)
print(results)
print(emmeans(model, pairwise ~ Age, infer=TRUE))
print(paste0('9m: M=',mean(df_exe$EndDur[df_exe$Age=='9m'],na.rm=TRUE),' SD=',sd(df_exe$EndDur[df_exe$Age=='9m'],na.rm = TRUE)))
print(paste0('12m: M=',mean(df_exe$EndDur[df_exe$Age=='12m'],na.rm = TRUE),' SD=',sd(df_exe$EndDur[df_exe$Age=='12m'],na.rm = TRUE)))
print(paste0('adult: M=',mean(df_exe$EndDur[df_exe$Age=='adult'],na.rm = TRUE),' SD=',sd(df_exe$EndDur[df_exe$Age=='adult'],na.rm = TRUE)))



print('')
print('')
print('')



df_obs<-rbind(df_9m_obs,df_12m_obs,df_adult_obs)

type3 <- list(Age = contr.sum)
model <- lmer(ReachDur ~ Age+(1|Subject),
              data = df_obs,
              contrasts=type3,
              control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
print(model)
results<-Anova(model, type = 2)
print(results)
print(emmeans(model, pairwise ~ Age, infer=TRUE))
print(paste0('9m: M=',mean(df_obs$ReachDur[df_obs$Age=='9m'],na.rm=TRUE),' SD=',sd(df_obs$ReachDur[df_obs$Age=='9m'],na.rm = TRUE)))
print(paste0('12m: M=',mean(df_obs$ReachDur[df_obs$Age=='12m'],na.rm = TRUE),' SD=',sd(df_obs$ReachDur[df_obs$Age=='12m'],na.rm = TRUE)))
print(paste0('adult: M=',mean(df_obs$ReachDur[df_obs$Age=='adult'],na.rm = TRUE),' SD=',sd(df_obs$ReachDur[df_obs$Age=='adult'],na.rm = TRUE)))


print('')
print('')
print('')


type3 <- list(Age = contr.sum)
model <- lmer(GraspDur ~ Age+(1|Subject),
              data = df_obs,
              contrasts=type3,
              control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
print(model)
results<-Anova(model, type = 2)
print(results)
print(emmeans(model, pairwise ~ Age, infer=TRUE))
print(paste0('9m: M=',mean(df_obs$GraspDur[df_obs$Age=='9m'],na.rm=TRUE),' SD=',sd(df_obs$GraspDur[df_obs$Age=='9m'],na.rm = TRUE)))
print(paste0('12m: M=',mean(df_obs$GraspDur[df_obs$Age=='12m'],na.rm = TRUE),' SD=',sd(df_obs$GraspDur[df_obs$Age=='12m'],na.rm = TRUE)))
print(paste0('adult: M=',mean(df_obs$GraspDur[df_obs$Age=='adult'],na.rm = TRUE),' SD=',sd(df_obs$GraspDur[df_obs$Age=='adult'],na.rm = TRUE)))


print('')
print('')
print('')


type3 <- list(Age = contr.sum)
model <- lmer(EndDur ~ Age+(1|Subject),
              data = df_obs,
              contrasts=type3,
              control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
print(model)
results<-Anova(model, type = 2)
print(results)
print(emmeans(model, pairwise ~ Age, infer=TRUE))
print(paste0('9m: M=',mean(df_obs$EndDur[df_obs$Age=='9m'],na.rm=TRUE),' SD=',sd(df_obs$EndDur[df_obs$Age=='9m'],na.rm = TRUE)))
print(paste0('12m: M=',mean(df_obs$EndDur[df_obs$Age=='12m'],na.rm = TRUE),' SD=',sd(df_obs$EndDur[df_obs$Age=='12m'],na.rm = TRUE)))
print(paste0('adult: M=',mean(df_obs$EndDur[df_obs$Age=='adult'],na.rm = TRUE),' SD=',sd(df_obs$EndDur[df_obs$Age=='adult'],na.rm = TRUE)))

sink()