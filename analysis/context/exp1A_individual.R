library(stringr)
library(dplyr)
library(ggplot2)
library(tidyr)

df = read.csv('exp1A.csv')

ggplot(filter(df,st > 1),
       aes(x=time,y=response,group=factor(st),color=factor(st))) +
  facet_grid(~factor(same_freq,
                     labels=c("Different Frequency","Same Frequency"))) +
  stat_summary(geom="line") +
  stat_summary(geom="line",color="black",linetype=2,
               data=rbind(mutate(filter(df,st == 0),same_freq=T),
                          mutate(filter(df,st == 0),same_freq=F)))

ggplot(filter(df,st > 1),
       aes(x=time,y=response,group=factor(st),color=factor(st))) +
  facet_grid(~factor(same_freq,
                     labels=c("Different Frequency","Same Frequency"))) +
  stat_summary(geom="ribbon",fun.data="mean_cl_boot",fill="lightgray",
               color="lightgray") +
  stat_summary(geom="line") +
  stat_summary(geom="line",color="black",linetype=2,
               data=rbind(mutate(filter(df,st == 0),same_freq=T),
                          mutate(filter(df,st == 0),same_freq=F)))

ggplot(filter(df,st > 1),
       aes(x=time,y=response,group=factor(st),color=factor(st))) +
  facet_grid(sid~factor(same_freq,
                     labels=c("Different Frequency","Same Frequency"))) +
  stat_summary(geom="line") +
  stat_summary(geom="line",color="black",linetype=2,
               data=rbind(mutate(filter(df,st == 0),same_freq=T),
                          mutate(filter(df,st == 0),same_freq=F)))

