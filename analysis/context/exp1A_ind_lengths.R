library(stringr)
library(dplyr)
library(ggplot2)
library(tidyr)

df = read.csv("exp1A.csv")

df = df %>%
  group_by(sid) %>%
  mutate(trial = cumsum(lag(sample,default=first(sample)) > sample))

# compute percept lens
lens = df %>%
  group_by(sid,trial) %>%
  mutate(percept = cumsum(lag(response,
                              default=first(response)) != response)) %>%
  group_by(sid,trial,percept) %>%
  summarize(context_A = first(context_A),context_B = first(context_B),
            test_A = first(test_A), test_B = first(test_B), st = first(st),
            len = last(time) - first(time),stimulus = first(response))

# plot some of the individual trials to verify the summary
ggplot(filter(df,sid == "EH",trial == 30),aes(x=time,y=response)) +
  geom_line()
filter(lens,sid == "EH",trial == 30)

ggplot(filter(df,sid == "RR",trial == 12),aes(x=time,y=response)) +
  geom_line()
filter(lens,sid == "RR",trial == 12)

# show some statistics across listeners

ggplot(lens,aes(x=len,fill=factor(stimulus))) +
  geom_density(bw=0.25*sd(lens$len),alpha=0.5) + facet_grid(st~sid)

ggplot(lens,aes(x=len)) +
  geom_density(fill="gray",bw=0.25*sd(lens$len)) + facet_grid(st~sid)

nb_lens = lens %>%
  group_by(sid,trial) %>%
  mutate(len = lead(len))

ggplot(nb_lens,aes(x=len,fill=factor(stimulus))) +
  geom_histogram() + facet_grid(st~sid) +
  coord_cartesian(ylim=c(0,10),xlim=c(0,10))

ggplot(nb_lens,aes(x=len,fill=factor(stimulus))) +
  geom_density(bw=0.25*sd(lens$len),alpha=0.5) + facet_grid(st~sid)

ggplot(nb_lens,aes(x=len)) +
  geom_density(fill="gray",bw=0.25*sd(lens$len)) + facet_grid(st~sid)

ggplot(nb_lens,aes(y=len,x=stimulus)) +
  geom_point(position=position_jitter(width=0.05),alpha=0.5,size=1) +
  facet_grid(st~sid) + xlim(-0.5,1.5)

# a few thoughts looking at the histograms (seems to be the most
# informative view)

# - 1. it looks like there is some consistentcy across
#      conditions for each participant (is this some
#      interaction with context?)
# - 2. there are a fair number of

# wait... what? there are a lot of 1-stream responses...
# for ALL conditions... that must be about the bias....???

