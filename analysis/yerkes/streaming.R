library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(zoo)

# load all of the data into memory

raw_df = NULL

dir = file.path("..","..","data","Yerkes_et_al_2019","Experiment 1 Data")
files = list.files(dir,"*.txt")
for(file in files){
  cur_df = read.table(file.path(dir,file),header=F)
  cur_df$experiment = 1
  raw_df = rbind(raw_df,cur_df)
}
dir = file.path("..","..","data","Yerkes_et_al_2019","Experiment 2 Data")
files = list.files(dir,"*.txt")
for(file in files){
  cur_df = read.table(file.path(dir,file),header=F)
  cur_df$experiment = 2
  raw_df = rbind(raw_df,cur_df)
}

# convert raw data file to a relatively readable format (but don't manipulate
# anything)
colnames(raw_df) = c("sid","pres_trial","code","sample","time","unknown",
  "experiment")
code_to_response = c("602" = "fuse_on","603" = "fuse_off",
                     "604" = "stream_on","605" = "stream_off")
code_to_task = c("101" = "AM","102" = "AM","103" = "AM",
                 "104" = "streaming","105" = "streaming","106" = "streaming",
                 "107" = "visual","108" = "visual","109" = "visual")
code_to_st = c("101" = 3, "102" = 6, "103" = 12,
               "104" = 3, "105" = 6, "106" = 12,
               "107" = 3, "108" = 6, "109" = 12)

all_events = raw_df %>%
  mutate(trial = cumsum(time == 0), # each trial starts at time 0
         task = code_to_task[as.character(code)],
         st = code_to_st[as.character(code)],
         response = code_to_response[as.character(code)],
         time = time / 10000) # from micro-seconds to seconds

# streaming events occur when the stream button is on non-streaming events
# occur when the fuse button is on all other events are missing values
all_with_streaming = all_events %>%
  mutate(streaming =
    ifelse(response == "stream_on",TRUE,
    ifelse(response == "fuse_on",FALSE,NA)))

# for all rows with missing values, use the most recent non-missing value.
# then remove all but the rows with a specific task and st code. this way each
# row represents a single moment in time, and contains all the information we
# need: the task listeners performed, the semitone separation of the stimuli,
# and the current response of the listeners.
combined = all_with_streaming %>%
  group_by(trial,sid,experiment) %>%
  mutate(streaming = na.locf(streaming,na.rm=FALSE)) %>%
  filter()

# compute, save, and plot the aggregate streaming values for each individual
stream_prop = combined %>%
  filter(time > 4,time < 6.72, task == "streaming") %>%
  group_by(experiment, sid, st) %>%
  summarize(streaming = mean(streaming))

write.csv(stream_prop,"stream_prop.csv")

pos = position_jitter(width=0.1)
ggplot(stream_prop,aes(x=factor(st),y=streaming)) +
  # geom_point(alpha=0.5,position=pos) +
  geom_line(alpha=0.5,#position=pos,
      aes(group=interaction(experiment,sid))) +
  stat_summary(geom="pointrange",fun.data=mean_cl_boot,
               fun.args=list(conf.int=0.95)) +
  stat_summary(geom="line",fun.data=mean_cl_boot,group=1) +
  coord_cartesian(ylim=c(0,1)) +
  theme_classic()
