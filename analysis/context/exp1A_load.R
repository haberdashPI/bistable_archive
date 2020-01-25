library(stringr)
library(dplyr)
library(tidyr)

dir = file.path("..","..","data","Snyder_et_al_2009","Exp1A")
files = list.files(dir,"*_Exp2.asc")

raw_df = NULL
for (file in files) {
  cur_df = read.table(file.path(dir,file),header=F)
  cur_df$sid = str_sub(file,1,2)
  raw_df = rbind(raw_df,cur_df)
}

column_names = sprintf("unknown%02d",1:20)
column_names[3] = "trial_index"
column_names[4] = "context_A"
column_names[5] = "test_A"
column_names[6] = "context_B"
column_names[7] = "test_B"
column_names[8] = "condition_index"
column_names[12] = "response"

colnames(raw_df) = c(column_names,"sid")

column_names = sprintf("unknown%02d",1:20)
SOA = 0.12
sr = SOA*2
t = (0:61)*sr

df = raw_df %>%
  mutate(sample = 2*(trial_index-1) + rep(1:2,length(trial_index)/2),
         time = sample*sr,
         same_freq = test_A == context_A,
         st = round(12*log(context_B / context_A)/log(2))) %>%
  select(sid,sample,time,context_A,context_B,st,test_A,test_B,response,
         same_freq)

write.csv(df,"exp1A.csv")
