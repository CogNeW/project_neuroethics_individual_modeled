#this is to check the response times and time taken during the both the survey and task

#load libs

#load data

clean.data$user_id = factor(clean.data$user_id)
RT = select_(clean.data,"user_id", "rt", "trial_index")

RT %>% 
  group_by(user_id) %>% 
  count(user_id) %>% 
  filter(! n == 75) 

RT %>% 
  group_by(user_id) %>% 
  count(user_id) %>% 
  filter(! n == 75) %>% 
  ggplot( aes(n))+geom_histogram(fill = "deepskyblue",
                                 color="black",
                                 bins = 60)+
  scale_x_continuous("number of trials", breaks=seq(0,75,5)) +
  theme_bw()

summary(RT$user_id)
##Plots for RT 
data_long = reshape2::melt(RT, 
                           id=c("user_id","trial_index"),
                           measure=c("rt"))

ggplot(data_long, aes(trial_index, value)) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun.data=mean_se, geom = "pointrange") +
  scale_y_continuous("RT (milliseconds)", breaks=seq(0,10000,300)) +
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) +
  labs(title="", 
       subtitle="",
       caption="",
       x="Trial Number",
       y="RT (milliseconds)")+
  theme_bw()

##we want to remove trials from subjecct with RT within the paramerts of 250 and 30000 Ms
##120000 would be two minutes!

datalist = list()
for (i in RT$user_id){
  x = filter(RT, user_id == i)
  dat = summary(x$rt) 
  datalist[[i]] = dat # add it to your list
}

user_id = factor(names(datalist)) #creates vector of user_id's 

for (i in user_id){
  x = filter(RT, user_id == i)
  sum = summary(x)
  print(sum)
}

rm(data_long, x, datalist, trial_index)




raw.data %>% 
  dplyr::select(user_id, time_elapsed, trial_index) %>% 
  filter(trial_index == 229) %>% 
  mutate(task_completed = time_elapsed/60000) %>% 
  dplyr::select(user_id, task_completed) -> task.time

task.time <- merge(clean.data, task.time, "user_id")
task.time %>% 
  dplyr::select(user_id, task_completed)  %>% distinct()-> sta
merge(clean.data, sta, "user_id")
task.time %>% 
  dplyr::select(user_id, task_completed) %>% 
  distinct() %>% 
  summarize(n = n(),
            mean = mean(task_completed),
            median =median(task_completed),
            sd = sd(task_completed),
            se = sd / sqrt(n),
            ci = qt(0.975, df = n - 1) * sd / sqrt(n))
psych::describe(task.time$task_completed)