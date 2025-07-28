library(tidyverse)
library(lme4)
library(lmerTest)
library(patchwork)

Neuroethics_Judgement <- read.csv(here::here("data", "Cleaned_Pilot", "clean_data_longV2.csv"), stringsAsFactors=TRUE) %>% 
  select(user_id, condition, risk, gain, experimental_treatment_selected, rt, Yes, No, age, sex, Income, education, race_ethnicity) %>% 
  mutate(experimental_treatment_selected = ifelse(experimental_treatment_selected == "Yes to treatment", 1, 0), trials = Yes + No) %>% 
  mutate(dominate = ifelse(Yes >= 75 *.90 | No >= 75 *.90, "dominate", "not_dominate")) %>% 
  filter(!dominate == "dominate") -> data

#Descriptive results
data %>% 
  select(user_id, No, trials) %>% 
  group_by(user_id) %>% 
  mutate(average_rejection_probability = (No/trials)*100,
         above_50 = factor(ifelse(average_rejection_probability > 50, 1, 0))) %>% 
  distinct() %>% 
  ungroup() -> Descriptive
Descriptive %>% 
  summarise(n = n(),
            mean = mean(average_rejection_probability),
            sd = sd(average_rejection_probability),
            se = sd / sqrt(n),
            ci = qt(0.975, df = n - 1) * sd / sqrt(n))
#calculate the participants being more likely to accept than reject the treament
(count(filter(Descriptive, above_50 == 1))/count(Descriptive))*100
#statisitcal test aganst 50
t.test(Descriptive$average_rejection_probability, mu=50)

#plot
data %>% 
  select(user_id, risk, gain, experimental_treatment_selected,Yes, No, trials) %>% 
  mutate(ratio = gain/risk) %>% 
  group_by(user_id) %>% 
  mutate(average_acceptance_probability = (Yes/trials))%>% 
  ungroup()%>% 
  ggplot(aes(ratio, average_acceptance_probability)) +
  stat_summary(fun.data = mean_se, geom="pointrange", color = "black", size = .25)+
  labs(x= "gain/risk",
       y= "P(Treatment)")+
  theme_bw(base_size = 12)+
  theme(legend.position = "none")+
  geom_hline(yintercept = .5, linetype = "dashed",color = "red")+
  geom_vline(xintercept = 3.65, linetype = "dashed",color = "red")+
  scale_y_continuous(labels = scales::percent)
  
data %>% 
  select(user_id, rt, experimental_treatment_selected) %>% 
  mutate(rt = rt/1000) %>%
  summarise(n = n(),
            mean = mean(rt),
            sd = sd(rt),
            se = sd / sqrt(n),
            ci = qt(0.975, df = n - 1) * sd / sqrt(n)) 

data %>% 
  filter(experimental_treatment_selected == 0) %>% 
  summarise(mean = mean(rt/1000))

data %>% 
  filter(experimental_treatment_selected == 1) %>% 
  summarise(mean = mean(rt/1000))


data %>% 
  select(user_id, rt, experimental_treatment_selected) %>% 
  mutate(rt = rt/1000) %>% 
  group_by(experimental_treatment_selected, user_id) %>%
  summarise(n = n(),
            rt_avg = mean(rt),
            sd = sd(rt),
            se = sd / sqrt(n),
            ci = qt(0.975, df = n - 1) * sd / sqrt(n)) -> test

# Sort by group then ID
test <- test[order(test$experimental_treatment_selected, test$user_id), ]
t.test(rt_avg ~ experimental_treatment_selected, test, paired=TRUE)


##how does this vary if I do not group by the subjects? 
ggplot(test, aes(rt_avg, color = factor(experimental_treatment_selected), fill = factor(experimental_treatment_selected)))+
  geom_density(color = "black", alpha = .55)+theme_classic2()

wilcox.test(rt_avg ~ experimental_treatment_selected, test, paired=TRUE, conf.int = TRUE)

test <- test %>% 
  select(user_id, rt_avg, experimental_treatment_selected) %>% 
  arrange(user_id) %>% 
  group_by(user_id) %>% 
  pivot_wider(names_from = experimental_treatment_selected, values_from = rt_avg) %>% 
  mutate(Z = factor(ifelse(`0` < `1`, 1, 0)))
summary(test)

#we want the correlation between average rejection probability and RT differences.

data %>% select(user_id, rt,experimental_treatment_selected) %>% 
correlation::correlation()

###looking at the response times more closely
Neuroethics_Judgement <- read.csv(here::here("data", "Cleaned_Pilot", "clean_data_longV2.csv"), stringsAsFactors=TRUE) %>% 
  select(user_id, condition, risk, gain, no_effect, experimental_treatment_selected, rt, Yes, No) %>% 
  mutate(experimental_treatment_selected = ifelse(experimental_treatment_selected == "Yes to treatment", 1, 0), trials = Yes + No) %>% 
  mutate(dominate = ifelse(Yes >= 75 *.90 | No >= 75 *.90, "dominate", "not_dominate")) %>% 
  filter(!dominate == "dominate")
#what if you do something about the correlation, 
#you just orthogonize the variables
pca_res <- prcomp(Neuroethics_Judgement[,c("risk", "gain", "no_effect")]) 
summary(pca_res)
biplot(pca_res)
Neuroethics_Judgement$pca1 <- pca_res$x[,1]
Neuroethics_Judgement$pca2 <-  pca_res$x[,2]

TEST <- Neuroethics_Judgement %>% 
  filter(!rt <= 250) %>% #before or after the intial cutoff
  filter(!rt >= 9000) %>% 
  group_by(user_id) %>% 
  mutate(Q1 = quantile(rt, .25),
         Q3 = quantile(rt, .75),
         IQR = IQR(rt)) %>% 
  filter(rt > (Q1 - 2*IQR) & rt < (Q3 + 2*IQR)) %>% 
  group_by(user_id) %>%
  mutate(
    risk = scale(risk),
    gain = scale(gain)
  ) %>%
  ungroup()

TEST$condition <- relevel(TEST$condition, "Motor Function")

psych::describeBy(TEST$rt, group = TEST$condition)
fivenum(TEST$rt)
plot(density(TEST$rt))

m_0 <- glmer(experimental_treatment_selected ~ (pca1 + pca2)*condition +
        (pca1 + pca2 | user_id), 
      data = TEST, family = binomial(link="logit"), 
      control=glmerControl(optimizer="bobyqa"))

summary(m_0)

TEST_2 <- TEST[TEST$experimental_treatment_selected == 1,]
m_1 <- lmer(log(rt) ~ (pca1 + pca2)*condition +
               (pca1 + pca2 | user_id), 
             data = TEST_2, 
             control=lmerControl(optimizer="bobyqa"))
summary(m_1)
#As the risk increases people take longer to select the experimental treatment,
#but the same is not true for gains

#so the choice is effected by the amount of gains, but
#only the are changes in risk effecting the RT
m_0 <- glmer(experimental_treatment_selected ~ (risk + gain)*condition + (risk + gain | user_id), 
             data = TEST, family = binomial(link="logit"), 
             control=glmerControl(optimizer="bobyqa"))
TEST_2 <- TEST[TEST$experimental_treatment_selected == 1,]
m_1 <- lmer(log(rt) ~ (risk + gain)*condition + (risk + gain | user_id), 
            data = TEST_2, control=lmerControl(optimizer="bobyqa"))
summary(m_1)
report::report(m_1)

TEST_3 <- TEST[TEST$experimental_treatment_selected == 0,]
m_1 <- lmer(log(rt) ~ (risk + gain)*condition + (risk + gain | user_id), 
            data = TEST_3, control=lmerControl(optimizer="bobyqa"))
summary(m_1)
report::report(m_1)

# summary(m_1)

#choice plot for manuscript
p1 <- TEST %>% 
  group_by(risk, condition) %>%
  mutate(n = n(),
         m_left = mean(experimental_treatment_selected),
         se = sqrt(var(experimental_treatment_selected) / length(experimental_treatment_selected))
  ) %>% 
  ungroup() %>% 
  ggplot(aes(x = risk, y = m_left, color = condition)) +
  geom_pointrange(aes(ymin = m_left - se, ymax = m_left + se)) +
  theme_light() +
  geom_line(size = .6) +
  geom_hline(yintercept = .5, linetype = "dashed",size = .2) +
  scale_color_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(.10, .8)) +
  scale_x_continuous(breaks =  seq(10, 70, 10), limits = c(5,70)) +
  labs(
    y = "Probability of Choosing Experimental Treatment",
    x = "Probablitiy of Worsening",
    color = "Domain"
  )+  theme(text = element_text(size = 15)) +
  theme(legend.position = "none")

p2 <- TEST %>% 
  group_by(gain, condition) %>%
  mutate(n = n(),
         m_left = mean(experimental_treatment_selected),
         se = sqrt(var(experimental_treatment_selected) / length(experimental_treatment_selected))
  ) %>% 
  ungroup() %>% 
  ggplot(aes(x = gain, y = m_left, color = condition)) +
  geom_pointrange(aes(ymin = m_left - se, ymax = m_left + se)) +
  theme_light() +
  geom_line(size = .6) +
  geom_hline(yintercept = .5, linetype = "dashed",size = .2) +
  scale_color_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(.10, .8)) +
  scale_x_continuous(breaks =  seq(10, 70, 10), limits = c(5,70)) +
  labs(
    y = "Probability of Choosing Experimental Treatment",
    x = "Probablitiy of Curing",
    color = "Domain"
  )+  theme(text = element_text(size = 15)) +
  theme(legend.position = "none")

p3 <- TEST %>% 
  group_by(no_effect, condition) %>%
  mutate(n = n(),
         m_left = mean(experimental_treatment_selected),
         se = sqrt(var(experimental_treatment_selected) / length(experimental_treatment_selected))
  ) %>% 
  ungroup() %>% 
  ggplot(aes(x = no_effect, y = m_left, color = condition)) +
  geom_pointrange(aes(ymin = m_left - se, ymax = m_left + se)) +
  theme_light() +
  geom_line(size = .6) +
  geom_hline(yintercept = .5, linetype = "dashed",size = .2) +
  scale_color_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(.10, .8)) +
  scale_x_continuous(breaks =  seq(10, 70, 10), limits = c(5,70)) +
  labs(
    y = "Probability of Choosing Experimental Treatment",
    x = "Probablitiy of No Effect",
    color = "Domain"
  )+  theme(text = element_text(size = 15)) 

p1 + p2 + p3 +  patchwork::plot_annotation(tag_levels = "A")


#response times plot for manuscript (just looking at the bin versions0)
p1 <- TEST %>% 
  group_by(risk,experimental_treatment_selected) %>%
  mutate(n = n(),
         m_left = mean(rt),
         se = sqrt(var(rt) / length(rt))
  ) %>% 
  ungroup() %>% 
  ggplot(aes(x = risk, y = m_left, color = experimental_treatment_selected)) +
  geom_pointrange(aes(ymin = m_left - se, ymax = m_left + se)) +
  theme_light() +
  geom_line(size = .6) +
  scale_color_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(1000, 1800)) +
  scale_x_continuous(breaks =  seq(10, 70, 10), limits = c(5,70)) +
  labs(title = "",
    y = "RT(ms)",
    x = "Probablitiy of Worsening",
    color = "Domain"
  )+  theme(text = element_text(size = 15))

p2 <- TEST %>% 
  group_by(gain,experimental_treatment_selected) %>%
  mutate(n = n(),
         m_left = mean(rt),
         se = sqrt(var(rt) / length(rt))
  ) %>% 
  ungroup() %>% 
  ggplot(aes(x = gain, y = m_left, color = experimental_treatment_selected)) +
  geom_pointrange(aes(ymin = m_left - se, ymax = m_left + se)) +
  theme_light() +
  geom_line(size = .6) +
  scale_color_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(1000, 1800)) +
  labs(
    y = "RT(ms)",
    x = "Probablitiy of Curing",
    color = "Domain"
  )+  theme(text = element_text(size = 15)) 

TEST %>% 
  group_by(no_effect, experimental_treatment_selected) %>%
  mutate(n = n(),
         m_left = mean(rt),
         se = sqrt(var(rt) / length(rt))
  ) %>% 
  ungroup() %>% 
  ggplot(aes(x = no_effect, y = m_left, color = experimental_treatment_selected)) +
  geom_pointrange(aes(ymin = m_left - se, ymax = m_left + se)) +
  theme_light() +
  geom_line(size = .6) +
  scale_color_brewer(palette = "Set1") +
  labs(
    y = "RT(ms)",
    x = "Probablitiy of No Effect",
    color = "Domain"
  )+  theme(text = element_text(size = 15)) 

pd <- position_dodge(0.2)

TEST %>% 
  group_by(condition, experimental_treatment_selected) %>%
  mutate(n = n(),
         m_left = mean(rt),
         se = sqrt(var(rt) / length(rt))
  ) %>% 
  ungroup() %>% 
  ggplot(aes(x = experimental_treatment_selected, y = m_left, color = condition)) +
  geom_pointrange(aes(ymin = m_left - se, ymax = m_left + se), position = pd) +
  theme_light() +
  geom_line(size = .6) +
  scale_color_brewer(palette = "Set1") +
  labs(
    y = "RT(ms)",
    x = "Probablitiy of No Effect",
    color = "Domain"
  )+  theme(text = element_text(size = 15)) 

TEST %>% 
  filter(experimental_treatment_selected == 0) %>% 
  group_by(condition, user_id) %>% 
  summarise(median = median(rt)) %>% 
  ggstatsplot::ggbetweenstats(
  x = condition,
  y = median,
  type = "nonparametric",
  title = "Response times across cognitive domains",
  ylab = "RT(ms)",
  xlab = "Cognitve Domain",
  bf.message = F,
) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(ggstatsplot)
## let's split the data frame and create a list by condition of interest
cue_list <- TEST %>% 
  group_by(condition, user_id, experimental_treatment_selected) %>% 
  summarise(median = median(rt)) %>% 
  split(f = .$condition, drop = TRUE)
## creating a list of plots by applying the same function for elements of the list
plot_list <- purrr::pmap(
  .l = list(
    data = cue_list,
    x = "experimental_treatment_selected",
    y = "median",
    type = "nonparametric",
    bf.message = F,
    ylab = "Median RT(ms)",
    xlab = "Experimental treatment selected",
    ggstatsplot.layer = T),
  .f = ggwithinstats
)

## combining all individual plots from the list into a single plot using combine_plots function
combine_plots(
  plotlist = plot_list,
  annotation.args = list(title = "Response times across cognitve domains"),
  plotgrid.args = list(ncol = 2)
)
