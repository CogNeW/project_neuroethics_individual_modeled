library(tidyverse)
raw <- read_csv("~/Dropbox/Medaglia_RiskBenefit_Neuroethics/DATA_ANALYSIS/data/raw_mturk_data/raw.csv")

dataset <- raw %>% 
  select(user_id, rt, experimental_treatment_selected,trial_index, condition, risk, gain, no_effect, Optimize)
  
exlusions_data <- dataset %>% 
  group_by(user_id) %>%
  mutate(
    Q1 = quantile(rt, .25, na.rm = T),
    Q3 = quantile(rt, .75, na.rm = T),
    IQR = IQR(rt, na.rm = T)
  ) %>%
  filter(rt > (Q1 - 2.5 * IQR) & rt < (Q3 + 2.5 * IQR)) %>%
  ungroup() %>%
  filter(!rt <= 300) %>% # response times cutoffs (think about)
  filter(!rt >= 10000) %>%
  # filter(!rt >= 5000) %>%
  mutate(rt = rt/1000,
         response = as.integer(experimental_treatment_selected),
         domain = as.integer(factor(condition)),
         subj_idx = as.integer(factor(user_id)))

write.csv(exlusions_data, here::here("data", "HDDM_data.csv"),row.names=FALSE)

exlusions_data %>% 
  group_by(risk, condition) %>%
  mutate(n = n(),
         m_left = mean(experimental_treatment_selected),
         se = sqrt(var(experimental_treatment_selected) / length(experimental_treatment_selected))
  ) %>% 
  ungroup() %>% 
  ggplot(aes(x = risk, y = m_left, color = condition)) +
  geom_pointrange(aes(ymin = m_left - se, ymax = m_left + se)) +
  theme_classic() +
  geom_line(size = .6) +
  geom_hline(yintercept = .5, linetype = "dashed",size = .2) +
  scale_color_brewer(palette = "Set1") +
  # scale_y_continuous(limits = c(.10, .8)) +
  # scale_x_continuous(breaks =  seq(10, 70, 10), limits = c(5,70)) +
  labs(
    y = "Probability of Choosing Experimental Treatment",
    x = "Probablitiy of Worsening",
    color = "Domain"
  )+  theme(text = element_text(size = 15)) +
  theme(legend.position = "none")

exlusions_data %>% 
  group_by(gain, condition) %>%
  mutate(n = n(),
         m_left = mean(experimental_treatment_selected),
         se = sqrt(var(experimental_treatment_selected) / length(experimental_treatment_selected))
  ) %>% 
  ungroup() %>% 
  ggplot(aes(x = gain, y = m_left, color = condition)) +
  geom_pointrange(aes(ymin = m_left - se, ymax = m_left + se)) +
  theme_classic() +
  geom_line(size = .6) +
  geom_hline(yintercept = .5, linetype = "dashed",size = .2) +
  scale_color_brewer(palette = "Set1") +
  # scale_y_continuous(limits = c(.10, .8)) +
  # scale_x_continuous(breaks =  seq(10, 70, 10), limits = c(5,70)) +
  labs(
    y = "Probability of Choosing Experimental Treatment",
    x = "Probablitiy of Curing",
    color = "Domain"
  )+  theme(text = element_text(size = 15)) +
  theme(legend.position = "none")

exlusions_data %>% 
  group_by(risk,experimental_treatment_selected) %>%
  mutate(n = n(),
         m_left = mean(rt),
         se = sqrt(var(rt) / length(rt))
  ) %>% 
  ungroup() %>% 
  ggplot(aes(x = risk, y = m_left, color = experimental_treatment_selected)) +
  geom_pointrange(aes(ymin = m_left - se, ymax = m_left + se)) +
  theme_classic() +
  geom_line(size = .6) +
  scale_color_brewer(palette = "Set1") +
  # scale_y_continuous(limits = c(1000, 1800)) +
  # scale_x_continuous(breaks =  seq(10, 70, 10), limits = c(5,70)) +
  labs(title = "",
       y = "RT(ms)",
       x = "Probablitiy of Worsening",
       color = "experimental treatment selected"
  )+  theme(text = element_text(size = 15))

exlusions_data %>% 
  group_by(gain,experimental_treatment_selected) %>%
  mutate(n = n(),
         m_left = mean(rt),
         se = sqrt(var(rt) / length(rt))
  ) %>% 
  ungroup() %>% 
  ggplot(aes(x = gain, y = m_left, color = experimental_treatment_selected)) +
  geom_pointrange(aes(ymin = m_left - se, ymax = m_left + se)) +
  theme_classic() +
  geom_line(size = .6) +
  scale_color_brewer(palette = "Set1") +
  # scale_y_continuous(limits = c(1000, 1800)) +
  labs(
    y = "RT(ms)",
    x = "Probablitiy of Curing",
    color = "experimental treatment selected"
  )+  theme(text = element_text(size = 15)) 


exlusions_data %>% 
  group_by(risk,condition,experimental_treatment_selected) %>%
  mutate(n = n(),
         m_left = mean(rt),
         se = sqrt(var(rt) / length(rt))
  ) %>% 
  ungroup() %>% 
  ggplot(aes(x = risk, y = m_left, color = condition)) +
  geom_pointrange(aes(ymin = m_left - se, ymax = m_left + se)) +
  theme_classic() +
  geom_line(size = .6) +
  scale_color_brewer(palette = "Set1") +
  # scale_y_continuous(limits = c(1000, 1800)) +
  # scale_x_continuous(breaks =  seq(10, 70, 10), limits = c(5,70)) +
  labs(title = "",
       y = "RT(ms)",
       x = "Probablitiy of Worsening",
       color = "domain"
  )+  theme(text = element_text(size = 15)) +
  facet_wrap(~experimental_treatment_selected)

exlusions_data %>% 
  group_by(gain,condition,experimental_treatment_selected) %>%
  mutate(n = n(),
         m_left = mean(rt),
         se = sqrt(var(rt) / length(rt))
  ) %>% 
  ungroup() %>% 
  ggplot(aes(x = gain, y = m_left, color = condition)) +
  geom_pointrange(aes(ymin = m_left - se, ymax = m_left + se)) +
  theme_classic() +
  geom_line(size = .6) +
  scale_color_brewer(palette = "Set1") +
  labs(
    y = "RT(ms)",
    x = "Probablitiy of Curing",
    color = "domain"
  )+  theme(text = element_text(size = 15))+
  facet_wrap(~experimental_treatment_selected)

exlusions_data %>% 
  group_by(risk,gain,condition,experimental_treatment_selected) %>%
  mutate(n = n(),
         m_left = mean(rt),
         se = sqrt(var(rt) / length(rt))
  ) %>% 
  ungroup() %>% 
  ggplot(aes(x = gain, y = risk, color = m_left)) +
  theme_classic() +
  geom_density_2d_filled(contour_var = "ndensity") +
  labs(
    y = "Worsening",
    x = "Curing",
    fill = "RT(ms)"
  )+  theme(text = element_text(size = 15))+
  facet_wrap(~condition+ experimental_treatment_selected)


#Scoring = Based on answers to 1-4 questions following
#Correct answers are as follows: 
#1 = 25; 2 = 30; 3 =20; 4 = 50.
Berlin <- raw %>% 
  select(user_id, Berlin1, Berlin2, Berlin3, Berlin4) %>% 
  distinct() %>%
  mutate(Berlin1 = Berlin1 == 25,
         Berlin2 = Berlin2 == 30, 
         Berlin3 = Berlin3 == 20, 
         Berlin4 = Berlin4 == 50, 
         BerlinScore = Berlin1 + Berlin2+ Berlin3 + Berlin4)
  
MC <- raw %>% 
  select(user_id, MC_Lying, MC_Assassin, MC_Torture, MC_Murder, MC_Stealing, MC_Steril) %>% 
  distinct() %>%
  group_by(user_id) %>% 
  mutate(count = sum(MC_Lying, MC_Assassin, MC_Torture, MC_Murder,MC_Stealing, MC_Steril))%>% 
  mutate(consequentialism_score = (count/6)) %>% 
  mutate(moral_commitment = if_else(consequentialism_score > .5,"Consequentialist", "Deontological"))
  # mutate(moral_commitment = case_when(
  #   consequentialism_score <  .5 ~ "Deontological",
  #   consequentialism_score >= 1 ~ "Strong Consequentialist",
  #   consequentialism_score >= .5 ~ "Consequentialist"
  # ))

MC$moral_commitment <- as.factor(MC$moral_commitment)
MC$moral_commitment <- relevel(MC$moral_commitment, "Deontological")

MC %>% 
  select(user_id, consequentialism_score, moral_commitment) %>% 
  distinct() %>% 
  group_by(moral_commitment) %>% 
  summarise(n = n())

MC %>% 
  select(user_id, consequentialism_score) %>% 
  distinct() %>% 
  ggplot(aes(consequentialism_score)) +
  geom_histogram(bins = 13)+
  geom_vline(xintercept = .5)

#look at the relationship across consequentialism score
test <- left_join(exlusions_data, MC, by = "user_id")
test <-  left_join(test, Berlin, by = "user_id")
test %>% 
  na.omit() %>% 
  group_by(risk, moral_commitment,experimental_treatment_selected) %>%
  mutate(n = n(),
         m_left = mean(rt),
         se = sqrt(var(rt) / length(rt))
  ) %>% 
  ungroup() %>%
  ggplot(aes(x = risk, y = m_left, color = moral_commitment)) +
  geom_pointrange(aes(ymin = m_left - se, ymax = m_left + se))+
  theme_classic() +
  geom_line(size = .6) +
  scale_color_brewer(palette = "Set2") +
  labs(
    y = "RT(ms)",
    x = "Probablitiy of Worsening",
    color = "moral commitment"
  )+  theme(text = element_text(size = 15))+
  facet_wrap(~experimental_treatment_selected)

test %>% 
  na.omit() %>% 
  group_by(gain, moral_commitment,experimental_treatment_selected) %>%
  mutate(n = n(),
         m_left = mean(rt),
         se = sqrt(var(rt) / length(rt))
  ) %>% 
  ungroup() %>%
  ggplot(aes(x = gain, y = m_left, color = moral_commitment)) +
  geom_pointrange(aes(ymin = m_left - se, ymax = m_left + se))+
  theme_classic() +
  geom_line(size = .6) +
  scale_color_brewer(palette = "Set2") +
  labs(
    y = "RT(ms)",
    x = "Probablitiy of Curing",
    color = "moral commitment"
  )+  theme(text = element_text(size = 15))+
  facet_wrap(~experimental_treatment_selected)

test %>% 
  na.omit() %>% 
  group_by(risk, moral_commitment) %>%
  mutate(n = n(),
         m_left = mean(experimental_treatment_selected),
         se = sqrt(var(experimental_treatment_selected) / length(experimental_treatment_selected))
  ) %>% 
  ungroup() %>% 
  ggplot(aes(x = risk, y = m_left, color = moral_commitment)) +
  geom_pointrange(aes(ymin = m_left - se, ymax = m_left + se)) +
  theme_classic() +
  geom_line(size = .6) +
  geom_hline(yintercept = .5, linetype = "dashed",size = .2) +
  scale_color_brewer(palette = "Set2") +
  labs(
    y = "Probability of Choosing Experimental Treatment",
    x = "Probablitiy of Worsening",
    color = "Moral Commitment"
  )+  theme(text = element_text(size = 15))

test %>% 
  na.omit() %>% 
  group_by(gain, moral_commitment) %>%
  mutate(n = n(),
         m_left = mean(experimental_treatment_selected),
         se = sqrt(var(experimental_treatment_selected) / length(experimental_treatment_selected))
  ) %>% 
  ungroup() %>% 
  ggplot(aes(x = gain, y = m_left, color = moral_commitment)) +
  geom_pointrange(aes(ymin = m_left - se, ymax = m_left + se)) +
  theme_classic() +
  geom_line(size = .6) +
  geom_hline(yintercept = .5, linetype = "dashed",size = .2) +
  scale_color_brewer(palette = "Set2") +
  labs(
    y = "Probability of Choosing Experimental Treatment",
    x = "Probablitiy of Curing",
    color = "Moral Commitment"
  )+  theme(text = element_text(size = 15))


test %>%
  na.omit() %>%
  group_by(risk,moral_commitment) %>%
  mutate(n = n(),
         m_left = mean(experimental_treatment_selected),
         se = sqrt(var(experimental_treatment_selected) / length(experimental_treatment_selected))
  ) %>%
  ungroup() %>%
  ggplot(aes(x = risk, y = m_left, color = moral_commitment)) +
  geom_pointrange(aes(ymin = m_left - se, ymax = m_left + se)) +
  theme_classic() +
  geom_line(size = .6) +
  geom_hline(yintercept = .5, linetype = "dashed",size = .2) +
  scale_color_brewer(palette = "Set2") +
  labs(
    y = "Probability of Choosing Experimental Treatment",
    x = "Probablitiy of Cure",
    color = "Moral Commitment"
  )+  theme(text = element_text(size = 15))


test %>%
  na.omit() %>%
  group_by(user_id,gain,experimental_treatment_selected) %>%
  mutate(
    q1 = quantile(rt, .1),
    q3 = quantile(rt, .3),
    q5 = quantile(rt, .5),
    q7 = quantile(rt, .7),
    q9 = quantile(rt, .9),
  ) %>%
  pivot_longer(cols = q1:q9,
               names_to = "quantiles",
               values_to = "rts"
  ) %>% ungroup() %>% 
  group_by(quantiles,risk,experimental_treatment_selected) %>%
  mutate(n = n(),
         m_rt = mean(rts),
         se = sqrt(var(rts) / length(rts))
  ) %>%
  ggplot(aes(x = risk, y = m_rt, color = factor(quantiles), shape =factor(experimental_treatment_selected), linetype = factor(experimental_treatment_selected))) +
  geom_pointrange(aes(ymin = m_rt - se, ymax = m_rt + se)) +
  theme_classic() +
  geom_line(size = 1) +
  scale_color_brewer(palette = "Set1") +
  labs(
    y = "RT(s)",
    x = "risk"
  ) + theme(legend.position="top")


# test %>% 
#   select(user_id, age) %>% 
#   distinct() %>% 
#   group_by(age) %>% 
#   summarise(n= n()) %>% 
#   ggplot(aes(age, n))+
#   geom_col()
