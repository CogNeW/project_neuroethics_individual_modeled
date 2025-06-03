get_treatment_counts <- function(data){
  data %>% 
    mutate(treatment = ifelse(experimental_treatment_selected == "Yes to treatment", 1, 0)) %>% 
  dplyr::select(user_id, treatment) %>% 
    group_by(user_id) %>% 
    summarise(N = n(),
              Yes = sum(treatment),
              No = N - Yes
    )
}
