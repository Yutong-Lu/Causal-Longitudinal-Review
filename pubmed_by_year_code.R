library(tidyverse)
library(readxl)
library(ggsci)

results <- read_excel("pubmed_by_year.xlsx")

results %>% 
  filter(year %in% 2000:2021) %>% 
  mutate(search = str_to_sentence(search)) %>% 
  ggplot(aes(x = year, y = count, fill = search))+
  geom_col(position = position_dodge2())+
  theme_bw(base_size = 12)+
  labs(y = "Number of publications in PubMed",
       x = "Year",
       fill = "Type of method")+
  scale_fill_aaas()+
  theme(legend.position = "bottom")

ggsave(filename = "pubmed_results.png",
       dpi = 300, units = "cm", height = 15,
       width = 20)