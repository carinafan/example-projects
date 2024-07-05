# Set up ----

# SAM bins based on region of significance

df.complete$sam_epi_group = NA

for (i in 1:nrow(df.complete)) {
  if (df.complete$sam_epi[i] < 78.1619) {
    df.complete$sam_epi_group[i] = "low"
  } 
  else if (df.complete$sam_epi[i] > 104.6184) {
    df.complete$sam_epi_group[i] = "high"
  }
  else {
    df.complete$sam_epi_group[i] = "mid"
  }
}

df.complete %>%
  group_by(sam_epi_group) %>%
  summarise(n = n())

# split data for plotting

df.samEnds = df.complete %>% 
  filter(sam_epi_group == "low" |
           sam_epi_group == "high")

df.samMid = df.complete %>% 
  filter(sam_epi_group == "mid")

nrow(df.samEnds) + nrow(df.samMid) # make sure n adds up correctly


#---- Figures in main text ----


# age -> face name

ggplot(df.fn.long, aes(x = age, y = fn_score, colour = fn_domain)) +
  geom_point(alpha = .3) +
  geom_smooth(method = "lm") + 
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 12)) +
  scale_colour_manual(name = "Recognition type", 
                      labels = c("associative", "item"),
                      values = c("darkcyan", "plum")) +
  labs(x = "Age", y = "Memory score")

ggsave("figures/fig1_age_fn.png", dpi = 600, width = 5, height = 5, units = "in")

# age x SAM-E -> CFQ

ggplot() +
  geom_point(data = df.samMid, 
             aes(x = age, y = cfq_total), 
             alpha = .2) +
  geom_point(data = df.samEnds, 
             aes(x = age, y = cfq_total, colour = sam_epi_group),
             alpha = .6) +
  geom_smooth(data = df.samEnds,
              aes(x = age, y = cfq_total, colour = sam_epi_group),
              method = "lm") +
  scale_colour_manual(name = "SAM-episodic\ngroup", 
                      values = c("indianred1", "blue3")) +
  geom_rect(aes(xmin = 57.76682, xmax = 73.04918, ymin = 0, ymax = Inf),
            fill = "whitesmoke",
            alpha = .6) +
  geom_vline(xintercept = 57.76682, linetype = "dashed") +
  geom_vline(xintercept = 73.04918, linetype = "dashed") +
  labs(x = "Age", y = "CFQ score (out of 100)") +
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 12))

ggsave("figures/fig2_age_samE_cfq.png", dpi = 600, width = 5, height = 5, units = "in")


#---- Supplemental figures ----

# figure S1: age x gender interaction in predicting CFQ

df.complete %>% 
  filter(gender != "Other/Prefer not to answer") %>% 
  ggplot(aes(x = age, y = cfq_total, colour = gender)) +
  geom_jitter(alpha = .2) +
  geom_smooth(method = "lm") +
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = 12)) +
  scale_colour_manual(name = "Gender", 
                      values = c("lightpink3", "palegreen3")) +
  labs(x = "Age", y = "CFQ score (out of 100)")

ggsave("figures/figS1_age_cfq_gender.png", 
       dpi = 300, width = 5, height = 5, units = "in")
