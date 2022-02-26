# Outplant Kelp Growth Code
# August 9, 2021
# By: Em Lim

#Load packages -----

library(tidyverse)
library(visreg)
library(ggplot2)
library(lme4)
library(lmerTest)
library(nlme)
library(PNWColors)
theme_set(theme_classic())
library(patchwork)

library(tidyverse)

source("Code/theme_black.R")
# why is this not working

# Load kelp growth data ------
kelp_pre <- read_csv("Data/Kelp_Outplant/2021_08_09_Alaria_pre_outplant.csv") %>%
  rename(pre_dist_hole_mm = dist_hole_mm,
         pre_total_length_cm = total_length_cm,
         pre_sporophylls = sporophylls,
         pre_repro = reproductive, 
         pre_notes = notes)
# piping to renames the data columns


kelp_mid <- read_csv("Data/Kelp_Outplant/2021_08_10_Alaria_mid_outplant.csv") %>%
  rename(mid_dist_hole_mm = dist_hole_mm) 

kelp_post <- read_csv("Data/Kelp_Outplant/2021_08_18_Alaria_post_outplant.csv") %>%
  rename(post_dist_hole_mm = dist_hole_mm,
         post_total_length_cm = total_length_cm,
         post_sporophylls = sporophylls,
         post_repro = reproductive, 
         post_notes = notes) 

kelp_nut_content <- read_csv("Data/Kelp_Outplant/2022_02_14_kelp_percent_N.csv")

#all above are just loading files with some renaming

# Join pre and mid and post
# then do a ton of manipulations
# standardize (i.e. center and divide by their standard deviation)
# standardizing means the coefficients will tell you about the "average" kelp ie the intercept lies at the mean kelp length and divorces effect sizes from units
# aka if you use mm instead of meters the effect size looks bigger
# Having our effect sizes listed in actual units (e.g. body size in cm) may be beneficial in some cases when you want to be able to make statements such as “for every cm increase in body size, we would expect to see an X increase in home range size). In this case, we don’t necessarily want to divide by the SD of variable, but we can still center using scale(). scale() has two arguments, center and scale, which by default are both set to true. If we use scale(scale = FALSE), we can now simply center our data around the mean value of the continuous variable(s).


#modifying the kelp1 data set. Taking kelp_pre and piping it through the following functions. Left_join adds columns from y to x, and includes all rows in the x
kelp1 <- kelp_pre %>%
  left_join(kelp_mid, by = c("site", "line", "kelp")) %>%
  left_join(kelp_post, by = c("site", "line", "kelp")) %>%
  mutate(first_change_mm = mid_dist_hole_mm - pre_dist_hole_mm, # growth in first week
         second_change_mm = post_dist_hole_mm - mid_dist_hole_mm,# growth in second week
         full_change_mm = post_dist_hole_mm - pre_dist_hole_mm, # total change
         post_total_length_mm = post_total_length_cm * 10, # change post total len to mm
         pre_total_length_mm = pre_total_length_cm * 10, # change pre total len to mm
         total_length_stand = scale(pre_total_length_mm), # standardize the pre length
         full_change_stand = scale(full_change_mm), # standardize the growth
         total_length_center = scale(pre_total_length_mm, scale = FALSE), # center the pre length around the mean
         
         full_change_center = scale(full_change_mm, scale = FALSE), # center the total growth around the mean
         pee_factor = ifelse(site == "Scotts", "high",
                             ifelse(site == "Ross", "high", "low")), # make high vs low pee
         pee_factor = factor(pee_factor, levels = c("low", "high")),
         full_growth_rate = full_change_mm / 14, # calculate a growth rate
         first_growth_rate = first_change_mm/ 7, # calculate growth in 1st week
         second_growth_rate = second_change_mm/ 7, # calculate growth in 2nd week
         change_length = post_total_length_mm - pre_total_length_mm, # how much did total length change
         percent_change = full_change_mm/pre_total_length_mm * 100) # percent change

# Add nutrients
nutrients <- read_csv("Data/Kelp_Outplant/2021_nutrients.csv") %>%
  mutate(NO3_uM = NO3_NO2_uM - NO2_uM)

nutrient_avg <- nutrients %>%
  group_by(site) %>%
  summarize(NO3_avg = mean(NO3_uM),
            NO2_avg = mean(NO2_uM),
            PO4_avg = mean(PO4_uM),
            SiO2_avg = mean(SiO2_uM)) 

# Create palettes ------
csee_pal <- pnw_palette("Starfish")
sailboat <- pnw_palette("Sailboat")

# Calculate NH4 for pre samples --------

# load data

#data with which bottles belong to which sites
ammonium_data <- read_csv("bottlesites_watersamples.csv")

# three samples per site, syringe 1, syringe 2, and a whirlpack (WP).
# The whirlpack provided a sample of water from which we drew two syringes from: one sample and one matrix. We added 200 uL of ammonium stock to the matrix sample and use the comparison between the WP sample and the WP matrix to calculate the matrix

# load fluorometry data
pre_fluorometry <- read_csv("Feb2022watersamples_fluorometry.csv") %>%
  mutate(mean_FLU = rowMeans(cbind(FLU1, FLU2, FLU3)))
# mean_FLU is the name of the new variable, the following is the modifying action youre doing. cbind means youre combining the columns. rowMeans?

# load standard curve data
standard_fluorometry <- read_csv("Feb2022watersamples_fluorometry_standards.csv") %>%
  mutate(mean_FLU = rowMeans(cbind(FLU1, FLU2, FLU3)))

#do the calculations for the standard curve
standard_fluorometry_f <- standard_fluorometry %>% 
  mutate(nh4_added_umol = nh4_vol_uL/1e6 * nh4_conc_og_umol, #amount of NH4
         total_vol_L = nh4_vol_uL/1e6 + og_vol_L, #new volume of sample + NH4
         nh4_conc_final_umol_L = nh4_added_umol / total_vol_L, #concentration
         #of NH4 in seawater sample
         mean_FLU = (FLU1 + FLU2 + FLU3)/3) #mean FLU reading

#graph standard curve. aes(x,y). dataframe, aes(x,y) then some more aestetic stuff
ggplot(standard_fluorometry_f, aes(mean_FLU, nh4_conc_final_umol_L)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE)

#model
# lm is used to fit linear models lm(y~x, data source)
sc_mod_aug_pre <- lm(nh4_conc_final_umol_L ~ mean_FLU, data = standard_aug_pre_f)
summary(sc_mod_aug_pre)

#save coefficients
# coef is a  function which extracts model coefficients from objects returned by modeling functions
int_aug_pre <- coef(sc_mod_aug_pre)[1]
slope_aug_pre <- coef(sc_mod_aug_pre)[2]

# calculate matrix as according to Taylor et al 2007
# % Matrix Effects = [(standard_spike - standard_zero) - (sample_spike - sample_zero)]/ (standard_spike - standard_zero) * 100

bottles_aug_pre2 <- bottles_aug_pre %>%
  left_join(glow_aug_pre, by = c("bottle")) %>%
  mutate(Fsm_zero = mean_FLU)

bottles_aug_pre_spike <- bottles_aug_pre2 %>%
  filter(sample == "matrix") %>%
  transmute(site_ID = site_ID,
            Fsm_spike = mean_FLU)

bottles_aug_pre_zero <- bottles_aug_pre2 %>%
  filter(syringe == "WP") %>%
  filter(sample != "matrix")

Fst_zero_aug_pre <- standard_aug_pre_f$mean_FLU[standard_aug_pre_f$nh4_vol_uL == "0"]
Fst_spike_aug_pre <- standard_aug_pre_f$mean_FLU[standard_aug_pre_f$nh4_vol_uL == "200"]

bottles_aug_pre_matrix <- bottles_aug_pre_zero %>% 
  left_join(bottles_aug_pre_spike, by = "site_ID") %>%
  mutate(
    Fst_zero = Fst_zero_aug_pre,
    Fst_spike = Fst_spike_aug_pre,
    matrix = 100 * ((Fst_spike_aug_pre - Fst_zero_aug_pre - (Fsm_spike - Fsm_zero))/
                      (Fst_spike_aug_pre - Fst_zero_aug_pre)),
    matrix_est = 20) %>%
  select(site_ID, Fsm_spike, Fst_zero, Fst_spike, matrix, matrix_est)

bottles_aug_pre3 <- bottles_aug_pre2 %>%
  left_join(bottles_aug_pre_matrix, by = "site_ID") %>%
  mutate(
    Fsm_cor_Holmes = Fsm_zero + (Fsm_zero * (matrix/100)),
    Fsm_cor_Taylor = (Fsm_zero / (1- matrix/100))
  ) %>%
  mutate(int = int_aug_pre, #include values for the int and slope in for every column
         slope = slope_aug_pre) %>% #to calculate the coversion to NH4 conc
  #those values come from our standard curve
  mutate(nh4_conc_pre = int + slope * Fsm_cor_Taylor) %>%
  filter(sample != "matrix") %>%
  select(-sample_matrix)
# you could stop here Kiara

nh4_pre <- bottles_aug_pre3 %>%
  select(c(site, depth, temp_est, sal_est, nh4_conc_pre)) %>%
  rename(depth_pre = depth,
         temp_pre = temp_est, 
         sal_pre = sal_est) %>%
  filter(site != "Taylor")

nh4_pre_summary <- nh4_pre %>%
  group_by(site) %>%
  summarize(nh4_conc_pre = mean(nh4_conc_pre)) %>%
  arrange(desc(nh4_conc_pre)) %>%
  filter(site != "Taylor")

# Calculate NH4 for post samples --------

# load data

#data with which bottles belong to which sites
bottles_aug_post <- read_csv("Data/Kelp_Outplant/2021_08_23_kelp_post_NH4_bottles.csv")

# load fluorometry data
glow_aug_post <- read_csv("Data/Kelp_Outplant/2021_08_23_kelp_post_fluorometry.csv") %>%
  mutate(mean_FLU = rowMeans(cbind(FLU1, FLU2, FLU3)))

# load standard curve data
standard_aug_post <- read_csv("Data/Kelp_Outplant/2021_08_23_kelp_post_standard.csv") %>%
  mutate(mean_FLU = rowMeans(cbind(FLU1, FLU2, FLU3)))

#do the calculations for the standard curve
standard_aug_post_f <- standard_aug_post %>% 
  mutate(nh4_added_umol = nh4_vol_uL/1e6 * nh4_conc_og_umol, #amount of NH4
         total_vol_L = nh4_vol_uL/1e6 + og_vol_L, #new volume of sample + NH4
         nh4_conc_final_umol_L = nh4_added_umol / total_vol_L, #concentration
         #of NH4 in seawater sample
         mean_FLU = (FLU1 + FLU2 + FLU3)/3) #mean FLU reading

#graph standard curve
ggplot(standard_aug_post_f, aes(mean_FLU, nh4_conc_final_umol_L)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE)

#model
sc_mod_aug_post <- lm(nh4_conc_final_umol_L ~ mean_FLU, data = standard_aug_post_f)
summary(sc_mod_aug_post)

#save coefficients
int_aug_post <- coef(sc_mod_aug_post)[1]
slope_aug_post <- coef(sc_mod_aug_post)[2]

# calculate matrix
bottles_aug_post2 <- bottles_aug_post %>%
  left_join(glow_aug_post, by = c("bottle")) %>%
  mutate(Fsm_zero = mean_FLU)

bottles_aug_post_spike <- bottles_aug_post2 %>%
  filter(sample == "matrix") %>%
  transmute(site_ID = site_ID,
            Fsm_spike = mean_FLU)

bottles_aug_post_zero <- bottles_aug_post2 %>%
  filter(syringe == "WP") %>%
  filter(sample != "matrix")

Fst_zero_aug_post <- standard_aug_post_f$mean_FLU[standard_aug_post_f$nh4_vol_uL == "0"]
Fst_spike_aug_post <- standard_aug_post_f$mean_FLU[standard_aug_post_f$nh4_vol_uL == "200"]

bottles_aug_post_matrix <- bottles_aug_post_zero %>% 
  left_join(bottles_aug_post_spike, by = "site_ID")%>%
  mutate(
    Fst_zero = Fst_zero_aug_post,
    Fst_spike = Fst_spike_aug_post,
    matrix = 100 * ((Fst_spike_aug_post - Fst_zero_aug_post - (Fsm_spike - Fsm_zero))/
                      (Fst_spike_aug_post - Fst_zero_aug_post)),
    matrix_est = 20) %>%
  select(site_ID, Fsm_spike, Fst_zero, Fst_spike, matrix, matrix_est)

bottles_aug_post3 <- bottles_aug_post2 %>%
  left_join(bottles_aug_post_matrix, by = "site_ID") %>%
  mutate(
    Fsm_cor_Holmes = Fsm_zero + (Fsm_zero * (matrix/100)),
    Fsm_cor_Taylor = (Fsm_zero / (1- matrix/100))
  ) %>%
  mutate(int = int_aug_post, #include values for the int and slope in for every column
         slope = slope_aug_post) %>% #to calculate the coversion to NH4 conc
  #those values come from our standard curve
  mutate(nh4_conc_post = int + slope * Fsm_cor_Taylor) %>%
  filter(sample != "matrix") %>%
  select(-sample_matrix)

nh4_post <- bottles_aug_post3 %>%
  select(c(site, depth, temp_est, sal_est, nh4_conc_post)) %>%
  rename(depth_post = depth,
         temp_post = temp_est, 
         sal_post = sal_est) %>%
  filter(site != "DixonInside")

nh4_post_summary <- nh4_post %>%
  group_by(site) %>%
  summarize(nh4_conc_post = mean(nh4_conc_post)) %>%
  arrange(desc(nh4_conc_post)) %>%
  filter(site != "DixonInside")

# look at unfiltered pee ------------------
post_pee <- bottles_aug_post3 %>% 
  rename(nh4_conc = nh4_conc_post)

unfiltered_pee <- bottles_aug_pre3 %>%
  rename(nh4_conc = nh4_conc_pre) %>%
  rbind(post_pee) %>% 
  select(c(site, site_ID, date, depth, temp_est, sal_est, matrix, nh4_conc)) %>% 
  mutate(period = ifelse(date == "2021-08-05", "one",
                         ifelse(date == "2021-08-06", "one",
                                "two")))

ggplot(unfiltered_pee, aes(period, nh4_conc, colour = date)) +
  geom_point() +
  geom_boxplot(aes(group = date)) +
  facet_grid(~site) +
  theme(legend.position = "none")

# Add the filtered NH4+ conc to the kelp data -------
kelp <- kelp1 %>%
  left_join(nh4_pre_summary, by = "site") %>%
  left_join(nh4_post_summary, by = "site") %>%
  mutate(nh4_avg = ifelse(site == "Taylor", nh4_conc_post,
                          ifelse(site == "DixonInside", nh4_conc_pre,
                                 (nh4_conc_post + nh4_conc_pre)/2))) %>%
  left_join(nutrient_avg, by = "site") %>%
  left_join(kelp_nut_content, by = c("site", "line",  "kelp")) %>%
  mutate(site = ifelse(site == "DixonInside", "Dixon", site))

# Just the nutrients
nutrients_all <- kelp %>%
  group_by(site) %>%
  select(site, nh4_avg) %>%
  unique() %>%
  left_join(nutrient_avg, by = "site")

#write_csv(nutrients_all, path = "Output/Output_data/nutrients_summary.csv")

# merge the nh4 data so we can look at it
nh4_pre1 <- nh4_pre %>%
  transmute(site = site,
            nh4_conc = nh4_conc_pre,
            period = "first")

nh4_post2 <- nh4_post %>%
  transmute(site = site,
            nh4_conc = nh4_conc_post,
            period= "second")

nh4_avg <- rbind(nh4_pre1, nh4_post2)

# Stats -------

# model with pee as factor and site as random factor
kelp_model <- lmer(full_change_stand ~ pee_factor * total_length_stand + (1|site), kelp)
summary(kelp_model)
visreg(kelp_model)

# model with pee as factor and site as random factor, but with percent change as response
# Worst with AIC

kelp_model2 <- lmer(percent_change ~ pee_factor + (1|site), kelp)
summary(kelp_model2)

# model with pee as a continuous variable
# I got rid of the abnormal nh4+ readings and averaged the rest

kelp_model3 <- lmer(full_change_stand ~ nh4_avg * total_length_stand + (1|site), kelp)
summary(kelp_model3) 
# singular fit is a prob, could fix with bayes baes
# could talk more to Hannah about Bayes
# could also try lme instead of lmer

# using lme instead of lmer gets rid of the singular fit error!
kelp_model4 <- lme(full_change_stand ~ nh4_avg * total_length_stand, random = ~ 1|site, kelp)
summary(kelp_model4) 
# its ok that i have a random effect of site even though there's only one avg per site
# we expect the sites to be affected by the ammonium, and we expect kelp to be affected by site
# AIC says this is the best model

# try with centered measurements (not divided by SD to keep original units)
kelp_model5 <- lme(full_change_center ~ nh4_avg * total_length_center, random = ~ 1|site, kelp)
summary(kelp_model5) 

# untransformed response, centred predictor
kelp_model6 <- lme(full_change_mm ~ nh4_avg * total_length_center, random = ~ 1|site, kelp)
summary(kelp_model6) 

AIC(kelp_model4, kelp_model5, kelp_model6)


# Ok I tried adding NO3 to the model and it wasn't significant
# And I tried swapping NO3 for nh4 and NO3 was still not significant!!!!!
# SO THE TRENDS ARE DRIVEN BY AMMONIUM! NOT NITRATE!

visreg(kelp_model4, "total_length_stand", by = "nh4_avg")

par(mfrow=c(2,2))
plot(kelp_model4)


AIC(kelp_model, kelp_model2, kelp_model3)


# Check effect of pee on %N !
kelp_na <- kelp %>% 
  drop_na(percent_N)

model_N_content <- lme(percent_N ~ nh4_avg + total_length_stand, random = ~ 1|site, kelp_na)
summary(model_N_content)
#ok so site ammonium concentration affected how much N there was in the kelps

model_N_content2 <- lme(full_change_stand ~ percent_N, random = ~ 1|site, kelp_na)
summary(model_N_content2) 
# And how much N there was in the kelps affected growth

# Graphing ---------

# pee vs site
ggplot(nh4_avg, aes(site, nh4_conc, colour = site)) +
  geom_boxplot() +
  geom_point()

# percent change vs pee
ggplot(kelp, aes(pee_factor, percent_change, colour = pee_factor)) +
  geom_boxplot()

# dot and line with continuous pee vs growth rate
ggplot(kelp, aes(nh4_avg, full_growth_rate)) +
  geom_boxplot(aes(fill = site), size = 1, colour = "white", alpha = 0.85) +
  geom_point(aes(colour = site), size = 2, alpha = 1) +
  geom_smooth(method = lm, colour = "white") +
  theme_black() +
  scale_fill_manual(values = c("#7bbcd5", "#d0e2af", "#f5db99", "#6e7cb9")) +
  scale_colour_manual(values = c("#7bbcd5", "#d0e2af", "#f5db99", "#6e7cb9")) +
  labs(x = "Ammonium Concentration (umol/L)", y = "Growth (mm/day)", fill = "Site", colour = "Site")

#ggsave("Output/Figures/Outplant_kelp_growth1.png", device = "png",
#       height = 9, width = 16, dpi = 400)

# boxplot with pee as a factor vs change standardized
ggplot(kelp, aes(pee_factor, full_change_stand, colour = pee_factor)) +
  geom_boxplot()

# dot and line wiht total length vs change standardized
ggplot(kelp, aes(total_length_stand, full_change_stand, colour = pee_factor)) +
  geom_point(size = 2) +
  geom_smooth(method = lm, size = 1.5) +
  labs(x = "Total length (standardized)", y = "Growth (standardized)", colour = "Ammonium") +
  theme_black() +
  scale_colour_manual(values = sailboat)

#ggsave("Output/Figures/Outplant_kelp_growth.png", device = "png",
#       height = 9, width = 16, dpi = 400)

# compare growth in first week to growth in second week
first_half <- ggplot(kelp, aes(pee_factor, first_growth_rate) ) +
  geom_boxplot() +
  ylim(c(-5, 20))

second_half <- ggplot(kelp, aes(pee_factor, second_growth_rate) ) +
  geom_boxplot()+
  ylim(c(-5, 20))

first_half + second_half


# How much much did the total length change?
ggplot(kelp, aes(full_change_mm, change_length, colour = pee_factor)) +
  geom_point() +
  geom_smooth(method = lm)
# Weird how many of them decreased in length

# percent change vs total length
ggplot(kelp, aes(post_total_length_mm, percent_change, colour = pee_factor)) +
  geom_point() +
  geom_smooth(method = lm)
# slopes are same but that's because the interaction was eaten by the percent change


# Nutrient content!
ggplot(kelp, aes(nh4_avg, percent_N)) +
  geom_boxplot(aes(fill = site), size = 1, colour = "white", alpha = 0.85) +
  geom_point(aes(colour = site), size = 2, alpha = 1) +
  geom_smooth(method = lm, colour = "white") +
  theme_black() +
  scale_fill_manual(values = c("#7bbcd5", "#d0e2af", "#f5db99", "#6e7cb9")) +
  scale_colour_manual(values = c("#7bbcd5", "#d0e2af", "#f5db99", "#6e7cb9")) +
  labs(x = "Ammonium Concentration (umol/L)", y = "% Total N", fill = "Site", colour = "Site")

#ggsave("Output/Figures/Outplant_nh4_vs_%N.png", device = "png",
#       height = 9, width = 16, dpi = 400)

# N content vs growth
ggplot(kelp, aes(percent_N, full_growth_rate)) +
  geom_point(aes(colour = site, pch = site), size = 4) +
  geom_smooth(method = lm) +
  scale_colour_manual(values = c("#7bbcd5", "#d0e2af", "#f5db99", "#6e7cb9")) +
  theme_black() +
  labs(x = "% Total Nitrogen", y = "Growth (mm/day)", fill = "Site", colour = "Site", pch = "Site")

#ggsave("Output/Figures/Outplant_%N_vs_growth.png", device = "png",
#       height = 9, width = 16, dpi = 400)

