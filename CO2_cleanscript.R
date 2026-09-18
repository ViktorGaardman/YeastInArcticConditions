library(tidyverse)
library(car)
library(emmeans)
library(DHARMa)
library(patchwork)
library(glmmTMB)

rm(list=ls())

options(contrasts = c("contr.sum", "contr.poly"))

set.seed(123)

#1. CO2 release rates
CO2_min <- read_csv("CO2_per_minute.csv")

#Add dummy 0s at 8:00 (start of yeast-kits)
CO2_df <- CO2_min |> 
  mutate(
    Year = str_sub(`EVENT DATE`, 7,8),
    Month = str_sub(`EVENT DATE`, 4,5),
    Day = str_sub(`EVENT DATE`, 1, 2)
  )

#Change faulty dates from raw data
CO2_df <- CO2_df |> 
  mutate(
    Day = fct_recode(
      Day,
      "03" = "13",
      "15" = "25",
      "31" = "10",
      
    ))

CO2_filt <- CO2_df |> 
  select(c(Hour, mlpermin, Day))

ReleaseRates <- rbind(CO2_filt, data.frame(Hour = 8, mlpermin = 0, Day = c("03", "15", "31")))

#Turn time of day to hours since start of fermentation
ReleaseRates <- ReleaseRates %>%
  mutate(
    HourStart = Hour - 8
  )

#Create a linear summary of release rates
CO2_mean <- ReleaseRates |>
  group_by(HourStart) |> 
  summarize(
    mean_rate = mean(mlpermin),
    sd_rate = sd(mlpermin)
  )

#Plot
ReleaseRate_plot <- ggplot(ReleaseRates, aes(x = HourStart, y = mlpermin)) + 
  theme_bw()+
  xlab("Hours of fermentation")+
  ylab("ml CO2 / min")+
  geom_point(color = "grey20", size = 2, alpha = 0.7)+
  geom_hline(yintercept = 185, col = "black", linetype = "dashed", lwd = 0.5) +
  geom_line(
    data = CO2_mean,
    aes(y = mean_rate,
        x = HourStart),
    linewidth = 0.8
  ) +
  geom_errorbar(
    data = CO2_mean,
    aes(
      x = HourStart,
      ymax = mean_rate + sd_rate,
      ymin = mean_rate - sd_rate),
    width = 0.2,
    linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  scale_x_continuous(limits = c(0,14), n.breaks = 8)+
  theme(axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text = element_text(size = 14),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)
  )

ReleaseRate_plot


#2. Comparison of attractants

df_attractant <- read.csv("TrappingCombinations.csv", sep =";", h = T)

#Turn total counts into capture rates
df_attractant <- df_attractant %>%
  filter(!is.na(Family)) %>%
  group_by(ID, Family, Species, TrapType, Date_Num) %>%
  summarize(
    CountPerHour = Count / Hours_sampling
  )

df_attractant$TrapType <- as_factor(df_attractant$TrapType)

#Rename species
df_attractant <- df_attractant %>% 
  mutate(
    Species = fct_recode(
      Species,
      "Oc. impiger" = "Oc_impiger",
      "Oc. nigripes" = "Oc_nigripes",
      "S. vittatum" = "S_vittatum",
      "S. rostratum" = "S_rostratum"
    )) %>% 
  mutate(Species = fct_relevel(
    Species,
    "Oc. impiger",
    "S. rostratum",
    "Oc. nigripes",
    "S. vittatum"))

#plot capture rate per species. Not enough species-specific data for analysis

species_boxplot <- ggplot(df_attractant, aes(x = TrapType, y = CountPerHour, color=
                                    Family)) +
  facet_wrap(~ Species, scales = "free_y") +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  geom_point(position = position_dodge(width = 0.75)) +
  theme_classic()+
  scale_color_manual(values = c("#4590ba", "#b87f00")) +
  ylab("Individuals per hour")+
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 14),
    legend.position = "none"
  )

species_boxplot


#Test effect of attractants on capture rates per family

#Sum the rate across all species per taxonomic family
df_attractant_fam <- df_attractant %>% 
  group_by(Family, TrapType, Date_Num) %>% 
  summarize(
    SumCountPerHour = sum(CountPerHour)
  )

df_attractant_fam$Family <- as.factor(df_attractant_fam$Family)
df_attractant_fam$TrapType <- as.factor(df_attractant_fam$TrapType)
df_attractant_fam$SumCountPerHour <- as.numeric(df_attractant_fam$SumCountPerHour)


attractant_mod <- glmmTMB(
  SumCountPerHour ~ Family * TrapType + (1|Date_Num),
  family = Gamma(link = "log"),
  data = df_attractant_fam
)

#Check model assumptions
simulationOutput <- simulateResiduals(fittedModel = attractant_mod)
plot(simulationOutput)

#Check model output
summary(attractant_mod)
Anova(attractant_mod, type = "3")

#Proceed with pairwise comparisons
pairwise <- emmeans(attractant_mod, pairwise ~ TrapType | Family, 
                type = "response", adjust = "none")


#plot predictions
plot_pred <- as.data.frame(emmeans(attractant_mod, ~ TrapType | Family, 
                                   type = "response", adjust = "none"))

predicted_attractant <- ggplot() +
  geom_jitter(data = df_attractant_fam, aes(x = TrapType, y = SumCountPerHour, color = Family), alpha = 0.5,
              width = 0.1) +
  geom_errorbar(
    data = plot_pred,
    aes(
      x = TrapType,
      ymin = asymp.LCL,
      ymax = asymp.UCL,
      color = Family
    ),
    width = 0.2,
    linewidth = 1
  ) +
  geom_point(data = plot_pred, aes(x = TrapType, y = response, color = Family), size = 4) +
  theme_classic() +
  ylab("Individuals per hour") +
  facet_wrap(~ Family, scales = "free_y") +
  scale_color_manual(values = c("#4590ba", "#b87f00")) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    legend.position = "none",
    strip.text = element_text(size = 14)
  )

predicted_attractant

attractant_plots <- (predicted_attractant / species_boxplot) + plot_annotation(tag_levels = "A")

attractant_plots <- attractant_plots + plot_layout(heights = c(1, 2))

#3. Yeast-canister comparison

comp_df <- read.csv("YeastVsCannister.csv", sep =";", h = T)

comp_df <- comp_df %>% 
  mutate_at(c("Site", "Type",
              "Date", "Species", "Group"), as.factor)


#Compare capture rates across species
comp_mod <- glmmTMB(
  CountPerHour ~ Type * Species + (1|Date) + (1|Site),
  family = Gamma(link = "log"),
  data = comp_df
)

#Check model diangostics
simulationOutput <- simulateResiduals(fittedModel = comp_mod)
plot(simulationOutput)

#Check model output
summary(comp_mod)
Anova(comp_mod, type = "3")

#Pairwise comparisons
emmeans(comp_mod, pairwise ~ Type | Species, 
                type = "response", adjust = "none")

plot_pred_comp_sp <- as.data.frame(emmeans(comp_mod, ~ Type | Species, 
                                   type = "response", adjust = "none"))

#Prepare data for plotting
plot_pred_comp_sp <- plot_pred_comp_sp %>% 
  mutate(
    Species = fct_recode(
      Species,
      "Oc. impiger" = "Oc_impiger",
      "Oc. nigripes" = "Oc_nigripes",
      "S. vittatum" = "S_vittatum",
      "S. rostratum" = "S_rostratum"
    )) %>% 
  mutate(Species = fct_relevel(
    Species,
    "Oc. impiger",
    "S. rostratum",
    "Oc. nigripes",
    "S. vittatum"))

comp_df <- comp_df %>% 
  mutate(
    Species = fct_recode(
      Species,
      "Oc. impiger" = "Oc_impiger",
      "Oc. nigripes" = "Oc_nigripes",
      "S. vittatum" = "S_vittatum",
      "S. rostratum" = "S_rostratum"
    )) %>%  mutate(Species = fct_relevel(
      Species,
      "Oc. impiger",
      "S. rostratum",
      "Oc. nigripes",
      "S. vittatum"))



predicted_comp_sp <- ggplot() +
  geom_jitter(data = comp_df, aes(x = Type, y = CountPerHour, color = Type), alpha = 0.5,
              width = 0.1) +
  geom_errorbar(
    data = plot_pred_comp_sp,
    aes(
      x = Type,
      ymin = asymp.LCL,
      ymax = asymp.UCL,
      color = Type
    ),
    width = 0.2,
    linewidth = 1
  ) +
  geom_point(data = plot_pred_comp_sp, aes(x = Type, y = response, color = Type), size = 4) +
  theme_classic() +
  ylab("Individuals per hour") +
  facet_wrap(~ Species, scales = "free_y") +
  scale_color_manual(values = c("#4590ba", "#b87f00")) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    legend.position = "none",
    strip.text = element_text(size = 14),
    axis.ticks.x = element_blank()
  )

#Repeat for taxonomic family

Fam_df <- comp_df %>% 
  group_by(Group, Date, Type, Site) %>% 
  summarize(
    Sum = sum(CountPerHour)
  )

Fam_comp_mod <- glmmTMB(
  Sum ~ Type * Group + (1|Date) + (1|Site),
  family = Gamma(link = "log"),
  data = Fam_df
)

#Check model diagnostics
simulationOutput <- simulateResiduals(fittedModel = Fam_comp_mod)
plot(simulationOutput)

#Check model output
summary(Fam_comp_mod)
Anova(Fam_comp_mod, type = "3")

#Pairwise comparisons
emmeans(Fam_comp_mod, pairwise ~ Type | Group, 
                type = "response", adjust = "none")

Fam_comp_pred <- as.data.frame(emmeans(Fam_comp_mod, ~ Type | Group, 
                            type = "response", adjust = "none"))

predicted_comp_fam <- ggplot() +
  geom_jitter(data = Fam_df, aes(x = Type, y = Sum, color = Type), alpha = 0.5,
              width = 0.1) +
  geom_errorbar(
    data = Fam_comp_pred,
    aes(
      x = Type,
      ymin = asymp.LCL,
      ymax = asymp.UCL,
      color = Type
    ),
    width = 0.2,
    linewidth = 1
  ) +
  geom_point(data = Fam_comp_pred, aes(x = Type, y = response, color = Type), size = 4) +
  theme_classic() +
  ylab("Individuals per hour") +
  facet_wrap(~ Group, scales = "free_y") +
  scale_color_manual(values = c("#4590ba", "#b87f00")) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 14),
    legend.position = c(0.95, 0.8)
  )

predicted_comp_fam

comparison_plots <- predicted_comp_fam / predicted_comp_sp

comparison_plots <- comparison_plots + plot_annotation(tag_levels = "A") + plot_layout(heights = c(1,2))
