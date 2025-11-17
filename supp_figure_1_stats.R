# Load libraries

library(easystats)
library(readxl)
library(dplyr)
library(ggplot2)
library(here)
library(gtsummary)
library(gt)
library(ggh4x)
library(ggsignif)

# Load data

base_dir = here()

data <- read_excel(file.path(base_dir, "LastScripts", "band_power_output.xlsx"))
data <- data %>% mutate(patient = as.factor(patient),
                       freq = as.factor(freq),
                       bpd = as.factor(bpd),
                       power = as.double(power))

data <- data %>% group_by(freq, patient) %>%
mutate(baseline = power[bpd == 0]) %>%
mutate(power_ref = power - baseline) 

data_no_baseline <- data %>%
  filter(bpd != 0)  

# Run stats


new_x_labs = c('4 mm', '8 mm', '12 mm', '16 mm', '20 mm')

freq_bands = c('Delta (2-4Hz)', 'Theta (4-8Hz)', 'Alpha (8-13Hz)', 'Beta (13-25Hz)', 'Gamma (25-50Hz)', 'High Gamma (50-200Hz)')

data_no_baseline <- data_no_baseline %>%
  mutate(bpd = factor(bpd, levels = c(1:5),
                      labels = new_x_labs))

data_no_baseline <- data_no_baseline %>%
mutate(freq = factor(freq, levels=c('delta', 'theta', 'alpha', 'beta', 'gamma', 'high_gamma'),
                     labels = freq_bands))


data_summary <- data_no_baseline %>%
  group_by(freq, bpd) %>%
  summarise(median_power_ref = median(power_ref),
            IQR_power_ref = IQR(power_ref),
            p_value = wilcox.test(power_ref, mu=0)$p.value,
            effect_size = rank_biserial(power_ref, mu = 0)$r)
# Adjust p-values for multiple comparisons

p_value_adjust <- p.adjust(data_summary$p_value, method="BH")

data_summary$p_value_correct <- p_value_adjust

data_no_baseline <- data_no_baseline %>%
  left_join(data_summary,by = c("freq", "bpd")) %>%
mutate(sig_star = case_when(
    p_value_correct < 0.001 & patient == "EC131" ~ "***",
    p_value_correct < 0.01 & patient == "EC131" ~ "**",
    p_value_correct < 0.05 & patient == "EC131" ~ "*",
    TRUE ~ ""))

data_stars <- data_no_baseline
data_stars <- data_stars %>% 
group_by(freq, bpd) %>%
  mutate(offset_y = max(power_ref)+0.5*max(abs(power_ref)))


adj_factor <- 1.6 # adjust this factor to give more space above the violins for the significance stars

position_scales <- list(
scale_y_continuous(
  limits = c(
    -max(abs(data_no_baseline[data_no_baseline$freq == "Delta (2-4Hz)", ]$power_ref)),
    adj_factor*max(abs(data_no_baseline[data_no_baseline$freq == "Delta (2-4Hz)", ]$power_ref))
  )
),
scale_y_continuous(
  limits = c(
    -max(abs(data_no_baseline[data_no_baseline$freq == "Theta (4-8Hz)", ]$power_ref)),
    adj_factor*max(abs(data_no_baseline[data_no_baseline$freq == "Theta (4-8Hz)", ]$power_ref))
  )
),
scale_y_continuous(
  limits = c(
    -max(abs(data_no_baseline[data_no_baseline$freq == "Alpha (8-13Hz)", ]$power_ref)),
    adj_factor*max(abs(data_no_baseline[data_no_baseline$freq == "Alpha (8-13Hz)", ]$power_ref))
  )
),
scale_y_continuous(
    limits = c(
        -max(abs(data_no_baseline[data_no_baseline$freq == "Beta (13-25Hz)", ]$power_ref)),
        adj_factor*max(abs(data_no_baseline[data_no_baseline$freq == "Beta (13-25Hz)", ]$power_ref))
    )
),
scale_y_continuous(
    limits = c(
        -max(abs(data_no_baseline[data_no_baseline$freq == "Gamma (25-50Hz)", ]$power_ref)),
        adj_factor*max(abs(data_no_baseline[data_no_baseline$freq == "Gamma (25-50Hz)", ]$power_ref))
    )
),
scale_y_continuous(
    limits = c(
        -max(abs(data_no_baseline[data_no_baseline$freq == "High Gamma (50-200Hz)", ]$power_ref)),
        adj_factor*max(abs(data_no_baseline[data_no_baseline$freq == "High Gamma (50-200Hz)", ]$power_ref))
    )
)
)

# plot
 g <- ggplot(data_no_baseline, aes(x=bpd, y=power_ref) ) +
    geom_violin(position=position_dodge(0.8), width=0.7,trim=FALSE, fill = "lightgray") +
    geom_boxplot(width=0.1) +
    scale_y_continuous() +
    theme_minimal() +
    labs(title="Difference in Re-referenced Band Power by Frequency Band and Re-Reference Distance",
        x="Bipolar Distance From Referential",
        y="Mean sqrt (power) relative to referential") +
    theme(legend.position="top") +
     geom_text(data=data_stars,aes(x=bpd,y=offset_y,label=sig_star), inherit.aes = FALSE,color = "black", size = 10)+
    facet_wrap2(~ freq, scales = "free_y") +
    facetted_pos_scales(y = position_scales) 



    
    # plot
    ggplot(data_summary, aes(x=bpd, y=effect_size)) +
    geom_point() +
    scale_y_continuous() +
    theme_minimal() +
    labs(title="Effect Size by Frequency Band and Re-Reference Distance",
        x="Bipolar Re-Reference Distance",
        y="Effect Size (rank biserial correlation coefficient)") +
    theme(legend.position="top") +
    facet_wrap(~ freq) 
    
readr::write_csv(data_summary, "output_table.csv")
