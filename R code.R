




library(zoo)
library(dplyr)
library(readr)
library(readxl)
library(DescTools)
library(Hmisc)
library(plyr)
library(psych)
library("fixest")
library("texreg")

PCI <- read_csv("import.csv")

# Calculate logarithmic transformation of sales
PCI <- PCI %>% mutate(lnsaleq = log(saleq))

# Winsorize variables at 5% bilateral level, new variables with _w suffix
winsor_cols <- c(
  "AIC_stock_it_2", "PCIC_stock_it_2", "AIPCI_stock_it_2",
  "AIC_stock_it_5", "PCIC_stock_it_5", "AIPCI_stock_it_5",
  "real_f", "cites_f", "patent_count_f",
  "ai_real_f", "ai_cites_f", "ai_patent_count_f",
  "real_i", "cites_i", "patent_count_i",
  "ai_real_i", "ai_cites_i", "ai_patent_count_i",
  "lnsaleq", "roa", "tobinQ", "xrdintensity",
  "cash_at", "leverage"
)

for (col in winsor_cols) {
  q_vals <- quantile(PCI[[col]], probs = c(0.05, 0.95), na.rm = TRUE)
  PCI[[paste0(col, "_w")]] <- DescTools::Winsorize(x = PCI[[col]], val = q_vals)
}

# Apply log(1+x) transformation to winsorized variables
log_cols_w <- c(
  "AIC_stock_it_2_w", "PCIC_stock_it_2_w", "AIPCI_stock_it_2_w",
  "AIC_stock_it_5_w", "PCIC_stock_it_5_w", "AIPCI_stock_it_5_w",
  "real_f_w", "cites_f_w", "patent_count_f_w",
  "ai_real_f_w", "ai_cites_f_w", "ai_patent_count_f_w",
  "real_i_w", "cites_i_w", "patent_count_i_w",
  "ai_real_i_w", "ai_cites_i_w", "ai_patent_count_i_w"
)

for (col_w in log_cols_w) {
  col_raw <- gsub("_w$", "", col_w)
  new_col <- paste0("ln.", col_raw)
  PCI[[new_col]] <- log(1 + PCI[[col_w]])
}

# Generate one-period lagged variables (fill missing values with 0)
lag_vars <- c(
  "lnsaleq_w", "roa_w", "tobinQ_w", "xrdintensity_w", "cash_at_w", 
  "leverage_w", "age", "ln.AIC_stock_it_2", "ln.PCIC_stock_it_2",
  "ln.AIPCI_stock_it_2", "ln.AIC_stock_it_5", "ln.PCIC_stock_it_5", 
  "ln.AIPCI_stock_it_5"
)

PCI <- dplyr::arrange(PCI, rcid, year_quarter)
PCI <- dplyr::group_by(PCI, rcid)

for (var in lag_vars) {
  new_var <- paste0("lag.", var)
  PCI[[new_var]] <- dplyr::lag(PCI[[var]], n = 1, default = 0)
}

PCI <- dplyr::ungroup(PCI)

write_csv(PCI, "output.csv")









# Z-score standardization for core variables and interaction term construction
PCI <- PCI %>%
  mutate(
    lag.ln.AIC_stock_it_2_z = (lag.ln.AIC_stock_it_2 - mean(lag.ln.AIC_stock_it_2, na.rm = TRUE)) / sd(lag.ln.AIC_stock_it_2, na.rm = TRUE),
    lag.ln.PCIC_stock_it_2_z = (lag.ln.PCIC_stock_it_2 - mean(lag.ln.PCIC_stock_it_2, na.rm = TRUE)) / sd(lag.ln.PCIC_stock_it_2, na.rm = TRUE),
    ai_pci_interaction_z = lag.ln.AIC_stock_it_2_z * lag.ln.PCIC_stock_it_2_z
  )

# Variable definition for regression models
controls <- "lag.lnsaleq_w + lag.roa_w + lag.tobinQ_w + lag.xrdintensity_w + lag.cash_at_w + lag.leverage_w + lag.age"
fe_term <- "| rcid + year_quarter"
model_names <- c("AI", "PCI", "AI+PCI+AI*PCI")

# Configuration for regression output
reg_args <- list(
  caption.above = FALSE,
  digits = 4,
  stars = c(0.1, 0.05, 0.01),
  star.symbols = c("*", "**", "***"),
  se.parentheses = TRUE,
  se.pos = "below",
  include.obs = TRUE,
  include.r2 = TRUE,
  include.adj.r2 = TRUE,
  include.fstat = FALSE,
  include.loglik = FALSE,
  include.aic = FALSE,
  add.lines = list(
    "Firm FE (rcid)" = rep("Yes", 3),
    "Quarter FE (year_quarter)" = rep("Yes", 3),
    "Controls" = rep("Yes", 3)
  ),
  custom.note = "* p<0.1, ** p<0.05, *** p<0.01",
  note.above = FALSE
)










# Regression analysis with specified dependent variables and 3 model specifications (AI only/PCI only/AI+PCI+interaction)
# Group 1: ln.patent_count_f and ln.cites_f
fe_patent_f1 <- feols(as.formula(paste("ln.patent_count_f", "~", "lag.ln.AIC_stock_it_2_z", "+", controls, fe_term)), data=PCI)
fe_patent_f2 <- feols(as.formula(paste("ln.patent_count_f", "~", "lag.ln.PCIC_stock_it_2_z", "+", controls, fe_term)), data=PCI)
fe_patent_f3 <- feols(as.formula(paste("ln.patent_count_f", "~", "lag.ln.AIC_stock_it_2_z + lag.ln.PCIC_stock_it_2_z + ai_pci_interaction_z", "+", controls, fe_term)), data=PCI)

fe_cites_f1 <- feols(as.formula(paste("ln.cites_f", "~", "lag.ln.AIC_stock_it_2_z", "+", controls, fe_term)), data=PCI)
fe_cites_f2 <- feols(as.formula(paste("ln.cites_f", "~", "lag.ln.PCIC_stock_it_2_z", "+", controls, fe_term)), data=PCI)
fe_cites_f3 <- feols(as.formula(paste("ln.cites_f", "~", "lag.ln.AIC_stock_it_2_z + lag.ln.PCIC_stock_it_2_z + ai_pci_interaction_z", "+", controls, fe_term)), data=PCI)

# Group 2: ln.ai_patent_count_f and ln.ai_cites_f
fe_ai_patent_f1 <- feols(as.formula(paste("ln.ai_patent_count_f", "~", "lag.ln.AIC_stock_it_2_z", "+", controls, fe_term)), data=PCI)
fe_ai_patent_f2 <- feols(as.formula(paste("ln.ai_patent_count_f", "~", "lag.ln.PCIC_stock_it_2_z", "+", controls, fe_term)), data=PCI)
fe_ai_patent_f3 <- feols(as.formula(paste("ln.ai_patent_count_f", "~", "lag.ln.AIC_stock_it_2_z + lag.ln.PCIC_stock_it_2_z + ai_pci_interaction_z", "+", controls, fe_term)), data=PCI)

fe_ai_cites_f1 <- feols(as.formula(paste("ln.ai_cites_f", "~", "lag.ln.AIC_stock_it_2_z", "+", controls, fe_term)), data=PCI)
fe_ai_cites_f2 <- feols(as.formula(paste("ln.ai_cites_f", "~", "lag.ln.PCIC_stock_it_2_z", "+", controls, fe_term)), data=PCI)
fe_ai_cites_f3 <- feols(as.formula(paste("ln.ai_cites_f", "~", "lag.ln.AIC_stock_it_2_z + lag.ln.PCIC_stock_it_2_z + ai_pci_interaction_z", "+", controls, fe_term)), data=PCI)

# Export results to two separate tables
do.call(htmlreg, c(list(
  list(fe_patent_f1, fe_patent_f2, fe_patent_f3, fe_cites_f1, fe_cites_f2, fe_cites_f3),
  file = "output_table1.doc",
  custom.model.names = c(model_names, model_names),
  caption = "Dependent Variables: ln.patent_count_f (Upper) & ln.cites_f (Lower)"
), reg_args))

do.call(htmlreg, c(list(
  list(fe_ai_patent_f1, fe_ai_patent_f2, fe_ai_patent_f3, fe_ai_cites_f1, fe_ai_cites_f2, fe_ai_cites_f3),
  file = "output_table2.doc",
  custom.model.names = c(model_names, model_names),
  caption = "Dependent Variables: ln.ai_patent_count_f (Upper) & ln.ai_cites_f (Lower)"
), reg_args))









# Figure 1

library(dplyr)
library(readr)
library(fixest)
library(sjPlot)
library(ggplot2)
library(patchwork)

PCI <- read_csv("input.csv") %>%
  mutate(
    AIC_std = scale(lag.ln.AIC_stock_it_2)[,1],
    PCIC_std = scale(lag.ln.PCIC_stock_it_2)[,1]
  ) %>%
  filter(!is.na(AIC_std) & !is.na(PCIC_std))

# Model specification and fitting with fixed effects
controls_str <- "lag.lnsaleq_w + lag.roa_w + lag.tobinQ_w + lag.xrdintensity_w + lag.cash_at_w + lag.leverage_w + lag.age"
fml1 <- as.formula(paste0("ln.patent_count_f ~ AIC_std * PCIC_std + ", controls_str, "| rcid + year_quarter"))
fml2 <- as.formula(paste0("ln.cites_f ~ AIC_std * PCIC_std + ", controls_str, "| rcid + year_quarter"))
fml3 <- as.formula(paste0("ln.ai_patent_count_f ~ AIC_std * PCIC_std + ", controls_str, "| rcid + year_quarter"))
fml4 <- as.formula(paste0("ln.ai_cites_f ~ AIC_std * PCIC_std + ", controls_str, "| rcid + year_quarter"))

model1 <- feols(fml1, data = PCI, cluster = "rcid")
model2 <- feols(fml2, data = PCI, cluster = "rcid")
model3 <- feols(fml3, data = PCI, cluster = "rcid")
model4 <- feols(fml4, data = PCI, cluster = "rcid")

# Core interaction plot function
plot_interact <- function(model, panel_label, panel_title, show_y_tick = TRUE, y_range = NULL) {
  p <- plot_model(model, 
                  type = "int", 
                  terms = c("AIC_std", "PCIC_std = c(-1,1)"),
                  ci.lvl = 0.95, 
                  show.data = FALSE,
                  colors = c("#B3D1FF", "#2E5984"),
                  line.size = 1.4,
                  ci.alpha = 0.2
  ) +
    labs(x = "AIC", y = "", title = paste0(panel_label, " ", panel_title)) +
    scale_x_continuous(expand = c(0.05, 0)) +
    {if (!is.null(y_range)) scale_y_continuous(expand = c(0.05, 0), limits = y_range)} +
    scale_color_manual(values = c("#B3D1FF", "#2E5984")) +
    scale_fill_manual(values = c("#B3D1FF", "#2E5984")) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold", color = "#333333", hjust = 0.5, margin = margin(b = 5)),
      axis.title.x = element_text(size = 10, face = "bold", color = "#333333"),
      axis.text.x = element_text(size = 8, color = "#333333"),
      axis.text.y = if (show_y_tick) element_text(size = 8, color = "#333333") else element_blank(),
      axis.ticks.y = if (show_y_tick) element_line(color = "#333333") else element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_blank(),
      legend.position = "none",
      aspect.ratio = 1,
      plot.margin = margin(8, 8, 8, 8, "mm")
    )
  
  return(p)
}

# Generate subplots
p1 <- plot_interact(model1, "(A)", "Patents", show_y_tick = TRUE)
p2 <- plot_interact(model2, "(B)", "Patent citations", show_y_tick = FALSE)

y_range_cd <- c(-0.4, 0.3)
p3 <- plot_interact(model3, "(C)", "AI patents", show_y_tick = TRUE, y_range = y_range_cd)
p4 <- plot_interact(model4, "(D)", "AI patent citations", show_y_tick = FALSE, y_range = y_range_cd)

# Combine subplots into final panel
final_plot <- (p1 + p2) / (p3 + p4) +
  plot_layout(guides = "collect", axes = "collect") &
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 9, face = "bold", color = "#333333"),
    legend.text = element_text(size = 8, color = "#333333"),
    legend.key.width = unit(2, "cm"),
    legend.margin = margin(t=10, b=5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) &
  scale_color_manual(
    name = "PCIC",
    values = c("#B3D1FF", "#2E5984"),
    labels = c("Low PCIC (−1 SD)", "High PCIC (+1 SD)")
  ) &
  scale_fill_manual(
    name = "PCIC",
    values = c("#B3D1FF", "#2E5984"),
    labels = c("Low PCIC (−1 SD)", "High PCIC (+1 SD)")
  )

out_path_prefix <- "output"


ggsave(
  filename = paste0(out_path_prefix, ".tiff"),
  plot = final_plot,
  width = 10,
  height = 10,
  dpi = 600,
  bg = "white",
  compression = "lzw"
)

ggsave(
  filename = paste0(out_path_prefix, ".pdf"),
  plot = final_plot,
  width = 10,
  height = 10,
  bg = "white",
  device = cairo_pdf
)

print(final_plot)





# Heterogeneity Analysis H1: SG&A intensity (stock-2, patent-count-f, lncites-f)
library(zoo)
library(dplyr)
library(readr)
library(readxl)
library(DescTools)
library(Hmisc)
library(plyr)
library(psych)
library("fixest")
library("texreg")

PCI <- read_csv("input.csv")

# Split sample by sgna_intensity (1=high, 0=low)
PCI_high <- PCI %>% filter(sgna_intensity == 1)
PCI_low <- PCI %>% filter(sgna_intensity == 0)

# Z-score standardization for AIC/PCIC and interaction term construction
PCI <- PCI %>%
  mutate(
    lag.ln.AIC_stock_it_2_z = (lag.ln.AIC_stock_it_2 - mean(lag.ln.AIC_stock_it_2, na.rm = TRUE)) / sd(lag.ln.AIC_stock_it_2, na.rm = TRUE),
    lag.ln.PCIC_stock_it_2_z = (lag.ln.PCIC_stock_it_2 - mean(lag.ln.PCIC_stock_it_2, na.rm = TRUE)) / sd(lag.ln.PCIC_stock_it_2, na.rm = TRUE),
    ai_pci_interaction_z = lag.ln.AIC_stock_it_2_z * lag.ln.PCIC_stock_it_2_z
  )

# Update grouped data with standardized variables
PCI_high <- PCI %>% filter(sgna_intensity == 1)
PCI_low <- PCI %>% filter(sgna_intensity == 0)

# Define control variables and fixed effects formula
controls <- "lag.lnsaleq_w + lag.roa_w + lag.tobinQ_w + lag.xrdintensity_w + lag.cash_at_w + lag.leverage_w + lag.age"
fe_term <- "| rcid + year_quarter"

# Define dependent variables in specified order
dv_target <- c(
  "ln.patent_count_f", 
  "ln.cites_f", 
  "ln.ai_patent_count_f", 
  "ln.ai_cites_f"
)

# Regression for high sgna_intensity group (AI+PCI+AI*PCI model)
fe_high1 <- feols(as.formula(paste(dv_target[1], "~", "lag.ln.AIC_stock_it_2_z + lag.ln.PCIC_stock_it_2_z + ai_pci_interaction_z", "+", controls, fe_term)), data=PCI_high)
fe_high2 <- feols(as.formula(paste(dv_target[2], "~", "lag.ln.AIC_stock_it_2_z + lag.ln.PCIC_stock_it_2_z + ai_pci_interaction_z", "+", controls, fe_term)), data=PCI_high)
fe_high3 <- feols(as.formula(paste(dv_target[3], "~", "lag.ln.AIC_stock_it_2_z + lag.ln.PCIC_stock_it_2_z + ai_pci_interaction_z", "+", controls, fe_term)), data=PCI_high)
fe_high4 <- feols(as.formula(paste(dv_target[4], "~", "lag.ln.AIC_stock_it_2_z + lag.ln.PCIC_stock_it_2_z + ai_pci_interaction_z", "+", controls, fe_term)), data=PCI_high)

# Regression for low sgna_intensity group (AI+PCI+AI*PCI model)
fe_low1 <- feols(as.formula(paste(dv_target[1], "~", "lag.ln.AIC_stock_it_2_z + lag.ln.PCIC_stock_it_2_z + ai_pci_interaction_z", "+", controls, fe_term)), data=PCI_low)
fe_low2 <- feols(as.formula(paste(dv_target[2], "~", "lag.ln.AIC_stock_it_2_z + lag.ln.PCIC_stock_it_2_z + ai_pci_interaction_z", "+", controls, fe_term)), data=PCI_low)
fe_low3 <- feols(as.formula(paste(dv_target[3], "~", "lag.ln.AIC_stock_it_2_z + lag.ln.PCIC_stock_it_2_z + ai_pci_interaction_z", "+", controls, fe_term)), data=PCI_low)
fe_low4 <- feols(as.formula(paste(dv_target[4], "~", "lag.ln.AIC_stock_it_2_z + lag.ln.PCIC_stock_it_2_z + ai_pci_interaction_z", "+", controls, fe_term)), data=PCI_low)

# Define model names for output table
model_names <- c(
  "ln.patent_count_f (high)", "ln.cites_f (high)", "ln.ai_patent_count_f (high)", "ln.ai_cites_f (high)",
  "ln.patent_count_f (low)",  "ln.cites_f (low)",  "ln.ai_patent_count_f (low)",  "ln.ai_cites_f (low)"
)

# Configuration for regression table output
reg_args <- list(
  caption.above = FALSE,
  digits = 4,
  stars = c(0.1, 0.05, 0.01),
  star.symbols = c("*", "**", "***"),
  se.parentheses = TRUE,
  se.pos = "below",
  include.obs = TRUE,
  include.r2 = TRUE,
  include.adj.r2 = TRUE,
  include.fstat = FALSE,
  include.loglik = FALSE,
  include.aic = FALSE,
  add.lines = list(
    "Firm FE (rcid)" = rep("Yes", 8), 
    "Quarter FE (year_quarter)" = rep("Yes", 8),
    "Controls" = rep("Yes", 8)
  ),
  custom.note = "* p<0.1, ** p<0.05, *** p<0.01. All independent variables are z-score standardized.",
  note.above = FALSE
)

# Combine all models and export to a single table
do.call(htmlreg, c(list(
  list(fe_high1, fe_high2, fe_high3, fe_high4, fe_low1, fe_low2, fe_low3, fe_low4),
  file = "output.doc",
  custom.model.names = model_names,
  caption = "Regression Results for AI+PCI+AI*PCI Model by SGNA Intensity (High vs. Low)"
), reg_args))









# Heterogeneity Analysis H2: Intangible assets intensity (stock-2, patent-count-f, lncites-f)
library(zoo)
library(dplyr)
library(readr)
library(readxl)
library(DescTools)
library(Hmisc)
library(plyr)
library(psych)
library("fixest")
library("texreg")

PCI <- read_csv("input.csv")

# Z-score standardization for core independent variables and interaction term construction
PCI <- PCI %>%
  mutate(
    lag.ln.AIC_stock_it_2_z = (lag.ln.AIC_stock_it_2 - mean(lag.ln.AIC_stock_it_2, na.rm = TRUE)) / sd(lag.ln.AIC_stock_it_2, na.rm = TRUE),
    lag.ln.PCIC_stock_it_2_z = (lag.ln.PCIC_stock_it_2 - mean(lag.ln.PCIC_stock_it_2, na.rm = TRUE)) / sd(lag.ln.PCIC_stock_it_2, na.rm = TRUE),
    ai_pci_interaction_z = lag.ln.AIC_stock_it_2_z * lag.ln.PCIC_stock_it_2_z
  )

# Split sample by intan_intensity (1=high, 0=low)
PCI_high <- PCI %>% filter(intan_intensity == 1)
PCI_low <- PCI %>% filter(intan_intensity == 0)

# Define control variables and fixed effects formula
controls <- "lag.lnsaleq_w + lag.roa_w + lag.tobinQ_w + lag.xrdintensity_w + lag.cash_at_w + lag.leverage_w + lag.age"
fe_term <- "| rcid + year_quarter"

# Regression for high intan_intensity group (AI+PCI+AI*PCI model)
fe_high1 <- feols(as.formula(paste("ln.patent_count_f", "~", "lag.ln.AIC_stock_it_2_z + lag.ln.PCIC_stock_it_2_z + ai_pci_interaction_z", "+", controls, fe_term)), data=PCI_high)
fe_high2 <- feols(as.formula(paste("ln.cites_f", "~", "lag.ln.AIC_stock_it_2_z + lag.ln.PCIC_stock_it_2_z + ai_pci_interaction_z", "+", controls, fe_term)), data=PCI_high)
fe_high3 <- feols(as.formula(paste("ln.ai_patent_count_f", "~", "lag.ln.AIC_stock_it_2_z + lag.ln.PCIC_stock_it_2_z + ai_pci_interaction_z", "+", controls, fe_term)), data=PCI_high)
fe_high4 <- feols(as.formula(paste("ln.ai_cites_f", "~", "lag.ln.AIC_stock_it_2_z + lag.ln.PCIC_stock_it_2_z + ai_pci_interaction_z", "+", controls, fe_term)), data=PCI_high)

# Regression for low intan_intensity group (AI+PCI+AI*PCI model)
fe_low1 <- feols(as.formula(paste("ln.patent_count_f", "~", "lag.ln.AIC_stock_it_2_z + lag.ln.PCIC_stock_it_2_z + ai_pci_interaction_z", "+", controls, fe_term)), data=PCI_low)
fe_low2 <- feols(as.formula(paste("ln.cites_f", "~", "lag.ln.AIC_stock_it_2_z + lag.ln.PCIC_stock_it_2_z + ai_pci_interaction_z", "+", controls, fe_term)), data=PCI_low)
fe_low3 <- feols(as.formula(paste("ln.ai_patent_count_f", "~", "lag.ln.AIC_stock_it_2_z + lag.ln.PCIC_stock_it_2_z + ai_pci_interaction_z", "+", controls, fe_term)), data=PCI_low)
fe_low4 <- feols(as.formula(paste("ln.ai_cites_f", "~", "lag.ln.AIC_stock_it_2_z + lag.ln.PCIC_stock_it_2_z + ai_pci_interaction_z", "+", controls, fe_term)), data=PCI_low)

# Configuration for top journal style regression table output
reg_args <- list(
  caption.above = FALSE,
  digits = 4,
  stars = c(0.1, 0.05, 0.01),
  star.symbols = c("*", "**", "***"),
  se.parentheses = TRUE,
  se.pos = "below",
  include.obs = TRUE,
  include.r2 = TRUE,
  include.adj.r2 = TRUE,
  include.fstat = FALSE,
  include.loglik = FALSE,
  include.aic = FALSE,
  add.lines = list(
    "Firm FE (rcid)" = rep("Yes", 8),          
    "Quarter FE (year_quarter)" = rep("Yes", 8),
    "Controls" = rep("Yes", 8)                
  ),
  custom.note = "* p<0.1, ** p<0.05, *** p<0.01. All core independent variables are standardized by z-score.",
  note.above = FALSE
)

# Combine all models and export to a single table
model_names <- c(
  "ln.patent_count_f (high)", "ln.cites_f (high)", "ln.ai_patent_count_f (high)", "ln.ai_cites_f (high)",
  "ln.patent_count_f (low)",  "ln.cites_f (low)",  "ln.ai_patent_count_f (low)",  "ln.ai_cites_f (low)"
)

do.call(htmlreg, c(list(
  list(fe_high1, fe_high2, fe_high3, fe_high4, fe_low1, fe_low2, fe_low3, fe_low4),
  file = "output.doc",
  custom.model.names = model_names,
  caption = "Regression Results (AI+PCI+AI*PCI) by Intan Intensity (High vs. Low)"
), reg_args))











# Placebo Test: Double Randomization + Lead T+8 + Firm-Year Clustering
library(zoo)
library(dplyr)
library(readr)
library(readxl)
library(DescTools)
library(Hmisc)
library(plyr)
library(psych)
library("fixest")
library("texreg")

PCI <- read_csv("input.csv")

# Set random seed for reproducibility and data preprocessing (sort + extract year for stratified sampling)
set.seed(9876)
PCI <- PCI %>%
  arrange(rcid, year_quarter) %>%
  mutate(year = substr(year_quarter, 1, 4)) %>%
  group_by(rcid) %>%
  mutate(year_quarter_cluster = year_quarter) %>%
  ungroup()

# Core placebo design: Lead T+8 + double random permutation + stratified sampling
PCI <- PCI %>%
  group_by(rcid) %>%
  mutate(
    lead.ln.AIC_stock_it_2 = lead(ln.AIC_stock_it_2, n = 8),
    lead.ln.PCIC_stock_it_2 = lead(ln.PCIC_stock_it_2, n = 8)
  ) %>%
  ungroup() %>%
  group_by(year) %>%
  mutate(
    lead.ln.AIC_stock_it_2_rand1 = sample(lead.ln.AIC_stock_it_2, replace = FALSE),
    lead.ln.PCIC_stock_it_2_rand1 = sample(lead.ln.PCIC_stock_it_2, replace = FALSE)
  ) %>%
  ungroup() %>%
  mutate(
    lead.ln.AIC_stock_it_2_rand = sample(lead.ln.AIC_stock_it_2_rand1, replace = FALSE),
    lead.ln.PCIC_stock_it_2_rand = sample(lead.ln.PCIC_stock_it_2_rand1, replace = FALSE),
    lag1.ln.AIC_stock_it_2 = lag(ln.AIC_stock_it_2, n = 1),
    lag1.ln.PCIC_stock_it_2 = lag(ln.PCIC_stock_it_2, n = 1),
    lag.lnsaleq_w = lag(lnsaleq_w, n = 1),
    lag.roa_w = lag(roa_w, n = 1),
    lag.tobinQ_w = lag(tobinQ_w, n = 1),
    lag.xrdintensity_w = lag(xrdintensity_w, n = 1),
    lag.cash_at_w = lag(cash_at_w, n = 1),
    lag.leverage_w = lag(leverage_w, n = 1),
    lag.age = lag(age, n = 1)
  ) %>%
  filter(
    !is.na(lead.ln.AIC_stock_it_2_rand),
    !is.na(lead.ln.PCIC_stock_it_2_rand),
    !is.infinite(lead.ln.AIC_stock_it_2_rand),
    !is.infinite(lead.ln.PCIC_stock_it_2_rand)
  )

# Z-score standardization for randomized lead variables and lag1 capability terms
PCI <- PCI %>%
  mutate(
    lead.ln.AIC_stock_it_2_z = (lead.ln.AIC_stock_it_2_rand - mean(lead.ln.AIC_stock_it_2_rand, na.rm = T)) / sd(lead.ln.AIC_stock_it_2_rand, na.rm = T),
    lead.ln.PCIC_stock_it_2_z = (lead.ln.PCIC_stock_it_2_rand - mean(lead.ln.PCIC_stock_it_2_rand, na.rm = T)) / sd(lead.ln.PCIC_stock_it_2_rand, na.rm = T),
    ai_pci_interaction_z = lead.ln.AIC_stock_it_2_z * lead.ln.PCIC_stock_it_2_z,
    lag1.ln.AIC_stock_it_2_z = (lag1.ln.AIC_stock_it_2 - mean(lag1.ln.AIC_stock_it_2, na.rm = T)) / sd(lag1.ln.AIC_stock_it_2, na.rm = T),
    lag1.ln.PCIC_stock_it_2_z = (lag1.ln.PCIC_stock_it_2 - mean(lag1.ln.PCIC_stock_it_2, na.rm = T)) / sd(lag1.ln.PCIC_stock_it_2, na.rm = T)
  )

# Define variables and fixed effects for regression
dv1 <- "ln.patent_count_f"
dv2 <- "ln.cites_f"
dv3 <- "ln.ai_patent_count_f"
dv4 <- "ln.ai_cites_f"
controls_base <- "lag.lnsaleq_w + lag.roa_w + lag.tobinQ_w + lag.xrdintensity_w + lag.cash_at_w + lag.leverage_w + lag.age + lag1.ln.AIC_stock_it_2_z + lag1.ln.PCIC_stock_it_2_z"
fe_term <- "| rcid + year_quarter"
model_names <- c("AI (Rand)", "PCI (Rand)", "AI+PCI+AI*PCI (Rand)")

# Regression with firm-year double clustering (placebo test models)
cluster_term <- ~rcid + year
fe1 <- feols(as.formula(paste(dv1, "~", "lead.ln.AIC_stock_it_2_z", "+", controls_base, fe_term)), data=PCI, cluster = cluster_term, na.rm = TRUE)
fe2 <- feols(as.formula(paste(dv1, "~", "lead.ln.PCIC_stock_it_2_z", "+", controls_base, fe_term)), data=PCI, cluster = cluster_term, na.rm = TRUE)
fe3 <- feols(as.formula(paste(dv1, "~", "lead.ln.AIC_stock_it_2_z + lead.ln.PCIC_stock_it_2_z + ai_pci_interaction_z", "+", controls_base, fe_term)), data=PCI, cluster = cluster_term, na.rm = TRUE)
fe4 <- feols(as.formula(paste(dv2, "~", "lead.ln.AIC_stock_it_2_z", "+", controls_base, fe_term)), data=PCI, cluster = cluster_term, na.rm = TRUE)
fe5 <- feols(as.formula(paste(dv2, "~", "lead.ln.PCIC_stock_it_2_z", "+", controls_base, fe_term)), data=PCI, cluster = cluster_term, na.rm = TRUE)
fe6 <- feols(as.formula(paste(dv2, "~", "lead.ln.AIC_stock_it_2_z + lead.ln.PCIC_stock_it_2_z + ai_pci_interaction_z", "+", controls_base, fe_term)), data=PCI, cluster = cluster_term, na.rm = TRUE)

fe7 <- feols(as.formula(paste(dv3, "~", "lead.ln.AIC_stock_it_2_z", "+", controls_base, fe_term)), data=PCI, cluster = cluster_term, na.rm = TRUE)
fe8 <- feols(as.formula(paste(dv3, "~", "lead.ln.PCIC_stock_it_2_z", "+", controls_base, fe_term)), data=PCI, cluster = cluster_term, na.rm = TRUE)
fe9 <- feols(as.formula(paste(dv3, "~", "lead.ln.AIC_stock_it_2_z + lead.ln.PCIC_stock_it_2_z + ai_pci_interaction_z", "+", controls_base, fe_term)), data=PCI, cluster = cluster_term, na.rm = TRUE)
fe10 <- feols(as.formula(paste(dv4, "~", "lead.ln.AIC_stock_it_2_z", "+", controls_base, fe_term)), data=PCI, cluster = cluster_term, na.rm = TRUE)
fe11 <- feols(as.formula(paste(dv4, "~", "lead.ln.PCIC_stock_it_2_z", "+", controls_base, fe_term)), data=PCI, cluster = cluster_term, na.rm = TRUE)
fe12 <- feols(as.formula(paste(dv4, "~", "lead.ln.AIC_stock_it_2_z + lead.ln.PCIC_stock_it_2_z + ai_pci_interaction_z", "+", controls_base, fe_term)), data=PCI, cluster = cluster_term, na.rm = TRUE)

# Configuration for placebo test table output
reg_args <- list(
  caption.above = FALSE,
  digits = 4,
  stars = c(0.1, 0.05, 0.01),
  star.symbols = c("*", "**", "***"),
  se.parentheses = TRUE,
  se.pos = "below",
  include.obs = TRUE,
  include.r2 = TRUE,
  include.adj.r2 = TRUE,
  include.fstat = FALSE,
  include.loglik = FALSE,
  include.aic = FALSE,
  add.lines = list(
    "Firm FE (rcid)" = rep("Yes", 6),
    "Quarter FE (year_quarter)" = rep("Yes", 6),
    "Controls (Lag T-1 + Lag1 Capability)" = rep("Yes", 6),
    "Dependent Var (T Period)" = rep("Yes", 6),
    "Double Randomized AI/PCI (Placebo)" = rep("Yes", 6),
    "Lead Period (T+8)" = rep("Yes", 6),
    "Clustering (Firm-Year)" = rep("Yes", 6)
  ),
  custom.note = "* p<0.1, ** p<0.05, *** p<0.01 (Placebo Test: Double Randomized AI/PCI, Lead T+8, Firm-Year Cluster)",
  note.above = FALSE
)

# Create output directory and export results
output_dir <- "./output/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

do.call(htmlreg, c(list(
  list(fe1, fe2, fe3, fe4, fe5, fe6),
  file = paste0(output_dir, "placebo_test_group1.doc"),
  custom.model.names = c(paste0(dv1, "_", model_names), paste0(dv2, "_", model_names)),
  caption = "Placebo Test Results: ln.patent_count_f & ln.cites_f (Double Randomized AI/PCI, Lead T+8)"
), reg_args))

do.call(htmlreg, c(list(
  list(fe7, fe8, fe9, fe10, fe11, fe12),
  file = paste0(output_dir, "placebo_test_group2.doc"),
  custom.model.names = c(paste0(dv3, "_", model_names), paste0(dv4, "_", model_names)),
  caption = "Placebo Test Results: ln.ai_patent_count_f & ln.ai_cites_f (Double Randomized AI/PCI, Lead T+8)"
), reg_args))










# GICS Industry
library(zoo)
library(dplyr)
library(readr)
library(readxl)
library(DescTools)
library(Hmisc)
library(plyr)
library(psych)
library("fixest")
library("texreg")

PCI <- read_csv("input.csv")

# Variable preprocessing: Industry-year detrending + firm-level standardization
PCI <- PCI %>%
  mutate(
    year = as.numeric(substr(year_quarter, 1, 4)),
    quarter = as.numeric(substr(year_quarter, 6, 6))
  ) %>%
  group_by(gics_code, year) %>%
  mutate(
    lag.ln.AIC_stock_it_2_detrend = lag.ln.AIC_stock_it_2 - mean(lag.ln.AIC_stock_it_2, na.rm=T),
    lag.ln.PCIC_stock_it_2_detrend = lag.ln.PCIC_stock_it_2 - mean(lag.ln.PCIC_stock_it_2, na.rm=T)
  ) %>%
  ungroup() %>%
  group_by(rcid) %>%
  mutate(
    lag.ln.AIC_stock_it_2_z = (lag.ln.AIC_stock_it_2_detrend - mean(lag.ln.AIC_stock_it_2_detrend, na.rm=T)) / sd(lag.ln.AIC_stock_it_2_detrend, na.rm=T),
    lag.ln.PCIC_stock_it_2_z = (lag.ln.PCIC_stock_it_2_detrend - mean(lag.ln.PCIC_stock_it_2_detrend, na.rm=T)) / sd(lag.ln.PCIC_stock_it_2_detrend, na.rm=T)
  ) %>%
  ungroup() %>%
  mutate(
    ai_pci_interaction_z = lag.ln.AIC_stock_it_2_z * lag.ln.PCIC_stock_it_2_z
  )

# Model specification: GICS×Year FE + Industry×Year clustered SE
controls <- "lag.lnsaleq_w + lag.roa_w + lag.tobinQ_w + lag.xrdintensity_w + lag.cash_at_w + lag.leverage_w + lag.age"
fe_term <- "| rcid + gics_code^year"  
cluster_term <- ~gics_code^year  
model_names <- c("AI", "PCI", "AI+PCI+AI*PCI")

reg_args <- list(
  caption.above = FALSE,
  digits = 4,
  stars = c(0.1, 0.05, 0.01),
  star.symbols = c("*", "**", "***"),
  se.parentheses = TRUE,
  se.pos = "below",
  include.obs = TRUE,
  include.r2 = TRUE,
  include.adj.r2 = TRUE,
  include.fstat = FALSE,
  include.loglik = FALSE,
  include.aic = FALSE,
  add.lines = list(
    "Firm FE" = rep("Yes", 3),
    "GICS × Year FE" = rep("Yes", 3),
    "Controls" = rep("Yes", 3),
    "Clustered SE (Industry×Year)" = rep("Yes", 3)
  ),
  custom.note = "* p<0.1, ** p<0.05, *** p<0.01",
  note.above = FALSE
)

# Regression analysis for ln.patent_count_f & ln.cites_f
fe_patent_f1 <- feols(
  as.formula(paste("ln.patent_count_f", "~", "lag.ln.AIC_stock_it_2_z", "+", controls, fe_term)),
  data = PCI,
  cluster = cluster_term
)
fe_patent_f2 <- feols(
  as.formula(paste("ln.patent_count_f", "~", "lag.ln.PCIC_stock_it_2_z", "+", controls, fe_term)),
  data = PCI,
  cluster = cluster_term
)
fe_patent_f3 <- feols(
  as.formula(paste("ln.patent_count_f", "~", "lag.ln.AIC_stock_it_2_z + lag.ln.PCIC_stock_it_2_z + ai_pci_interaction_z", "+", controls, fe_term)),
  data = PCI,
  cluster = cluster_term
)

fe_cites_f1 <- feols(
  as.formula(paste("ln.cites_f", "~", "lag.ln.AIC_stock_it_2_z", "+", controls, fe_term)),
  data = PCI,
  cluster = cluster_term
)
fe_cites_f2 <- feols(
  as.formula(paste("ln.cites_f", "~", "lag.ln.PCIC_stock_it_2_z", "+", controls, fe_term)),
  data = PCI,
  cluster = cluster_term
)
fe_cites_f3 <- feols(
  as.formula(paste("ln.cites_f", "~", "lag.ln.AIC_stock_it_2_z + lag.ln.PCIC_stock_it_2_z + ai_pci_interaction_z", "+", controls, fe_term)),
  data = PCI,
  cluster = cluster_term
)

# Regression analysis for ln.ai_patent_count_f & ln.ai_cites_f
fe_ai_patent_f1 <- feols(
  as.formula(paste("ln.ai_patent_count_f", "~", "lag.ln.AIC_stock_it_2_z", "+", controls, fe_term)),
  data = PCI,
  cluster = cluster_term
)
fe_ai_patent_f2 <- feols(
  as.formula(paste("ln.ai_patent_count_f", "~", "lag.ln.PCIC_stock_it_2_z", "+", controls, fe_term)),
  data = PCI,
  cluster = cluster_term
)
fe_ai_patent_f3 <- feols(
  as.formula(paste("ln.ai_patent_count_f", "~", "lag.ln.AIC_stock_it_2_z + lag.ln.PCIC_stock_it_2_z + ai_pci_interaction_z", "+", controls, fe_term)),
  data = PCI,
  cluster = cluster_term
)

fe_ai_cites_f1 <- feols(
  as.formula(paste("ln.ai_cites_f", "~", "lag.ln.AIC_stock_it_2_z", "+", controls, fe_term)),
  data = PCI,
  cluster = cluster_term
)
fe_ai_cites_f2 <- feols(
  as.formula(paste("ln.ai_cites_f", "~", "lag.ln.PCIC_stock_it_2_z", "+", controls, fe_term)),
  data = PCI,
  cluster = cluster_term
)
fe_ai_cites_f3 <- feols(
  as.formula(paste("ln.ai_cites_f", "~", "lag.ln.AIC_stock_it_2_z + lag.ln.PCIC_stock_it_2_z + ai_pci_interaction_z", "+", controls, fe_term)),
  data = PCI,
  cluster = cluster_term
)

# Export regression results
reg_args_patent <- reg_args
reg_args_patent$add.lines <- lapply(reg_args_patent$add.lines, function(x) rep(x, 6))
do.call(htmlreg, c(list(
  list(fe_patent_f1, fe_patent_f2, fe_patent_f3, fe_cites_f1, fe_cites_f2, fe_cites_f3),
  file = "output_table1.doc",
  custom.model.names = c(model_names, model_names),
  caption = "Dependent Variables: ln.patent_count_f (Upper) & ln.cites_f (Lower) - Final Optimized Model"
), reg_args_patent))

reg_args_ai <- reg_args
reg_args_ai$add.lines <- lapply(reg_args_ai$add.lines, function(x) rep(x, 6))
do.call(htmlreg, c(list(
  list(fe_ai_patent_f1, fe_ai_patent_f2, fe_ai_patent_f3, fe_ai_cites_f1, fe_ai_cites_f2, fe_ai_cites_f3),
  file = "output_table2.doc",
  custom.model.names = c(model_names, model_names),
  caption = "Dependent Variables: ln.ai_patent_count_f (Upper) & ln.ai_cites_f (Lower) - Final Optimized Model"
), reg_args_ai))








# 4-Quarter Hurdle Model (Two-Part Model: FE-Logit + FE-PPML)
library(zoo)
library(dplyr)
library(tidyr)
library(readr)
library(fixest)
library(texreg)

PCI <- read_csv("input.csv", show_col_types = FALSE)

# Data preprocessing: Standardization + 4-quarter cumulative variables + dummy variables
PCI <- PCI %>%
  mutate(
    AIC_std = scale(lag.ln.AIC_stock_it_2)[,1],
    PCIC_std = scale(lag.ln.PCIC_stock_it_2)[,1],
    AIPCI_std = scale(lag.ln.AIPCI_stock_it_2)[,1],
    ai_pci_interaction = AIC_std * PCIC_std
  ) %>%
  arrange(rcid, year_quarter) %>% 
  group_by(rcid) %>%
  mutate(
    patent_count_f_4q = rollsumr(patent_count_f, k = 4, fill = NA),
    cites_f_4q = rollsumr(cites_f, k = 4, fill = NA),
    ai_patent_count_f_4q = rollsumr(ai_patent_count_f, k = 4, fill = NA),
    ai_cites_f_4q = rollsumr(ai_cites_f, k = 4, fill = NA),
    patent_dummy_f_4q = as.integer(patent_count_f_4q > 0),
    cites_dummy_f_4q = as.integer(cites_f_4q > 0),
    ai_patent_dummy_f_4q = as.integer(ai_patent_count_f_4q > 0),
    ai_cites_dummy_f_4q = as.integer(ai_cites_f_4q > 0)
  ) %>%
  ungroup() %>%
  drop_na(contains("_4q"), lag.lnsaleq_w, lag.roa_w)

# Regression specification: 3 independent variable combinations (AI/PCI/AI*PCI)
controls <- "lag.lnsaleq_w + lag.roa_w + lag.tobinQ_w + lag.xrdintensity_w + lag.cash_at_w + lag.leverage_w + lag.age"
fe_term <- "| rcid + year_quarter"
vcov_cluster <- ~rcid

iv_list <- list(
  iv1 = "AIC_std",          
  iv2 = "PCIC_std",         
  iv3 = "AIC_std + PCIC_std + ai_pci_interaction"
)

model_labels <- c(
  "AI (Extensive: FE-Logit)", "PCI (Extensive: FE-Logit)", "AI*PCI (Extensive: FE-Logit)",
  "AI (Intensive: FE-PPML)", "PCI (Intensive: FE-PPML)", "AI*PCI (Intensive: FE-PPML)"
)

# Core function: Run two-part model and export merged table
run_merged_two_part_model <- function(dv_dummy, dv_count, dv_name, output_file) {
  # Extensive margin: FE-Logit (Pr>0)
  extensive_models <- lapply(iv_list, function(iv) {
    fml <- as.formula(paste(dv_dummy, "~", iv, "+", controls, fe_term))
    feglm(fml, data = PCI, family = binomial(link = "logit"), vcov = vcov_cluster)
  })
  
  # Intensive margin: FE-PPML (Cond>0)
  data_intensive <- filter(PCI, .data[[dv_count]] > 0)
  intensive_models <- lapply(iv_list, function(iv) {
    fml <- as.formula(paste(dv_count, "~", iv, "+", controls, fe_term))
    fepois(fml, data = data_intensive, vcov = vcov_cluster)
  })
  
  # Merge and export results to DOC table
  merged_models <- c(extensive_models, intensive_models)
  
  col_headers <- c(
    rep(paste0("Extensive Margin (FE-Logit)\nPr(", dv_name, ">0)"), 3),
    rep(paste0("Intensive Margin (FE-PPML)\nCond(", dv_name, ">0)"), 3)
  )
  
  htmlreg(
    l = merged_models,
    file = output_file,
    custom.model.names = model_labels,
    custom.col.headers = col_headers,
    caption = paste0("Two-Part Model: Extensive (FE-Logit) and Intensive (FE-PPML) Margin for ", dv_name, " (4-Quarter Cumulative)"),
    caption.above = TRUE,
    digits = 3,
    stars = c(0.1, 0.05, 0.01),
    star.symbols = c("*", "**", "***"),
    include.obs = TRUE,
    include.r2 = TRUE,
    include.adj.r2 = FALSE,
    include.loglik = FALSE,
    include.aic = FALSE,
    include.bic = FALSE,
    add.lines = list(
      "Firm Fixed Effects" = rep("Yes", 6),
      "Quarter Fixed Effects" = rep("Yes", 6),
      "Control Variables" = rep("Yes", 6),
      "Clustered SE (Firm)" = rep("Yes", 6)
    ),
    doctype = FALSE,
    html.tag = FALSE,
    body.tag = FALSE,
    table.tag = TRUE,
    css = "",
    custom.note = "* p<0.1, ** p<0.05, *** p<0.01. Standard errors are clustered at the firm level."
  )
  
  return(list(extensive = extensive_models, intensive = intensive_models))
}

# Run models and export results
output_path <- "./output/"

# Dependent variable 1: patent_count_f_4q
model1 <- run_merged_two_part_model(
  dv_dummy = "patent_dummy_f_4q",
  dv_count = "patent_count_f_4q",
  dv_name = "patent_count_f_4q",
  output_file = paste0(output_path, "merged_table_patent_count_f_4q.doc")
)

# Dependent variable 2: cites_f_4q
model2 <- run_merged_two_part_model(
  dv_dummy = "cites_dummy_f_4q",
  dv_count = "cites_f_4q",
  dv_name = "cites_f_4q",
  output_file = paste0(output_path, "merged_table_cites_f_4q.doc")
)

# Dependent variable 3: ai_patent_count_f_4q
model3 <- run_merged_two_part_model(
  dv_dummy = "ai_patent_dummy_f_4q",
  dv_count = "ai_patent_count_f_4q",
  dv_name = "ai_patent_count_f_4q",
  output_file = paste0(output_path, "merged_table_ai_patent_count_f_4q.doc")
)

# Dependent variable 4: ai_cites_f_4q
model4 <- run_merged_two_part_model(
  dv_dummy = "ai_cites_dummy_f_4q",
  dv_count = "ai_cites_f_4q",
  dv_name = "ai_cites_f_4q",
  output_file = paste0(output_path, "merged_table_ai_cites_f_4q.doc")
)








# Appendix Dataset Construction
library(zoo)
library(dplyr)
library(readr)
library(readxl)
library(DescTools)
library(Hmisc)
library(plyr)
library(psych)
library("fixest")
library("texreg")

PCI <- read_csv("input.csv")

# Remove "_ρk" suffix from column names (header modification only)
colnames(PCI) <- gsub("_ρk$", "", colnames(PCI))

# Verify column name modification results (optional)
cat("Modified relevant column names:\n")
grep("(AIC|PCIC|AIPCI)_stock_it_(2|5)", colnames(PCI), value = TRUE)
invisible(sapply(colnames(PCI), cat, "\n")) 

# Define 18 target variables for 5% winsorization
winsor_cols <- c(
  "AI_workforce_share_it_2", "PCI_workforce_share_it_2", "AIPCI_workforce_share_it_2",
  "AI_workforce_share_it_5", "PCI_workforce_share_it_5", "AIPCI_workforce_share_it_5",
  "AI_post_share_it_2", "PCI_post_share_it_2", "AIPCI_post_share_it_2",
  "AI_post_share_it_5", "PCI_post_share_it_5", "AIPCI_post_share_it_5",
  "AIC_stock_it_2", "PCIC_stock_it_2", "AIPCI_stock_it_2",
  "AIC_stock_it_5", "PCIC_stock_it_5", "AIPCI_stock_it_5"
)

# Batch winsorization (5th/95th percentiles) with "_w" suffix
for (col in winsor_cols) {
  q_vals <- quantile(PCI[[col]], probs = c(0.05, 0.95), na.rm = TRUE)
  PCI[[paste0(col, "_w")]] <- DescTools::Winsorize(x = PCI[[col]], val = q_vals)
}

write_csv(PCI, "output.csv")








# Tbele B1
# B1 Partial Correlation Analysis (Pooled + Within-Firm)
library(readr)
library(lmtest) 

PCI <- read_csv("input.csv", show_col_types = FALSE)

# Select core variables and rename (avoid select conflict)
var_check <- c("rcid", "AI_workforce_share_it_2_w", "AI_post_share_it_2_w", 
               "PCI_workforce_share_it_2_w", "PCI_post_share_it_2_w")
if(all(var_check %in% colnames(PCI))){
  data_all <- PCI[, var_check]
  colnames(data_all) <- c("rcid", "AI_work", "AI_post", "PCI_work", "PCI_post")
}else{
  stop("Core variables not found, check variable names!")
}
data_all <- data_all[complete.cases(data_all), ]

# Step 1: Pooled partial correlations (consistent with B2.2 format)
# PartialCorr(AI_work, AI_post | PCI_work)
r1 <- residuals(lm(AI_work ~ PCI_work, data = data_all))
r2 <- residuals(lm(AI_post ~ PCI_work, data = data_all))
partial_corr_AI_pooled <- cor(r1, r2)

# PartialCorr(PCI_work, PCI_post | AI_work)
r3 <- residuals(lm(PCI_work ~ AI_work, data = data_all))
r4 <- residuals(lm(PCI_post ~ AI_work, data = data_all))
partial_corr_PCI_pooled <- cor(r3, r4)

N_pooled <- nrow(data_all)

# Output Pooled results
cat("=== B1.1a Pooled Partial Correlation Results ===\n")
cat("PartialCorr(AI_work, AI_post | PCI_work) =", round(partial_corr_AI_pooled, 4), "\n")
cat("PartialCorr(PCI_work, PCI_post | AI_work) =", round(partial_corr_PCI_pooled, 4), "\n")
cat("Sample Size N =", N_pooled, "\n\n")

# Step 2: Within-Firm partial correlations (consistent with B2.2 format)
# Firm-level demeaning (core logic from B3)
AI_work_raw <- data_all$AI_work
AI_post_raw <- data_all$AI_post
PCI_work_raw <- data_all$PCI_work
PCI_post_raw <- data_all$PCI_post

AI_work_mean <- ave(AI_work_raw, data_all$rcid, FUN=function(x) mean(x, na.rm=T))
AI_post_mean <- ave(AI_post_raw, data_all$rcid, FUN=function(x) mean(x, na.rm=T))
PCI_work_mean <- ave(PCI_work_raw, data_all$rcid, FUN=function(x) mean(x, na.rm=T))
PCI_post_mean <- ave(PCI_post_raw, data_all$rcid, FUN=function(x) mean(x, na.rm=T))

data_all$AI_work_W <- AI_work_raw - AI_work_mean
data_all$AI_post_W <- AI_post_raw - AI_post_mean
data_all$PCI_work_W <- PCI_work_raw - PCI_work_mean
data_all$PCI_post_W <- PCI_post_raw - PCI_post_mean

# Within-Firm partial correlations calculation
r1_W <- residuals(lm(AI_work_W ~ PCI_work_W, data = data_all))
r2_W <- residuals(lm(AI_post_W ~ PCI_work_W, data = data_all))
partial_corr_AI_within <- cor(r1_W, r2_W)

r3_W <- residuals(lm(PCI_work_W ~ AI_work_W, data = data_all))
r4_W <- residuals(lm(PCI_post_W ~ AI_work_W, data = data_all))
partial_corr_PCI_within <- cor(r3_W, r4_W)

N_within <- nrow(data_all)
N_firms <- length(unique(data_all$rcid))

# Output Within-Firm results
cat("=== B1.1b Within-Firm Partial Correlation Results ===\n")
cat("PartialCorr(AI_work_W, AI_post_W | PCI_work_W) =", round(partial_corr_AI_within, 4), "\n")
cat("PartialCorr(PCI_work_W, PCI_post_W | AI_work_W) =", round(partial_corr_PCI_within, 4), "\n")
cat("Sample Size N =", N_within, "\n")
cat("Number of Firms =", N_firms, "\n")







# Tbele B2
# B2.1a: Pearson Correlation Matrix
library(dplyr)
library(readr)

PCI <- read_csv("input.csv")  

# Calculate 3×3 Pearson correlation matrix
PCI %>% 
  select(AIC_stock_it_2_w, PCIC_stock_it_2_w, AIPCI_stock_it_2_w) %>%
  na.omit() %>%
  cor(method = "pearson") %>%
  round(4) %>%
  print()







# B2.1b: Within-Firm Correlation Matrix
library(dplyr)
library(readr)

PCI <- read_csv("input.csv")  

# Filter core variables and remove NA values
PCI_filtered <- PCI %>%
  select(rcid, AIC_stock_it_2_w, PCIC_stock_it_2_w, AIPCI_stock_it_2_w) %>%
  na.omit()

# Firm-level demeaning
AIC_mean <- ave(PCI_filtered$AIC_stock_it_2_w, PCI_filtered$rcid, FUN = function(x) mean(x, na.rm = TRUE))
PCIC_mean <- ave(PCI_filtered$PCIC_stock_it_2_w, PCI_filtered$rcid, FUN = function(x) mean(x, na.rm = TRUE))
AIPCI_mean <- ave(PCI_filtered$AIPCI_stock_it_2_w, PCI_filtered$rcid, FUN = function(x) mean(x, na.rm = TRUE))

PCI_filtered <- PCI_filtered %>%
  mutate(
    AIC_dm = AIC_stock_it_2_w - AIC_mean,
    PCIC_dm = PCIC_stock_it_2_w - PCIC_mean,
    AIPCI_dm = AIPCI_stock_it_2_w - AIPCI_mean
  )

# Calculate and print within-firm correlation matrix
PCI_filtered %>% 
  select(AIC_dm, PCIC_dm, AIPCI_dm) %>% 
  na.omit() %>% 
  cor() %>% 
  round(4) %>% 
  print()







# B2.2a: Pooled Partial Correlations
library(dplyr)
library(readr)
library(lmtest)

PCI <- read_csv("input.csv")

# Define core winsorized variables
data <- PCI %>%
  select(
    AIPCI = AIPCI_stock_it_2_w,
    AIC = AIC_stock_it_2_w,
    PCIC = PCIC_stock_it_2_w
  ) %>%
  na.omit()

# Calculate pooled partial correlations
# PartialCorr(AIPCI, AIC | PCIC)
r1 <- residuals(lm(AIPCI ~ PCIC, data = data))
r2 <- residuals(lm(AIC ~ PCIC, data = data))
partial_corr1 <- cor(r1, r2)

# PartialCorr(AIPCI, PCIC | AIC)
r3 <- residuals(lm(AIPCI ~ AIC, data = data))
r4 <- residuals(lm(PCIC ~ AIC, data = data))
partial_corr2 <- cor(r3, r4)

N <- nrow(data)

cat("=== B2.2a Pooled Partial Correlation Results ===\n")
cat("PartialCorr(AIPCI, AIC | PCIC) =", round(partial_corr1, 4), "\n")
cat("PartialCorr(AIPCI, PCIC | AIC) =", round(partial_corr2, 4), "\n")
cat("Sample Size N =", N, "\n")







# B2.2b: Within-Firm Partial Correlations
library(dplyr)
library(readr)

PCI <- read_csv("input.csv")

# Filter core variables and remove NA values
PCI_filtered <- PCI %>%
  select(rcid, AIPCI_stock_it_2_w, AIC_stock_it_2_w, PCIC_stock_it_2_w) %>%
  na.omit()

# Rename core variables
data <- PCI_filtered %>%
  mutate(
    AIPCI = AIPCI_stock_it_2_w,
    AIC = AIC_stock_it_2_w,
    PCIC = PCIC_stock_it_2_w
  ) %>%
  select(rcid, AIPCI, AIC, PCIC)

# Firm-level demeaning
AIPCI_mean <- ave(data$AIPCI, data$rcid, FUN = function(x) mean(x, na.rm = TRUE))
AIC_mean <- ave(data$AIC, data$rcid, FUN = function(x) mean(x, na.rm = TRUE))
PCIC_mean <- ave(data$PCIC, data$rcid, FUN = function(x) mean(x, na.rm = TRUE))

data_demeaned <- data %>%
  mutate(
    AIPCI_dm = AIPCI - AIPCI_mean,
    AIC_dm = AIC - AIC_mean,
    PCIC_dm = PCIC - PCIC_mean
  ) %>%
  na.omit()

# Calculate within-firm partial correlations
# PartialCorr(AIPCI_dm, AIC_dm | PCIC_dm)
r1_dm <- residuals(lm(AIPCI_dm ~ PCIC_dm, data = data_demeaned))
r2_dm <- residuals(lm(AIC_dm ~ PCIC_dm, data = data_demeaned))
partial_corr1_dm <- cor(r1_dm, r2_dm)

# PartialCorr(AIPCI_dm, PCIC_dm | AIC_dm)
r3_dm <- residuals(lm(AIPCI_dm ~ AIC_dm, data = data_demeaned))
r4_dm <- residuals(lm(PCIC_dm ~ AIC_dm, data = data_demeaned))
partial_corr2_dm <- cor(r3_dm, r4_dm)

N_filtered <- nrow(data_demeaned)
cat("=== B2.2b Within-Firm Partial Correlation Results (Filtered) ===\n")
cat("Filtered Sample Size N =", N_filtered, "\n")
cat("PartialCorr(AIPCI_dm, AIC_dm | PCIC_dm) =", round(partial_corr1_dm, 4), "\n")
cat("PartialCorr(AIPCI_dm, PCIC_dm | AIC_dm) =", round(partial_corr2_dm, 4), "\n")









# B2.3 Residualized AIPCI Check
library(dplyr)
library(readr)
library(lfe)

# Load data and extract core variables
PCI <- read_csv("input.csv") %>%
  select(rcid, year_quarter,
         AIPCI = AIPCI_stock_it_2_w,
         AIC = AIC_stock_it_2_w,
         PCIC = PCIC_stock_it_2_w) %>%
  na.omit()

# B2.3a: Two-way FE regression + residual extraction + correlation calculation
fe_reg <- felm(AIPCI ~ AIC | rcid + year_quarter, data = PCI)
PCI$AIPCI_netAI <- residuals(fe_reg)

# Pooled correlation
cor_pooled <- cor(PCI$AIPCI_netAI, PCI$PCIC)

# Within-firm correlation
AIPCI_netAI_mean <- ave(PCI$AIPCI_netAI, PCI$rcid, FUN = function(x) mean(x, na.rm = TRUE))
PCIC_mean <- ave(PCI$PCIC, PCI$rcid, FUN = function(x) mean(x, na.rm = TRUE))

PCI <- PCI %>%
  mutate(
    AIPCI_netAI_dm = AIPCI_netAI - AIPCI_netAI_mean,
    PCIC_dm = PCIC - PCIC_mean
  )

cor_within <- cor(PCI$AIPCI_netAI_dm, PCI$PCIC_dm)

# B2.3b: Extract R²
r2 <- as.numeric(summary(fe_reg)$r2)

# Output results
cat("=== B2.3 Results ===\n")
cat("B2.3a Pooled Corr(AIPCI_netAI, PCIC) =", round(cor_pooled,4), "\n")
cat("B2.3a Within-firm Corr =", round(cor_within,4), "\n")
cat("B2.3b Within R² =", round(r2,4), "\n")
cat("Sample Size N =", nrow(PCI), "\n")






# Table B3: Variance Decomposition (AIC/PCIC/AIPCI) - SD & Economic Significance
library(dplyr)
library(readr)

# Load data and extract core winsorized variables
PCI <- read_csv("input.csv") %>%
  select(rcid,
         AIC   = AIC_stock_it_2_w,
         PCIC  = PCIC_stock_it_2_w,
         AIPCI = AIPCI_stock_it_2_w) %>%
  na.omit()

# Variance decomposition function
decomp <- function(df, varname) {
  # Calculate firm mean and within-firm deviation
  firm_mean <- tapply(df[[varname]], df$rcid, mean, na.rm = TRUE)
  within_deviation <- df[[varname]] - firm_mean[df$rcid]
  
  # Variance and standard deviation calculation
  var_total  <- var(df[[varname]], na.rm = TRUE)
  var_within <- var(within_deviation, na.rm = TRUE)
  var_between <- var(firm_mean, na.rm = TRUE)
  
  sd_total  <- sd(df[[varname]], na.rm = TRUE)
  sd_within <- sd(within_deviation, na.rm = TRUE)
  sd_between <- sd(firm_mean, na.rm = TRUE)
  
  # Variance share percentage
  between_share_pct <- 100 * var_between / var_total
  within_share_pct  <- 100 * var_within / var_total
  
  # Return results dataframe
  return(data.frame(
    Capability = varname,
    Var_total = var_total,
    Var_between = var_between,
    Var_within = var_within,
    Between_share_pct = between_share_pct,
    Within_share_pct = within_share_pct,
    Total_SD = sd_total,
    Between_SD = sd_between,
    Within_SD = sd_within,
    N_firm_quarters = nrow(df),
    N_firms = length(unique(df$rcid))
  ))
}

# Calculate variance decomposition for AIC/PCIC/AIPCI
results_AIC <- decomp(PCI, "AIC")
results_PCIC <- decomp(PCI, "PCIC")
results_AIPCI <- decomp(PCI, "AIPCI")

# Combine and print final results
final_results <- rbind(results_AIC, results_PCIC, results_AIPCI)
print(final_results)




# Demo

library(zoo)
library(dplyr)
library(readr)
library(readxl)
library(DescTools)
library(Hmisc)
library(plyr)
library(psych)
library("fixest")
library("texreg")

PCI <- read_csv("input.csv")

PCI %>%
  distinct(rcid) %>% 
  sample_n(size = 100) %>%  
  inner_join(PCI, by = "rcid") %>% 
  write_csv("output.csv") 


