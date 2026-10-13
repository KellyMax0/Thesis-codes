suppressPackageStartupMessages({
  library(survival); library(ggplot2); library(dplyr); library(readr); library(patchwork)
})

# Load data
surv_dat <- read_csv("patient_groups.csv", show_col_types = FALSE) %>%
  filter(survival_available == "yes", !is.na(os_months), !is.na(os_event), os_months > 0) %>%
  mutate(os_event = as.integer(os_event)) %>%
  select(patient_id, group, os_months, os_event)
expr_sub <- read_csv("smarca4_split_normalised/expression_subgroups_deseq2.csv", show_col_types = FALSE) %>%
  filter(gene == "ARID1A") %>%
  transmute(patient_id, group = paste(gene, expression_group))
sub_dat <- surv_dat %>% select(-group) %>% inner_join(expr_sub, by = "patient_id")
surv_all <- bind_rows(surv_dat, sub_dat)

# Groups and colours
CBIO_COLORS <- c("SWI/SNF WT" = "#FF9900", "ARID1A Higher" = "#DC3912", "ARID1A Lower" = "#109618")
GROUP_NAMES <- c("SWI/SNF WT" = "SWI WT",
                 "ARID1A Higher" = "ARID1A HIGHER EXPR", "ARID1A Lower" = "ARID1A LOWER EXPR")

# Plot settings
BASE_FONT  <- if (Sys.info()[["sysname"]] == "Darwin") "Helvetica" else "Liberation Sans"
PANEL_W    <- 6.8
PANEL_H    <- 4.0
RISK_ROW_H <- 0.183
TXT_PT     <- 9.0
Y_TOP      <- 0.046
RISK_LAB_IN <- 1.47
MAX_GROUPS <- 5
FIG_W      <- 12.5
FIG_H      <- 6.3
TXT_COL    <- "grey20"

# Kaplan-Meier plot function
km_cbioportal <- function(data, groups_keep, comparison_label, outfile,
                          names = GROUP_NAMES, colours = CBIO_COLORS) {
  d <- data %>% filter(group %in% groups_keep) %>%
    mutate(group = factor(group, levels = groups_keep))
  fit <- survfit(Surv(os_months, os_event) ~ group, data = d)
  lr  <- survdiff(Surv(os_months, os_event) ~ group, data = d)
  p   <- pchisq(lr$chisq, df = length(lr$n) - 1, lower.tail = FALSE)
  p_lab <- if (p < 1e-4) "Logrank Test P-Value: < 0.0001" else sprintf("Logrank Test P-Value: %.3g", p)
  leg_labels <- setNames(paste0("(", LETTERS[seq_along(groups_keep)], ") ",
                                names[groups_keep]), groups_keep)

  s  <- summary(fit, censored = TRUE)
  km <- data.frame(time = s$time, surv = s$surv, n.censor = s$n.censor,
                   group = sub("^group=", "", s$strata))
  km <- bind_rows(km, data.frame(time = 0, surv = 1, n.censor = 0, group = groups_keep)) %>%
    mutate(group = factor(group, levels = groups_keep))
  cens <- km %>% filter(n.censor > 0)

  xmax      <- max(d$os_months)
  x_breaks  <- seq(0, floor(xmax / 10) * 10, by = 10)
  last_tick <- max(x_breaks)
  risk_times <- seq(0, last_tick, by = 20)
  if (last_tick %% 20 != 0)
    risk_times <- c(risk_times[risk_times <= last_tick - 20], last_tick)
  x_scale <- scale_x_continuous(limits = c(0, xmax), breaks = x_breaks,
                                expand = expansion(0, 0))

  p_main <- ggplot(km, aes(x = time, y = surv, color = group)) +
    geom_step(linewidth = 0.35) +
    geom_point(data = cens, shape = 3, size = 1.5, stroke = 0.9) +
    scale_color_manual(values = colours[groups_keep], labels = leg_labels, name = p_lab) +
    x_scale +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1),
                       labels = function(x) paste0(round(x * 100), "%"),
                       expand = expansion(mult = c(0, Y_TOP))) +
    coord_cartesian(clip = "off") +
    labs(x = "Overall Survival (Months)", y = "Probability of Overall Survival") +
    guides(color = guide_legend(override.aes = list(shape = 15, size = 3.2,
                                                    linewidth = 0, linetype = "blank"))) +
    theme_classic(base_size = TXT_PT, base_family = BASE_FONT) +
    theme(axis.line   = element_line(linewidth = 0.42, colour = "black"),
          axis.ticks  = element_line(linewidth = 0.42, colour = "black"),
          axis.ticks.length = unit(0.042, "in"),
          axis.text   = element_text(size = TXT_PT, colour = "black"),
          axis.text.x = element_text(margin = margin(t = 2.7)),
          axis.text.y = element_text(margin = margin(r = 2.7)),
          axis.title  = element_text(size = TXT_PT, colour = TXT_COL),
          axis.title.x = element_text(margin = margin(t = 12.5)),
          axis.title.y = element_text(margin = margin(r = 26)),
          legend.position = "right", legend.justification = "top",
          legend.title = element_text(size = TXT_PT, colour = TXT_COL,
                                      margin = margin(b = 8)),
          legend.text  = element_text(size = TXT_PT, colour = TXT_COL,
                                      margin = margin(l = 7)),
          legend.key = element_blank(),
          legend.key.height = unit(0.26, "in"), legend.key.width = unit(0.12, "in"),
          legend.box.spacing = unit(0.06, "in"),
          legend.margin = margin(0, 0, 0, 0),
          plot.margin = margin(5.5, 5.5, 0, 5.5))

  sr <- summary(fit, times = risk_times, extend = TRUE)
  n_groups <- length(groups_keep)
  risk_df <- data.frame(time = sr$time, n = sr$n.risk,
                        group = sub("^group=", "", sr$strata)) %>%
    mutate(row = n_groups + 1 - match(group, groups_keep))

  lab_w <- suppressWarnings(max(sapply(leg_labels, function(l) grid::convertWidth(grid::grobWidth(
    grid::textGrob(l, gp = grid::gpar(fontsize = TXT_PT, fontfamily = BASE_FONT))),
    "inches", valueOnly = TRUE))))
  lab_in <- max(RISK_LAB_IN, lab_w + 0.5)
  x_lab <- -lab_in / PANEL_W * xmax
  lab_df <- data.frame(row = n_groups:1, lab = leg_labels[groups_keep])

  p_risk <- ggplot(risk_df, aes(x = time, y = row, label = n)) +
    geom_text(size = TXT_PT / .pt, family = BASE_FONT, colour = "black") +
    geom_text(data = lab_df, aes(x = x_lab, y = row, label = lab), inherit.aes = FALSE,
              hjust = 0, size = TXT_PT / .pt, family = BASE_FONT, colour = "black") +
    annotate("text", x = x_lab, y = n_groups + 1, label = "Number at risk (n)",
             hjust = 0, fontface = "bold", size = TXT_PT / .pt, family = BASE_FONT) +
    scale_x_continuous(expand = expansion(0, 0)) +
    scale_y_continuous(expand = expansion(0, 0)) +
    coord_cartesian(xlim = c(0, xmax), ylim = c(0.5, n_groups + 1.5), clip = "off") +
    theme_void(base_family = BASE_FONT) +
    theme(plot.margin = margin(4, 5.5, 5.5, lab_in * 72 + 6))

  p_out <- (p_main / p_risk / plot_spacer()) +
    plot_layout(guides = "collect",
                widths  = unit(PANEL_W, "in"),
                heights = unit(c(PANEL_H * (1 + Y_TOP), RISK_ROW_H * (n_groups + 1),
                                 RISK_ROW_H * (MAX_GROUPS - n_groups)), "in")) &
    theme(legend.position = "right", legend.justification = "top")
  ggsave(outfile, p_out, width = FIG_W, height = FIG_H, dpi = 300, bg = "white")

  tab <- as.data.frame(summary(fit)$table)
  data.frame(
    comparison = comparison_label,
    group = names[sub("^group=", "", rownames(tab))],
    n = tab$records, events = tab$events,
    median_os_months = round(tab$median, 1),
    median_95lcl = round(tab$`0.95LCL`, 1), median_95ucl = round(tab$`0.95UCL`, 1),
    logrank_p = signif(p, 3), row.names = NULL)
}

# Comparisons
WT <- "SWI/SNF WT"
summary_tab <- bind_rows(
  km_cbioportal(surv_all, c(WT, "ARID1A Lower"), "ARID1A lower expr vs SWI WT",
                "figS10A_km_arid1a_lower_vs_swi_wt.png"),
  km_cbioportal(surv_all, c(WT, "ARID1A Higher"), "ARID1A higher expr vs SWI WT",
                "figS10B_km_arid1a_higher_vs_swi_wt.png"))
print(summary_tab)
write_csv(summary_tab, "arid1a_km_summary.csv")
