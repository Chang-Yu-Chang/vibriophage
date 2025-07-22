#'

library(tidyverse)
library(ggh4x)
library(cowplot)
library(ggrepel)
library(ggsci)
source(here::here("tools.R"))
load(file = here::here("data/mains.rdata"))

tb_sc <- tibble(
    scenario = factor(1:4),
    sc_name = c("Scenario 1: SIRB", "Scenario 2: Phage", "Scenario 3: Plankton", "Scenario 4: Full model")
)

## Infection peak ----
tb <- tb_results %>%
    select(scenario, tb) %>%
    unnest(tb) %>%
    filter(type == "human", name == "I") %>%
    group_by(scenario, name) %>%
    left_join(tb_sc)
tb_peaks <- tb %>%
    filter(value == max(value))
p1 <- tb %>%
    ggplot() +
    geom_line(aes(x = time, y = value, color = sc_name), linewidth = 1) +
    geom_text_repel(data = tb_peaks, aes(x = time, y = value, label = sc_name, color = sc_name), nudge_x = 15, nudge_y = 1000, vjust = 1) +
    scale_color_npg() +
    scale_x_continuous(limits = c(0, 50)) +
    scale_y_continuous(limits = c(0, 10500)) +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        panel.grid = element_blank()
    ) +
    guides(color = "none") +
    labs(x = "Time", y = "Number of infected (I)")
p2 <- tb_peaks %>%
    mutate(sc_name = factor(sc_name, rev(tb_sc$sc_name))) %>%
    ggplot() +
    geom_col(aes(x = sc_name, y = time, fill = sc_name)) +
    scale_fill_npg() +
    coord_flip(clip = "off") +
    theme_bw() +
    theme(
        axis.text.y = element_text(color = pal_npg()(4)),
        panel.grid = element_blank(),
        plot.background = element_blank()
    ) +
    guides(fill = "none") +
    labs(x = "", y = "Time to peak")

p <- ggdraw(p1) +
    draw_plot(p2, .45, .4, .5, .3)
ggsave(here::here("plots/infected.png"), p, width = 5, height = 5)
ggsave(here::here("plots/infected.pdf"), p, width = 5, height = 5)

## Full scenarios ----
tb <- tb_results %>%
    select(scenario, tb) %>%
    unnest(tb) %>%
    group_by(scenario, name) %>%
    left_join(tb_sc)

p <- tb %>%
    ggplot() +
    geom_line(aes(x = time, y = value, color = name, linetype = name), linewidth = .5) +
    scale_color_manual(breaks = c("S", "I", "R", "V", "P", "X", "Z_F", "Z_V"), values = c(S = 1, I = 1, R = 1, V = "royalblue4", P = "maroon", X = "seagreen", Z_F = "grey60", Z_V = "grey20")) +
    scale_linetype_manual(breaks = c("S", "I", "R", "V", "P", "X", "Z_F", "Z_V"), values = c(S = 1, I = 2, R = 3, V = 1, P = 1, X = 1, Z_F = 1, Z_V = 1, X = 1)) +
    facet_grid2(sc_name ~ type, independent = "y", scales = "free_y", switch = "y",
                strip = strip_themed(
                    background_y = list(element_rect(fill = pal_npg()(4)[1]), element_rect(fill = pal_npg()(4)[2]), element_rect(fill = pal_npg()(4)[3]), element_rect(fill = pal_npg()(4)[4])),
                    background_x = element_blank()
    )) +
    scale_y_facet(type == "microbe", type = "log10") +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        legend.title = element_blank(),
        legend.position = "top",
        panel.grid = element_blank(),
        strip.placement = "outside",
        axis.title.y = element_blank()
    ) +
    guides(
        color = guide_legend(nrow = 1)
    ) +
    labs()

ggsave(here::here("plots/mains_full.png"), p, width = 6, height = 8)
ggsave(here::here("plots/mains_full.pdf"), p, width = 6, height = 8)

## vibrio decay rate ----
tb <- tb_results %>%
    select(scenario, tb) %>%
    unnest(tb) %>%
    group_by(scenario, name) %>%
    left_join(tb_sc) %>%
    filter(name == "V")
tb_fin <- tb %>%
    filter(time == max(time)) %>%
    bind_cols(vjust = c(1, -7, -1, -1))

p1 <- tb %>%
    ggplot() +
    geom_line(aes(x = time, y = value, color = sc_name), linewidth = 1) +
    geom_text_repel(data = tb_fin, aes(x = time, y = value, label = sc_name, color = sc_name, vjust = vjust), hjust = 0) +

    scale_color_npg() +
    scale_x_continuous(limits = c(0, 100)) +
    scale_y_log10(limits = c(1e-6, 1e6)) +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        panel.grid = element_blank()
    ) +
    guides(color = "none") +
    labs(x = "Time", y = "Number of Vibrio (V)")


p2 <- tb %>%
    filter(value == max(value)) %>%
    ggplot() +
    geom_col(aes(x = sc_name, y = value, fill = sc_name)) +
    scale_y_log10() +
    scale_fill_npg() +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
        axis.text.x = element_text(color = pal_npg()(4), angle = 25, hjust = 1),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_blank()
    ) +
    guides(fill = "none") +
    labs(y = "Peak value")

p <- ggdraw(p1) +
    draw_plot(p2, .25, .08, .3, .35)
ggsave(here::here("plots/vibrio.png"), p, width = 6, height = 6)
ggsave(here::here("plots/vibrio.pdf"), p, width = 6, height = 6)

