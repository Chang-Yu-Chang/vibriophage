


clean_ode <- function (out) {
    as_tibble(out) %>%
        pivot_longer(-time) %>%
        mutate(time = as.numeric(time), value = as.numeric(value))
}


tb_types <- tibble(
    type = c(rep("human", 3), rep("microbe", 2), rep("plankton", 3)),
    name = c("S", "I", "R", "V", "P", "Z_F", "Z_V", "X")
)


plot_one_out <- function (out, p_title = "") {

    p1 <- out %>%
        filter(type == "human") %>%
        mutate(name = factor(name, c("S", "I", "R"))) %>%
        ggplot() +
        geom_line(aes(x = time, y = value, linetype = name)) +
        scale_linetype_manual(values = c(S = 1, I = 2, R = 3)) +
        #facet_nested_wrap(~type + name, scales = "free_y") +
        coord_cartesian(clip = "off") +
        theme_bw() +
        theme(
            legend.position = "top",
            legend.title = element_blank()
        ) +
        guides() +
        labs()

    p2 <- out %>%
        filter(type == "microbe") %>%
        mutate(name = factor(name, c("V", "P"))) %>%
        ggplot() +
        geom_line(aes(x = time, y = value, color = name)) +
        scale_color_manual(values = c(V = "royalblue4", P = "maroon")) +
        scale_y_log10() +
        coord_cartesian(clip = "off") +
        theme_bw() +
        theme(
            legend.position = "top",
            legend.title = element_blank()
        ) +
        guides() +
        labs()

    p3 <-  out %>%
        filter(type == "plankton") %>%
        mutate(name = factor(name, c("Z_F", "Z_V", "X"))) %>%
        ggplot() +
        geom_line(aes(x = time, y = value, color = name)) +
        scale_color_manual(values = c(X = "seagreen", Z_F = "grey60", Z_V = "grey20")) +
        scale_y_log10() +
        coord_cartesian(clip = "off") +
        theme_bw() +
        theme(
            legend.position = "top",
            legend.title = element_blank()
        ) +
        guides() +
        labs()

    p_lower <- plot_grid(p1, p2, p3, nrow = 1, align = "h", axis = "tb")
    p <- plot_grid(ggdraw() + draw_text(p_title, hjust = 0, x = 0.05), p_lower, ncol = 1, rel_heights = c(.1, 1)) +
        theme(plot.background = element_rect(fill = "white", color = NA))
    return(p)
}
