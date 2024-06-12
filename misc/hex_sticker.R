

library(tidyverse)

hex_plot_data <-
  readRDS('../publication/06_data/wang_panel_B_dat.rds') %>%
  filter(str_detect(sample, '3792')) %>%
  select(state, exp_baf, baf)

plot <-
  hex_plot_data %>%
  ggplot(aes(state, baf)) +
  geom_hline(
    data =
      hex_plot_data %>%
      select(state, exp_baf) %>%
      distinct(),
    aes(yintercept = exp_baf),
    col = 'black',
    lty = 2
  ) +
  geom_violin(
    scale = 'width',
    col = "#ECE953",
    fill = "#ECE953",
    linewidth = 0.5,
    width = 1.2
  ) +
  scale_y_continuous(
    trans = scales::asn_trans(),
    limits = c(0,1),
    breaks = c(seq(0.1, 0.9, by = 0.2))
  ) +
  theme_void()

# Create a hex sticker

hexSticker::sticker(
  subplot = plot,
  package = "TBtypeR",
  p_color = "#ECE953",
  p_family = 'sans',
  p_fontface = 'bold',
  p_x = 1.02,
  p_y = 1.5,
  p_size = 23,
  s_x = 1,
  s_y = 0.875,
  s_width = 1.4,
  s_height = 0.9,
  h_fill = "#2372B9",
  h_color = "#49A942",
  spotlight = TRUE,
  l_x = 0.95,
  l_y = 0.2,
  filename = "TBtypeR_hex_logo.png"
)

