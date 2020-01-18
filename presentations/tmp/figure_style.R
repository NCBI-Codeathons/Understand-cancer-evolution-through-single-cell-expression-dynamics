# figure styling
customPlot = list(
    theme_light(),
    theme(plot.title = element_text(size = 14, 
                                    face = "bold", 
                                    margin = margin(10, 0, 10, 0)),
          axis.text.x = element_text(size = 10, face = "bold"), 
          axis.text.y = element_text(size = 10, face = "bold"), 
          axis.title.x = element_text(size = 12, face = "bold"), 
          axis.title.y = element_text(size = 12, face = "bold"),
          axis.line = element_line(colour = "black", size = 0.5),
          panel.grid.major = element_line(), 
          panel.grid.minor = element_line(),
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 10, face = "bold"),
          # plot.margin = unit(c(0.5,0.9,0.3,0.7), "cm"))
          plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
)

# plotting
res_plot = ggplot(df) +
    geom_point(aes(x = df[,x], y = df[,y]),
               col = "black",
               alpha = 1/3) + 
    customLabs +
    # geom_abline(col = "grey", lwd = 0.5) + 
    customPlot

# add concordance
library(cowplot)
if (cor == TRUE) {
    cor = round(cor(df[,x], df[,y], use = "complete.obs"), 2)
    final_plot = ggdraw(add_sub(res_plot, paste("Concordance\ncorrelation =", cor), 
                                vpadding = grid::unit(0, "lines"),
                                y = 11, x = 0.05,
                                hjust = 0, 
                                fontface = "italic",
                                size = 12))
} else {final_plot = res_plot}


print(final_plot)