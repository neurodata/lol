require(ggplot2)
require(lol)
require(reshape2)
require(Rmisc)
require(randomForest)
require(gridExtra)
require(latex2exp)

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

performance <- readRDS('simulations.rds')

plot_example <- function(data, sim) {
  p1 <- ggplot(data[data$classification=='lda',], aes(x=dimensions, y=lhat, color=algorithm, group=algorithm)) +
    stat_summary(geom="line", fun.y="mean", size=2) +
    xlab("Dimensions") +
    ylab(TeX("$\\hat{L}$")) +
    ggtitle(paste(sim, "Simulation, LDA Classifier")) +
    scale_color_discrete(name="Algorithm") +
    theme_bw()
  p2 <- ggplot(data[data$classification=='rf',], aes(x=dimensions, y=lhat, color=algorithm, group=algorithm)) +
    stat_summary(geom="line", fun.y="mean", size=2) +
    xlab("Dimensions") +
    ylab(TeX("$\\hat{L}$")) +
    ggtitle(paste(sim, "Simulation, RF Classifier")) +
    theme_bw() +
    scale_color_discrete(name="Algorithm")

  my_legend <- g_legend(p1)
  p3 <- grid.arrange(arrangeGrob(p1 + theme(legend.position=NaN), p2 + theme(legend.position=NaN), nrow=2), my_legend, nrow=1, widths=c(.88, .12))

}

for (sim in unique(performance$simulation)) {
  subset <- performance[performance$simulation == sim,]
  plotlist <- list()
  undims <- unique(subset$dimensions)
  plot_example(subset, sim)
}
