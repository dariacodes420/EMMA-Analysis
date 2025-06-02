#plotting clustering of different water sources (minus the GMWL because I couldn't figure that out lol)
EMMA <- read.csv("/Users/dariasmac/Desktop/Watershed Biogeochemistry/EMMA.csv") #read in dataset
EMMA_lm <- lm(EMMA$dD.H....VSMOW. ~ EMMA$d18O....VSMOW.)
plot(EMMA_lm)
summary(EMMA_lm)
?EMMA
plot(EMMA$dD.H....VSMOW., EMMA$d18O....VSMOW.)
ggplot(gc, aes(x=Clades, y=GC, group=Genes, colour=Genes))
?ggplot
ggplot(EMMA, aes(x = dD.H....VSMOW., y = d18O....VSMOW., group = X.2, colour = X.2)) +
  labs(x = "dD.H", y = "d18O")  +
  geom_point(size=3) +
  geom_abline(aes(intercept = 10, slope = 8, color = "black", linetype = "dashed", size = 1.5))
library(ggplot2)


#plotting elevation against isotopes
EMMA <- read.csv("/Users/dariasmac/Desktop/Watershed Biogeochemistry/EMMA.csv")
plot(EMMA$X.7, EMMA$d18O....VSMOW., type = 'p', lwd = 1.25, xlab = "Elevation (m asl)", ylab = "δ18O Concentration (% VSMOW)", col = "lightsteelblue1", pch = 16)
par(new = TRUE)
par(mar = c(5, 4, 4, 5))
plot(EMMA$X.7, EMMA$dD.H....VSMOW., type = 'p', lwd = 1.25, col = "lightsteelblue4", axes = FALSE, xlab = "", ylab = "", pch = 17)
axis(4)
legend(x = 300, y = -10, c("δ18O (‰ VSMOW)", "δ2H (‰ VSMOW)"), col = c("lightsteelblue1", "lightsteelblue4"), lty = 1, lwd = 2)
axis(4)
mtext(text =, "δ2H Concentration (% VSMOW)", side = 4, line = 3)

