#PERFORMING ANOVA ON EMMA ANALYSIS VARIABLES

#read in data
isotopes <- read.csv("/Users/dariasmac/Desktop/Watershed Biogeochemistry/ESCI 431_Water_Isotopes-2 (2).csv")

isotopes$d18O....VSMOW.<-as.numeric(isotopes$d18O....VSMOW.) # Make sure your variables are numeric and not character strings.

boxplot(CDBR) # View the dataset as a box plot. 

#You are already familiar with one visualization that will help you assess normality, the histogram. The only trick here is that you need to make a separate histogram for each group and it is nice to have them on the same page, and with the same scale on the horizontal axis, so you can compare them.

#categorize data by end-member
d18O <- isotopes$d18O....VSMOW.
snow <- d18O[288:311]
ice <- d18O[313:322]
glacier_melt <- d18O[327:334]
groundwater <- d18O[337:344]
rain_a <- d18O[346:348]
rain_b <- d18O[349:356]

#plot histograms to check for normality  
par(mfrow = c(2,3))# par sets graphical parameters. (2,3) creates a 2x3 array of graphs so up to 6 histograms can be viewed in one plot window

hist(snow)
hist(ice)
hist(glacier_melt)
hist(groundwater)
hist(rain_a)
hist(rain_b)


#use shapiro.test to test for normality Here we are choosing an alpha = 0.05. The null hypothesis for the Shapiro Wilk test is that the data are normally distributed. 
#If the p value is low, there is evidence to reject the null hypothesis of the Shapiro-Wilk normality test. If that is the case, then we should try transforming the data. 
shapiro.test(snow)
shapiro.test(ice)
shapiro.test(glacier_melt)
shapiro.test(groundwater)
shapiro.test(rain_a)
shapiro.test(rain_b)

#make all variables same length to turn into a matrix then into a data frame
isotope_list <- list(snow, ice, glacier_melt, groundwater, rain_a, rain_b)
maxl <- max(sapply(isotope_list, length))
sapply(isotope_list, FUN = function(x, ml) {
       difference <- ml - length(x)
       c(x, rep(NA, difference))
   }, ml = maxl, simplify = FALSE)

ice_1 <- c(ice, rep(NA, maxl - length(ice)))
glacier_melt_1 <- c(glacier_melt, rep(NA, maxl - length(glacier_melt)))
groundwater_1 <- c(groundwater, rep(NA, maxl - length(groundwater)))
rain_a_1 <- c(rain_a, rep(NA, maxl - length(rain_a)))
rain_b_1 <- c(rain_b, rep(NA, maxl - length(rain_b)))

#turn into matrix
iso_mat <- cbind(snow, ice_1, glacier_melt_1, groundwater_1, rain_a_1, rain_b_1)

#turn into data frame
iso_frame <- data.frame(iso_mat)

#turn into long form data so you can use aov()
library(tidyr)
data_long <- gather(iso_frame, location, d18O, snow:rain_b_1, factor_key=TRUE)

#perform anova
iso_aov <- aov(d18O ~ location, data = data_long)

##             Df Sum Sq Mean Sq F value   Pr(>F)    
## location     5 178.16   35.63   33.49 1.57e-15 ***
## Residuals   55  58.51    1.06                     
## ---
##  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
## 83 observations deleted due to missingness

#perform kruskall-wallis test
kruskal.test(d18O ~ location, data = data_long)

##	Kruskal-Wallis rank sum test
##
## data:  d18O by location
## Kruskal-Wallis chi-squared = 30.899, df = 5, p-value = 9.808e-06

#DO THIS AS NEEDED, DIDNT WORK IN THIS CASE: log transform non-normal data (in this case rain_b_1 is not normal, can't use log transformation because values are negative)
iso_frame$rain_b_sqrt <- sqrt(iso_frame$rain_b_1)
#perform shapiro.test again to test for normality after transforming


#DO AS NEEDED: exclude rain_b_1:
iso_mat_2 <- cbind(snow, ice_1, glacier_melt_1, groundwater_1, rain_a_1)

iso_frame_2 <- data.frame(iso_mat_2)

library(tidyr)
data_long_2 <- gather(iso_frame_2, location, d18O, snow:rain_a_1, factor_key=TRUE)

iso_aov_2 <- aov(d18O ~ location, data = data_long_2)

summary(iso_aov_2)

##             Df Sum Sq Mean Sq F value  Pr(>F)    
## location     4  39.89   9.972   8.625 2.5e-05 ***
## Residuals   48  55.50   1.156                    
## ---
##   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
## 67 observations deleted due to missingness
