#EMMA for diff months
#Nooksack River Watershed, WA EMMA
#Samples collected 2022 - 2024

############################################################
###                   Read in the raw data               ###
############################################################
###Comments:Cannot have NA values for any of these tests, so make sure no NA, <, or > values. 
###Original code provided by: Taylor Joe Mills (University of Colorado Boulder), 
###Modified by Alana Wilson (University of Colorado Boulder)
###Modifed again by Alia Khan (University of Colorado Boulder)
####Modifed again by Alia Khan (Western Washington University)
############################################################

wd<-setwd("/Users/dariasmac/Desktop/Watershed Biogeochemistry") #change the working directory to the location of your files
#Read in the raw data
#raw data should be in .csv format and have column headers (e.g. names of solutes)
#RawData is CSV format of stream water chemistry data, each row is a different sample, with chemistry data in columns
RawData <- read.csv("River_Nooksack_condensed_All_v2.csv", header=T, na.strings="NA") #Cannot have NA values for any of these tests, so make sure no NA, <, or > values.  Make sure your data has column headers
#Reduce raw data to just tracers, removing columns with sample name, date, elevation, etc.
Tracers<-na.omit(RawData[,c(3:6)]) 
#Columns with potential chemical tracers, here they are pulled out of RawData to form a new dataset of streamwater called 'Tracers'.  Should only include chemistry columns, not date, not time, etc.
#The 'Tracers' file now has 15 columns, one for each solute

#Read in raw data for End Members.  
#End members with multiple samples should have 5 rows (5th, 25th percentile, median, and 75th percentile, 95th)
#Use CSV format for end member data (e.g. precipitation, groundwater, snow melt, etc.), where each row is a different sample and chemistry data is in the same columns as RawData file
#the 'EMall' input file should include all potential end members that you have chemistry data for
EMall<-read.csv("End_Members_Nooksack_avg_v2.csv", header=T, na.strings="NA")

#Reduce End Member file to tracers, removing columns with sample name, date, elevation, etc.
EMTracers<-na.omit(EMall[,c(2:5)])
#Columns with your chemical tracers, here they are pulled out of EMall to form a new dataset called EMTracers. Should only include chemistry data, not date, not time, etc.
#The 'EMTracers' file now has 15 columns, one for each solute and should be the same solutes as 'Tracers' file

NumTracers <- length(Tracers[1,]) #creates a value that is equivalent to the number of tracers, to be used as a count for running loops.

#Standardize the stream tracer data
Tracers.std<-scale(Tracers) #creates a standardized matrix of 'Tracers' data

#Standardize End Member tracer data based on stream chemistry, i changes automatically based on number of variables
EMTracers.std <- scale(EMTracers) #creates a standardized matrix of 'EMTracers' data
for (i in 1:NumTracers) {
  EMTracers.std[,i]<-(EMTracers[,i]-mean(Tracers[,i])) / sd(Tracers[,i])
}
#no need to modify the above loop, because 'i' is defined by 1:NumTracers & NumTracers is already modified to reflect the data set you are working with.  

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
###The following large block of code runs the diagnostic tests on tracers described in Hooper 2003. 
#Reference: Hooper, Richard P. "Diagnostic tools for mixing models of stream water chemistry." Water resources research 39.3 (2003)###
#For explanation of the math, see Hooper, 2003 reference above.
###The purpose of these tests is to evaluate conservative behavior of tracers and determine dimensionality of the mixing space     ###
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

############################################################
###                   Manual PCA                         ###
############################################################
###Comments:This was done manually to ensure no strange behavior in built-in PCA functions
############################################################
Cormatrix <- cor(Tracers.std) #constructs correlation matrix of standardized tracer data
Evectors <- eigen(Cormatrix) #perform eigen decomposition on the correlation matrix
Evectors <- Evectors$vectors #extract eigen vectors

############################################################
###       Evaluate tracers in 1-Dimensional space        ###
############################################################
###Comments: the variable names were used to reflect those used in Hooper, 2003 equation(2). Equation is solved below in steps.  To add dimensions, change the number of vectors retained in V#DT variable
############################################################
Xstar <- Tracers.std  #just renames the Tracers matrix
V <- (Evectors) #this is actually VT (V transposed) b/c R outputs VT not V
V1DT <- as.matrix(V[,1]) #There is no transpose operation because V is actually VT; this is where you change the number of dimensions - 2 dimensions would be as.matrix(V[,1:2])
V1D <- t(V1DT) #transposes V1DT so that we now have V (not V Transposed)
V1DtV1D <- V1D%*%V1DT
InvV1DtV1D <- solve(V1DtV1D)
XstarV1DT <- Xstar%*%V1DT
Proj1Dstd <- XstarV1DT%*%InvV1DtV1D%*%V1D #This solves the Hooper 2003 EQ and gives you the standardized projection
##############################
#Destandardize the projection#
##############################
Proj1D <- Proj1Dstd
for (i in 1:NumTracers) {
  Proj1D[,i]<-(Proj1Dstd[,i]*sd(Tracers[,i])+mean(Tracers[,i])) ###Destandardizes the projection by multiplying by standard deviation of Tracers + mean of Tracers
}
##############################
###Calculate residuals       #
##############################
Resids1D <-Proj1D - Tracers #Calculate residuals where Proj1D is river data projected and de-standardized, and Tracers is original river solute values for each sample 
write.table((cbind(Resids1D,Tracers)), "full1DResids.txt", sep="\t") #Output table of residuals
##############################
#Plot the 1D residuals          #
##############################
par(mfrow=c(3,3)) #creates a window for multiple plots (3x3 = 3 rows by 3 columns), change based on number of tracers
titlevector <- c("dD","NO3", "Al", "K")#label all chemical tracers in their column order.  This is a vector of chemistry data names, typically the same as the header of the Tracer datadet.  MAKE SURE the title vector assigns the proper names to the proper columns!
for (i in 1:NumTracers) {
  lm1D<-(lm(Resids1D[,i]~Tracers[,i]))
  values <- (paste("R2 =",format(summary(lm1D)$adj.r.squared, digits=2),"P-value =",format(summary(lm1D)$coefficients[2,4], digits=2)))
  plot(Tracers[,i],Resids1D[,i],xlab="Measured concentration", ylab="Residuals", main=titlevector[i], sub=values)
  abline(lm1D$coefficients[1],lm1D$coefficients[2])
}

#The residuals should be random if the 1-D space correctly represents the data.
#If the residuals have correlation, they are not conservative.
#There is no strict definition of correlation.  
#A possible threshold here is R^2>0.5, which would eliminate dO18, dD, Cl, NO3, SO42, Na,Mg, Al, K, Ca, Fe57, Zn, Ba from the tracer set

##############################
###Calculate Residual Root Mean Square Error (RRMSE)   #
##############################
RRMSE1D<-1
for (i in 1:NumTracers) { 
  RRMSE1D[i] <- (sqrt(sum(Resids1D[,i]^2)))/(length(Tracers[,i]*mean(Tracers[,i])))
}  #Where NumTracers is original number of tracers, but they are not all conservative
par(mfrow=c(1,1)) #sets plotting space for one graph
barplot(RRMSE1D,names.arg=titlevector, ylab = "Relative Root Square Error (RRMSE, %)", main="RRMSE for all tracers in 1D mixing space")

#The RRMSE should be minimized if possible.  
#When this is not possible, the RRMSE should not change from one dimension to the next (e.g. from 1-D to 2-D mixing space)

############################################################
###       Evaluate tracers in 2-Dimensional space        ###
############################################################
###Comments: the variable names were used to reflect those used in Hooper, 2003 equation(2). Equation is solved below in steps
############################################################
Xstar <- Tracers.std #just renames the Tracers matrix
V <- (Evectors)#this is actually VT b/c R outputs VT not V
V2DT <- as.matrix(V[,1:2])#There is no transpose operation because V is actually VT (V Transposed); this is where you change the number of dimensions
V2D <- t(V2DT) #Now you actually have V
V2DtV2D <- V2D%*%V2DT
InvV2DtV2D <- solve(V2DtV2D)
XstarV2DT <- Xstar%*%V2DT
Proj2Dstd <- XstarV2DT%*%InvV2DtV2D%*%V2D  #This solves the Hooper 2003 EQ and gives you the standardized projection
##############################
#Destandardize the projection#
##############################
Proj2D <- Proj2Dstd
for (i in 1:NumTracers) {
  Proj2D[,i]<-(Proj2Dstd[,i]*sd(Tracers[,i])+mean(Tracers[,i]))###Destandardize the projection
}
##############################
###Calculate residuals       #
##############################
Resids2D<-Proj2D - Tracers #Calculate residuals by subtracting original river values
write.table((cbind(Resids2D,Tracers)), "full2DResids.txt", sep="\t") #Output data table of residuals

##############################
#Plot the residuals          #
##############################
par(mfrow=c(3,3)) #change with # of variables
titlevector <- c("dD","NO3", "Al", "K")
for (i in 1:NumTracers) {       #changes the extent of loop based on the number of tracers
  lm2D<-(lm(Resids2D[,i]~Tracers[,i]))
  values <- (paste("R2 =",format(summary(lm2D)$adj.r.squared, digits=2),"P-value =",format(summary(lm2D)$coefficients[2,4], digits=2)))
  plot(Tracers[,i],Resids2D[,i],xlab="Measured concentration", ylab="Residuals", main=titlevector[i], sub=values)
  abline(lm2D$coefficients[1],lm2D$coefficients[2])
} 

#The residuals should be random if the 2-D space correctly represents the data.
#If the residuals have correlation, they are not conservative.
#There is no strict definition of correlation.  
#For 2D space the change in RRMSE from 1D space is used to define conservative tracers. 

##############################
###Calculate RRMSE           #
##############################
RRMSE2D<-1
for (i in 1:NumTracers) {   #changes the extent of loop based on the number of tracers
  RRMSE2D[i] <- (sqrt(sum(Resids2D[,i]^2)))/(length(Tracers[,i]*mean(Tracers[,i])))
}
par(mfrow=c(1,1)) #changes plotting window dimensions to 1 row x 1 column
barplot(RRMSE2D,names.arg=titlevector, ylab = "Relative Root Square Error (RRMSE, %)", main="RRMSE for all tracers in 2D mixing space")

#The RRMSE should be minimized if possible.  
#When this is not possible, the RRMSE should not change from one dimension to the next (e.g. from 1-D to 2-D mixing space)
#Comparing the RRMSE bar plots for 1D and 2D, Na, K, and Cl all change noticably, so they are excluded from the PCA.

######################################################################################################################################
######################################################################################################################################
######################################################################################################################################
###The following large block of code runs Princiapl Component Analysis on the stream dataset and performs EMMA based on the tracers selected above.         ###
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

############################################################
###           Run PCA for selected tracers               ###
############################################################
###Comments: Column #s 1-(tracer#1), 2-(tracer#2), 3-(tracer#3), 4-(tracer#4), 5-(tracer#5), etc. PCA run with tracers defined as conservative from diagnostics run above
############################################################
#In 2-D space, if the RRMSE changes noticably, they are excluded from the PCA.

#Run final PCA with 4 solutes: dD (1), NO3(2), Al(3), K(4)
PCA=prcomp(Tracers.std[,c(1,2,3,4)])#Must use standardized data for PCA WHERE Tracers.std is original stream water sample chemistry.  Columns selected are for the conservative tracers you want to use for EMMA
summary(PCA) #Check Proportion of Variance, xx% of variance explained after how many PCs?
#Here 99% of variance explained by first 2 Principal Components

par(mfrow=c(1,1))


write.table(PCA$sdev, "PCAsummaryfull.txt", sep="\t") #Outputs square root of eigenvalue. A SDEV of <1 means the eigenvalue is also <1#!!!!!!!!!!!!!!!!!!!!!!!!!SDEV!!!!!!!!!!!!!!!!$%^%@#%#$%#$%@#$@^@$%@^@$%^%^%^*&^%$%(*&^%$
write.table(PCA$rotation, "PCAloadingsfull.txt", sep="\t") #Outputs the weighting for each Principal Component (what percentage of the variance each one explains)
############################################################
###     Run EMMA using selected conservative tracers     ###
############################################################
###Comments: Column #s 1-(tracer#1), 2-(tracer#2), 3-(tracer#3), 4-(tracer#4), 5-(tracer#5), etc. PCA run with tracers defined as conservative from diagnostics
############################################################
#Run EMMA with only conservative tracers
#Project end members to U-space. Coordinates of projected end members are the Principal Component scores for each EM in the stream PC space.
EMScores= EMTracers.std[,c(1,2,3, 4)] %*% PCA$rotation #Make sure you subset the EM dataset to the same tracers as the stream dataset used above in the PCA

####################################
#Plot the mixing diagram in U-space#
####################################
#Here we use 2-D space (3 End Members)
#Set up the plotting space
par(mfrow=c(1,1))
#par(mar=c(6.1,4.1,4,2)) #can define the margins if you want to, but not necessary
#Plot stream values (scores from PCA).  
plot(PCA$x[1,1], PCA$x[1,2],xlab="PC1", ylab="PC2", cex.axis=1.8, cex.lab=1.4, lwd=1, main="Nooksack River Watershed End Member Mixing Analysis (EMMA)", xlim=range(-5:5), ylim=range(-5:5))  
points(PCA$x[1:119,1], PCA$x[1:119,2],xlab="PC1", pch=20,cex=2, col="cadetblue")  #rain (sep-mar)
points(PCA$x[120:225,1], PCA$x[120:225,2],xlab="PC1", pch=20, cex=2, col="coral4") #snowmelt (april-jun)
points(PCA$x[226:298,1], PCA$x[226:298,2],xlab="PC1", pch=20,cex=2, col="pink") #glacier melt (jul-)


#PCA$x[1:10,1], PCA$x[1:10,2], PCA$x[1:10,3] would be used for 10 stream samples and 3 dimensions (4 end members)
#with(PCA[1:22,], text(sr~dpi, labels = row.names(RawData[1:22,]), pos = 4))
#text(PCA$x[,1],PCA$x[,2], pos=3, cex=1, col="blue") #This labels the streamwater points but is not necessary
#Plot the end members in U-space on the same figure.  Also plot the 25th and 75th percentiles as dashed #lines for each end member.  

#Rain, End Member #1
points(EMScores[1,1], EMScores[1,2] , col="red", cex=3, pch=13)
#lines(c(EMScores[3,1], EMScores[3,1]),c(EMScores[2,2], EMScores[4,2]), col="red", lty=2) #median, median, 25th percentile, 75th percentile
#lines(c(EMScores[2,1], EMScores[4,1]),c(EMScores[3,2], EMScores[3,2]), col="red", lty=2 )#25th percentile, 75th percentile, median, median

#Rain2, End Member #2
points(EMScores[2,1], EMScores[2,2] , col="purple", cex=3, pch=13)
#lines(c(EMScores[8,1], EMScores[8,1]),c(EMScores[7,2], EMScores[9,2]), col="purple", lty=2)#median, median, 25th percentile, 75th percentile
#lines(c(EMScores[7,1], EMScores[9,1]),c(EMScores[8,2], EMScores[8,2]), col="purple", lty=2)#25th percentile, 75th percentil, median, median

#Groundwater, End Member#3
points(EMScores[3,1], EMScores[3,2] , col="brown", cex=3, pch=3)
#lines(c(EMScores[13,1], EMScores[13,1]),c(EMScores[12,2], EMScores[14,2]), col="magenta", lty=2)
#ines(c(EMScores[12,1], EMScores[14,1]),c(EMScores[13,2], EMScores[13,2]), col="magenta", lty=2)

#Snow, End Member#4
points(EMScores[4,1], EMScores[4,2] , col="chartreuse", cex=3, pch=8)
#lines(c(EMScores[18,1], EMScores[18,1]),c(EMScores[17,2], EMScores[19,2]), col="green", lty=2)
#lines(c(EMScores[17,1], EMScores[19,1]),c(EMScores[18,2], EMScores[18,2]), col="green", lty=2)

#Ice, End Member #5
points(EMScores[5,1], EMScores[5,2] , col="black", cex=3, pch=9)
#lines(c(EMScores[23,1], EMScores[23,1]View),c(EMScores[22,2], EMScores[24,2]), col="orange2", lty=2)
#lines(c(EMScores[22,1], EMScores[24,1]),c(EMScores[23,2], EMScores[23,2]), col="orange2", lty=2)

#Glacial Melt, End Member #6
points(EMScores[6,1], EMScores[6,2] , col="royalblue", cex=3, pch=10)

#TRIANGLE of chosen end members that best surround the stream data points (blue dots)
lines(c(EMScores[1,1],EMScores[5,1]),c(EMScores[1,2],EMScores[5,2]),lwd=2, lty=3)
lines(c(EMScores[5,1],EMScores[4,1]),c(EMScores[5,2],EMScores[4,2]),lwd=2, lty=3) 
lines(c(EMScores[4,1],EMScores[3,1]),c(EMScores[4,2],EMScores[3,2]),lwd=2, lty=3)
lines(c(EMScores[3,1],EMScores[2,1]),c(EMScores[3,2],EMScores[2,2]),lwd=2, lty=3)  
lines(c(EMScores[2,1],EMScores[6,1]),c(EMScores[2,2],EMScores[6,2]),lwd=2, lty=3)
lines(c(EMScores[6,1],EMScores[1,1]),c(EMScores[6,2],EMScores[1,2]),lwd=2, lty=3) 

#Make a legend
#legend("topleft",title="End Members", c("Rain1", "Rain2", "Groundwater", "Snow", "Ice"),
#col=c("red","purple","brown","grey","black"),pch=c(13,13, 3,8,9),lty=c(-1,-1,-1), pt.cex=1.5, pt.lwd=.8, cex=0.8, ncol=1, bg="white")
#legend("topright",title="Samples", c("NF","HP","MD","MF","RLP","SF","VDY"), 
#col=c("cadetblue","coral4", "pink","aquamarine", "darkblue",
#"blanchedalmond", "cornflowerblue"), pch=c(20, 20, 20,20, 20, 20, 20), lty=c(-1,-1,-1), pt.cex=1.5, pt.lwd=.8, cex=0.8, ncol=1, bg="white")
legend("topleft",title="End Members", c("Rain1", "Rain2", "Groundwater", "Snow", "Ice","Glacier_melt"),
       col=c("red","purple","brown","chartreuse","black", "royalblue"), pch=c(13, 13, 3,8,9, 10), lty=c(-1,-1,-1), pt.cex=1.5, pt.lwd=.8, cex=0.8, ncol=1, bg="white")
legend("topright",title="Time of Year", c("Sep-Mar","Apr-Jun","Jul-Aug"), 
       col=c("cadetblue","coral4", "pink"),
           pch=c(20, 20, 20), lty=c(-1,-1,-1), pt.cex=1.5, pt.lwd=.8, cex=0.8, ncol=1, bg="white")
#legend("bottomleft",c("Samples with E.coli"), col=c("black"),pch=c(24),lty=c(-1,-1,-1), pt.cex=1.5, pt.lwd=.8, cex=0.8, ncol=1, bg="white")
######################################################################################################################################

######################################################################################################################################
######################################################################################################################################
###The following large block of code calculates the relative contributions of endmembers. This is currently setup for 2D (3 End Member) mixing ###
###space and 3 end members, but can be expanded or reduced as needed by modifying the coefficient and solution matrices               ###
######################################################################################################################################
######################################################################################################################################
######################################################################################################################################

############################################################
###   Calculate relative contributions of endmembers     ###
############################################################
###Comments: This converts the system of equations in Liu 2008 eq(6-8) to matrix-vector form and solves for  the relative contributions (f1, f2, and f3). Tracers A1 and B1 are the two PC scores
############################################################
####################################
#Construct coefficient matrix. Row1=1, Row2=PC1 score, Row3=PC2 score. #Col1=median of snow,End Member 1 (????), Col2=median of Langtang R., EM2 (????), #Col3 = EM3 ***********CLARIFY WHICH EMS, Median of groundwaterEM3 (????)
####################################
CoefMatrix <- matrix(nrow=5,ncol=6)
CoefMatrix[1,c(1:5)] <- 1
# choose appropriate EM here to "bound" the stream samples***subjective
# Here row 1 = Rain1, row 2 = Rain2, row 3 = Groundwater, row 4= Snow, row 5 = Ice
CoefMatrix[2,c(1:4)] <- c(EMScores[1,1],EMScores[1,2],EMScores [1,3],EMScores[1,4])#median values of each end member, here 2 EMs are used, only 1 PC
CoefMatrix[3,c(1:4)] <- c(EMScores[2,1], EMScores[2,2], EMScores[2,3], EMScores[2,4])#add this line if there is a 2nd Principal Component (3 end members)
CoefMatrix[4,c(1:4)]<- c(EMScores [3,1], EMScores[3,2], EMScores[3,3], EMScores[3,4])#add this line if there is a 3rd Principal Component (4 end members)
CoefMatrix[5,c(1:4)]<- c(EMScores [4,1], EMScores[4,2], EMScores[4,3], EMScores[4,4])
CoefMatrix[6,c(1:4)]<- c(EMScores [5,1], EMScores[5,2], EMScores[5,3], EMScores[5,4])

#add this line if there is a 3rd Principal Component (4 end members)

####################################
#Construct Solution vector for each sample, 1st entry=1, 2nd entry=PC1 score of stream sample, 3rd entry=PC2 score of stream sample
####################################

Solvector <- matrix(nrow=length(PCA$x[,1]),ncol=6) #ncol = number of end members
for (i in 1:length(PCA$x[,1])) 
{Solvector[i,c(1:6)] <- cbind(1,PCA$x[i,1],PCA$x[i,2],PCA$x[i,3],PCA$x[i,4])} #i,c(i:n) where n is number of end members and where cbind(1,PCA$x[i,1], PCA#x[i,2] is used for 2 PCs, etc.)

#############################################################
#Calculate relative contribution in % for each stream sample#
#############################################################
RelContr <- matrix(nrow=length(PCA$x[,1]),ncol=5) #5x3 matrix where ncol=# of end members and five rows are number of samples
for (i in 1:length(PCA$x[,1])) {
  RelContr[i,c(1:5)] <- (solve(CoefMatrix,Solvector[i,c(1:5)]))*100  #for both instances of 'i,c(1:n)', n is number of end members
}
write.table(RelContr, "RelContributions.txt", sep="\t")

#Outputs file of relative contribution in % for each endmember, be sure to change name
#Here column 1=groundwater, column 2 = snow, column 3 = Upper Langtang River (inferred to be glacial melt)

#COmpare the Relative Contributions table to the elevation and date of sample collection from full_Langtang_river spreadsheet
#Note that the date of sample collection has an influence on the results and must considered in the interpretation

#Outputs file of relative contribution in % for each endmember, be sure to change name
#Here column 1=groundwater, column 2 = snow, column 3 = Upper Langtang River (inferred to be glacial melt)

#COmpare the Relative Contributions table to the elevation and date of sample collection from full_Langtang_river spreadsheet
#Note that the date of sample collection has an influence on the results and must considered in the interpretation

