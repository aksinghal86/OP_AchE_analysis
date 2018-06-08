library(readxl)
library(data.table)
library(dplyr)
library(car)


setwd("/Users/asinghal/Documents/Projects/OP AchE analysis/")


read_excel_dt <- function(..., col_names=FALSE, na="NA"){
  as.data.table(read_excel(..., col_names=col_names, na=na))
}

# Comparison between species for each chemical 
interspecies <- function(dat) {
  #' If both the rat and human distributions are normally distributed, 
  #' then do the t test 
  #' If both the rat and human distributions are lognormally distributed, 
  #' then do the t test on the log-transformed values
  #' If neither of the above, then do the wilcoxon rank sum test 
  pv <- by(dat,    
     dat[,Chemical], # Comparison within each chemical
     function (x) {
       ifelse(all(x[,Distr]=="Normal"), 
              # Condition: Both rat and human normally distributed?
              t.test(ki~Species, x)$p.value, 
              # Extract p value of the t test
              ifelse(all(x[,Distr]=="Lognormal"), 
                     # Condition: Both rat and human lognormally distributed
                     t.test(log(ki)~Species, x)$p.value,
                     # Extract p value of the t test (log transformed ki)
                     wilcox.test(ki~Species, x)$p.value
                     # Wilcoxon rank sum test; extract p values
              )
       )
     })
  
  #' Trick to convert the output from above into a data table 
  #' for easy merging later
  pv <- suppressWarnings(melt(data.table(t(sapply(pv, I)))))
  names(pv) <- c("Chemical", "P_Interspecies")
  
  return(pv)
}



# Comparison by Age, Sex and Ethnicity per chemical
intraspecies <- function(dat) {
  #' Create a new column Transformedki, whick keeps normally distributed and 
  #' "Neither" ki values as-is, log-transforms lognormally distributed values
  #' Do an ANOVA on normally and lognormally distributed values; otherwise
  #' Kruskal-Wallis (KW) test
  dat[, Transformedki:=ki];dat[Distr=="Lognormal", dat:=log(ki)]
  #' This is for KW test, which only takes a group
  dat[, Group:=paste(Sex, AgeGroup, Ethnicity)]
  
  #' ANOVA 
  t1 <- by(dat[Distr!="Neither",], 
           dat[Distr!="Neither",Chemical],
           function(x) summary(aov(Transformedki~AgeGroup+Sex+Ethnicity, x)))
  
  #' KW  
  t2 <- by(dat[Distr=="Neither",], 
                    dat[Distr=="Neither",Chemical], 
                    function(x) kruskal.test(Transformedki~factor(Group), x))

  return(c(t1, t2))
}

#### Load and format data ####
options(stringsAsFactors = FALSE)



wb <- "preliminary inhibition data.xlsx"
ws <- excel_sheets(wb)
# Select all sheets except for Summary sheet
ws <- ws[!grepl("Summary", ws)]

# Read and convert to data table
dt1 <- sapply(ws, function(s) {
  read_excel_dt(wb, sheet=s)
}, simplify=FALSE, USE.NAMES=TRUE)

#' Assign column names (they are all in row 2) 
#' Remove the first two rows 
#' Select relevant cols, i.e., "Sample", "Species", "Sex", "Age", "Ethnicity", "ki"
dt2 <-  lapply(dt1, function(x) { 
  names(x) <- as.character(x[2,])
  x[,AgeGroup:="Child"]; x[Age==0, AgeGroup:="Fetal"]; x[Age>=18, AgeGroup:="Adult"]
  x <- x[!is.na(Sample),]
  x[3:nrow(x), .(Sample, Species, Sex, Age, AgeGroup, Ethnicity, ki)]})


# Combine all the different tables into one with the Chemical as the first column
dt3 <- dplyr::bind_rows(dt2, .id='Chemical')
dt3 <- dt3[!is.na(Sample),]

# Convert Age and ki into numeric columns
dt3[,c("Age", "ki"):=lapply(.SD, as.numeric), .SDcols=c("Age", "ki")]

#### Statistical Comparisons ####
dt3[, .N, by=c("Chemical", "Species")]
par(mfrow=c(4,4))

## Normality Q-Q plots
dt3[, qqPlot(ki), by=c("Chemical", "Species")]
dt3[, plot(density(ki, na.rm=TRUE)), by=c("Chemical", "Species")]

## Shapiro test for normality
# Shapiro-Wilk test to check for normality
SW <- dt3[, shapiro.test(ki), by=c("Chemical", "Species")]
#' If p-value >= 0.05, normally distributed
#' Assign to new column Distr; all others will get NA (useful later)
dt3[SW[p.value>=0.05,], Distr:="Normal", on=c("Chemical", "Species")]

#' If not normally distributed (NA values of Distr), check for lognormal
SWlog <- dt3[is.na(Distr), shapiro.test(log(ki)), by=c("Chemical", "Species")]
#' Assign to "Distr"
dt3[SWlog[p.value>=0.05,], Distr:="Lognormal", on=c("Chemical", "Species")]
#' If neither normal or lognormal, assign "Neither" to Distr
dt3[is.na(Distr), Distr:="Neither"] 

# Check 
dt3[Distr=="Neither", qqPlot(ki), by=c("Chemical", "Species")]



#### Inter- and intraspecies comparison ####
# Data table of p values
InterPvalues <- interspecies(dt3)

# Data frame of ANOVA and KW tests 
IntraComp <- intraspecies(dt3[Species=="Human",])

#### Stats and UF calculations ####

#' Calculate mean, sd, 50th and 90th percentile for the human data
#' and merge with the InterPvalues to create a new data table 
Stats <- InterPvalues[dt3[Species=="Human", 
                     .(mean(ki, na.rm=TRUE), 
                       sd(ki, na.rm=TRUE), 
                       quantile(ki, .50, na.rm=TRUE),
                       quantile(ki, .95, na.rm=TRUE)),
                     by="Chemical"], 
                 on="Chemical"]

#' Calculate mean and sd and merge
Stats <- Stats[dt3[Species=="Rat", 
                     .(mean(ki, na.rm=TRUE), 
                       sd(ki, na.rm=TRUE)),
                     by="Chemical"],
                 on="Chemical"]

names(Stats) <- c("Chemical","P_Interspecies", "Mean_Human", "SD_Human", 
                   "P50_Human","P95_Human", "Mean_Rat", "SD_Rat")

# Calculate Coefficient of Variation (CV) for humans
Stats[,CV_Human:=SD_Human/Mean_Human*100]
# Calculate interspecies factor as described in EPA DDEF guidance
Stats[, InterFactor:=1]
Stats[P_Interspecies<0.05, InterFactor:=ceiling(Mean_Human/Mean_Rat)]
# Calculate intraspecies factor as described in EPA DDEF guidance
Stats[, IntraFactor:=ceiling(P95_Human/P50_Human)]
# Calculate uncertainty factor (UF) by multiplying inter- and intraspecies factors
Stats[, UF:=InterFactor*IntraFactor]


#### Output to files #####
# Text ouput of intraspecies comparison, which contains info from models
sink("Intraspecies comparison.txt") # Attach file
IntraComp # Send output to text file
sink() # back to console

# .csv output of uncertainty factors, etc. 
write.csv(Stats, "Stats and uncertainty factors.csv", row.names = FALSE)

##### FINISH ##### 



