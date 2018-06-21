#' -----------------------------------------------------------------------------
#' -----------------------------------------------------------------------------
#' Analysis of experimental data collected on 16 organopesticide (OP) in 
#' rats and humans to understand differences in acetylcholine esterase (AchE)
#' inhibition between rats and humans to ultimately estimate uncertainty factors
#' contingent upon intra and interspecies differences.  Uncertainty factors are 
#' calculated based on U.S. EPA DDEF (data-derived extrapolation factors) 
#' guidance 
#' ----------------------- Ankur Singhal ---------------------------------------
#' --------------------- Managing Scientist ------------------------------------
#' ----------------------- Exponent, Inc. --------------------------------------
#' ----------------------- June 20, 2018 ---------------------------------------

library(readxl)
library(data.table)
library(dplyr)
library(car)
library(ggplot2)


setwd("/Users/asinghal/Documents/Projects/OP AchE analysis/")


ReadExcelDt <- function(..., col_names = FALSE, na = "NA"){
  as.data.table(read_excel(..., col_names = col_names, na = na))
}


# Comparison between species for each chemical 
Interspecies <- function(dat) {
  #' If both the rat and human distributions are normally distributed, 
  #' then do the t test 
  #' If both the rat and human distributions are lognormally distributed, 
  #' then do the t test on the log-transformed values
  #' If neither of the above, then do the wilcoxon rank sum test 
  pv <- by(dat,    
     dat[, Chemical], # Comparison within each chemical
     function (x) {
       ifelse (all(x[, Distr] == "Normal"), 
              # Condition: Both rat and human normally distributed?
              t.test(ki ~ Species, x)$p.value, 
              # Extract p value of the t test
              ifelse (all(x[, Distr] == "Lognormal"), 
                     # Condition: Both rat and human lognormally distributed
                     t.test(log(ki) ~ Species, x)$p.value,
                     # Extract p value of the t test (log transformed ki)
                     wilcox.test(ki ~ Species, x)$p.value
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
Intraspecies <- function(dat) {
  #' Create a new column Transformedki, whick keeps normally distributed and 
  #' "Neither" ki values as-is, log-transforms lognormally distributed values
  #' Do an ANOVA on normally and lognormally distributed values; otherwise
  #' Kruskal-Wallis (KW) test
  dat[, Transformedki := ki]; dat[Distr == "Lognormal", dat := log(ki)]
  #' This is for KW test, which only takes a group
  dat[, Group := paste(Sex, AgeGroup, Ethnicity)]
  
  #' ANOVA 
  t1 <- by(dat[Distr != "Neither",], 
           dat[Distr != "Neither", Chemical],
           function(x) {
             summary(aov(Transformedki ~ AgeGroup + Sex + Ethnicity, x))
           })
  
  #' KW  
  t2 <- by(dat[Distr == "Neither",], 
                    dat[Distr == "Neither", Chemical], 
                    function(x) {
                      kruskal.test(Transformedki ~ factor(Group), x)
                    })

  return(c(t1, t2))
}


# Load and format all inhibition data ------------------------------------------
options(stringsAsFactors = FALSE)

wb <- "preliminary inhibition data.xlsx"
ws <- excel_sheets(wb)
# Select all sheets except for Summary sheet
ws <- ws[!grepl("Summary", ws)]


# Read and convert to data table
dt1 <- sapply(ws, function(s) {
  ReadExcelDt(wb, sheet=s)
}, simplify = FALSE, USE.NAMES = TRUE)

#' Assign column names (they are all in row 2) 
#' Remove the first row 
#' Select relevant cols -- "Sample", "Species", "Sex", "Age", "Ethnicity", "ki"
dt2 <-  lapply(dt1, function(x) { 
  names(x) <- as.character(x[2, ])
  x[, AgeGroup := "Child"]
  x[Age == 0, AgeGroup := "Fetal"]
  x[Age >= 18, AgeGroup := "Adult"]
  x <- x[!is.na(Sample), ]
  x[2:nrow(x), .(Sample, Species, Sex, Age, AgeGroup, Ethnicity, ki)]
  })


# Combine all the different tables into one with Chemical as the first column
dt3 <- dplyr::bind_rows(dt2, .id = 'Chemical')

# Convert Age and ki into numeric columns
dt3[, c("Age", "ki") := lapply(.SD, as.numeric), .SDcols = c("Age", "ki")]


# Load and format replicate data -----------------------------------------------
#' Replicate data only collected in humans for samples 3, 14, 16 and 17 and only
#' for chemicals Omethoate, Phosmet Oxon and Naled
options(stringsAsFactors = FALSE)

rm(wb, ws)
wb <- "replicate data.xlsx"
ws <- excel_sheets(wb)

# Read and convert to data table
rep <- sapply(ws, function(s) {
  ReadExcelDt(wb, sheet=s)
}, simplify = FALSE, USE.NAMES = TRUE)

#' Assign column names (they are all in row 2) 
#' Remove the first two rows 
#' Select relevant cols -- "Sample", "Species", "Sex", "Age", "Ethnicity", "ki"
rep2 <- lapply(rep, function(x) { 
  names(x) <- as.character(x[1, ])
  x[, AgeGroup := "Child"]
  x[Age == 0, AgeGroup := "Fetal"]
  x[Age >= 18, AgeGroup := "Adult"]
  x <- x[!is.na(Sample), ]
  x[2:nrow(x), .(Sample, Species, Sex, Age, AgeGroup, Ethnicity, ki)]
  })

# Combine all the different tables into one with Chemical as the first column
rep3 <- dplyr::bind_rows(rep2, .id = 'Chemical')

# Convert Age and ki into numeric columns
rep3[, c("Age", "ki") := lapply(.SD, as.numeric), .SDcols = c("Age", "ki")]

# Basic stats on the replicate data --------------------------------------------
#' Only three replicates per sample, so cannot do any intrasample comparisons, 
#' though it looks some values are outliers
rep.stats <- rep3[, .(mean(ki, na.rm = TRUE), sd(ki, na.rm = TRUE)), 
                  by = c("Chemical", "Sample")]
names(rep.stats) <- c("Chemical", "Sample", "Mean", "SD")


# Replace all inhibition data with replicate data, where applicable ------------
all.dt <- rep.stats[dt3, on=c("Chemical", "Sample")]
all.dt[!is.na(Mean), ki := Mean]

# Statistical comparisons of the new data set with replicate data in it---------
all.dt[, .N, by = c("Chemical", "Species")]

## Normality Q-Q plots
par(mfrow = c(4,4))
all.dt[, qqPlot(ki), by = c("Chemical", "Species")]
all.dt[, plot(density(ki, na.rm = TRUE)), by = c("Chemical", "Species")]

# Shapiro test for normality
sw <- all.dt[, shapiro.test(ki), by = c("Chemical", "Species")]
#' If p-value >= 0.05, normally distributed
#' Assign to new column Distr; all others will get NA (useful later)
all.dt[sw[p.value >= 0.05, ], Distr := "Normal", on = c("Chemical", "Species")]

#' If not normally distributed (NA values of Distr), check for lognormality
sw.log <- all.dt[is.na(Distr), shapiro.test(log(ki)), by = c("Chemical", "Species")]
#' Assign to "Distr"
all.dt[sw.log[p.value >= 0.05,], Distr := "Lognormal", on = c("Chemical", "Species")]
#' If neither normal or lognormal, assign "Neither" to Distr
all.dt[is.na(Distr), Distr := "Neither"] 

# Check 
all.dt[Distr == "Neither", qqPlot(ki), by = c("Chemical", "Species")]


# Inter- and intraspecies comparison -------------------------------------------
# Data table of p values
inter.p.vals <- Interspecies(all.dt)

# Data frame of ANOVA and KW tests 
intra.comp <- Intraspecies(all.dt[Species == "Human", ])


# Stats and UF calculations ----------------------------------------------------
#' Calculate mean, sd, 50th and 90th percentile for the human data
#' and merge with the InterPvalues to create a new data table 
stats <- inter.p.vals[all.dt[Species == "Human", 
                     .(mean(ki, na.rm = TRUE), 
                       sd(ki, na.rm = TRUE), 
                       quantile(ki, .50, na.rm = TRUE),
                       quantile(ki, .95, na.rm = TRUE)),
                     by = "Chemical"], 
                 on = "Chemical"]

#' Calculate mean and sd and merge
stats <- stats[all.dt[Species == "Rat", 
                     .(mean(ki, na.rm = TRUE), 
                       sd(ki, na.rm = TRUE)),
                     by = "Chemical"],
                 on = "Chemical"]

names(stats) <- c("Chemical","P_Interspecies", "Mean_Human", "SD_Human", 
                   "P50_Human","P95_Human", "Mean_Rat", "SD_Rat")

# Calculate Coefficient of Variation (CV) for humans
stats[, CV_Human := SD_Human/Mean_Human * 100]
# Calculate interspecies factor as described in EPA DDEF guidance
stats[, InterFactor := 1]
stats[P_Interspecies < 0.05, InterFactor := round(Mean_Human / Mean_Rat, 1)]
# Calculate intraspecies factor as described in EPA DDEF guidance
stats[, IntraFactor := round(P95_Human/P50_Human, 1)]
# Calculate total uncertainty factor (UF) 
stats[, UF := ceiling(InterFactor * IntraFactor)]


# Evaluation of Ki Variability -------------------------------------------------
# In an effort to tease apart measurement variability from true variability
ggplot(stats, aes(x=log(Mean_Human), y=log(CV_Human)))+
  geom_point(aes(color=Chemical), shape=21, size=5)+
  stat_smooth(method="lm", se=F, color="black", size=1) +
  labs(x="Log of Ki Mean", y="Log of Ki Coefficient of Variation") + 
  theme_minimal() +
  theme(axis.title=element_text(size=24), 
        axis.text=element_text(size=18),
        legend.text=element_text(size=14),
        legend.title=element_blank())


# Output to files --------------------------------------------------------------
# Text ouput of intraspecies comparison, which contains info from models
sink("Intraspecies comparison.txt") # Attach file
intra.comp # Send output to text file
sink() # back to console

# .csv output of uncertainty factors, etc. 
write.csv(stats, "Stats and uncertainty factors.csv", row.names = FALSE)




