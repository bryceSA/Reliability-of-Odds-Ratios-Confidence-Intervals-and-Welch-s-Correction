rm(list=ls())

#The OR_CI() function computes the odds ratio confidence interval.

#Input argument:
#n11: A numeric specifying the cell count for (1,1).
#n1d: A numeric specifying the total cell count for the first row.
#n21: A numeric specifying the cell count for (2,1).
#n2d: A numeric specifying the total cell count for the second row.
#dist: A character string specifying the distribution to be used. 
#"z" for standard normal, and "t" for t with Welch's adjustment
#adjust: A character string specifying the cell count adjustment. 
#"Agresti" for the independent-smoothed adjustment of Agresti.
#"Gart" for the Gart adjustment.
#"Woolf" for no adjustment, corresponding to the original Woolf logit interval
#alpha: A numeric related to confidence level. Note that the confidence level
#is given by 1-alpha. It is set to 0.05 by default.

#Change alpha to .1 for 90% CI
OR_CI <- function(n11, n1d, n21, n2d, dist=c("z","t"), 
                  adjust=c("Agresti","Gart","Woolf"), alpha=0.05)
{
  dist <- match.arg(dist)
  adjust <- match.arg(adjust)
  #The second column counts
  n12 <- n1d - n11
  n22 <- n2d - n21
  #Total column counts
  nd1 <- n11 + n21
  nd2 <- n12 + n22
  #Grand total count
  N <- n11 + n12 + n21 + n22
  
  #Cell count adjustments
  if(adjust == "Gart")
  {
    n11a <- n11 + 0.5
    n12a <- n12 + 0.5
    n21a <- n21 + 0.5
    n22a <- n22 + 0.5
  }
  else if(adjust == "Agresti")
  {
    c11 <- 2*n1d*nd1/(N^2)
    c12 <- 2*n1d*nd2/(N^2)
    c21 <- 2*n2d*nd1/(N^2)
    c22 <- 2*n2d*nd2/(N^2)
    
    n11a <- n11 + c11
    n12a <- n12 + c12
    n21a <- n21 + c21
    n22a <- n22 + c22
  }
  else #adjust == "Woolf"
  {
    n11a <- n11
    n12a <- n12
    n21a <- n21
    n22a <- n22
  }
  
  #theta is the odds ratio (OR) estimate
  theta  <- (n11a*n22a)/(n12a*n21a)
  
  #Standard error
  se <- sqrt(1/n11a + 1/n12a + 1/n21a + 1/n22a)
  
  #Critical value calculation
  if(dist=="z")
  {
    crit <- qnorm(alpha/2, lower.tail=FALSE)
  }
  if(dist=="t")
  {
    Na <- n11a + n12a + n21a + n22a
    p11a <- n11a/Na
    p12a <- n12a/Na
    p21a <- n21a/Na
    p22a <- n22a/Na
    
    f1 <- 1/p11a + 1/p12a + 1/p21a + 1/p22a
    f3 <- 1/(p11a^3) + 1/(p12a^3) + 1/(p21a^3) + 1/(p22a^3)
    
    #Welch's degrees of freedom
    nu <- (2*Na*f1^2)/(f3 - f1^2)
    
    crit <- qt(alpha/2, df=nu, lower.tail=FALSE)
  }
  
  #Lower and upper bound of the confidence interval on the log scale
  loglower <- log(theta) - crit*se
  logupper <- log(theta) + crit*se
  
  #Lower and upper bound of the confidence inerval on the original scale
  lower <- exp(loglower)
  upper <- exp(logupper)
  
  #Handling NAs by assigning 0 and infinity for the lower and upper bound.
  lowerNAs <- which(is.na(lower))
  upperNAs <- which(is.na(upper))
  
  if(length(lowerNAs) > 0) lower[lowerNAs] <- 0
  if(length(upperNAs) > 0) upper[upperNAs] <- Inf
  
  results <- rbind(lower, upper)
  rownames(results) <- c("lower", "upper")
  
  return(results)
}

#The simulations matrix below gives all the different simulation scenarios to consider.
#n1d and n2d gives the total number of trials for the first and second population.
#p1 is the probability of success for the first population.
#OR is the odds ratio.
#alpha is the value so that 100(1-alpha)% represents the confidence level.

simulations <- expand.grid(n1d=c(5,10,50,100,500,1000,2000), n2d=c(5,10,50,100,500,1000,2000), 
                           p1=c(.05, .1, .3, .5, .7, .9), OR=c(.001, .01, .3, .5, .7, .95, .99, .999), alpha=0.05)
simulations

#numsimcases gives the number of simulation cases.
numsimcases <- dim(simulations)[1]
numsimcases

#Number of Monte Carlo simulations 
#(10000 or more is recommended, but use 1000 for preliminary study)
MC <- 10000

#Initializing some vectors
empirical_CR_za <- c()
empirical_CR_ta <- c()
empirical_CR_zg <- c()
empirical_CR_tg <- c()
empirical_CR_zw <- c()
empirical_CR_tw <- c()
p2 <- c()
for(i in 1:numsimcases)
{
  print(paste("Running simulation ", i, sep=""))
  #Retrieving information about the i-th simulation scenario.
  curr_sim <- simulations[i,]
  curr_n1d <- curr_sim$n1d
  curr_n2d <- curr_sim$n2d
  curr_p1 <- curr_sim$p1
  curr_OR <- curr_sim$OR
  curr_p2 <- curr_p1/(curr_OR*(1-curr_p1)+curr_p1)
  p2[i] <- curr_p2
  curr_alpha <- curr_sim$alpha
  
  #Generating random n11's and n21's.
  n11s <- rbinom(n=MC, size=curr_n1d, prob=curr_p1)
  n21s <- rbinom(n=MC, size=curr_n2d, prob=curr_p2)
  
  #dist = "z" and adjust = "Agresti"
  case_za <- OR_CI(n11=n11s, n1d=curr_n1d, n21=n21s, n2d=curr_n2d, 
                   dist="z", adjust="Agresti", alpha=curr_alpha)
  
  #dist = "t" and adjust = "Agresti"
  case_ta <- OR_CI(n11=n11s, n1d=curr_n1d, n21=n21s, n2d=curr_n2d, 
                   dist="t", adjust="Agresti", alpha=curr_alpha)
  
  #dist = "z" and adjust = "Gart"
  case_zg <- OR_CI(n11=n11s, n1d=curr_n1d, n21=n21s, n2d=curr_n2d, 
                   dist="z", adjust="Gart", alpha=curr_alpha)
  
  #dist = "t" and adjust = "Gart"
  case_tg <- OR_CI(n11=n11s, n1d=curr_n1d, n21=n21s, n2d=curr_n2d, 
                   dist="t", adjust="Gart", alpha=curr_alpha)  
  
  #dist = "z" and adjust = "Woolf"
  case_zw <- OR_CI(n11=n11s, n1d=curr_n1d, n21=n21s, n2d=curr_n2d, 
                   dist="z", adjust="Woolf", alpha=curr_alpha)
  
  #dist = "t" and adjust = "Woolf"
  case_tw <- OR_CI(n11=n11s, n1d=curr_n1d, n21=n21s, n2d=curr_n2d, 
                   dist="t", adjust="Woolf", alpha=curr_alpha)    
  
  #Empirical coverage rate calculations
  empirical_CR_za[i] <- mean((case_za[1,] <= curr_OR)*(curr_OR <= case_za[2,])==1)
  empirical_CR_ta[i] <- mean((case_ta[1,] <= curr_OR)*(curr_OR <= case_ta[2,])==1)  
  empirical_CR_zg[i] <- mean((case_zg[1,] <= curr_OR)*(curr_OR <= case_zg[2,])==1)
  empirical_CR_tg[i] <- mean((case_tg[1,] <= curr_OR)*(curr_OR <= case_tg[2,])==1)   
  empirical_CR_zw[i] <- mean((case_zw[1,] <= curr_OR)*(curr_OR <= case_zw[2,])==1)
  empirical_CR_tw[i] <- mean((case_tw[1,] <= curr_OR)*(curr_OR <= case_tw[2,])==1)   
}

#Empirical coverage rates
empirical_CR_za
empirical_CR_ta
empirical_CR_zg
empirical_CR_tg
empirical_CR_zw
empirical_CR_tw

#All the results can be found here.
all_results <- cbind(simulations, 
                     p2,
                     empirical_CR_za,
                     empirical_CR_ta,
                     empirical_CR_zg,
                     empirical_CR_tg,
                     empirical_CR_zw,
                     empirical_CR_tw                 
)
all_results

#Save simulation results as a CSV file.
write.csv(all_results, file="all_results.csv", quote=FALSE, row.names=FALSE)

results <- read.csv("all_results.csv")
#install.packages("plotly")
library(plotly)


##########################
#OR and Sample Size vs CR#
#     GRAPHS             #
##########################

#: za_CR vs OR, n1d, and n2d #jordon
fig <- plot_ly(results, x = ~OR, y = ~n1d, z = ~n2d, color = ~empirical_CR_za, colors = c('red', 'green', 'blue'), 
               type = "scatter3d", mode = "markers") %>%
  layout(scene = list(
    xaxis = list(title = 'OR'),
    yaxis = list(title = 'n1d'),
    zaxis = list(title = 'n2d')
  ))

fig

#: ta_CR vs OR, n1d, and n2d #jordon
fig <- plot_ly(results, x = ~OR, y = ~n1d, z = ~n2d, color = ~empirical_CR_ta, colors = c('red', 'green', 'blue'), 
               type = "scatter3d", mode = "markers") %>%
  layout(scene = list(
    xaxis = list(title = 'OR'),
    yaxis = list(title = 'n1d'),
    zaxis = list(title = 'n2d')
  ))

fig

#: zg_CR vs OR, n1d, and n2d #lormel
fig <- plot_ly(results, x = ~OR, y = ~n1d, z = ~n2d, color = ~empirical_CR_zg, colors = c('red', 'green', 'blue'), 
               type = "scatter3d", mode = "markers") %>%
  layout(scene = list(
    xaxis = list(title = 'OR'),
    yaxis = list(title = 'n1d'),
    zaxis = list(title = 'n2d')
  ))

fig


#: tg_CR vs OR, n1d, and n2d #lormel
fig <- plot_ly(results, x = ~OR, y = ~n1d, z = ~n2d, color = ~empirical_CR_tg, colors = c('red', 'green', 'blue'), 
               type = "scatter3d", mode = "markers") %>%
  layout(scene = list(
    xaxis = list(title = 'OR'),
    yaxis = list(title = 'n1d'),
    zaxis = list(title = 'n2d')
  ))
fig

#:zw_CR vs OR, n1d, and n2d
fig <- plot_ly(results, x = ~OR, y = ~n1d, z = ~n2d, color = ~empirical_CR_zw, colors = c('red', 'green', 'blue'), 
               type = "scatter3d", mode = "markers") %>%
  layout(scene = list(
    xaxis = list(title = 'OR'),
    yaxis = list(title = 'n1d'),
    zaxis = list(title = 'n2d')
  ))
fig

#:tw_CR vs OR, n1d, and n2d
fig <- plot_ly(results, x = ~OR, y = ~n1d, z = ~n2d, color = ~empirical_CR_tw, colors = c('red', 'green', 'blue'), 
               type = "scatter3d", mode = "markers") %>%
  layout(scene = list(
    xaxis = list(title = 'OR'),
    yaxis = list(title = 'n1d'),
    zaxis = list(title = 'n2d')
  ))
fig

##########################
#   p1 and Sample Size   #
#    vs. CR Graphs       #
##########################

#: za_CR vs p1, n1d, and n2d #jordon
fig <- plot_ly(results, x = ~p1, y = ~n1d, z = ~n2d, color = ~empirical_CR_za, colors = c('red', 'green', 'blue'), 
               type = "scatter3d", mode = "markers") %>%
  layout(scene = list(
    xaxis = list(title = 'p1'),
    yaxis = list(title = 'n1d'),
    zaxis = list(title = 'n2d')
  ))

fig

#: ta_CR vs p1, n1d, and n2d #jordon
fig <- plot_ly(results, x = ~p1, y = ~n1d, z = ~n2d, color = ~empirical_CR_ta, colors = c('red', 'green', 'blue'), 
               type = "scatter3d", mode = "markers") %>%
  layout(scene = list(
    xaxis = list(title = 'p1'),
    yaxis = list(title = 'n1d'),
    zaxis = list(title = 'n2d')
  ))

fig

#: zg_CR vs p1, n1d, and n2d #lormel
fig <- plot_ly(results, x = ~p1, y = ~n1d, z = ~n2d, color = ~empirical_CR_zg, colors = c('red', 'green', 'blue'), 
               type = "scatter3d", mode = "markers") %>%
  layout(scene = list(
    xaxis = list(title = 'p1'),
    yaxis = list(title = 'n1d'),
    zaxis = list(title = 'n2d')
  ))

fig


#: tg_CR vs p1, n1d, and n2d #lormel
fig <- plot_ly(results, x = ~p1, y = ~n1d, z = ~n2d, color = ~empirical_CR_tg, colors = c('red', 'green', 'blue'), 
               type = "scatter3d", mode = "markers") %>%
  layout(scene = list(
    xaxis = list(title = 'p1'),
    yaxis = list(title = 'n1d'),
    zaxis = list(title = 'n2d')
  ))
fig

#:zw_CR vs p1, n1d, and n2d
fig <- plot_ly(results, x = ~p1, y = ~n1d, z = ~n2d, color = ~empirical_CR_zw, colors = c('red', 'green', 'blue'), 
               type = "scatter3d", mode = "markers") %>%
  layout(scene = list(
    xaxis = list(title = 'p1'),
    yaxis = list(title = 'n1d'),
    zaxis = list(title = 'n2d')
  ))
fig

#:tw_CR vs p1, n1d, and n2d
fig <- plot_ly(results, x = ~p1, y = ~n1d, z = ~n2d, color = ~empirical_CR_tw, colors = c('red', 'green', 'blue'), 
               type = "scatter3d", mode = "markers") %>%
  layout(scene = list(
    xaxis = list(title = 'p1'),
    yaxis = list(title = 'n1d'),
    zaxis = list(title = 'n2d')
  ))
fig

##########################
#   Histograms of CR     #
##########################
par(mfrow = c(2, 3))
hist(results$empirical_CR_ta, main = "Histogram of ta_CR", xlab = "empirical_CR_ta")
hist(results$empirical_CR_za, main = "Histogram of za_CR", xlab = "empirical_CR_za")
hist(results$empirical_CR_tg, main = "Histogram of tg_CR", xlab = "empirical_CR_tg")
hist(results$empirical_CR_zg, main = "Histogram of zg_CR", xlab = "empirical_CR_zg")
hist(results$empirical_CR_zw, main = "Histogram of zw_CR", xlab = "empirical_CR_zw")
hist(results$empirical_CR_tw, main = "Histogram of tw_CR", xlab = "empirical_CR_tw")

###################################################
#      2D Cross-Section (Line Plot)               #   
#    Fix n2d and p1, then plot coverage vs. OR    #
#    for multiple values of n1d.                  #
###################################################

#Fix n2d = 10, p1 = 0.3
df_slice <- subset(results, n2d == 1000 & p1 == 0.3)

#Plot coverage vs. OR grouping by n1d
fig_slice <- plot_ly(
  data = df_slice,
  x = ~OR,
  y = ~empirical_CR_zg,
  color = ~factor(n1d),  #color lines by n1d
  type = 'scatter',
  mode = 'lines+markers'
) %>%
  layout(
    title = "ZG Coverage vs. OR (n2d=1000, p1=0.3)",
    xaxis = list(title = 'OR'),
    yaxis = list(title = 'Empirical Coverage (ZG)'),
    legend = list(title = list(text = 'n1d'))
  )

fig_slice  


#Fix n2d = 10, p1 = 0.3
df_slice <- subset(results, n2d == 1000 & p1 == 0.3)

# Plot coverage vs. OR
fig_slice <- plot_ly(
  data = df_slice,
  x = ~OR,
  y = ~empirical_CR_tg,
  color = ~factor(n1d),  #color lines by n1d
  type = 'scatter',
  mode = 'lines+markers'
) %>%
  layout(
    title = "TG Coverage vs. OR (n2d=1000, p1=0.3)",
    xaxis = list(title = 'OR'),
    yaxis = list(title = 'Empirical Coverage (TG)'),
    legend = list(title = list(text = 'n1d'))
  )

fig_slice 


#Fix n2d = 10, p1 = 0.3
df_slice <- subset(results, n2d == 1000 & p1 == 0.3)

#Plot coverage vs. OR 
fig_slice <- plot_ly(
  data = df_slice,
  x = ~OR,
  y = ~empirical_CR_zg,
  color = ~factor(n1d),  
  type = 'scatter',
  mode = 'lines+markers'
) %>%
  layout(
    title = "ZA Coverage vs. OR (n2d=1000, p1=0.3)",
    xaxis = list(title = 'OR'),
    yaxis = list(title = 'Empirical Coverage (ZA)'),
    legend = list(title = list(text = 'n1d'))
  )

fig_slice  


#Fix n2d = 10, p1 = 0.3
df_slice <- subset(results, n2d == 1000 & p1 == 0.3)

# Plot coverage vs. OR
fig_slice <- plot_ly(
  data = df_slice,
  x = ~OR,
  y = ~empirical_CR_ta,
  color = ~factor(n1d),  
  type = 'scatter',
  mode = 'lines+markers'
) %>%
  layout(
    title = "TA Coverage vs. OR (n2d=1000, p1=0.3)",
    xaxis = list(title = 'OR'),
    yaxis = list(title = 'Empirical Coverage (TA)'),
    legend = list(title = list(text = 'n1d'))
  )

fig_slice 


#Fix n2d = 10, p1 = 0.3
df_slice <- subset(results, n2d == 1000 & p1 == 0.3)

#Plot coverage vs. OR
fig_slice <- plot_ly(
  data = df_slice,
  x = ~OR,
  y = ~empirical_CR_zw,
  color = ~factor(n1d),  
  type = 'scatter',
  mode = 'lines+markers'
) %>%
  layout(
    title = "ZW Coverage vs. OR (n2d=1000, p1=0.3)",
    xaxis = list(title = 'OR'),
    yaxis = list(title = 'Empirical Coverage (ZW)'),
    legend = list(title = list(text = 'n1d'))
  )

fig_slice  


#Fix n2d = 10, p1 = 0.3
df_slice <- subset(results, n2d == 1000 & p1 == 0.3)

# Plot coverage vs. OR
fig_slice <- plot_ly(
  data = df_slice,
  x = ~OR,
  y = ~empirical_CR_tw,
  color = ~factor(n1d), 
  type = 'scatter',
  mode = 'lines+markers'
) %>%
  layout(
    title = "TW Coverage vs. OR (n2d=1000, p1=0.3)",
    xaxis = list(title = 'OR'),
    yaxis = list(title = 'Empirical Coverage (TW)'),
    legend = list(title = list(text = 'n1d'))
  )

fig_slice 

###################
#   ZG ANALYSIS   #
###################
df_slice <- subset(results, n2d == 5 & n1d == 5)

# Plot coverage vs. p1
fig_slice <- plot_ly(
  data = df_slice,
  x = ~p1,
  y = ~empirical_CR_tg,
  color = ~factor(OR),  
  type = 'scatter',
  mode = 'lines+markers'
) %>%
  layout(
    title = "ZG Coverage vs. p1 (n2d=5, n1d=5)",
    xaxis = list(title = 'p1'),
    yaxis = list(title = 'Empirical Coverage (ZG)'),
    legend = list(title = list(text = 'OR'))
  )

fig_slice 

# Fix n2d = 10, p1 = 0.3
df_slice <- subset(results, n2d == 1000 & n1d == 5)

# Plot coverage vs. p1
fig_slice <- plot_ly(
  data = df_slice,
  x = ~p1,
  y = ~empirical_CR_zg,
  color = ~factor(OR),  
  type = 'scatter',
  mode = 'lines+markers'
) %>%
  layout(
    title = "ZG Coverage vs. p1 (n2d=1000, n1d=5)",
    xaxis = list(title = 'p1'),
    yaxis = list(title = 'Empirical Coverage (ZG)'),
    legend = list(title = list(text = 'OR'))
  )

fig_slice 


# Fix n2d = 10, p1 = 0.3
df_slice <- subset(results, n2d == 1000 & n1d == 1000)

# Plot coverage vs. p1
fig_slice <- plot_ly(
  data = df_slice,
  x = ~p1,
  y = ~empirical_CR_zg,
  color = ~factor(OR),  
  type = 'scatter',
  mode = 'lines+markers'
) %>%
  layout(
    title = "ZG Coverage vs. p1 (n2d=1000, n1d=1000)",
    xaxis = list(title = 'p1'),
    yaxis = list(title = 'Empirical Coverage (ZG)'),
    legend = list(title = list(text = 'OR'))
  )

fig_slice

df_slice <- subset(results, n2d == 5 & n1d == 5)

# Plot coverage vs. p1
fig_slice <- plot_ly(
  data = df_slice,
  x = ~p1,
  y = ~empirical_CR_tg,
  color = ~factor(OR),
  type = 'scatter',
  mode = 'lines+markers'
) %>%
  layout(
    title = "TG Coverage vs. p1 (n2d=5, n1d=5)",
    xaxis = list(title = 'p1'),
    yaxis = list(title = 'Empirical Coverage (TG)'),
    legend = list(title = list(text = 'OR'))
  )

fig_slice 

# Fix n2d = 10, p1 = 0.3
df_slice <- subset(results, n2d == 1000 & n1d == 5)

# Plot coverage vs. p1
fig_slice <- plot_ly(
  data = df_slice,
  x = ~p1,
  y = ~empirical_CR_tg,
  color = ~factor(OR),  
  type = 'scatter',
  mode = 'lines+markers'
) %>%
  layout(
    title = "TG Coverage vs. p1 (n2d=1000, n1d=5)",
    xaxis = list(title = 'p1'),
    yaxis = list(title = 'Empirical Coverage (TG)'),
    legend = list(title = list(text = 'OR'))
  )

fig_slice 


# Fix n2d = 10, p1 = 0.3
df_slice <- subset(results, n2d == 1000 & n1d == 1000)

# Plot coverage vs. p1, grouping by OR
fig_slice <- plot_ly(
  data = df_slice,
  x = ~p1,
  y = ~empirical_CR_tg,
  color = ~factor(OR),  
  type = 'scatter',
  mode = 'lines+markers'
) %>%
  layout(
    title = "TG Coverage vs. p1 (n2d=1000, n1d=1000)",
    xaxis = list(title = 'p1'),
    yaxis = list(title = 'Empirical Coverage (TG)'),
    legend = list(title = list(text = 'OR'))
  )

fig_slice

##############################
# WOOLF TEST PERFORMS POORLY #
#   UNLESS O.R IS approx 1   #
##############################
# Fix n2d = 10, p1 = 0.3, show ZW coverage
df_slice <- subset(results, n2d == 5 & n1d == 5)

# Plot coverage vs. p1, grouping by OR
fig_slice <- plot_ly(
  data = df_slice,
  x = ~OR,
  y = ~empirical_CR_zw,
  color = ~factor(p1),  # color lines by n1d
  type = 'scatter',
  mode = 'lines+markers'
) %>%
  layout(
    title = "ZW Coverage vs. OR (n2d=5,  n1d=5)",
    xaxis = list(title = 'OR'),
    yaxis = list(title = 'Empirical Coverage (ZW)'),
    legend = list(title = list(text = 'p1'))
  )

fig_slice  

######################
# Real Data Analysis #
######################

#install.packages("vcd")
library(vcd)
data("Arthritis")
data
#If a patient’s Improved status is "Some" or "Marked" group into "Better".
Arthritis$Better <- ifelse(Arthritis$Improved %in% c("Some","Marked"),
                           "Better","None")

Arthritis$Better <- factor(Arthritis$Better, levels=c("None","Better"))

contab <- with(Arthritis, table(Treatment, Better))
contab



OR_CI <- function(n11, n1d, n21, n2d, dist=c("z","t"), 
                  adjust=c("Agresti","Gart","Woolf"), alpha=0.05)
{
  dist <- match.arg(dist)
  adjust <- match.arg(adjust)
  #The second column counts
  n12 <- n1d - n11
  n22 <- n2d - n21
  #Total column counts
  nd1 <- n11 + n21
  nd2 <- n12 + n22
  #Grand total count
  N <- n11 + n12 + n21 + n22
  
  #Cell count adjustments
  if(adjust == "Gart")
  {
    n11a <- n11 + 0.5
    n12a <- n12 + 0.5
    n21a <- n21 + 0.5
    n22a <- n22 + 0.5
  }
  else if(adjust == "Agresti")
  {
    c11 <- 2*n1d*nd1/(N^2)
    c12 <- 2*n1d*nd2/(N^2)
    c21 <- 2*n2d*nd1/(N^2)
    c22 <- 2*n2d*nd2/(N^2)
    
    n11a <- n11 + c11
    n12a <- n12 + c12
    n21a <- n21 + c21
    n22a <- n22 + c22
  }
  else #adjust == "Woolf"
  {
    n11a <- n11
    n12a <- n12
    n21a <- n21
    n22a <- n22
  }
  
  #theta is the odds ratio (OR) estimate
  theta  <- (n11a*n22a)/(n12a*n21a)
  
  #Standard error
  se <- sqrt(1/n11a + 1/n12a + 1/n21a + 1/n22a)
  
  #Critical value calculation
  if(dist=="z")
  {
    crit <- qnorm(alpha/2, lower.tail=FALSE)
  }
  if(dist=="t")
  {
    Na <- n11a + n12a + n21a + n22a
    p11a <- n11a/Na
    p12a <- n12a/Na
    p21a <- n21a/Na
    p22a <- n22a/Na
    
    f1 <- 1/p11a + 1/p12a + 1/p21a + 1/p22a
    f3 <- 1/(p11a^3) + 1/(p12a^3) + 1/(p21a^3) + 1/(p22a^3)
    
    #Welch's degrees of freedom
    nu <- (2*Na*f1^2)/(f3 - f1^2)
    
    crit <- qt(alpha/2, df=nu, lower.tail=FALSE)
  }
  
  #Lower and upper bound of the confidence interval on the log scale
  loglower <- log(theta) - crit*se
  logupper <- log(theta) + crit*se
  
  #Lower and upper bound of the confidence inerval on the original scale
  lower <- exp(loglower)
  upper <- exp(logupper)
  
  #Handling NAs by assigning 0 and infinity for the lower and upper bound.
  lowerNAs <- which(is.na(lower))
  upperNAs <- which(is.na(upper))
  
  if(length(lowerNAs) > 0) lower[lowerNAs] <- 0
  if(length(upperNAs) > 0) upper[upperNAs] <- Inf
  
  results <- rbind(lower, upper)
  rownames(results) <- c("lower", "upper")
  
  return(results)
}

n11 <- 29   #Placebo, None
n12 <- 14   #Placebo, Better
n21 <- 13   #Treated, None
n22 <- 28   #Treated, Better
n1d <- n11 + n12  #total placebo = 43
n2d <- n21 + n22  #total treated = 41

#C.I's
ci_zA <- OR_CI(n11, n1d, n21, n2d, dist="z", adjust="Agresti", alpha=0.05)
ci_tA <- OR_CI(n11, n1d, n21, n2d, dist="t", adjust="Agresti", alpha=0.05)
ci_zG <- OR_CI(n11, n1d, n21, n2d, dist="z", adjust="Gart", alpha=0.05)
ci_tG <- OR_CI(n11, n1d, n21, n2d, dist="t", adjust="Gart", alpha=0.05)
ci_zW <- OR_CI(n11, n1d, n21, n2d, dist="z", adjust="Woolf", alpha=0.05)
ci_tW <- OR_CI(n11, n1d, n21, n2d, dist="t", adjust="Woolf", alpha=0.05)

results <- data.frame(
  method = c("z_Agresti","t_Agresti","z_Gart","t_Gart","z_Woolf","t_Woolf"),
  lower = c(ci_zA["lower",1], ci_tA["lower",1],
            ci_zG["lower",1], ci_tG["lower",1],
            ci_zW["lower",1], ci_tW["lower",1]),
  upper = c(ci_zA["upper",1], ci_tA["upper",1],
            ci_zG["upper",1], ci_tG["upper",1],
            ci_zW["upper",1], ci_tW["upper",1])
)
results
