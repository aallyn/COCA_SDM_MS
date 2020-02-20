#####
## SDM prediction validation plots
#####

## Going to try to replicate the Taylor Diagram, which combines information on the root mean square error (standard deviation of the residuals -- how far off our predictions are on the measurement scale), coefficient of determination (correlation between predicted and observed value, as well as deviation off of the 1:1 line as you could have a correlation of 1 and be systematically biased away from 1:1 line), and bias (ratio of standard deviation in the model to the data to assess if the model is under or over predicting the amount of variability observed in the real world, making things too smooth or too dynamic).

library(plotrix) # Includes a taylor.diagram function, with arguments for reference and model values. Going to modify this a bit I think...
library(tidyverse)

# Model results
mod.stats<- read_csv("~/Github/COCA/Results/NormalVoting_BiomassIncPresNoExposure_03152019/mod.results.csv")

taylor_diagram_func<- function(mydata, obs = "obs", mod = "mod", group = NULL, type = "default", 
 normalise = FALSE, cols = "brewer1", rms.col = "darkgoldenrod", cor.col = "black", arrow.lwd = 3, annotate = "centred\nRMS error", key = TRUE, key.title = group, key.columns = 1, key.pos = "right", strip = TRUE, auto.text = TRUE, ...) {
  
  ## Details
  # This function plots a Taylor Diagram of model prediction accuracy, sumamrizing the root mean square error, the coefficient of determination, and the ration of standard deviations. 
  
  # Args:
  # mod.stats = Path to the model results csv folder, which needs at least the three values (RMSE, Coeff Determination, SD ratio)
  # group = Grouping variable -- defaults to species
  # out.path = Path to save the Taylor Diagram
  
  # Returns: NULL; saves plot to output directory
  
  ## Start function
  # Install libraries
  library_check(c("tidyverse"))
  
  # Set arguments for debugging -- this will NOT run when you call the function. Though, you can run each line inside the {} and then you will have everything you need to walk through the rest of the function.
  if(FALSE){
    dat<- "~/Github/COCA/Results/NormalVoting_BiomassIncPresNoExposure_03152019/mod.results.csv"
    group<- "COMNAME"
    out.path<- "~/Desktop/"
  }

  if(FALSE){
    mod.res
    grad.corr.lines = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1)
    pcex = 1
    cex.axis = 1
    normalize = TRUE
    mar = c(5, 4, 6, 6)  
  }
  
  # Real data
  mod.res<- read_csv(mod.stats) %>%
    drop_na(RMSE.SDM.B, CorrCoeff.SDM.B, Bias.SDM.B, RMSE.NEVA.B, CorrCoeff.NEVA.B, Bias.NEVA.B) %>%
    filter(., AUC.SDM > 0.7)
  
  # Getting maxSD for plotting
  maxsd<- 1.5 * max(mod.res$Bias.SDM.B, 1)
  
  # Empty plot first
  # Creating empty plot first
  plot.base<- ggplot() + 
    scale_x_continuous(name = "Standard deviation (normalized)", limits = c(0, maxsd), breaks = seq(from = 0, to = maxsd, by = 0.5), expand = c(0, 0, 0.05, 0)) +
    scale_y_continuous(name = "Standard deviation (normalized)", limits = c(0, maxsd), breaks = seq(from = 0, to = maxsd, by = 0.5), expand = c(0, 0, 0.05, 0)) +
    theme_classic()
  
  # Coeff D rays 
  for(i in 1:length(grad.corr.lines)){
    x.vec<- c(0, maxsd*grad.corr.lines[i])
    y.vec<- c(0, maxsd*sqrt(1 - grad.corr.lines[i]^2))
    
    if(i ==1){
      coeffd.rays.df<- data.frame("Ray" = rep(1, length(x.vec)), "x" = x.vec, "y" = y.vec)
    } else {
      temp<- data.frame("Ray" = rep(i, length(x.vec)), "x" = x.vec, "y" = y.vec)
      coeffd.rays.df<- bind_rows(coeffd.rays.df, temp)
    }
  }
  
  # Add rays
  plot.coeffd<- plot.base +
    geom_line(data = coeffd.rays.df, aes(x = x, y = y, group = Ray), lty = "longdash")
  
  coeffd.labs<- coeffd.rays.df %>%
    group_by(Ray) %>%
    summarize(., 
              "x" = max(x, na.rm = TRUE), 
              "y" = max(y, na.rm = TRUE)) %>%
    data.frame()
  
  coeffd.labs$Label<- grad.corr.lines
  
  plot.coeffd<- plot.coeffd +
    geom_label(data = coeffd.labs, aes(x = x, y = y, label = Label), fill = "white", label.size = NA)
  
  # SD arcs
  # Need to add in SD arcs
  sd.arcs<- seq(from = 0, to = maxsd, by = 1)
  
  for(i in 1:length(sd.arcs)){
    x.vec<- sd.arcs[i]*cos(seq(0, pi/2, by = 0.03))
    y.vec<- sd.arcs[i]*sin(seq(0, pi/2, by = 0.03))
    
    if(i ==1){
      sd.arcs.df<- data.frame("Arc" = rep(1, length(x.vec)), "x" = x.vec, "y" = y.vec)
    } else {
      temp<- data.frame("Arc" = rep(i, length(x.vec)), "x" = x.vec, "y" = y.vec)
      sd.arcs.df<- bind_rows(sd.arcs.df, temp)
    }
  }
  
  # Add arcs to plot.base
  plot.sd<- plot.coeffd +
    geom_line(data = sd.arcs.df, aes(x = x, y = y, group = Arc), lty = "dotted") 
  
  # Now gamma? -- Standard deviation arcs around the reference point
  gamma<- pretty(c(0, maxsd), n = 4)[-1]
  gamma<- gamma[-length(gamma)]
  labelpos<- seq(45, 70, length.out = length(gamma))
  
  for(gindex in 1:length(gamma)) {
    xcurve <- cos(seq(0, pi, by = 0.03)) * gamma[gindex] + sd.r
    endcurve <- which(xcurve < 0)
    endcurve <- ifelse(length(endcurve), min(endcurve) - 1, 105)
    ycurve <- sin(seq(0, pi, by = 0.03)) * gamma[gindex]
    maxcurve <- xcurve * xcurve + ycurve * ycurve
    startcurve <- which(maxcurve > maxsd * maxsd)
    startcurve <- ifelse(length(startcurve), max(startcurve) + 1, 0)
    x.vec<- xcurve[startcurve:endcurve]
    y.vec<- ycurve[startcurve:endcurve]
    
    if(gindex ==1){
      gamma.df<- data.frame("Gamma" = rep(1, length(x.vec)), "x" = x.vec, "y" = y.vec)
    } else {
      temp<- data.frame("Gamma" = rep(gindex, length(x.vec)), "x" = x.vec, "y" = y.vec)
      gamma.df<- bind_rows(gamma.df, temp)
    }
  }
  
  gamma.df$Gamma<- factor(gamma.df$Gamma, levels = unique(gamma.df$Gamma))
  
  # Add em
  plot.gamma<- plot.sd +
    geom_line(data = gamma.df, aes(x = x, y = y, group = Gamma), lty = "solid", col = "lightgray")
  
  # Label...
  gamma.labs<- gamma.df %>%
    group_by(Gamma) %>%
    summarize("x" = mean(x, na.rm = TRUE), 
              "y" = median(y, na.rm = TRUE))
  
  plot.gamma<- plot.gamma +
    geom_label(data = gamma.labs, aes(x = x, y = y, label = Gamma), fill = "white", label.size = NA)
  
  # Add in reference point
  plot.all<- plot.gamma +
    geom_point(aes(x = sd.r, y = 0), color = "red")
  
  # Add in species points
  mod.results.td<- mod.results %>%
    filter(., CorrCoeff.SDM.B >= 0) %>%
    mutate(., "TD.X" = Bias.SDM.B * CorrCoeff.SDM.B,
           "TD.Y" = Bias.SDM.B * sin(acos(CorrCoeff.SDM.B)))
  
  plot.all.f<- plot.all +
    geom_point(data = subset(mod.results.td, SEASON == "FALL"), aes(x = TD.X, y = TD.Y, color = COMNAME), show.legend = FALSE)
  
  plot.all.s<- plot.all +
    geom_point(data = subset(mod.results.td, SEASON == "SPRING"), aes(x = TD.X, y = TD.Y, color = COMNAME), show.legend = FALSE)
}
