########################################################
### Function for plotting
########################################################

# This fucntion was extracted from "10_Analysis clade niche.R" for Acaena sister clades

library(ggplot2)

plotAnalysis <- function(data, 
                         m, # linear model object
                         xv, yv, # column names of responce and explanatory variable
                         xlabname, ylabname, # axes names for plot
                         nodeNumber,
                         showStats = T # TRUE; Show p value and slope of linear model and colour points, FALSE; No stat values and black points
){
  
  if(showStats == T){
    myplot <- ggplot(data, aes_string(x = xv, y = yv, label = nodeNumber)) +
      geom_point() +
      # text label for points
      geom_text(size=5) +
      # change xy labels
      labs(x = xlabname, y = ylabname) +
      # change text size
      theme(text = element_text(size = 20),
            axis.text.x = element_text(size = 20)) +
      # drow LM line & confident intervals 
      stat_smooth(method = "lm", col = "red") +
      # show stats result as title
      labs(title = paste("Adj R2 =", signif(summary(m)$adj.r.squared, digits = 2),
                         "Intercept =", signif(m$coef[[1]], 2),
                         " Slope =", signif(m$coef[[2]], 2),
                         " P =", signif(summary(m)$coef[2, 4], 2))) +
      theme(panel.background = element_rect(fill = "gray95"))
  } else {
    myplot <- ggplot(data, aes_string(x = xv, y = yv, label = nodeNumber)) +
      geom_point() +
      # change xy labels
      labs(x = xlabname, y = ylabname) +
      # change text size
      theme(text = element_text(size = 20),
            axis.text.x = element_text(size = 20)) +
      # drow LM line & confident intervals 
      stat_smooth(method = "lm", col = "red") +
      theme(panel.background = element_rect(fill = "gray95"), legend.position="none")
  }
  return(myplot)
}
