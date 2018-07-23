########################################################
### Function for plotting
########################################################

plotAnalysis <- function(data,
                         genus_name, # Name of genus as plot title
                         xv, yv, # column names of responce and explanatory variable
                         xlabname, ylabname, # axes names for plot
                         nodeNumbercol = NULL,
                         showStats = FALSE, # TRUE; Show p value and slope of linear model and colour points, FALSE; No stat values and black points
                         label.point = FALSE, # TRUE; Show labels above points
                         cex = 10 # title and axis label size
                         ){
  
  
  myplot <- ggplot(data, aes_string(x = xv, y = yv)) +
    geom_point() +
    ggtitle(genus_name) +
    # change xy labels
    labs(x = xlabname, y = ylabname) +
    # change text size
    theme(text = element_text(size = cex),
          axis.text.x = element_text(size = cex)) +
    # drow LM line & confident intervals 
    stat_smooth(method = "lm", col = "red") +
    theme(panel.background = element_rect(fill = "gray95"), legend.position="none") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  
  if(showStats == TRUE){

    # linear model object
    m <- lm(formula(paste(yv, "~" ,xv)), data)
    
    myplot <- myplot +
      # show stats result as title
      labs(title = paste(genus_name, "Adj R2 =", signif(summary(m)$adj.r.squared, digits = 2),
                         "Intercept =", signif(m$coef[[1]], 2),
                         " Slope =", signif(m$coef[[2]], 2),
                         " P =", signif(summary(m)$coef[2, 4], 2)))
  }

  if(label.point == TRUE){
    myplot <- myplot +
      # text label for points
      geom_text(aes_string(label = nodeNumbercol),
                vjust = -0.5,
                size = 5)
  }
  
  return(myplot)
}


