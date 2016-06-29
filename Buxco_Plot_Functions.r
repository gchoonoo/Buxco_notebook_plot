library(colorRamps)
library(RColorBrewer)
library(lattice)

# Error bar function
error.bar = function(x, y, upper, lower=upper, arrow.length=0.05, ...) {
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper)) {
    #print('x')
    #print(length(x))
    #print('y')
    #print(length(y))
    #print('upper')
    #print(length(upper))
    stop("vectors must be same length")
  }
  arrows(x, y+upper, x, y-lower, angle=90, code=3, length=arrow.length, ...)
}

# Standard error function
std_error <- function(x) sd(x)/sqrt(length(x))

# Time series plot function
buxco_plot = function(data1=full_data, batch_date=unique(full_data[which(full_data[,"Mating"] == unique(full_data[,"Mating"])[1]),"Date"])[1], mating=as.character(unique(full_data[,"Mating"]))[1], var_data=vars[1], virus=virus){
  
  #print(head(data1))
  
  # Subset date, mating, and variable
  if(sum(data1[,"Mating"] %in% mating & data1[,"Date"] %in% batch_date & data1[,"Variable_Name"] == var_data & data1[,"Virus"] %in% virus) >0 & sum(is.na(data1[data1[,"Mating"] %in% mating & data1[,"Date"] %in% batch_date & data1[,"Variable_Name"] == var_data & data1[,"Virus"] %in% virus, "mean_per_day"])) < length(data1[data1[,"Mating"] %in% mating & data1[,"Date"] %in% batch_date & data1[,"Variable_Name"] == var_data & data1[,"Virus"] %in% virus, "mean_per_day"])){
    
    data1[data1[,"Mating"] %in% mating,] -> data1
    
    data1[data1[,"Date"] %in% batch_date & data1[,"Variable_Name"] == var_data & data1[,"Virus"] %in% virus,] -> data1
    
    #print(head(data1))
    
    data1$line_virus_date = paste(data1$Mating,data1$Virus,data1$Date,sep="_")
    
    #print(head(data1))
    
    # subset columns
    data1[, c("line_virus_date","line_virus_mean_per_day","Days_PI","mean_per_day")] -> data1_v2
    
    #print(head(data1_v2))
    
    table(data1_v2[,"Days_PI"],data1_v2[,"line_virus_date"]) -> completeness
    
    colnames(completeness)[which(sapply(1:ncol(completeness),function(x)sum(completeness[,x] ==0) > 0))] -> id_missing_data1
    
    # Fill in NAs
    if(length(id_missing_data1) >0){
      for(i in 1:length(id_missing_data1)){
        
        missing_days = unique(data1_v2$Days_PI)[which(!as.numeric(as.character(unique(data1_v2$Days_PI))) %in% data1_v2[which(data1_v2[,"line_virus_date"] == id_missing_data1[i]),"Days_PI"])]
        
        dat = vector("list",length(missing_days))
        
        for(j in 1:length(missing_days)){
          
          dat[[j]]=c(id_missing_data1[i], NA, missing_days[j], NA)
          
          data1_v2 <- rbind(data1_v2, dat[[j]])
          
        }
        
      }
    }
    
    data1_v2[order(data1_v2[,"line_virus_date"], as.numeric(as.character(data1_v2[,"Days_PI"]))),] -> data1_v3
    
    unique(data1_v3[,c("line_virus_date","line_virus_mean_per_day","Days_PI")]) -> data1_v3
    
    matrix(ncol=length(unique(data1_v3[,"Days_PI"])),nrow=length(unique(data1_v3[,"line_virus_date"])), data=as.numeric(as.character(data1_v3[,"line_virus_mean_per_day"])), byrow=T) -> data1_v4
    
    row.names(data1_v4) <- unique(data1_v3[,"line_virus_date"])
    colnames(data1_v4) <- unique(data1_v3[,"Days_PI"])
    
    ###############################
    # aggregate standard error
    ###############################
    aggregate(data1_v2$mean_per_day ~ data1_v2$line_virus_date + data1_v2$Days_PI, data=data1_v2, std_error) -> data1.error
    
    names(data1.error) <- c("line_virus_date","Days_PI","Std_error")    
    
    table(data1.error[,"Days_PI"],data1.error[,"line_virus_date"]) -> completeness2
    
    colnames(completeness2)[which(sapply(1:ncol(completeness2),function(x)sum(completeness2[,x] ==0) >0))] -> id_missing_data12
    
    if(length(id_missing_data12) >0){
      for(i in 1:length(id_missing_data12)){
        missing_days_e = unique(data1.error$Days_PI)[which(!as.numeric(as.character(unique(data1.error$Days_PI))) %in% data1.error[which(data1.error[,"line_virus_date"] == id_missing_data12[i]),"Days_PI"])]
        dat.e = vector("list",length(missing_days_e))
        for(j in 1:length(missing_days_e)){
          dat.e[[j]]=c(id_missing_data12[i], missing_days_e[j], NA)
          data1.error <- rbind(data1.error, dat.e[[j]])
        }
      }
    }
    
    data1.error[order(data1.error[,"line_virus_date"], as.numeric(as.character(data1.error[,"Days_PI"]))),] -> data1.error_v2
    
    matrix(ncol=length(unique(data1.error_v2[,"Days_PI"])),nrow=length(unique(data1.error_v2[,"line_virus_date"])), data=as.numeric(as.character(data1.error_v2[,"Std_error"])), byrow=T) -> data1.error_v3
    
    row.names(data1.error_v3) <- unique(data1.error_v2[,"line_virus_date"])
    colnames(data1.error_v3) <- unique(data1.error_v2[,"Days_PI"])
    
    ################################
    
    # Plots
    xrange <- c(1,length(as.numeric(as.character(colnames(data1_v4)))))
    
    # number of lines to draw
    ntrees <- nrow(data1_v4)
    
    # ncols for legend
    n_cols = floor(ntrees/4) + ceiling((ntrees%%4)/4)
    
    yrange <- range(data1_v4+data1.error_v3,data1_v4-data1.error_v3,data1_v4, na.rm=T)
    
    # set up the plot 
    plot(xrange, yrange, type="n", xlab="Day Post Infection", ylab=var_data, main=var_data, xaxt='n')
    axis(side=1,1:length(as.numeric(as.character(colnames(data1_v4)))),labels=colnames(data1_v4))
    colors1 <- rainbow(ntrees) 
    linetype1 <- 1:ntrees
    
    # add lines 
    for (i in 1:ntrees) { 
      #print(i)
      lines(as.numeric(as.character(data1_v4[i,])), type="l", lwd=4, col=colors1[i], pch="")
      error.bar(1:xrange[2], as.numeric(as.character(data1_v4[i,])), as.numeric(as.character(data1.error_v3[i,])), lwd=4, col=colors1[i])
    } 
    
    # legend
    legend('topright', row.names(data1_v4), cex=0.8, col=colors1, title="Mating + Virus + Batch Date",lwd=4)
    
  }else{
    return(NULL)
  }
}

# Dotplot function

dot_plot_data = function(var, virus, lines, xlab=NULL, day=NULL, day_summary, vert_line=0) {
  ## Get the appropriate dataframe
  
  # Subsetting data frame read in "Buxco_Plots.r" (GC)
  # df = eval(parse(text=var))
  df=buxco_annot[buxco_annot[,"Variable_Name"] %in% var,]
  
  ## Summarize across all days or on a single day
  if (is.null(day) & day_summary %in% c(2, 3)) {
    df_sub = df[df$Virus %in% c(virus, 'Mock'),]
  } else {
    df_sub = df[df$Virus %in% c(virus, 'Mock') & df$Days_PI==day,]
  }
  
  if (nrow(df_sub) == 0) {return(NULL)}
  
  ## Subset selected lines
  df_sub = df_sub[df_sub$Mating %in% lines, ]
  if (nrow(df_sub) == 0) {return(NULL)}
  
  ## Aggregate to create a single value for each animal
  if (day_summary == 1 | day_summary == 2) {
    buxco_col = 'mean_diff_infected_mock_per_day'
  } else {
    buxco_col = 'AUC_diff_infected_mock'
  }
  pheno = aggregate(df_sub[, buxco_col], list(df_sub$Mating, df_sub$ID, df_sub$Date, df_sub$Virus), mean, na.rm=T)
  colnames(pheno) = c('Mating', 'ID', 'Date', 'Virus', 'value')
  
  ## Set label; a combination of animal Mating and Date
  # Added as.character(pheno$Date) (GC)
  pheno$label = paste0(substring(pheno$Date, 1, 3), "-", sapply(strsplit(as.character(pheno$Date), '_'), '[', 2), "__", pheno$Mating)
  
  ## Create dotplot
  if (dim(pheno)[1] > 0) {
    plot_obj = dotplot(reorder(pheno[,'label'], pheno[,'value'], mean, na.rm=T) ~
                         pheno[,'value'] | pheno[,'Virus'],
                       panel = function(x,y,...) {panel.dotplot(x,y,...); panel.abline(v=vert_line, col.line="red")},
                       pch=19, ylab='Date + Mating', xlab=xlab)
  } else {
    return(NULL)
  }
  
  ## Return formatted data
  return(list(p_obj=plot_obj, var=var, virus=virus, day=day, xlab=xlab))
}


dot_plot = function(dp_data) {
  ## Get heritability
  # Updated this to all_buxco_herit$variable instead of tolower(variable) (GC)
  herit = all_buxco_herit[with(all_buxco_herit, virus==dp_data$virus & all_buxco_herit$variable==dp_data$var), 'icc_gelman_hill']
  
  ## Create plot title
  if (dp_data$virus == 'FLU') {
    virus_label = 'Influenza'
  } else {
    virus_label = 'SARS'
  }
  
  if (!is.null(dp_data$day)) {
    dp_data$p_obj$main = paste(virus_label, ' D', dp_data$day, ' ', dp_data$xlab, ';  ICC = ', herit, sep='')
  } else {
    dp_data$p_obj$main = paste(virus_label, ' ', dp_data$xlab, ';  ICC = ', herit, sep='')
  }
  plot(dp_data$p_obj)
  
  return(1)
}