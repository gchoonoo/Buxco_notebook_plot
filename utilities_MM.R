library(plethy)
library(R.utils)
library(plyr)
library(RColorBrewer)
library(reshape2)
library(IRanges)
library(ggplot2)

# Function for annotating virus column
virus.query <- function(obj) {
  query = "SELECT Break_Chunk_ID, substr(substr(Sample_Name, instr(Sample_Name, ' ') + 1), instr(substr(Sample_Name, instr(Sample_Name, ' ') + 1), ' ') + 1) as Virus FROM Sample JOIN Chunk_Time USING (Sample_ID)" 
  return(query)
}

#outer.cols=c(Flu="black", SARS="brown", Mock="blue")
mvtsplot.data.frame <- function(use.dta, plot.value="Penh",main=plot.value, outer.group.name=NULL, inner.group.name=NULL, outer.cols=NULL, colorbrewer.pal="PRGn")
{
    
    if (is.data.frame(use.dta) == FALSE || ncol(use.dta) < 3 || all(c('Sample_Name', 'Days') %in% colnames(use.dta)) == FALSE)
    {
        stop("ERROR: use.dta needs to be a data.frame with at least 3 columns, two of which need to be named Sample_Name and Days")
    }
    
    if (length(plot.value) != 1 || is.character(plot.value) == FALSE || plot.value %in% colnames(use.dta) == FALSE || is.numeric(use.dta[,plot.value]) == FALSE)
    {
        stop("ERROR: plot.value needs to be a character vector of length 1 corresponding to a numeric column in use.dta")
    }
    
    if (missing(outer.group.name) || is.null(outer.group.name))
    {
        use.dta$temp_outer <- factor("temp")
        outer.group.name <- "temp_outer"
    }
    else if (length(outer.group.name) != 1 || is.character(outer.group.name) == FALSE || outer.group.name %in% colnames(use.dta) == FALSE || (is.character(use.dta[,outer.group.name]) == FALSE && is.factor(use.dta[,outer.group.name]) == FALSE))
    {
        stop("ERROR: outer.group.name needs to a character vector of length one corresponding to a character or factor column in use.dta")
    }
    else if (is.character(use.dta[,outer.group.name]))
    {
        use.dta[,outer.group.name] <- factor(use.dta[,outer.group.name])
    }
    
    if (missing(inner.group.name) || is.null(inner.group.name))
    {
        use.dta$temp_inner <- factor("temp")
        inner.group.name <- "Sample_Name"
    }
    if (length(inner.group.name) != 1 || is.character(inner.group.name) == FALSE || inner.group.name %in% colnames(use.dta) == FALSE || (is.character(use.dta[,inner.group.name]) == FALSE && is.factor(use.dta[,inner.group.name]) == FALSE))
    {
        stop("ERROR: inner.group.name needs to a character vector of length one corresponding to a character or factor column in use.dta")
    }
    else if (is.character(use.dta[,inner.group.name]))
    {
        use.dta[,inner.group.name] <- factor(use.dta[,inner.group.name])
    }
    
    if (outer.group.name == "temp_outer")
    {
        outer.cols <- c(temp="black")
    }
    else if ((length(outer.cols) != nlevels(use.dta[,outer.group.name]) || is.null(names(outer.cols)) || all(names(outer.cols) %in% levels(use.dta[,outer.group.name])) == FALSE || all(outer.cols %in% colors() == FALSE)))
    {
        stop("ERROR: outer.cols needs to be a named character vector corresponding to levels in the outer.group.name column containing the names of colors")
    }
    
    if (is.character(colorbrewer.pal) == FALSE || length(colorbrewer.pal) != 1 || colorbrewer.pal %in% rownames(brewer.pal.info) == FALSE)
    {
        stop("ERROR: colorbrewer.pal needs to be a single valid RColorBrewer palette.")
    }
    
    #pal.list <- rep(colorbrewer.pal, nlevels(use.dta[,outer.group.name]))
    #names(pal.list) <- levels(use.dta[,outer.group.name])
    
    mean.mat <- acast(formula=Days~Sample_Name, data=use.dta, value.var=plot.value)
	#print(colnames(mean.mat))
    
    sample.cross <- use.dta[,c("Sample_Name", inner.group.name, outer.group.name)]
    sample.cross <- sample.cross[!duplicated(sample.cross),]
    rownames(sample.cross) <- as.character(sample.cross$Sample_Name)
    
	#print(rownames(sample.cross))
    sample.cross <- sample.cross[colnames(mean.mat),]
    
	#print(sample.cross)
    clear.ord <- do.call("order", list(sample.cross[,outer.group.name], sample.cross[,inner.group.name]))
	#print(clear.ord)
    
    mean.mat <- mean.mat[,clear.ord]
    
    sample.cross <- sample.cross[clear.ord,]
    
    outer.group <- sample.cross[,outer.group.name]
    inner.group <- sample.cross[,inner.group.name]
    
    cx <- t(scale(t(mean.mat)))
    
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    par(las = 1, cex.axis = 0.6)
    cn <- colnames(cx)
    nc <- ncol(cx)
    
    utime <- 1:nrow(cx)
    
    use.layout.list <- lapply(seq(from=1, length.out=nlevels(outer.group), by=2), function(x)
           {
                return (c(rep(x, 6), x+1))
           })
    
    use.layout <- do.call("rbind", use.layout.list)
    
    if (nrow(use.layout) == 1)
    {
        use.layout <- rbind(use.layout, use.layout, use.layout)
    }
    
    use.layout <- rbind(use.layout, c(rep((nlevels(outer.group)*2)+1, 6), (nlevels(outer.group)*2)+2))
  
    layout(use.layout)
    
    colm.list <- lapply(levels(outer.group), function(x)
                        {
                            apply(mean.mat[,outer.group == x], 2, boxplot, plot=F)
                        })
    
    names(colm.list) <- levels(outer.group)
    
    #right.xlim <- range(do.call("cbind", lapply(unlist(colm.list, recursive=F), "[", "stats")), na.rm = TRUE)
	right.xlim = c(-2.5, 3.5)
   
    #side2 <- max(0.55, max(strwidth(cn, "inches")) + 0.1)
    side2 <- max(strwidth(cn, "inches")) + .8
    
    use.heat.colors <- brewer.pal(4, colorbrewer.pal)
    
    for (i in levels(outer.group))
    {
        #c(bottom, left, top, right)
        
        bottom <- .05
        
        if (i == levels(outer.group)[1])
        {
            top <- .1
        }else
        {
            top <- 0 
        }
        
        par(mai = c(bottom, side2, top, 0.05))
        cur.cx <- cx[,outer.group == i]
        image(utime, seq_len(sum(outer.group == i)), cur.cx, col = use.heat.colors, xlim = range(utime),
            xaxt = "n", yaxt = "n", ylab = "", xlab = "")
        #axis(2, at = seq_len(nc), cn)
        
        usrpar <- par("usr")
        par(usr = c(usrpar[1:2], 0, 1))
        
        cur.inner <- inner.group[outer.group == i, drop=T]
        
        group.bounds <- lapply(levels(cur.inner), function(x) end(IRanges(Rle(cur.inner == x))))
        names(group.bounds) <- unique(cur.inner)
        group.lines <- sort(as.numeric(unlist(group.bounds)))
        last.line <- group.lines[length(group.lines)]
        group.lines <- group.lines[-length(group.lines)]/length(cur.inner)
        
        if (length(group.lines) > 0)
        {
            abline(h = group.lines, lwd = 1, col = 1)
        
            for (j in names(group.bounds)){
                for (k in group.bounds[[j]]){
                    dist.val <- k/(length(cur.inner))
                    half.dist <- dist.val - min(abs(dist.val - group.lines[group.lines != dist.val]))/2
                    mtext(text=j, side=2, at=half.dist, las=1, cex=.5, line=.5, col=outer.cols[i])
                }
                
            }
            
            mtext(text=i, side=2, las=0, line=.5+(ceiling(max(nchar(cn))/2)), col=outer.cols[i])
        }
        
        colm <- colm.list[[i]]
        
        par(mai = c(bottom, .05, top, .1))
        
        if (i == levels(outer.group)[nlevels(outer.group)])
        {
            plot(1:length(colm), type = "n", ylab = "", yaxt = "n",
                xlab = "", xlim = right.xlim)
        }else
        {
            plot(1:length(colm), type = "n", ylab = "", yaxt = "n", xaxt = "n",
                xlab = "", xlim = right.xlim)
        }
       
        usrpar <- par("usr")
        par(usr = c(usrpar[1:2], 0, 1))
        
        ypos <- (1:length(colm) - 1/2)/length(colm)
        
        for (cur_el in 1:length(colm))
        {
            bxp(colm[[cur_el]], horizontal=T, at=ypos[cur_el],add=T, boxwex=.1, xaxt="n", yaxt="n", ylim=c(-2.5, 3.5))
            #temporary...
            #mtext(text=names(colm)[cur_el], side=2, at=ypos[cur_el], cex=.5)
        }
        
        abline(h = group.lines, lwd = 1, col = 1)
    }
    
    all.meds <- sapply(levels(outer.group), function(x)
                       {
                            apply(mean.mat[,outer.group == x], 1, median, na.rm = TRUE)
                       })
    
    #bottom.ylim <- range(all.meds, na.rm = TRUE)
    
    bottom.ylim <- c(-3,3)
    
    par(mai = c(0.6, side2, 0.05, 0.05))
    plot(utime, all.meds[,1], type = "l", ylim = bottom.ylim, xaxt = "n", 
        xlab = "", ylab = "Median", col=outer.cols[colnames(all.meds)[1]])
    
    if (ncol(all.meds) > 1)
    {
        for (i in 2:ncol(all.meds))
        {
            lines(all.meds[,i], type="l", col=outer.cols[colnames(all.meds)[i]])
        }

    }
    
    par(usr = c(c(0,nrow(all.meds)) , par("usr")[3:4]))
    Axis(at=1:nrow(all.meds), labels=rownames(all.meds), side = 1)
    mtext("Days", side=1, outer=F, padj=3, cex=.7)
    
    #as the above margins are too large, mainly as a placeholder...
    par(mai = c(0.05, 0.05, 0.1, 0.1))
    plot(0, 0, xlab = "", ylab = "", axes = FALSE, type = "n")
    #maybe make a legend here
    #text(0, 0, main)
    legend("center", legend=rev(c("Low", "", "", "High")), fill=rev(use.heat.colors), title=main,y.intersp = .75)
    
}



