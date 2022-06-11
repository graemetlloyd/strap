#' Plots a phylogeny against the geological time scale
#'
#' @description
#'
#' Plots a time-scaled phylogeny against the international geological time scale.
#'
#' @param tree A tree as a phylo object.
#' @param ages A dataset containing the first and last appearence datums,"FAD" and "LAD" respectively, of all taxa in the phylogeny. See the object $ages in utils::data(Dipnoi) for an example.
#' @param direction The direction the tree is to be plotted in, options include "rightwards" and "upwards", see help(plot.phylo).
#' @param units The temporal unit(s) to be included in the timescale, options include: "Eon", "Era", "Period", "Epoch", "Age" and "User". The option "User" is required when including a user-defined timescale. This also requires an object to be assigned to user.scale (see Details).
#' @param boxes Option for including grey boxes at a certain temporal resolution, options are the same as for units.
#' @param tick.scale The resolution of the tick marks at the base of the timescale, the default is the same as units. The resolution of the scale can also be chosen by specifiying a value or removed entirely by using "no".
#' @param user.scale The data object to be used when including a user-defined time scale, requires the option "User" to be included in units. See utils::data(UKzones) as an example of the required data format.
#' @param cex.age Size of the text on the scale bar.
#' @param cex.ts Size of the text on the geological time scale.
#' @param cex.tip Size of the tip labels on the phylogeny
#' @param width Width of the edges of the phylogeny.
#' @param label.offset A value for the distance between the nodes and tip labels, see help(plot.phylo).
#' @param ts.col Option for using standard ICS colours on the time scale.
#' @param vers The version of the geological time scale to use. Options include: "ICS2013", "ICS2012", "ICS2010", "ICS2009" or "ICS2008".
#' @param x.lim A two item statement for the geological range, in millions of years, of the plot i.e. (0,65). If only one value is used it will be used as the upper limit, see help(plot.phylo).
#' @param quat.rm Option to remove the names from Quaternary time bins, useful when plotting clades with long durations that range through to the recent.
#' @param erotate A numerical value for the rotation for the Epoch/Series temporal units, default values are 0 when direction = "upwards" and 90 when direction = "rightwards".
#' @param arotate A numerical value for the rotation for the Age/Stage temporal units, default values are 0 when direction = "upwards" and 90 when direction = "rightwards".
#' @param urotate A numerical value for the rotation for the User temporal units, default values are 0 when direction = "upwards" and 90 when direction = "rightwards".
#' @param ... All other arguments passed to plot.phylo
#'
#' @details
#'
#' Often palaeontologists wish to display phylogenetic trees against a geological time scale for easy visualization of the temporal position of the taxa the generation of which can be time-consuming. geoscalePhylo fills this need and works in ths same way as geoscale.plot from the package geoscale and allows users to plot a time-scaled phylogeny against the International Chronostratigraphic Chart (Gradstein, 2014). This function accepts any tree time-scaled through either the function DatePhylo in this package or the timePaleoPhy function in the library paleotree.
#'
#' Built-in options allows the user control over which direction the tree is plotted in (either horizonally or vertically) as well as deciding which temporal units are included in the time scale (see below for example).
#'
#' Temporal units
#'
#' The function geoscalePhylo allows for a time-scaled phylogeny to be plotted against geologic time using either the current geologic time scale of Gradstein et al., 2012 or previously published time scales by the International Commisioin on Stratigraphy. The time scale that is plotted is comprised of a number of temporal components representing the different units that the geological time scale is divided into. There are five main temporal units that can be included, each of which have two alternative names and are as follows: Eon (Eonothem), Era (Erathem), System (Period), Series (Epoch), and Stage (Age). These alternative names can be used interchangably i.e. both Eon and Erathem are accepted, however should both these alternative names be included then that temporal unit will only be included once. In addition, the order in which they are included into units does not affect the order in which they appear in the chart so units=c("Period","Epoch","Age") will produce the same results as units=c("Age","Epoch","Period") with the default order as they were listed previously with Eons plotted at the base and Stages at the top.
#'
#' Including a user-defined time scale
#'
#' There is a sixth option that can be included into the units argument. "User" allows for an additional temporal unit to be plotted i.e. biozonal or terrane-specific time scales. This requires a matrix of three columns named "Start", "End" and "Name" representing the bottom, top of each temporal bin (in millions of years) and the name to be plotted respectively. An example dataset called UKzones representing Stages of the UK Ordovician System is included in the package. See below for an example of how to implement this option.
#'
#' Stratigraphic ranges
#'
#' geoscalePhylo allows for the stratigraphic ranges to be included in the plot. This requires an matrix with the first appearance and last appearance dates in millions of years (FAD and LAD respectively) with the row names containing all the tip labels of the taxa in the tree, exactly as they appear in tree$tip.label and the column names should be "FAD" and "LAD". In order to add the stratigraphic ranges to the plot this matrix should be attached to the argument ages. See below for an example of this option.
#'
#' Apparent appearance of polytomies
#'
#' It should be noted that using certain methods for time-scaling a tree, such as the "basic" method (the default), it can create the appearance of polytomies in a tree is otherwise fully resolved due to the presence of a large number of zero length branches. This can be solved by using another timescaling method such as the "equal" method which will enforce all the branches to have a positive length.
#'
#' @return Nothing (simply produces a plot of the tree against geologic time).
#'
#' @author Mark A. Bell \email{mark.bell521@@gmail.com}
#'
#' @references
#'
#' Gradstein, F. M., Ogg, J. M., and Schmitz, M. 2012. A Geologic Time Scale. Elsevier, Boston, USA.
#'
#' @examples
#'
#' ### Example lungfish data
#' utils::data(Dipnoi)
#'
#' tree_l <- DatePhylo(Dipnoi$tree, Dipnoi$ages, method = "equal", rlen = 1)
#'
#' geoscalePhylo(tree = tree_l, boxes = "Age", cex.tip = 0.4)
#'
#' # Plotting the tree with the stratigraphical ranges included
#' geoscalePhylo(tree = tree_l, ages = Dipnoi$ages, boxes = "Age", cex.tip = 0.4)
#'
#' # Including all temporal units into the stratigraphic column
#' geoscalePhylo(tree_l, Dipnoi$ages, units = c("Eon", "Era", "Period", "Epoch", "Age"),
#'   boxes = "Age", cex.tip = 0.4)
#'
#' # Plotting the numerical values on the time scale at Age resolution
#' geoscalePhylo(tree_l, Dipnoi$ages, units = c("Eon", "Era", "Period", "Epoch", "Age"),
#'   boxes="Age", cex.tip = 0.4, tick.scale = "Age")
#'
#' ### Example trilobite data
#' utils::data(Asaphidae)
#'
#' tree_a <- DatePhylo(Asaphidae$trees[[1]], Asaphidae$ages, method = "equal", rlen = 1)
#'
#' geoscalePhylo(ladderize(tree_a), Asaphidae$ages, boxes = "Age", x.lim = c(504, 435),
#'   cex.tip = 0.5, cex.ts = 0.5, vers = "ICS2009")
#'
#' # Plotting the tree vertically
#' geoscalePhylo(ladderize(tree_a), Asaphidae$ages, boxes = "Age", x.lim = c(504, 435),
#'   cex.tip = 0.5, cex.ts = 0.5, direction = "upwards", vers = "ICS2009")
#'
#' # Including a user-defined time scale
#' utils::data(UKzones)
#' utils::data(Asaphidae)
#'
#' tree_a <- DatePhylo(Asaphidae$trees[[1]], Asaphidae$ages, method = "equal", rlen = 1)
#'
#' geoscalePhylo(ladderize(tree_a), Asaphidae$ages, units = c("Eon", "Era", "Period",
#'   "Epoch", "User"), boxes = "Age", cex.tip = 0.4, user.scale = UKzones,
#'   vers = "ICS2009", cex.ts = 0.5, x.lim = c(520, 440), direction = "upwards")
#'
#' # Rotating the text on the time scale
#' tree_a <- DatePhylo(Asaphidae$trees[[1]], Asaphidae$ages, method = "equal", rlen = 1)
#'
#' #geoscalePhylo(ladderize(tree_a), Asaphidae$ages, units = c("Period",
#' #  "Epoch", "Age", "User"), boxes = "Age", cex.tip = 0.4, user.scale = UKzones,
#' #  vers = "ICS2009", cex.ts = 0.5, x.lim = c(520, 440), arotate = 0, erotate = 0, urotate = 0)
#'
#' @export geoscalePhylo
geoscalePhylo <- function(tree, ages, direction = "rightwards", units = c("Period", "Epoch", "Age"), boxes = "Age", tick.scale = "myr", user.scale, cex.age = 0.3, cex.ts = 0.3, cex.tip = 0.3, width = 1, label.offset,ts.col = TRUE, vers = "ICS2013", x.lim, quat.rm = FALSE, erotate, arotate, urotate, ...) {
  
  options <- as.list(match.call())
   if(any(names(options) == "type")){
     if(all(options$type != c("phylogram", "cladogram", "p", "c"))){
       return(cat("type must be either 'phylogram' or 'cladogram'."))
     }
   }

  if(all(direction != c("rightwards","upwards"))){
      return(cat("direction must be either 'rightwards' or 'upwards', here set to 'rightwards'."))
      direction <- "rightwards"
    }
  
  if(is.null(tree$root.time)){     
    return(cat("\n tree$root.time is missing, check tree is time scaled."))
  } else {root.age <- tree$root.time}
  
  if(boxes == "User" && any(units != "User")){
    boxes <- "no"
  }
  
  if(all(boxes != units)){
    boxes <- "no"
  }
  
  if(tick.scale == "User" && all(units != "User")){
    tick.scale <- "myr"
  }
   
  if(missing(ages) == FALSE){
      ranges <- TRUE
  } else{
      ranges <- FALSE
  }
  
  if(missing(user.scale) & any(units == "User")){
    units <- units[units != "User"]
    cat("\n user.scale not provided, 'Other' removed from units.")
  }
  
  if(missing(ages) == FALSE){
    ages<-ages[tree$tip.label,]    
  }
  
  if(any(units == "User") & !missing(user.scale)){   
    Midpoint <- matrix(ncol=1,nrow=length(user.scale[,1]))
      Midpoint[,1] <- (user.scale[,"Start"] + user.scale[,"End"])/2
        user.scale <- cbind(user.scale,Midpoint)  
  }
  
  if(all(units != "Age") && boxes == "Age"){
    boxes <- "no"
  }
  
  units <- paste(toupper(substring(units,1,1)),substring(units,2),sep="")
    
  # Standardizing the names of temporal units
  
  units[units == "Eonothem"] <- "Eon"
  units[units == "Erathem"] <- "Era"
  units[units == "Series"] <- "Epoch"
  units[units == "System"] <- "Period"  
  units[units == "Stage"] <- "Age"
    units <- unique(units)
  
  boxes[boxes == "Eonothem"] <- "Eon"
  boxes[boxes == "Erathem"] <- "Era"
  boxes[boxes == "Series"] <- "Epoch"
  boxes[boxes == "System"] <- "Period"  
  boxes[boxes == "Stage"] <- "Age"
      
  if(length(units) == 1){
    ts.width=0.15
  } else if(length(units) == 2){
    ts.width=0.2
  } else if(length(units) >= 3){
    ts.width=0.25
  }

  if(ranges == TRUE && missing(ages) == FALSE){
     missing.tip.names <- setdiff(tree$tip.label,row.names(ages))
      if(length(missing.tip.names) > 0){            
        cat(paste("\n",missing.tip.names,"not present in ages file, ranges set to FALSE"))
        cat("\n ranges set to FALSE")
          ranges <- FALSE
      }    
  }
   
  timescales <- NULL
   utils::data(timescales,envir=environment())
    timescale <- timescales[[vers]]
      if(quat.rm == TRUE){
       timescale[(timescale[,"Midpoint"] < 3),"Name"] <- NA
      }
    
  tscale.data<-matrix(ncol=3,nrow=6)
    colnames(tscale.data) <-c("srt","Depth","size")
    rownames(tscale.data) <-c("Eon","Era","Period","Epoch","Age","User")
      if(direction == "upwards"){
        tscale.data[,"srt"] <- c(90,90,90,0,0,0)
      } else tscale.data[,"srt"] <- c(0,0,0,90,90,90)
      
      tscale.data[,"Depth"] <- c(1,1,1,2,3.5,3.5)
      tscale.data[,"size"] <- c(1,1,1,0.8,0.8,0.8)

  ## ROTATING THE NAMES
  
  if(!missing(erotate) && !is.numeric(erotate)){
    return(cat("\n value for protate must be numeric."))
  }
  if(!missing(arotate) && !is.numeric(arotate)){
    return(cat("\n value for arotate must be numeric."))
  }
  if(!missing(urotate) && !is.numeric(urotate)){
    return(cat("\n value for urotate must be numeric."))
  }
  if(!missing(erotate)){
    tscale.data["Epoch","srt"] <- erotate
  }
  if(!missing(arotate)){
    tscale.data["Age","srt"] <- arotate
  }
  if(!missing(urotate)){
    tscale.data["User","srt"] <- urotate
  }

  ## GEOLOGICAL RANGES OF TAXA

  units<-rownames(tscale.data)[sort(match(units,rownames(tscale.data)),decreasing=T)] 
    
  if(!missing(x.lim) ){
    x.lim <- sort(root.age - x.lim)
  } else if(ranges == TRUE && !missing(ages) && missing(x.lim)){
    x.lim <- (root.age - min(ages)) + diff(range(ages))*0.05
  } else {
    x.lim <- NULL
  }
    
    timescale<-timescale[order(timescale[,1],decreasing=T),]
    	timescale.rescaled <- timescale
  		timescale.rescaled[,c("Start","End","Midpoint")] <- root.age - timescale[,c("Start","End","Midpoint")]

  #if(root.age > 5){
   # timescale[timescale[,"Name"] == "Holocene" && !is.na(timescale[,"Name"] == "Holocene"),"Name"] <- NA
  #}
  
  first_la <- tree$root.time - dist.nodes(tree)[1,Ntip(tree)+1]
	
    if(ranges==TRUE && missing(ages) == FALSE){
      offset<-array(dim=length(tree$tip.label),data=1)
       offset.correction <- diff(range(ages)) * 0.01 
        taxon.ranges <- root.age - ages[,c("FAD","LAD")]
         if(first_la != ages[1,"LAD"]){
          if(!missing(label.offset)){
           offset <- array(dim=length(ages[,"FAD"]),data=(ages[,"FAD"] - ages[,"LAD"])+label.offset)
          } else {
           offset <- array(dim=length(ages[,"FAD"]),data=(ages[,"FAD"] - ages[,"LAD"])+offset.correction)
          }
         }
    } else if(!missing(label.offset)){
      offset <- label.offset 
    } else {
      offset = 0
    }

  ### ADDING A TIMESCALE

  if(tick.scale != "n" | tick.scale != "no"){
    if(tick.scale == "myr" | is.numeric(tick.scale)){    
      scale.ticks=1
      
      if(is.numeric(tick.scale)){
       scale.ages <- tick.scale
      } else {scale.ages=10}
        
        tick.position <- root.age - seq(0,4600,scale.ticks)
          age.name <- seq(0,4600,scale.ages)
          age.position <- root.age - age.name
            lwd<-c(1,0.5,0.5,0.5,0.5,0.7,0.5,0.5,0.5,0.5)
            col<-c("black","grey","grey","grey","grey")
    }

    if(tick.scale != "myr" & is.numeric(tick.scale) == FALSE){
      age.name<-subset(timescale,timescale[,"Type"] == tick.scale & timescale[,"Source"] == "ICS")
        age.name<-sort(unique(c(age.name[,"Start"],age.name[,"End"])))
        	age.position<- tick.position <- root.age - age.name
            lwd=1
              col="black" 
    } 
    
    if(tick.scale == "User"){     
      age.name <- sort(unique(c(user.scale[,"Start"],user.scale[,"End"])))
       age.position <- tick.position <- root.age - age.name
          lwd=1
            col="black"        
      }
   }
  
  # Plotting the tree vertically
  if(direction == "upwards"){
    
    # ADDING THE TIME SCALE
    graphics::par(fig=c(0,ts.width,0,1))
     graphics::par(mar=c(3,1,2,5))
    
      graphics::par(lend=1); plot.phylo(tree,plot=FALSE,label.offset=offset,no.margin=T,y.lim=x.lim,direction="upwards",cex=cex.tip,...)
          
        timescale.names <- timescale.rescaled
         timescale.names <- timescale.names[timescale.names[,"End"] > min(c(graphics::par()$usr[3],graphics::par()$usr[4])) & timescale.names[,"Start"] < max(c(graphics::par()$usr[3],graphics::par()$usr[4])),]
          timescale.names[timescale.names[,"End"] > max(c(graphics::par()$usr[3],graphics::par()$usr[4])),"End"] <- max(c(graphics::par()$usr[3],graphics::par()$usr[4]))
          timescale.names[timescale.names[,"Start"] < min(c(graphics::par()$usr[3],graphics::par()$usr[4])),"Start"] <- min(c(graphics::par()$usr[3],graphics::par()$usr[4]))
    
          timescale.names[,"Midpoint"] <- (timescale.names[,"Start"] + timescale.names[,"End"]) /2     
          timescale.names[,"Range"] <- abs(timescale.names[,"Start"] - timescale.names[,"End"])    
    
            for(t in 1:length(timescale.names[,1])){
             if(timescale.names[t,"Range"] < timescale[rownames(timescale.names)[t],"Range"]*0.25){
              timescale.names[t,"Name"] <- NA 
              }
            }
    
          unit.depths <- tscale.data[units,"Depth"]
            if(tick.scale == "n" | tick.scale == "no"){
              unit.depths <- c(unit.depths,0.5)          
            } else if(length(units) <= 3){
              unit.depths <- c(unit.depths,2)
            } else if(length(units) > 3) {
              unit.depths <- c(unit.depths,2)
            }

           unit.depths <- cumsum(unit.depths/sum(unit.depths))
            unit.depths<-c(graphics::par()$usr[2],graphics::par()$usr[2]-(unit.depths*(graphics::par()$usr[2]-graphics::par()$usr[1])))
             
         depth <- unit.depths[length(unit.depths) - 1] - unit.depths[length(unit.depths)]
    
        if(tick.scale != "n" && tick.scale != "no"){
          graphics::text((unit.depths[length(unit.depths)]+depth*0.3),age.position,age.name,cex=cex.age,srt=0)           
           graphics::segments((unit.depths[length(unit.depths)-1]),tick.position,(unit.depths[length(unit.depths)]+depth*0.75),tick.position,lwd=lwd,col=col)
        }
    
    for(t in 1:length(units)){
      
      if(units[t] == "User"){
        tscale<-user.scale
         tscale[,c("Start","End","Midpoint")] <- root.age - tscale[,c("Start","End","Midpoint")]
         tscale.names <- tscale
      } else {
        tscale<-subset(timescale.rescaled,timescale.rescaled[,"Type"] == units[t])
        tscale.names<-subset(timescale.names,timescale.names[,"Type"] == units[t])
      }
          
      if(ts.col == TRUE & units[t] != "User"){
        graphics::rect(tscale[,"Start"],unit.depths[t],tscale[,"End"],unit.depths[t+1], col=grDevices::rgb(tscale[,"Col_R"], tscale[,"Col_G"], tscale[,"Col_B"], maxColorValue=255))
      } else graphics::rect(tscale[,"Start"],unit.depths[t],tscale[,"End"],unit.depths[t+1],col="white")
      graphics::text((unit.depths[t] + unit.depths[t+1])/2,tscale.names[,"Midpoint"],tscale.names[,"Name"],cex=cex.ts*tscale.data[match(units[t],rownames(tscale.data)),"size"],srt=tscale.data[match(units[t],rownames(tscale.data)),"srt"])
    }
    
    # ADDING THE PHYLOGENY
    graphics::par(fig=c(ts.width,1,0,1),new=T)
     graphics::par(mar=c(3,0,2,2))
    
      graphics::par(lend=1); plot.phylo(tree,plot=FALSE,label.offset=offset,no.margin=T,y.lim=x.lim,direction="upwards",cex=cex.tip,...)
        lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    
      if (!missing(boxes) && boxes != "no" && boxes !="n"){
       if(boxes == "User"){
         tscale <- root.age - user.scale[,c("Start","End")]
       } else {
         tscale<-subset(timescale.rescaled,timescale.rescaled[,"Type"] == boxes)}
       
           graphics::rect(graphics::par()$usr[3],tscale[,"Start"],graphics::par()$usr[4],tscale[,"End"],col=c("grey90","white"),border=NA)
       }
        
      graphics::par(fig=c(ts.width,1,0,1),new=T)
        graphics::par(mar=c(3,0,2,2))
    
         graphics::par(lend=1); plot.phylo(tree,label.offset=offset,edge.width=width,no.margin=T,y.lim=x.lim,cex=cex.tip,direction="upwards",...)
          if (ranges == TRUE){
            graphics::par(lend=1); graphics::segments(lastPP$xx[c(1:length(tree$tip.label))],taxon.ranges[,"FAD"],lastPP$xx[c(1:length(tree$tip.label))],taxon.ranges[,"LAD"],col="black",lwd=width*2)
          }
    
    if(units[1] == "User"){
      graphics::segments(graphics::par()$usr[1],min(root.age - user.scale[,"Start"]),graphics::par()$usr[1],max(root.age - user.scale[,"End"]))
    } else {      
      graphics::segments(graphics::par()$usr[1],min(timescale.rescaled[timescale.rescaled[,"Type"] == units[1],"Start"]),graphics::par()$usr[1],max(timescale.rescaled[timescale.rescaled[,"Type"] == units[1],"End"]))
    }
       
  } else{
    # Plotting the tree horizontally
    
  # PLOT 1 - TIMESCALE
    graphics::par(fig=c(0,1,0,ts.width))
    graphics::par(mar=c(1,3,0,2))

      graphics::par(lend=1); plot.phylo(tree,plot=FALSE,label.offset=offset,no.margin=T,x.lim=x.lim,direction="rightwards",cex=cex.tip,...)
            
        timescale.names <- timescale.rescaled
         timescale.names <- timescale.names[timescale.names[,"End"] > min(c(graphics::par()$usr[1],graphics::par()$usr[2])) & timescale.names[,"Start"] < max(c(graphics::par()$usr[1],graphics::par()$usr[2])),]
          timescale.names[timescale.names[,"End"] > max(c(graphics::par()$usr[1],graphics::par()$usr[2])),"End"] <- max(c(graphics::par()$usr[1],graphics::par()$usr[2]))
          timescale.names[timescale.names[,"Start"] < min(c(graphics::par()$usr[1],graphics::par()$usr[2])),"Start"] <- min(c(graphics::par()$usr[1],graphics::par()$usr[2]))
    
            timescale.names[,"Midpoint"] <- (timescale.names[,"Start"] + timescale.names[,"End"]) /2     
            timescale.names[,"Range"] <- abs(timescale.names[,"Start"] - timescale.names[,"End"])    
    
              for(t in 1:length(timescale.names[,1])){
               if(timescale.names[t,"Range"] < timescale[rownames(timescale.names)[t],"Range"]*0.25){
                timescale.names[t,"Name"] <- NA 
               }
              }
        
       unit.depths <- tscale.data[units,"Depth"]
        if(tick.scale == "n" & tick.scale == "no"){
          unit.depths <- c(unit.depths,0.5)          
        } else if(length(units) <= 3){
          unit.depths <- c(unit.depths,2)
        } else if(length(units) > 3) {
          unit.depths <- c(unit.depths,2)
        }
     
         unit.depths <- cumsum(unit.depths/sum(unit.depths))
          unit.depths<-c(graphics::par()$usr[4],graphics::par()$usr[4]-(unit.depths*(graphics::par()$usr[4]-graphics::par()$usr[3])))
    
       depth <- unit.depths[length(unit.depths) - 1] - unit.depths[length(unit.depths)]
    
    if(tick.scale != "n" && tick.scale != "no"){
      graphics::text(age.position,(unit.depths[length(unit.depths)]+depth*0.3),age.name, cex=cex.age,srt=90)
  	  graphics::segments(tick.position,(unit.depths[length(unit.depths)-1]),tick.position,(unit.depths[length(unit.depths)]+depth*0.6),lwd=lwd,col=col)
    }

    for(t in 1:length(units)){
 	
 	    if(units[t] == "User"){
 	    	tscale<-user.scale
  	    	tscale[,c("Start","End","Midpoint")] <- root.age - tscale[,c("Start","End","Midpoint")]
     	    	tscale.names <- tscale
      } else {
 	      tscale<-subset(timescale.rescaled,timescale.rescaled[,"Type"] == units[t])
 	      tscale.names<-subset(timescale.names,timescale.names[,"Type"] == units[t])
 	    }
 	 	
 	    if(ts.col == TRUE & units[t] != "User"){graphics::rect(tscale[,"Start"],unit.depths[t],tscale[,"End"],unit.depths[t+1],col=grDevices::rgb(tscale[,"Col_R"],tscale[,"Col_G"],tscale[,"Col_B"],maxColorValue=255))
 		    } else graphics::rect(tscale[,"Start"],unit.depths[t],tscale[,"End"],unit.depths[t+1],col="white")
          graphics::text(tscale.names[tscale.names[,"Type"] == units[t],"Midpoint"],(unit.depths[t] + unit.depths[t+1])/2,tscale.names[tscale.names[,"Type"] == units[t],"Name"],cex=cex.ts*tscale.data[match(units[t],rownames(tscale.data)),"size"],srt=tscale.data[match(units[t],rownames(tscale.data)),"srt"])
    }
 
  ## PLOT 2: PHYLOGENY

    graphics::par(fig=c(0,1,ts.width,1),new=T)
    graphics::par(mar=c(0,3,2,2))
  
      graphics::par(lend=1); plot.phylo(tree,plot=FALSE,label.offset=offset,no.margin=T,x.lim=x.lim,cex=cex.tip,...)
        lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  
      if (!missing(boxes) && boxes != "no" && boxes !="n"){
        if(boxes == "User"){
          tscale <- root.age - user.scale[,c("Start","End")]
        } else {tscale<-subset(timescale.rescaled,timescale.rescaled[,"Type"] == boxes)} 
  	        graphics::rect(tscale[,"Start"],graphics::par()$usr[3],tscale[,"End"],graphics::par()$usr[4],col=c("grey90","white"),border=NA)}

    graphics::par(fig=c(0,1,ts.width,1),new=T)
    graphics::par(mar=c(0,3,2,2))
  
      graphics::par(lend=1); plot.phylo(tree,label.offset=offset,edge.width=width,no.margin=T,x.lim=x.lim,cex=cex.tip,...)
        if (ranges == TRUE){
          graphics::par(lend=1); graphics::segments(taxon.ranges[,"FAD"],lastPP$yy[c(1:length(tree$tip.label))],taxon.ranges[,"LAD"],lastPP$yy[c(1:length(tree$tip.label))],col="black",lwd=width*2)
        }
  
  if(units[1] == "User"){
    graphics::segments(min(root.age - user.scale[,"Start"]),graphics::par()$usr[3],max(root.age - user.scale[,"End"]),graphics::par()$usr[3])
  } else {      
    graphics::segments(min(timescale.rescaled[timescale.rescaled[,"Type"] == units[1],"Start"]),graphics::par()$usr[3],max(timescale.rescaled[timescale.rescaled[,"Type"] == units[1],"End"]),graphics::par()$usr[3])
  }
  
  }
}
