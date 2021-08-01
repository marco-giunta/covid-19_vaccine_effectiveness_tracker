DbdaAcfPlot = function( codaObject , parName=varnames(codaObject)[1] , plColors=NULL ) {
  if ( all( parName != varnames(codaObject) ) ) { 
    stop("parName must be a column name of coda object")
  }
  nChain = length(codaObject)
  if ( is.null(plColors) ) plColors=1:nChain
  xMat = NULL
  yMat = NULL
  for ( cIdx in 1:nChain ) {
    acfInfo = acf(codaObject[,c(parName)][[cIdx]],plot=FALSE) 
    xMat = cbind(xMat,acfInfo$lag)
    yMat = cbind(yMat,acfInfo$acf)
  }
  matplot( xMat , yMat , type="o" , pch=20 , col=plColors , ylim=c(0,1) ,
           main="" , xlab="Lag" , ylab="Autocorrelation" )
  abline(h=0,lty="dashed")
  EffChnLngth = effectiveSize(codaObject[,c(parName)])
  text( x=max(xMat) , y=max(yMat) , adj=c(1.0,1.0) , cex=1.25 ,
        labels=paste("ESS =",round(EffChnLngth,1)) )
}


DbdaDensPlot = function( codaObject , parName=varnames(codaObject)[1] , plColors=NULL ) {
  if ( all( parName != varnames(codaObject) ) ) { 
    stop("parName must be a column name of coda object")
  }
  nChain = length(codaObject) # or nchain(codaObject)
  if ( is.null(plColors) ) plColors=1:nChain
  xMat = NULL
  yMat = NULL
  hdiLims = NULL
  for ( cIdx in 1:nChain ) {
    densInfo = density(codaObject[,c(parName)][[cIdx]]) 
    xMat = cbind(xMat,densInfo$x)
    yMat = cbind(yMat,densInfo$y)
    hdiLims = cbind(hdiLims,HDIofMCMC(codaObject[,c(parName)][[cIdx]]))
  }
  matplot( xMat , yMat , type="l" , col=plColors , 
           main="" , xlab="Param. Value" , ylab="Density" )
  abline(h=0)
  points( hdiLims[1,] , rep(0,nChain) , col=plColors , pch="|" )
  points( hdiLims[2,] , rep(0,nChain) , col=plColors , pch="|" )
  text( mean(hdiLims) , 0 , "95% HDI" , adj=c(0.5,-0.2) )
  EffChnLngth = effectiveSize(codaObject[,c(parName)])
  MCSE = sd(as.matrix(codaObject[,c(parName)]))/sqrt(EffChnLngth) 
  text( max(xMat) , max(yMat) , adj=c(1.0,1.0) , cex=1.25 ,
        paste("MCSE =\n",signif(MCSE,3)) )
}



HDIofMCMC = function( sampleVec , credMass=0.95 ) {
  # Computes highest density interval from a sample of representative values,
  #   estimated as shortest credible interval.
  # Arguments:
  #   sampleVec
  #     is a vector of representative values from a probability distribution.
  #   credMass
  #     is a scalar between 0 and 1, indicating the mass within the credible
  #     interval that is to be estimated.
  # Value:
  #   HDIlim is a vector containing the limits of the HDI
  sortedPts = sort( sampleVec )
  ciIdxInc = ceiling( credMass * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
  }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
  HDIlim = c( HDImin , HDImax )
  return( HDIlim )
}



plotPost = function( paramSampleVec , cenTend=c("mode","median","mean")[1] , 
                     credMass=0.95, HDItextPlace=0.7, hdicol = 'dodgerblue',
                     xlab=NULL , xlim=NULL , yaxt=NULL , ylab=NULL , 
                     main=NULL , cex=NULL , cex.lab=NULL ,
                     col=NULL , border=NULL , showCurve=TRUE , breaks=NULL  , 
                     ... ) { #
  # Override defaults of hist function, if not specified by user:
  # (additional arguments "..." are passed to the hist function)
  if ( is.null(xlab) ) xlab="Param. Val."
  if ( is.null(cex.lab) ) cex.lab=1.5
  if ( is.null(cex) ) cex=1.4
  if ( is.null(xlim) ) xlim=range( c(paramSampleVec ) )
  if ( is.null(main) ) main=""
  if ( is.null(yaxt) ) yaxt="n"
  if ( is.null(ylab) ) ylab=""
  if ( is.null(col) ) col="grey"
  if ( is.null(border) ) border="white"
  
  # convert coda object to matrix:
  if (any(class(paramSampleVec) == "mcmc.list")) { # modifica messa qui
    paramSampleVec = as.matrix(paramSampleVec)
  }
  
  summaryColNames = c("ESS","mean","median","mode",
                      "hdiMass","hdiLow","hdiHigh")
  postSummary = matrix( NA , nrow=1 , ncol=length(summaryColNames) , 
                        dimnames=list( c( xlab ) , summaryColNames ) )
  
  # require(coda) # for effectiveSize function
  postSummary[,"ESS"] = effectiveSize(paramSampleVec)
  
  postSummary[,"mean"] = mean(paramSampleVec)
  postSummary[,"median"] = median(paramSampleVec)
  mcmcDensity = density(paramSampleVec)
  postSummary[,"mode"] = mcmcDensity$x[which.max(mcmcDensity$y)]
  
  HDI = HDIofMCMC( paramSampleVec , credMass )
  postSummary[,"hdiMass"]=credMass
  postSummary[,"hdiLow"]=HDI[1]
  postSummary[,"hdiHigh"]=HDI[2]
  
  # Plot histogram.
  cvCol = "darkgreen"
  ropeCol = "darkred"
  if ( is.null(breaks) ) {
    if ( max(paramSampleVec) > min(paramSampleVec) ) {
      breaks = c( seq( from=min(paramSampleVec) , to=max(paramSampleVec) ,
                       by=(HDI[2]-HDI[1])/18 ) , max(paramSampleVec) )
    } else {
      breaks=c(min(paramSampleVec)-1.0E-6,max(paramSampleVec)+1.0E-6)
      border="skyblue"
    }
  }
  if ( !showCurve ) {
    par(xpd=NA)
    histinfo = hist( paramSampleVec , xlab=xlab , yaxt=yaxt , ylab=ylab ,
                     freq=F , border=border , col=col ,
                     xlim=xlim , main=main , cex=cex , cex.lab=cex.lab ,
                     breaks=breaks  , ...) #
  }
  if ( showCurve ) {
    par(xpd=NA)
    histinfo = hist( paramSampleVec , xlab=xlab , yaxt=yaxt , ylab=ylab ,
                     freq=F , border=border , col=col ,
                     xlim=xlim , main=main , cex=cex , cex.lab=cex.lab ,
                     breaks=breaks  , ... ) #
    densCurve = density( paramSampleVec , adjust=2 )
    lines( densCurve$x , densCurve$y , type="l" , lwd=3 , col="black" , bty="n" ,
          xlim=xlim , xlab=xlab , yaxt=yaxt , ylab=ylab ,
          main=main , cex=cex , cex.lab=cex.lab  , ... ) #
    abline(v = postSummary[,"mode"], lwd = 3)
    l <- min(which(densCurve$x >= HDI[1]))
    h <- max(which(densCurve$x < HDI[2]))

    polygon(c(densCurve$x[c(l, l:h, h)]),
            c(0, densCurve$y[l:h], 0),
            col = hdicol, density = 60)
    
    
  }
  cenTendHt = 0.9*max(histinfo$density)
  cvHt = 0.7*max(histinfo$density)
  ROPEtextHt = 0.55*max(histinfo$density)
  # Display central tendency:
  mn = mean(paramSampleVec)
  med = median(paramSampleVec)
  mcmcDensity = density(paramSampleVec)
  mo = mcmcDensity$x[which.max(mcmcDensity$y)]
  if ( cenTend=="mode" ){ 
    text( mo , cenTendHt ,
          bquote(mode==.(signif(mo,3))) , adj=c(.5,0) , cex=cex )
  }
  if ( cenTend=="median" ){ 
    text( med , cenTendHt ,
          bquote(median==.(signif(med,3))) , adj=c(.5,0) , cex=cex , col=cvCol )
  }
  if ( cenTend=="mean" ){ 
    text( mn , cenTendHt ,
          bquote(mean==.(signif(mn,3))) , adj=c(.5,0) , cex=cex )
  }
    
  # Display the HDI.
  lines( HDI , c(0,0) , lwd=4 , lend=1 )
  text( mean(HDI) , 0 , bquote(.(100*credMass) * "% HDI" ) ,
        adj=c(.5,-1.7) , cex=cex )
  text( HDI[1] , 0 , bquote(.(signif(HDI[1],3))) ,
        adj=c(HDItextPlace,-0.5) , cex=cex )
  text( HDI[2] , 0 , bquote(.(signif(HDI[2],3))) ,
        adj=c(1.0-HDItextPlace,-0.5) , cex=cex )
  par(xpd=F)
  #
  return( postSummary )
}



diagMCMC = function( codaObject , parName=varnames(codaObject)[1], clr = TRUE, title = 'theta') {
  if(isTRUE(clr)){
    DBDAplColors = c("skyblue","black","royalblue","steelblue")}
    else { 
    DBDAplColors = c("green","darkgreen","seagreen","greenyellow")
    }
  par( mar=0.5+c(3,4,1,0) , oma=0.1+c(0,0,2,0) , mgp=c(2.25,0.7,0) , 
       cex.lab=1.5 )
  layout(matrix(1:4,nrow=2))
  # traceplot and gelman.plot are from CODA package:
  require(coda)
  coda::traceplot( codaObject[,c(parName)] , main="" , ylab="Param. Value" ,
                   col=DBDAplColors) 
  tryVal = try(
    coda::gelman.plot( codaObject[,c(parName)] , main="" , auto.layout=FALSE , 
                       col=DBDAplColors )
  )  
  # if it runs, gelman.plot returns a list with finite shrink values:
  if ( class(tryVal)=="try-error" ) {
    plot.new() 
    print(paste0("Warning: coda::gelman.plot fails for ",parName))
  } else { 
    if ( class(tryVal)=="list" & !is.finite(tryVal$shrink[1]) ) {
      plot.new() 
      print(paste0("Warning: coda::gelman.plot fails for ",parName))
    }
  }
  DbdaAcfPlot(codaObject,parName,plColors=DBDAplColors)
  DbdaDensPlot(codaObject,parName,plColors=DBDAplColors)
  mtext( text=title , outer=TRUE , adj=c(0.5,0.5) , cex=2.0 )
}


barcomparator <- function(lista){
    medie <- c()
    up <- c()
    down <- c()
    title <- c('Pfizer old', 'Pfizer new', 'Moderna', 'Astrazeneca', 'Janssen')
    for (i in 1:length(lista)) {
        medie[i] <- mean(lista[[i]])
        up[i] <- HDIofMCMC(lista[[i]], credMass = 0.95)[1]
        down[i] <- HDIofMCMC(lista[[i]], credMass = 0.95)[2]
    }

    bar <- barplot(medie*100, xpd = F, names=title, ylim=c(50,100), density=60, axis.lty=1, las=1, 
                   col=c('yellow','greenyellow','red','blue','green'),
                   main='Vaccine effectiveness\nand their credibility interval',
                   ylab='VE')
    arrows(bar, up*100, bar, down*100, angle=90, code=3, length=0.2,
               col=c('yellow','greenyellow','red','blue','green'), lwd=4)
    points(bar,medie*100, lwd=6, col=c('yellow','greenyellow','red','blue','green'))
}



ageinference <- function(df, pv, pp){
    indici <- grepl("^.*[0-9]", rownames(df))
    paramSampleVec <- list()
    for (i in 1:nrow(df[indici,])) {
        paramSampleVec[[i]] <- (1 - pv[[i]]/pp[[i]])*100
    }
    summaryColNames = c("ESS","mean","median","mode",
                          "hdiMass","hdiLow","hdiHigh")
    postSummary = matrix( NA , nrow=nrow(df[indici,]) , ncol=length(summaryColNames) , 
                            dimnames=list(c(rownames(df[indici,])), summaryColNames ) )

    for (i in 1:nrow(df[indici,])) {
    # require(coda) # for effectiveSize function
        postSummary[i,"ESS"] = effectiveSize(paramSampleVec[[i]])

        postSummary[i,"mean"] = mean(paramSampleVec[[i]])
        postSummary[i,"median"] = median(paramSampleVec[[i]])
        mcmcDensity = density(paramSampleVec[[i]])
        postSummary[i,"mode"] = mcmcDensity$x[which.max(mcmcDensity$y)]

        HDI = HDIofMCMC( paramSampleVec[[i]], credMass = 0.95)
        postSummary[i,"hdiMass"]=0.95
        postSummary[i,"hdiLow"]=HDI[1]
        postSummary[i,"hdiHigh"]=HDI[2]
    }
    postSummary

}

load('./mcmc_data/mcmc_ai.RData')

pfizer.age.fun <- function(df = pfizer.df, pv = pfizer.pv, pp = pfizer.pp, col = c("#FFCC00","#FFCC66","#CC9900","#996600")) {
    pfizer.age <- ai_p #read.csv('./age_data/pfizer_age.csv') # ageinference(df, pv, pp)
    plot(c(8,100), c(10,105), type= "n", xaxt ="n", xlab = "Age range", ylab = "Effectiveness", main = "Pfizer", cex.main=2, cex.axis=1.5, cex.lab=1.5)
    rect(12, pfizer.age[1,][[6]], 15, pfizer.age[1,][[7]], density = 80, border = col[1], col = col[1])
    rect(16, pfizer.age[2,][[6]], 64, pfizer.age[2,][[7]], density = 80, border = col[2], col = col[2])
    rect(65, pfizer.age[3,][[6]], 74, pfizer.age[3,][[7]], density = 80, border = col[3], col = col[3])
    rect(75, pfizer.age[4,][[6]], 100,pfizer.age[3,][[7]], density = 80, border = col[4], col = col[4])
    points((12+15)/2,pfizer.age[1,][[2]], lwd=5, col=col[1])
    points((16+64)/2,pfizer.age[2,][[2]], lwd=5, col=col[2])
    points((65+74)/2,pfizer.age[3,][[2]], lwd=5, col=col[3])
    points((75+100)/2,pfizer.age[4,][[2]], lwd=5, col=col[4])
    axis(side = 1, at = c(12,16,65,75,100), las=2)
}

moderna.age.fun <- function(df = moderna.df, pv = moderna.pv, pp = moderna.pp, col = c("#FF3300","#CC3300","#990000")) {
    moderna.age <- ai_m #read.csv('./age_data/moderna_age.csv') # ageinference(df, pv, pp)
    plot(c(8,100), c(10,105), type= "n", xaxt ="n", xlab = "Age range", ylab = "Effectiveness", main = "Moderna", cex.main=2, cex.axis=1.5, cex.lab=1.5)
    rect(12, moderna.age[1,][[6]], 17, moderna.age[1,][[7]], density = 80, border = col[1], col = col[1])
    rect(18, moderna.age[2,][[6]], 64, moderna.age[2,][[7]], density = 80, border = col[2], col = col[2])
    rect(65, moderna.age[3,][[6]], 100, moderna.age[3,][[7]], density = 80,border = col[3], col = col[3])
    points((12+17)/2,moderna.age[1,][[2]], lwd=5, col=col[1])
    points((18+64)/2,moderna.age[2,][[2]], lwd=5, col=col[2])
    points((65+100)/2,moderna.age[3,][[2]], lwd=5,col=col[3])
    axis(side = 1, at = c(12,18,65,100), las=2)
}

janssen.age.fun <- function(df = janssen.df, pv = janssen.pv, pp = janssen.pp, col = c("#99FF33","#339900")) {
    janssen.age <- ai_j #read.csv('./age_data/janssen_age.csv') # ageinference(df, pv, pp)
    plot(c(8,100), c(10,105), type= "n", xaxt ="n", xlab = "Age range", ylab = "Effectiveness", main = "Janssen", cex.main=2, cex.axis=1.5, cex.lab=1.5)
    rect(18, janssen.age[1,][[6]], 59, janssen.age[1,][[7]], density = 80, border = col[1], col = col[1])
    rect(60, janssen.age[2,][[6]], 100, janssen.age[2,][[7]], density = 80, border = col[2], col = col[2])
    points((18+59)/2,janssen.age[1,][[2]], lwd=5, col=col[1])
    points((60+100)/2,janssen.age[2,][[2]], lwd=5, col=col[2])
    axis(side = 1, at = c(18,60,100), las=2)
}



plotPost1 <- function(paramSampleVec, postdensity ,cenTend=c("mode","median","mean")[1] , 
                     credMass=0.95, HDItextPlace=0.7, hdicol1 = 'dodgerblue',hdicol2 = 'red',
                     xlab=NULL , xlim=NULL , yaxt=NULL , ylab=NULL , 
                     main=NULL , cex=NULL , cex.lab=NULL ,
                     col=NULL , border=NULL , breaks=NULL , 
                     ... ) {
  # Override defaults of hist function, if not specified by user:
  # (additional arguments "..." are passed to the hist function)
  if ( is.null(xlab) ) xlab="Param. Val."
  if ( is.null(cex.lab) ) cex.lab=1.5
  if ( is.null(cex) ) cex=1.4
  if ( is.null(xlim) ) xlim=range( c(paramSampleVec ) )
  if ( is.null(main) ) main=""
  if ( is.null(yaxt) ) yaxt="n"
  if ( is.null(ylab) ) ylab=""
  if ( is.null(col) ) col="grey"
  if ( is.null(border) ) border="white"
  
  # convert coda object to matrix:
  if (any(class(paramSampleVec) == "mcmc.list")) {
    paramSampleVec = as.matrix(paramSampleVec)
  }
  
  summaryColNames <- c("ESS","mean","median","mode",
                      "hdiMass","hdiLow","hdiHigh")
  postSummary <- matrix( NA , nrow=1 , ncol=length(summaryColNames) , 
                        dimnames=list( c( xlab ) , summaryColNames ) )
  
  # require(coda) # for effectiveSize function
  postSummary[,"ESS"] <- effectiveSize(paramSampleVec)
  
  postSummary[,"mean"] <- mean(paramSampleVec)
  postSummary[,"median"] <- median(paramSampleVec)
  mcmcDensity <- density(paramSampleVec)
  postSummary[,"mode"] <- mcmcDensity$x[which.max(mcmcDensity$y)]
  
  HDI1 <- HDIofMCMC( paramSampleVec , credMass )
  postSummary[,"hdiMass"] <- credMass
  postSummary[,"hdiLow"] <- HDI1[1]
  postSummary[,"hdiHigh"] <- HDI1[2]
  
    densCurve = density(paramSampleVec , adjust=2 )
    plot( densCurve$x , densCurve$y , type="l" , lwd=4 , col="black" , bty="n" ,
          xlim=xlim , xlab=xlab , yaxt=yaxt , ylab=ylab ,
          main=main , cex=cex , cex.lab=cex.lab , ... )
#     abline(v = postSummary[,"mode"], lwd = 3)
    l1 <- min(which(densCurve$x >= HDI1[1]))
    h1 <- max(which(densCurve$x < HDI1[2]))
    
    polygon(c(densCurve$x[c(l1, l1:h1, h1)]),
            c(0, densCurve$y[l1:h1], 0),
            col = hdicol1, density = 60)
    
    lines(postdensity, col="red",lwd=3)
    l2 <- min(which(postdensity$x >= HDI1[1]))
    h2 <- max(which(postdensity$x < HDI1[2]))
    polygon(c(postdensity$x[c(l2, l2:h2, h2)]),
            c(0, postdensity$y[l2:h2], 0),
            col = hdicol2, density = 60)
    
    legend("topleft", legend = c("Posterior density with conjugate prior", "Posterior density with mcmc"), col = c(hdicol1,hdicol2),
           fill = c(hdicol1,hdicol2))
    
    lines( HDI1 , c(0,0) , lwd=4 , lend=1 )
    text( mean(HDI1) , 0 , bquote(.(100*credMass) * "% HDI" ) ,
          adj=c(.5,-1.7) , cex=cex )
  }




priorcomparison <- function(posterior1, posterior2, 
                     credMass=0.95, HDItextPlace=0.7, hdicol1 = 'dodgerblue',hdicol2 = 'green',
                     xlab=NULL , xlim=NULL , ylim=NULL, yaxt=NULL , ylab=NULL , 
                     main=NULL) 
{
  if ( is.null(xlab) ) xlab = "Param. Val."
  if ( is.null(xlim) ) xlim = range( c(posterior1$x ) )
  if ( is.null(ylim) ) ylim = range( c(posterior2$y ) )
  if ( is.null(main) ) main = ""
  if ( is.null(ylab) ) ylab = "Posterior distribution"
  
    plot(posterior1$x , posterior1$y , type="l" , bty="n", lwd=3 , col=hdicol1,
          xlim=xlim , ylim = ylim, xlab=xlab , yaxt=yaxt , ylab=ylab , main=main)
    abline(v = posterior1$x[which.max(posterior1$y)], lwd = 3, col="darkblue")

    lines(posterior2, col=hdicol2,lwd=3)
    abline(v = posterior2$x[which.max(posterior2$y)], lwd = 3, col="darkgreen")
    
    legend("topleft", legend = c("Posterior density with beta prior", "Posterior density with uniform prior", 
           sprintf("mode = %.2f",posterior1$x[which.max(posterior1$y)]), sprintf("mode = %.2f",posterior2$x[which.max(posterior2$y)])),
           col = c(hdicol1,hdicol2,"darkblue","darkgreen"),fill = c(hdicol1,hdicol2,"darkblue","darkgreen"))
    
    # HDI1 = c(89,97)
    # l1 <- min(which(posterior1$x >= HDI1[1]))
    # h1 <- max(which(posterior1$x < HDI1[2]))
    # polygon(c(posterior1$x[c(l1, l1:h1, h1)]),
    #         c(0, posterior1$y[l1:h1], 0),
    #         col = hdicol1, density = 60)
    
    # HDI2 = c(92,99)
    # l2 <- min(which(posterior2$x >= HDI2[1]))
    # h2 <- max(which(posterior2$x < HDI2[2]))
    # polygon(c(posterior2$x[c(l2, l2:h2, h2)]),
    #         c(0, posterior2$y[l2:h2], 0),
    #         col = hdicol2, density = 60)
    
    
    # text( mean(HDI1) , 0 , bquote(.(100*credMass) * "% HDI" ) ,
    #       adj=c(.5,-1.7) , cex=1.4 )

  }

posterior <- function(df.pv, df.pp) {
    yv <- sum(df.pv)
    nv <- length(df.pv)
    yp <- sum(df.pp)
    np <- length(df.pp)
    #     set.seed(1234)
    pv <- rbeta(10**5, shape1 = 2.88 + yv, shape2 = 93.12 + nv - yv)
    #     set.seed(1234)
    pp <- rbeta(10**5, shape1 = 1.41 + yp, shape2 = 66.11 + np - yp)
    eff <- ((pp - pv)/pp)*100
    den <- density(eff)
    return(den)
}

posterior.uniform <- function(df.pv, df.pp) {
    yv <- sum(df.pv)
    nv <- length(df.pv)
    yp <- sum(df.pp)
    np <- length(df.pp)
    
    pv <- rbinom(10^5,size = nv, prob = yv/nv)
    pp <- rbinom(10^5,size = np, prob = yp/np)
    eff <- ((pp - pv)/pp)*100
    den <- density(eff)
    return(den)
}