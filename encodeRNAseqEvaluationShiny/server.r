library(shiny)

## calculate standardizing coefficients based on house keeping genes
hkctr <- function(quant){
    library(matrixStats,quietly=T)
    load("housekeep_index.rda")
    colMeans(log2(sapply(quant,function(x) colMedians(x[housekeep_index,]))))
}

## filter quantification table by gene types
typefilter <- function(quant,genetype){
    load("genetypes.rda")
    for(i in 1:length(quant)){
        quant[[i]] <- quant[[i]][genetypes==genetype,]
    }
    quant
}

## set signals less than cutoffs to zero
cutfilter <- function(quant,cuts){
    for(i in 1:length(quant)){
        quant[[i]][quant[[i]]<cuts[i]] <- 0
    }
    quant
}

## plot repliates variation stratified by D 
sdaplot_scale <- function(quant,medians,colors,text,xlim,ylim,xlab,ylab,main){
    for(i in 1:length(quant)){
        quant[[i]] <- log2(quant[[i]])
        rep11 <- quant[[i]][quant[[i]][,1]!=-Inf & quant[[i]][,2]!=-Inf,1]
        rep12 <- quant[[i]][quant[[i]][,1]!=-Inf & quant[[i]][,2]!=-Inf,2]
        rep21 <- quant[[i]][quant[[i]][,3]!=-Inf & quant[[i]][,4]!=-Inf,3]
        rep22 <- quant[[i]][quant[[i]][,3]!=-Inf & quant[[i]][,4]!=-Inf,4]
        A1 <- (rep11+rep12)/2 - medians[i]
        A2 <- (rep21+rep22)/2 - medians[i]
        SD1 <- abs(rep11-rep12)
        SD2 <- abs(rep21-rep22)
        tmp1 <- loess.smooth(A1, SD1, span = 2/3, degree = 1,
                            family = "symmetric", evaluation = 1000)
        tmp2 <- loess.smooth(A2, SD2, span = 2/3, degree = 1,
                             family = "symmetric", evaluation = 1000)
        if(i==1){
            plot(tmp1$x,tmp1$y,type='l',col=colors[i],lwd=2,xlim=xlim,ylim=ylim,
                 xlab=xlab,ylab=ylab,main=main)
        }else{
            lines(tmp1$x,tmp1$y,type='l',col=colors[i],lwd=2)
        }
        lines(tmp2$x,tmp2$y,type='l',col=colors[i],lty=2,lwd=2)
        lines(quantile(xlim,c(0.6,0.7)),quantile(ylim,rep(1-i*0.04,2)),
              type='l',col=colors[i],lwd=2)
        text(quantile(xlim,0.73),quantile(ylim,1-i*0.04),pos=4,labels=text[i])
    }
    lines(quantile(xlim,c(0.2,0.3)),quantile(ylim,rep(0.92,2)),
          type='l',lwd=2)
    text(quantile(xlim,0.33),quantile(ylim,0.92),pos=4,labels="GM12878")
    lines(quantile(xlim,c(0.2,0.3)),quantile(ylim,rep(0.88,2)),
          type='l',lwd=2,lty=2)
    text(quantile(xlim,0.33),quantile(ylim,0.88),pos=4,labels="K562")
}

## plot 0s' information stratified by D 
p0plot_scale0 <- function(quant,medians,step,colors,text,xlim,ylim,xlab,ylab,main){
    for(i in 1:length(quant)){
        index1_0 <- (quant[[i]][,1]==0 & quant[[i]][,2]!=0) |
                    (quant[[i]][,1]!=0 & quant[[i]][,2]==0)
        index2_0 <- (quant[[i]][,3]==0 & quant[[i]][,4]!=0) |
                    (quant[[i]][,3]!=0 & quant[[i]][,4]==0)
        index1_00 <- quant[[i]][,1]==0 & quant[[i]][,2]==0
        index2_00 <- quant[[i]][,3]==0 & quant[[i]][,4]==0
        rep11_0 <- quant[[i]][index1_0,1]
        rep12_0 <- quant[[i]][index1_0,2]
        rep21_0 <- quant[[i]][index2_0,3]
        rep22_0 <- quant[[i]][index2_0,4]
        allunits1 <- length(index1_00)
        both0prop1 <- round(sum(index1_00) / allunits1,3)
        one0prop1 <- round(sum(index1_0) / allunits1, 3)
        bothnon0prop1 <- 1-both0prop1-one0prop1
        allunits2 <- length(index2_00)
        both0prop2 <- round(sum(index2_00) / allunits2,3)
        one0prop2 <- round(sum(index2_0) / allunits2, 3)
        bothnon0prop2 <- 1-both0prop2-one0prop2
        k <- seq(xlim[1],xlim[2],step)
        p1 <- sapply(k,function(x)
                     sum(log2(rep11_0+rep12_0)-medians[i]>x)) / allunits1
        p2 <- sapply(k,function(x)
                     sum(log2(rep21_0+rep22_0)-medians[i]>x)) / allunits2
        if(i==1){
            plot(k,p1,type='l',col=colors[i],lwd=2,xlim=xlim,ylim=ylim,
                 xlab=xlab,ylab=ylab,main=main)
        }else{
            lines(k,p1,type='l',col=colors[i],lwd=2)
        }
        lines(quantile(xlim,c(0.4,0.45)),quantile(ylim,rep(1.04-i*0.08,2)),
              type='l',col=colors[i],lwd=2)
        text(quantile(xlim,0.47),quantile(ylim,1.04-i*0.08),pos=4,
             labels=paste0(text[i],"_GM12878"))
        text(quantile(xlim,seq(0.7,0.9,0.1)),quantile(ylim,1.04-i*0.08),pos=4,
             labels=c(both0prop1,one0prop1,bothnon0prop1))
        lines(k,p2,type='l',col=colors[i],lty=2,lwd=2)
        lines(quantile(xlim,c(0.4,0.45)),quantile(ylim,rep(1-i*0.08,2)),
              type='l',col=colors[i],lwd=2,lty=2)
        text(quantile(xlim,0.47),quantile(ylim,1-i*0.08),pos=4,
             labels=paste0(text[i],"_K562"))
        text(quantile(xlim,seq(0.7,0.9,0.1)),quantile(ylim,1-i*0.08),pos=4,
             labels=c(both0prop2,one0prop2,bothnon0prop2))
    }
    text(quantile(xlim,seq(0.7,0.9,0.1)),quantile(ylim,1),pos=4,
         labels=c("0,0","0,1","1,1"))    
}


p0plot_scale1 <- function(quant,medians,step,colors,text,xlim,ylim,xlab,ylab,main){
    for(i in 1:length(quant)){
        index1_0 <- (quant[[i]][,1]==0 & quant[[i]][,2]!=0) |
                    (quant[[i]][,1]!=0 & quant[[i]][,2]==0)
        index2_0 <- (quant[[i]][,3]==0 & quant[[i]][,4]!=0) |
                    (quant[[i]][,3]!=0 & quant[[i]][,4]==0)
        rep11_0 <- quant[[i]][index1_0,1]
        rep12_0 <- quant[[i]][index1_0,2]
        rep21_0 <- quant[[i]][index2_0,3]
        rep22_0 <- quant[[i]][index2_0,4]
        allunits1 <- length(index1_0)
        allunits2 <- length(index2_0)
        k <- seq(xlim[1],xlim[2],step)
        p1 <- sapply(k,function(x)
                     sum(log2(rep11_0+rep12_0)-medians[i]>x)) / allunits1
        p2 <- sapply(k,function(x)
                     sum(log2(rep21_0+rep22_0)-medians[i]>x)) / allunits2
        if(i==1){
            plot(k,p1,type='l',col=colors[i],lwd=2,xlim=xlim,ylim=ylim,
                 xlab=xlab,ylab=ylab,main=main)
        }else{
            lines(k,p1,type='l',col=colors[i],lwd=2)
        }
        lines(k,p2,type='l',col=colors[i],lty=2,lwd=2)
        lines(quantile(xlim,c(0.6,0.7)),quantile(ylim,rep(1-i*0.04,2)),
              type='l',col=colors[i],lwd=2)
        text(quantile(xlim,0.73),quantile(ylim,1-i*0.04),pos=4,labels=text[i])
    }
    lines(quantile(xlim,c(0.2,0.3)),quantile(ylim,rep(0.92,2)),
          type='l',lwd=2)
    text(quantile(xlim,0.33),quantile(ylim,0.92),pos=4,labels="GM12878")
    lines(quantile(xlim,c(0.2,0.3)),quantile(ylim,rep(0.88,2)),
          type='l',lwd=2,lty=2)
    text(quantile(xlim,0.33),quantile(ylim,0.88),pos=4,labels="K562")
}

p0plot_scale2 <- function(quant,text){
  p0s <- matrix(0,length(quant),6)
    for(i in 1:length(quant)){
        index1_0 <- (quant[[i]][,1]==0 & quant[[i]][,2]!=0) |
                    (quant[[i]][,1]!=0 & quant[[i]][,2]==0)
        index2_0 <- (quant[[i]][,3]==0 & quant[[i]][,4]!=0) |
                    (quant[[i]][,3]!=0 & quant[[i]][,4]==0)
        index1_00 <- quant[[i]][,1]==0 & quant[[i]][,2]==0
        index2_00 <- quant[[i]][,3]==0 & quant[[i]][,4]==0
        rep11_0 <- quant[[i]][index1_0,1]
        rep12_0 <- quant[[i]][index1_0,2]
        rep21_0 <- quant[[i]][index2_0,3]
        rep22_0 <- quant[[i]][index2_0,4]
        allunits1 <- length(index1_00)
        both0prop1 <- round(sum(index1_00) / allunits1,3)
        one0prop1 <- round(sum(index1_0) / allunits1, 3)
        bothnon0prop1 <- 1-both0prop1-one0prop1
        allunits2 <- length(index2_00)
        both0prop2 <- round(sum(index2_00) / allunits2,3)
        one0prop2 <- round(sum(index2_0) / allunits2, 3)
        bothnon0prop2 <- 1-both0prop2-one0prop2
        p0s[i,1:3] <- c(both0prop1,one0prop1,bothnon0prop1)
        p0s[i,4:6] <- c(both0prop2,one0prop2,bothnon0prop2)
    }
  nacol <- paste(rep(c("GM12878","K562"),each=3),rep(c("0,0","0,1","1,1"),2),sep=" ")
  rownames(p0s) <- text
  colnames(p0s) <- nacol
    p0s
}

## CAT plot between replicates
catplot <- function(quant,medians,fpkmmedian,cut,colors,text,xlim,ylim,xlab,ylab,
                    main,stringent=TRUE){
    cutD <- log2(cut)-fpkmmedian
    cuts <- 2^(cutD+medians)
    for(i in 1:length(quant)){
        if(!stringent)  quant[[i]][quant[[i]]==0] <- cuts[i]
        quant[[i]] <- log2(quant[[i]])
        fc1 <- quant[[i]][,1] - quant[[i]][,3]
        fc2 <- quant[[i]][,2] - quant[[i]][,4]
        if(stringent){
            #fc1tmp <- fc1[!is.na(fc1) & !is.na(fc2)]
            #fc2tmp <- fc2[!is.na(fc1) & !is.na(fc2)]
            #fc1 <- fc1tmp[fc1tmp!=Inf & fc1tmp!=-Inf & fc2tmp!=Inf & fc2tmp!=-Inf]
            #fc2 <- fc2tmp[fc1tmp!=Inf & fc1tmp!=-Inf & fc2tmp!=Inf & fc2tmp!=-Inf]
            fc1[is.na(fc1) | fc1==-Inf | fc1==Inf] <- 0
            fc2[is.na(fc2) | fc2==-Inf | fc2==Inf] <- 0
        }
        names1 <- names(fc1)[sort.list(abs(fc1),decreasing=T)]
        names2 <- names(fc2)[sort.list(abs(fc2),decreasing=T)]
        endrank <- min(xlim[2],length(names1),length(names2))
        prop <- sapply(xlim[1]:endrank, function(x)
                      round(length(intersect(names1[1:x],names2[1:x]))/x,4))
        if(i==1){
            plot(xlim[1]:endrank,prop,type='l',col=colors[i],lwd=2,
                 ylim=ylim,xlim=xlim,xlab=xlab,ylab=ylab,main=main)
        }else{
            lines(xlim[1]:endrank,prop,type='l',col=colors[i],lwd=2)
        }
        lines(quantile(xlim,c(0.6,0.7)),quantile(ylim,rep(0.42-i*0.04,2)),
              type='l',col=colors[i],lwd=2)
        text(quantile(xlim,0.73),quantile(ylim,0.42-i*0.04),
             pos=4,labels=text[i])
    }
}

## CAT plot between average of two replicates and microarray
catplot_microarray <- function(quant,medians,fpkmmedian,cut,array,colors,text,xlim,ylim,
                               xlab,ylab,main,stringent=TRUE){
    load(array)
    seq_names <- rownames(quant[[1]])[!is.na(match(rownames(quant[[1]]),array_names))]
    cutD <- log2(cut)-fpkmmedian
    cuts <- 2^(cutD+medians)
    for(i in 1:length(quant)){
        if(!stringent)  quant[[i]][quant[[i]]==0] <- cuts[i]
        quant[[i]] <- log2(quant[[i]][seq_names,])
        fc1 <- (quant[[i]][,1] + quant[[i]][,2] - quant[[i]][,3] - quant[[i]][,4])/2
        if(stringent){
            fc1 <- fc1[!is.na(fc1) & fc1!=Inf & fc1!=-Inf]
        }
        names_seq <- names(fc1)[sort.list(abs(fc1),decreasing=T)]
        names_array <- array_names[!is.na(match(array_names,names_seq))]
        endrank <- min(xlim[2],length(names_seq),length(names_array))
        prop <- sapply(xlim[1]:endrank, function(x)
               round(length(intersect(names_seq[1:x],unlist(names_array[1:x])))/x,4))
        if(i==1){
            plot(xlim[1]:endrank,prop,type='l',col=colors[i],lwd=2,xlim=xlim,
                 ylim=ylim,xlab=xlab,ylab=ylab,main=main)
        }else{
            lines(xlim[1]:endrank,prop,type='l',col=colors[i],lwd=2)
        }
        lines(quantile(xlim,c(0.6,0.7)),quantile(ylim,rep(1-i*0.04,2)),
              type='l',col=colors[i],lwd=2)
        text(quantile(xlim,0.73),quantile(ylim,1-i*0.04),
             pos=4,labels=text[i])
    }
}

## shiny server
colors <- c("brown","red","royalblue","seagreen","olivedrab1","purple",
            "maroon1","black","orange","yellow")
labels <- c("RSEM tpm","RSEM tpmpme","Flux Capacitor","Cufflinks with STAR","Cufflinks with Tophat","Sailfish","eXpress","Naive")
shinyServer(function(input, output) {
    packs <- reactive({
      cat(input$protocol,'\t',input$genetype,"\n")
        load(paste0(input$protocol,"_g.rda"))
        medians <- hkctr(quant)
      FPKMmedian <- sum(medians[4:5])/2
      cutD <- log2(input$cutFPKM)-FPKMmedian
        thresholds <- 2^(cutD+medians)
        quant <- typefilter(quant,input$genetype)
        quant <- cutfilter(quant,thresholds)
        cat(thresholds,"\n")
        all0kick <- sapply(quant,max)>0
        quant <- quant[all0kick]
        labels <- labels[all0kick]
        colors2 <- colors[which(all0kick)]
        medians <- medians[all0kick]
        thresholds <- thresholds[all0kick]
        list(quant=quant,medians=medians,labels=labels,colors=colors2,FPKMmedian=FPKMmedian)
    })
    output$caption <- renderText({
        input$protocol
    })
    output$sdplot <- renderPlot({
        pack <- packs()
        sdaplot_scale(pack$quant,pack$medians,pack$colors,pack$labels,
                      xlim=c(input$xstart1,input$xend1),ylim=c(input$ystart1,input$yend1),
                      xlab="Detrended log signal",ylab="SD",main=NULL)
    })
    output$p0plot <- renderPlot({
        pack <- packs()
        p0plot_scale1(pack$quant,pack$medians,step=0.1,pack$colors,pack$labels,
                     xlim=c(input$xstart2,input$xend2),ylim=c(input$ystart2,input$yend2),
                     xlab="Detrended log signal",ylab="Proportion",main=NULL)
    })
    output$p0stbl <- renderTable({
        pack <- packs()
        p0plot_scale2(pack$quant,pack$labels)
    })
    output$catplot1 <- renderPlot({
        pack <- packs()
        catplot(pack$quant,pack$medians,pack$FPKMmedian,input$constant1,pack$colors,pack$labels,
                xlim=c(input$xstart3,input$xend3),ylim=c(input$ystart3,input$yend3),
                xlab="Size of list",ylab="Proportion in common",main=NULL)
    })
    output$catplot2 <- renderPlot({
        pack <- packs()
        catplot(pack$quant,pack$medians,pack$FPKMmedian,input$constant1,pack$colors,pack$labels,
                xlim=c(input$xstart3,input$xend3),ylim=c(input$ystart3,input$yend3),
                xlab="Size of list",ylab="Proportion in common",main=NULL,stringent=F)
    })
    output$catplotarray1 <- renderPlot({
        pack <- packs()
        catplot_microarray(pack$quant,pack$medians,pack$FPKMmedian,input$constant2,"array.rda",pack$colors,pack$labels,
                           xlim=c(input$xstart4,input$xend4),ylim=c(input$ystart4,input$yend4),
                           xlab="Size of list",ylab="Proportion in common",main=NULL)
    })
    output$catplotarray2 <- renderPlot({
        pack <- packs()
        catplot_microarray(pack$quant,pack$medians,pack$FPKMmedian,input$constant2,"array.rda",pack$colors,pack$labels,
                           xlim=c(input$xstart4,input$xend4),ylim=c(input$ystart4,input$yend4),
                           xlab="Size of list",ylab="Proportion in common",main=NULL,stringent=F)
    })
})
