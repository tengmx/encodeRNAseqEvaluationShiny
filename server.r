library(shiny)
hkctr <- function(quant,hkpfile,genefile){
    library(matrixStats,quietly=T)
    gene <- read.delim(genefile,skip=5,header=FALSE,stringsAsFactors=FALSE,
            colClasses=c("character","character","character","integer","integer",
              "character","character","character"),comment.char="",quote="")
    housekeep <- read.delim(hkpfile,header=FALSE,stringsAsFactors=FALSE,
                            colClasses=c("character","character"),sep=" ")
    feature <- strsplit(gene[,9],"; ")
    gene_name <- sapply(sapply(feature,function(x) x[substr(x,1,9)=="gene_name"]),
                        function(x) gsub("\"","",gsub("gene_name ","",x)))
    housekeep_index <- !is.na(match(gene_name,housekeep[,1])) #56564:56668 spike-in
    colMeans(log2(sapply(quant,function(x) colMedians(x[housekeep_index,]))))
}

typefilter <- function(quant,genetype,genefile){
    gene <- read.delim(genefile,skip=5,header=FALSE,stringsAsFactors=FALSE,
            colClasses=c("character","character","character","integer","integer",
                  "character","character","character"),comment.char="",quote="")
    feature <- strsplit(gene[,9],"; ")
    genetypes <- sapply(sapply(feature,function(x) x[substr(x,1,9)=="gene_type"]),
                        function(x) gsub("\"","",gsub("gene_type ","",x)))
    names(genetypes) <- NULL
    quantable <- list()
    for(i in 1:length(quant)){
        quantable[[i]] <- quant[[i]][genetypes==genetype,]
    }
    quantable
}

cutfilter <- function(quant,cuts){
    for(i in 1:length(quant)){
        quant[[i]][quant[[i]]<cuts[i]] <- 0
    }
    quant
}

sdaplot_scale <- function(quant,medians,colors,text,xlim,ylim,xlab,ylab,main){
    for(i in 1:length(quant)){
        quant[[i]] <- log2(quant[[i]])
        rep11 <- quant[[i]][quant[[i]][,1]!=-Inf & quant[[i]][,2]!=-Inf,1]
        rep12 <- quant[[i]][quant[[i]][,1]!=-Inf & quant[[i]][,2]!=-Inf,2]
        rep21 <- quant[[i]][quant[[i]][,3]!=-Inf & quant[[i]][,4]!=-Inf,3]
        rep22 <- quant[[i]][quant[[i]][,3]!=-Inf & quant[[i]][,4]!=-Inf,4]
        A1 <- (rep11+rep12)/2 -medians[i]
        A2 <- (rep21+rep22)/2 -medians[i]
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
        lines(quantile(xlim,c(0.5,0.6)),quantile(ylim,rep(1.04-i*0.08,2)),
              type='l',col=colors[i],lwd=2)
        text(quantile(xlim,0.65),quantile(ylim,1.04-i*0.08),pos=4,
             labels=paste0(text[i],"_GM12878"))
        lines(quantile(xlim,c(0.5,0.6)),quantile(ylim,rep(1-i*0.08,2)),
              type='l',col=colors[i],lwd=2,lty=2)
        text(quantile(xlim,0.65),quantile(ylim,1-i*0.08),pos=4,
             labels=paste0(text[i],"_K562"))
    }
    text(quantile(xlim,0.1),quantile(ylim,0.95),pos=4,
                      labels=paste0("FPKM:1 ~ D:",round(-mean(medians[4:5]),2)))
}

p0plot_scale <- function(quant,medians,step,colors,text,xlim,ylim,xlab,ylab,main){
    for(i in 1:length(quant)){
        index1 <- quant[[i]][,1]!=0 & quant[[i]][,2]!=0
        index2 <- quant[[i]][,3]!=0 & quant[[i]][,4]!=0
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
        allunits1 <- length(index1)
        both0prop1 <- round(sum(index1_00) / allunits1,3)
        one0prop1 <- round(sum(index1_0) / allunits1, 3)
        bothnon0prop1 <- round(sum(index1) / allunits1, 3)
        allunits2 <- length(index2)
        both0prop2 <- round(sum(index2_00) / allunits2,3)
        one0prop2 <- round(sum(index2_0) / allunits2, 3)
        bothnon0prop2 <- round(sum(index2) / allunits2, 3)
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
             labels=paste0(text[i],"_G"))
        text(quantile(xlim,seq(0.7,0.9,0.1)),quantile(ylim,1.04-i*0.08),pos=4,
             labels=c(both0prop1,one0prop1,bothnon0prop1))
        lines(k,p2,type='l',col=colors[i],lty=2,lwd=2)
        lines(quantile(xlim,c(0.4,0.45)),quantile(ylim,rep(1-i*0.08,2)),
              type='l',col=colors[i],lwd=2,lty=2)
        text(quantile(xlim,0.47),quantile(ylim,1-i*0.08),pos=4,
             labels=paste0(text[i],"_K"))
        text(quantile(xlim,seq(0.7,0.9,0.1)),quantile(ylim,1-i*0.08),pos=4,
             labels=c(both0prop2,one0prop2,bothnon0prop2))
    }
    text(quantile(xlim,seq(0.7,0.9,0.1)),quantile(ylim,1),pos=4,
         labels=c("0,0","0,1","1,1"))
}

## CAT plot on 4-column table list
catplot <- function(quant,colors,text,xlim,ylim,xlab,ylab,main,stringent=TRUE){
    for(i in 1:length(quant)){
        quant[[i]] <- log2(quant[[i]])
        fc1 <- quant[[i]][,1] - quant[[i]][,3]
        fc2 <- quant[[i]][,2] - quant[[i]][,4]
        if(stringent){
            fc1tmp <- fc1[!is.na(fc1) & !is.na(fc2)]
            fc2tmp <- fc2[!is.na(fc1) & !is.na(fc2)]
            fc1<- fc1tmp[fc1tmp!=Inf & fc1tmp!=-Inf & fc2tmp!=Inf & fc2tmp!=-Inf]
            fc2<- fc2tmp[fc1tmp!=Inf & fc1tmp!=-Inf & fc2tmp!=Inf & fc2tmp!=-Inf]
        }else{
            fc1 <- fc1[!is.na(fc1) & fc1!=Inf & fc1!=-Inf]
            fc2 <- fc2[!is.na(fc2) & fc2!=Inf & fc2!=-Inf]
        }
        names1 <- names(fc1)[sort.list(abs(fc1),decreasing=T)]
        names2 <- names(fc2)[sort.list(abs(fc2),decreasing=T)]
        endrank <- min(xlim[2],length(names1),length(names2))
        prop<- sapply(xlim[1]:endrank, function(x)
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

## CAT plot on 4-column table list (microarray, the same name set)
catplot_microarray <- function(quant,array,colors,text,xlim,ylim,
                               xlab,ylab,main,stringent=TRUE){
    load(array)
    fc2 <- fc[sort.list(abs(fc),decreasing=T)]
    array_names <- strsplit(names(fc2),split="_")
    seq_names <- sapply(rownames(quant[[1]]),function(x) strsplit(x,"\\.")[[1]][1])
    idx1 <- sapply(array_names,function(x) sum(!is.na(match(x,seq_names)))>0)
    array_names <- array_names[idx1]
    seq_names <- rownames(quant[[1]])[!is.na(match(seq_names,unlist(array_names)))]
    for(i in 1:length(quant)){
        quant[[i]] <- log2(quant[[i]][seq_names,])
        fc1 <- (quant[[i]][,1] + quant[[i]][,2] - quant[[i]][,3] - quant[[i]][,4])/2
        fc1 <- fc1[!is.na(fc1) & fc1!=Inf & fc1!=-Inf]
        names_seq <- names(fc1)[sort.list(abs(fc1),decreasing=T)]
        names_seq <- sapply(names_seq,function(x) strsplit(x,"\\.")[[1]][1])
        if(stringent){
            names_array <- array_names[sapply(array_names,function(x)
                               sum(!is.na(match(x,names_seq)))>0)]
        }else{
            names_array <- array_names
        }
        endrank <- min(xlim[2],length(names_seq),length(names_array))
        prop<- sapply(xlim[1]:endrank, function(x)
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

colors <- c("brown","red","royalblue","seagreen","olivedrab1","purple",
            "maroon1","black","orange","yellow")
labels <- c("rsem","rsem_pme","flux","cuff_s","cuff_t","sailfish","express","naive")
shinyServer(function(input, output) {
  packs <- reactive({
    load(paste0(input$protocol,"_g.rda"))
    medians <- hkctr(quant,"HK_genes.txt","gene.gtf")
    thresholds <- 2^(input$cutD+medians)
    quant <- typefilter(quant,input$genetype,"gene.gtf")
    quant <- cutfilter(quant,thresholds)
    cat(thresholds,"\n")
    all0kick <- sapply(quant,max)>0
    quant <- quant[all0kick]
    labels <- labels[all0kick]
    colors2 <- colors[which(all0kick)]
    medians <- medians[all0kick]
    list(quant=quant,medians=medians,labels=labels,colors=colors2)
  })
  output$caption <- renderText({
    input$protocol
  })
  output$sdplot <- renderPlot({
    pack <- packs()
    sdaplot_scale(pack$quant,pack$medians,pack$colors,pack$labels,
                  xlim=c(input$xstart1,input$xend1),ylim=c(input$ystart1,input$yend1),
                  xlab="D: A-log2(median(ObsSgl_ctr))",ylab="SD",main=input$genetype)
  })
  output$p0plot <- renderPlot({
      pack <- packs()
      p0plot_scale(pack$quant,pack$medians,step=0.1,pack$colors,pack$labels,
                   xlim=c(input$xstart2,input$xend2),ylim=c(input$ystart2,input$yend2),
                   xlab="D: log2(k)-log2(median(ObsSgl_ctr))",ylab="Proportion",main=input$genetype)
  })
  output$catplot1 <- renderPlot({
      pack <- packs()
      catplot(pack$quant,pack$colors,pack$labels,
              xlim=c(input$xstart3,input$xend3),ylim=c(input$ystart3,input$yend3),
              xlab="Size of list",ylab="proportion in common",main=input$genetype)
  })
  output$catplot2 <- renderPlot({
    pack <- packs()
    catplot(pack$quant,pack$colors,pack$labels,
            xlim=c(input$xstart3,input$xend3),ylim=c(input$ystart3,input$yend3),
            xlab="Size of list",ylab="proportion in common",main=input$genetype,stringent=F)
  })
  output$catplotarray1 <- renderPlot({
      pack <- packs()
      catplot_microarray(pack$quant,"cel.rda",pack$colors,pack$labels,
                         xlim=c(input$xstart4,input$xend4),ylim=c(input$ystart4,input$yend4),
                         xlab="Size of list",ylab="proportion in common",main=input$genetype)
  })
  output$catplotarray2 <- renderPlot({
      pack <- packs()
      catplot_microarray(pack$quant,"cel.rda",pack$colors,pack$labels,
                         xlim=c(input$xstart4,input$xend4),ylim=c(input$ystart4,input$yend4),
                         xlab="Size of list",ylab="proportion in common",main=input$genetype,stringent=F)
  })
})
