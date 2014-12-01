#Description:
#   This function clusters, then splits by clusters a dataframe of statistics
#   **Note: with a large number of genes (rows) and/or statistics (columns), this function can be very slow
#     Plotting as a pdf file speeds things up 2x. 
#input format:
#stats.in-- a dataframe with your differential expression statistics
#   **Note: the first two columns need to be named "main.model.categories" and "PACID"
#     These columns should have the categories for each gene... e.g. significant in trt only, not sig. etc.
#     PACID is just the name of the gene
#     The remaining columns should be the statistics or values (e.g. lsmeans) you want to cluster and plot

#categories-- character vector: the factors in stats.in$main.model.categories you want to focus on
#tstat.type-- character vector: the columns to use in stats.in
#n.clusters-- the number of groups to make in hclust
#plot.title-- the title of the plot
#pdf.file--the name of the pdf file to write to... to write into the R window, specific "to.window"
summary.exp1<-function(stats.in, 
                       categories,
                       tstat.type,
                       n.clusters,
                       plot.title,
                       pdf.file){
  sub.data<-stats.in[stats.in[,1] %in% categories,c("main.model.categories","PACID",tstat.type)]
  
  dist.t<-dist(sub.data[,tstat.type])
  hc<-hclust(dist.t)
  d.out<-colour_clusters(hc,n.clusters, col=rainbow(n.clusters))
  #split clusters into 6 groups, do line plots for each group
  cut.in<-slice(hc,k=n.clusters)[order(hc$order)]
  sub.data$hcgroup<-data.frame(cut.in)[,1]
  
  t.melt<-melt(sub.data,id.vars=c("main.model.categories","PACID","hcgroup"))
  colnames(t.melt)[4:5]<-c("comparison","t.stat")
  t.melt$hcgroup<-as.factor(t.melt$hcgroup)
  
  cat.counts<-data.frame(table(t.melt$hcgroup))
  colnames(cat.counts)<-c("hcgroup","count")
  cat.counts$trt.time<-"dry.eve"
  cat.counts$count<-cat.counts$count/length(tstat.type)
  cat.counts$PACID<-10
  
  line.plot<-ggplot(t.melt, aes(x=comparison, y=t.stat, group=PACID, color=hcgroup, fill=hcgroup))+
    scale_color_manual(values=rainbow(n.clusters))+
    geom_line(y=0, size= 1, col="black")+
    geom_line(size=.5)+
    facet_wrap(~hcgroup)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90),legend.position="none")+
    geom_rect(col="lightgrey",fill="lightgrey", xmax=1.5, xmin=0, ymax=1,ymin=-1)+
    geom_text(aes(x=trt.time, label=count), color="black", x=1, y=0,  data=cat.counts)+
    ggtitle(plot.title)  
  #make the plot window
  if(pdf.file!="in.window"){
    pdf(pdf.file) 
  }
  par(mar=c(0,0,0,0), no.readonly=TRUE)
  plot.new() 
  gl <- grid.layout(nrow=2, ncol=1, heights=c(64,256))
  vp.1 <- viewport(layout.pos.col=1, layout.pos.row=1) # dend + profiles
  vp.2 <- viewport(layout.pos.col=1, layout.pos.row=2) # depth functions
  pushViewport(viewport(layout=gl))
  #plot 1
  pushViewport(vp.1)
  par(new=TRUE, fig=gridFIG())
  plot(d.out)
  popViewport()
  #plot2
  par(new=TRUE, fig=gridFIG())
  pushViewport(vp.2)
  print(line.plot, newpage=FALSE)
  popViewport(1)
  if(pdf.file!="in.window"){
    dev.off()
  }
}