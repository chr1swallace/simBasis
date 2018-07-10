#!/home/cew54/localc/bin/Rscript

                                        #!/usr/bin/env Rscript
## simGWAS helper functions

library(optparse)
library(magrittr)


## atm set to true untill fully validated
TEST<-FALSE

option_list = list(
        make_option(c("-s", "--scenario_file"), type="character", default=NULL,
            help="scenario file to use", metavar="character"),
        make_option(c("-o", "--out_dir"), type="character", default=NULL,
                help="where to put plots and source data generated", metavar="character"),
        make_option(c("-b", "--ld_source"), type="character", default='jp_ld',
              help="ld block source to use one of jp_ld or hm_ld",metavar="character"),
        make_option(c("-n", "--nsim"), type="numeric", default=100,
              help="Number of simulations to run",metavar="numeric"),
        make_option(c("-p", "--prefix"), action="store_true", default=FALSE,
              help="prefix to add to output files, automatically turns off plotting",metavar="numeric")
)

library(randomFunctions)
## if(!TEST){
##   opt_parser = OptionParser(option_list=option_list);
##   args = parse_args(opt_parser)
## }else{
  args <- getArgs(default=list(
      out_dir="/home/cew54/Projects/simBasis/output",
      #ldblock='77',
      #scenario_file="/home/ob219/rds/hpc-work/simBasis/support/scenarios/scenario1.yml",
      scenario_file="/rds/project/cew54/rds-cew54-wallace-share/Projects/simBasis/scenarios/jp_ld.backbone_share.yml",
      #scenario_file="/rds/project/cew54/rds-cew54-wallace-share/Projects/simBasis/scenarios/scenario2.yml",
      #out_dir="/home/ob219/rds/hpc-work/simBasis/simulations/",
      ld_source='hm_ld',
      nsim=200,
      prefix=FALSE
  ))
## }

print(args)

PLOT <- FALSE
ofile <- gsub(".yml","",basename(args$scenario_file))
if(args$prefix){
  PLOT <- FALSE
  PREFIX <- sample(letters,6) %>% paste(.,collapse="")
  ofile <- paste(PREFIX,ofile,sep='_')
}

if(file.exists(args$out_dir)){
  ofile <- file.path(args$out_dir,ofile)
}else{
  stop(sprintf("out_dir %s does not exist",args$out_dir))
}


library(data.table)
library(devtools)
## install_github("chr1swallace/simGWAS")
## install_github("ollyburren/cupcake")
library(simGWAS)
library(cupcake)
library(mvtnorm)
library(yaml)
library(ggplot2)
library(cowplot)
library(latex2exp)



loader <- function(f) {
    (load(f))
    ## (load(file=paste0(ofile,".RData"))) #"~/Projects/simBasis/100sims-v3.RData")
table(all.pd$metric)

all.pd[,mtype:=sub("_.*","",label1)]
all.pd[,qtype:=sub("_.*","",label2)]
all.pd[,grp:="non"]
with(all.pd, table(grp))
all.pd[sub("R","",mtype) < sub("R","",qtype), c("mtype","qtype") := list(qtype,mtype)]
all.pd[,grp:=gsub("R","",paste(mtype,qtype,sep="-"))]
all.pd[label1==label2 | label1==paste0("R",label2) | label1==paste0("R",label2),grp:=paste0(grp,"-exact")]
all.pd[qtype==paste0("R",mtype) | mtype==paste0("R",qtype) | mtype==qtype, grp:="b_similar"]
all.pd[label1==paste0("R",label2) | label1==paste0("R",label2), grp:="a_exact"]
all.pd[(grepl("GWAS",label1) & grepl("share",label2)) |
       (grepl("GWAS",label2) & grepl("share",label1)), grp:="50%"]
all.pd[grepl("random",label1) & grepl("GWAS",label2), grp:="random-GWAS"]
all.pd[grepl("random",label2) & grepl("GWAS",label1), grp:="random-GWAS"]
all.pd[grepl("random",label1) & grepl("share",label2), grp:="random-share"]
all.pd[grepl("random",label2) & grepl("share",label1), grp:="random-share"]
tmp <- all.pd[,.(n=.N),by=c("mtype","qtype","grp")]
tmp[order(grp,mtype,qtype),]

## all.pd <- all.pd[grep("beta_ws",metric),]
all.pd[,sd:=sd(d),by=c("metric")]
 all.pd[,stat:="old"]
all.pd[grepl("ws",metric),stat:="new"]
all.pd[,gamma:="gamma"]
all.pd[grepl("nog",metric),gamma:="nog"]
    all.pd[,size:=paste(sub(".*_","",label1),sub(".*_","",label2),sep="/")]
    copy(all.pd)
}

plot1 <- function(all.pd,facets=stat ~ size) {
s <- all.pd[label1!=label2 & grp=="a_exact" & grepl("GWAS10",label1) & grepl("GWAS10",label2),
            .(emn=mean(d),esd=sd(d)),by=c("metric","cond")]

x <- merge(all.pd[label1!=label2 & grepl("shrinkage",metric),],s,by=c("metric","cond"))
x[,size:=paste(sub(".*_","",label1),sub(".*_","",label2),sep="/")]
x[,Z:=(d)/esd] #/esd]

unique(x[,.(mtype,qtype,grp)])[grp=="a_exact",]

x[grp=="a_exact" ,mean(Z),by=c("mtype","qtype")]
x[grp=="a_exact" ,sd(Z),by=c("mtype","qtype")]

by <- unique(c("grp","stat","gamma","size","cond",as.character(facets)[-1]))
xm <- x[,.(mnZ=median(Z)),by=by]

ggplot(x[size=="5000/5000",], aes(x=grp,y=Z,col=gamma)) +
  geom_violin(draw_quantiles=c(0.25,.5,.75)) + #geom_point() +
  geom_hline(aes(yintercept=mnZ,col=gamma), data=xm[grp=="a_exact" & size=="5000/5000",],linetype="dashed") +
  theme(axis.text.x=element_text(angle=90)) + facet_grid(facets) + #,scales="free") + 
  background_grid() 
}

## limit to a few pairs of scenarios to see more sample sizes
plot2 <- function(all.pd,facets=stat~size) {
s <- all.pd[label1!=label2 & grp=="a_exact" & grepl("GWAS10",label1) & grepl("GWAS10",label2),
            .(emn=mean(d),esd=sd(d)),by=c("metric","cond")]
x <- merge(all.pd[label1!=label2 & grepl("shrinkage",metric),],s,by=c("metric","cond"))
x[,Z:=(d)/esd] #/esd]
x <- x[grp %in% c("50%","a_exact","b_similar","GWAS10-back","random-GWAS") &
       sub(".*_","",label1)==sub(".*_","",label2)]
by <- unique(c("grp","stat","gamma","size","cond",as.character(facets)[-1]))
xm <- x[,.(mnZ=median(Z)),by=by]
ggplot(x, aes(x=grp,y=Z,col=gamma)) +
  geom_violin(draw_quantiles=c(0.25,.5,.75)) + #geom_point() +
  geom_hline(aes(yintercept=mnZ,col=gamma), data=xm[grp=="a_exact" ,],linetype="dashed") +
  theme(axis.text.x=element_text(angle=90)) + facet_grid(facets) + #,scales="free") + 
  background_grid() 
}

d <- dirname(args$scenario_file)
files <- list.files(d,pattern="yml") %>% sub(".yml",".RData",.)
setwd("~/Projects/simBasis")
files
(labs <- gsub("jp_ld.|_.*","",files))

X12 <- lapply(files, function(f) loader(file.path("output_1.2",f)))
for(i in seq_along(X12))
    X12[[i]][,cond:=paste0(labs[i],"-1.2")]
X14 <- lapply(files, function(f) loader(file.path("output_1.44",f)))
for(i in seq_along(X14))
    X14[[i]][,cond:=paste0(labs[i],"-1.44")]




X <- X12
plot1(X[[1]]) # this plot shows new > old, if we want 50% > exact
plot1(X[[2]]) # this plot shows gamma can matter, if we want 50% > exact in new shrinkage
plot1(X[[3]])# this plot shows gamma can matter, if we want random-GWAS > exact in new shrinkage

X <- X14
plot1(X[[1]]) # this plot shows new > old, if we want 50% > exact
plot1(X[[2]]) # this plot shows gamma can matter, if we want 50% > exact in new shrinkage
plot1(X[[3]])# this plot shows gamma can matter, if we want random-GWAS > exact in new shrinkage

Z <- rbindlist(c(X12,X14))
Z <- Z[stat=="new",]
Z[,lab:=sub("-.*","",cond)]
Z[,or:=sub(".*-","",cond)]
plot1(Z,facets=or ~ lab)

plot12 <- lapply(X12,plot1)
plot14 <- lapply(X14,plot1)

library(cowplot)
plot_grid(plotlist=c(plot12,plot14))

################################################################################

## junkyard

if(!interactive())
    q()

## d <- all.pd[xor(grepl("^R",label2),grepl("^R",label1)),.(mn=mean(d),sd=mean(sd)),
##             by=c("label1","label2","metric","grp","qtype","mtype")]
s <- all.pd[label1!=label2 & grp=="exact",.(emn=mean(d),esd=sd(d)),by="metric"]
d <- all.pd[label1!=label2,.(mn=mean(d)),
            by=c("label1","label2","metric","grp","qtype","mtype")]
d <- merge(d,s,by="metric")
d[,Z:=(mn - emn)/esd]
with(d, table(grp))
with(all.pd[label1!=label2,], table(grp))
d[,size:=paste(sub(".*_","",label1),sub(".*_","",label2),sep="/")]
head(d)

d[,stat:="old"]
d[grepl("ws",metric),stat:="new"]
d[,gamma:="gamma"]
d[grepl("nog",metric),gamma:="nog"]

## don't use z - simrep can be large
## ggplot(d[grepl("^z",metric),], aes(x=metric,y=Z,col=grp)) + geom_violin() + geom_point() + background_grid()

## beta gets at biology better, not misled by sample size - simrep and rep v v similar
ggplot(d[grepl("^beta.*shrinkage",metric),], aes(x=grp,y=Z,col=gamma)) + geom_violin(draw_quantiles=c(0.25,.5,.75)) + geom_point() +
  theme(axis.text.x=element_text(angle=90)) + facet_grid(stat~size) + 
  background_grid() 

ggplot(d[grepl("^beta.*shrinkage",metric),], aes(x=grp,y=Z,col=grp)) + geom_violin(draw_quantiles=c(0.25,.5,.75)) + #geom_point() +
  theme(axis.text.x=element_blank()) +
  background_grid() + facet_grid(.~metric)

dc <- dcast(d, metric + label1 ~ label2, value.var="Z")


toplot <- all.pd[qtype==paste0("R",mtype) | mtype==paste0("R",qtype),]
toplot <- all.pd[qcase==5000 & mcase==5000, ] # restrict to 5000_5000 to plot more easily

library(ggridges)
ggplot(toplot,aes(x=d,y=label1)) + geom_density_ridges() + facet_wrap( ~ label2)
ggplot(toplot,aes(x=d,y=label1)) + geom_density_ridges() + facet_grid(metric ~ label2)
ggplot(toplot,aes(x=d,y=qtype)) + geom_density_ridges() + facet_grid(mtype ~ metric,scales="free")

ggplot(toplot[metric=="z_emp_shrinkage",],
       aes(x=d,y=label1)) + geom_density_ridges() + facet_wrap( ~ label2,scales="free")

if(PLOT){
sum.all.pd <- all.pd[,list(mean.d=median(d)),by=c('label1','label2','metric')]

lab.lev <- c('GWAS10_500_3000','GWAS10_2000_2000',
'GWAS10_5000_5000','share_500_3000','share_2000_2000',
'share_5000_5000','random_500_3000','random_2000_2000','random_5000_5000')

lab.lev2 <- c('ID1','ID2',
'ID3','S1','S2',
'S3','R1','R2','R3')

sum.all.pd[,label1:=factor(label1,levels=lab.lev)]
sum.all.pd[,label2:=factor(label2,levels=lab.lev)]
levels(sum.all.pd$label1)<-lab.lev2
levels(sum.all.pd$label2)<-lab.lev2


all.plots <- lapply(split(sum.all.pd,sum.all.pd$metric),function(dat){
  ggplot(data = dat, aes(x=label1, y=label2, fill=mean.d)) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0,vjust=0.5)) +
    geom_tile() + ggtitle(unique(dat$metric)) +
    xlab("Scenario") + ylab("Scenario")
})

ppd <- plot_grid(plotlist=all.plots[c('z','beta','beta_r_emp_maf_se','beta_emp_shrinkage','beta_shrinkage_nog','beta_ws_emp_shrinkage','beta_ws_shrinkage_nog')])
save_plot(paste(out.file.stub,'pw_dist_beta','pdf',sep="."),ppd,base_height=10,base_width=15,base_aspect=1)

ppd <- plot_grid(plotlist=all.plots[c('z','beta','beta_r_emp_maf_se','z_emp_shrinkage','z_shrinkage_nog','z_ws_emp_shrinkage','z_ws_shrinkage_nog')])
save_plot(paste(out.file.stub,'pw_dist_z','pdf',sep="."),ppd,base_height=10,base_width=15,base_aspect=1)
}
saveRDS(all.pd,file=paste(out.file.stub,'pw_dist','RDS',sep="."))

getDist <- function(S,ref='control'){
  R <- S$basis$x[ref,]
  tDT <- data.table(reference=ref,t(apply(S$proj,1,function(x) sqrt(sum((x - R)^2)))))
}

## compute dist matrix for all simulations considered
dist.res <- lapply(DT.sims,function(db){
  ts <- createBasisAndProj(db,basis.sims,proj.sims)
  ctres <- lapply(ts,getDist) %>% rbindlist
  ctres[,metric:=names(ts)]
  ares <- lapply(ts,getDist,'GWAS10') %>% rbindlist
  ares[,metric:=names(ts)]
  rbind(ctres,ares)
})


all.distances <- rbindlist(dist.res)

mall <- melt(all.distances,id.vars=c('metric','reference'))

#mall[,scale.value:=scale(value),by=c('reference','metric')]

mall[,variable:=factor(variable,levels=c('GWAS10_500_3000','GWAS10_2000_2000',
'GWAS10_5000_5000','share_500_3000','share_2000_2000',
'share_5000_5000','random_500_3000','random_2000_2000','random_5000_5000'))]



#mall[metric=='beta',metric:=latex2exp("$\beta")]

metric_names <- c(
TeX('$\\hat{\\beta}$'),
  TeX('$\\hat{\\gamma}$'),
  TeX('$Z$'),
  #TeX('$\\hat{\\gamma}_{SS}$'),
  #TeX('$\\hat{\\gamma}_{MAF}$'),
  TeX('Method 1'),
  #TeX('Method1$_{SS}$'),
  #TeX('Method1$_{MAF}$'),
  TeX('Method 2')
  #TeX('Method2$_{SS}$'),
  #TeX('Method2$_{MAF}$')
)

mall[,metric:=factor(metric,levels=c('beta','r_emp_maf_se','z','emp_shrinkage','ws_emp_shrinkage'),labels=metric_names)]
mall[,variable:=factor(gsub("_"," ",as.character(variable)),levels=gsub("_"," ",levels(variable)))]

if(PLOT){
  pp <- ggplot(mall,aes(x=variable,y=value,col=reference)) + geom_boxplot() +
  facet_wrap( ~ metric,scale="free",ncol=3,labeller = label_parsed) +
  xlab("Scenario") + ylab("Distance") +
  theme(axis.text.x=element_text(angle = -90, hjust = 0,vjust=0.5),
  strip.background =element_rect(fill="grey95")) +
  scale_colour_manual("Origin",values=c("control"="steelblue1","GWAS10"="firebrick1"),labels=c('Control','GWAS10')) +
  background_grid(major = "xy", minor = "none")
  save_plot(paste(out.file.stub,'pdf',sep="."),pp,base_height=10,base_width=10)
  #save_plot(paste(out.file.stub,'scale','pdf',sep="."),pb,base_height=10,base_width=10)
}
saveRDS(mall,file=paste(out.file.stub,'RDS',sep="."))

## code for running on the Q
if(FALSE){
  SCEN.DIR <- '/rds/project/cew54/rds-cew54-wallace-share/Projects/simBasis/scenarios'
  TOTAL <- 1000
  BATCH <- 100
  OUT_DIR <- '/home/ob219/tmp/test_scen/'

  n <- TOTAL/BATCH

  cmd <- sapply(list.files(path=SCEN.DIR,pattern="*.yml",full.names=TRUE),function(f){
    sapply(1:n,function(i){
      sprintf("Rscript /home/ob219/git/simBasis/R/simFullGWAS.R -s %s -o %s -n %d -p",f,OUT_DIR,BATCH)
    })
  })

  write(cmd,file="/home/ob219/git/simBasis/sh/scen.txt")
}



## BELOW HERE IS EXPERIMENTAL CODE FOR malhalanobis distance stuff.


## I think that this is equivalent to the malhalanobis distance
## here we multiply the loading for a given axis by the variance explained by that

if(FALSE){

mahalanobis_pairwise<-function(pc,proj){
  ## add in zero for mahalanobis distance
  load <- pc$x
  den <- pc$sdev^2
  if(exists("proj"))
    load<-rbind(load,proj)
  ## pairwise distance function
  #dist<-apply(load,1,function(p){
  #  apply(pc$x,1,function(x) (x-p)^2/pc$sdev^2) %>% colSums %>% sqrt
  #})
  ## equivalent and faster
  ## note we can't use scale if we add in the projections - we assume centred
  ## to standardise divide through by the sdev
  stand <- apply(load,1,function(x) x/pc$sdev) %>% t
  as.matrix(dist(stand))
  #mahalanobis
  #man <- apply(pc$x,1,function(x) x^2/den) %>% colSums
}


## not sure that this is what we want as it massively inflates distances for
## components that explain a tiny amount of variance but I need to read more about it
## see code in mahalanobis.R for an explanation
mahalanobis_pairwise(sim1$emp_shrinkage$basis,sim1$emp_shrinkage$proj)

}
