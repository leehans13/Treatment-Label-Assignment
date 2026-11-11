suppressPackageStartupMessages({library(data.table)})
a<-commandArgs(trailingOnly=TRUE);stopifnot(length(a)==1,a[1]%in%c('yazar','reichart','bca'));name<-a[1];root<-file.path('data/nullpower-reproduction-20261104',name);meta<-readRDS(file.path(root,'prepared.rds'));p<-file.path(root,'full');man<-readRDS(file.path(p,'manifest.rds'));stopifnot(identical(unname(man$script_md5),unname(tools::md5sum(names(man$script_md5)))),identical(man$prepared_md5,unname(tools::md5sum(file.path(root,'prepared.rds')))))
expected<-unlist(lapply(1:6,function(t)unlist(lapply(c(10L,50L),function(n)unlist(lapply(c('cell','subject'),function(s)sprintf('ct%d_n%d_%s_iter%02d.rds',t,n,s,1:10)))))));paths<-file.path(p,expected);if(!all(file.exists(paths)))stop('Missing checkpoints: ',sum(!file.exists(paths)))
labels<-c('glmGamPoi (Pb)','edgeR (Pb)','limma (Pb)','glmGamPoi (cell)','robust GLM (cell)','NB GLMM')
for(t in 1:6)for(n in c(10L,50L))for(s in c('cell','subject')){
 files<-file.path(p,sprintf('ct%d_n%d_%s_iter%02d.rds',t,n,s,1:10));items<-lapply(files,readRDS)
 for(x in items) {
  stopifnot(identical(dim(x$p),c(100L,6L)),identical(colnames(x$p),labels),identical(rownames(x$p),meta$feature[x$genes]),length(unique(x$genes))==100L,identical(x$genes,items[[1]]$genes),length(x$donors)==n,length(x$cell_indices)==length(x$group),all(x$group %in% c(0,1)),identical(x$scenario,s),identical(x$cell_type,meta$top[t]),length(x$rng_after)>0L)
  finite<-x$p[is.finite(x$p)];stopifnot(all(finite>=0 & finite<=1+1e-12))
 }
}
out<-file.path(root,'outputs');if(dir.exists(out)&&length(list.files(out,all.files=FALSE)))stop('Output directory not empty; refusing to overwrite.');dir.create(out,recursive=TRUE,showWarnings=FALSE)
stats<-list();errors<-list();bounds<-list()
for(t in 1:6)for(n in c(10L,50L))for(s in c('cell','subject')){
 files<-file.path(p,sprintf('ct%d_n%d_%s_iter%02d.rds',t,n,s,1:10));items<-lapply(files,readRDS)
 for(i in seq_along(items)){e<-items[[i]]$errors;for(k in which(nzchar(e)))errors[[length(errors)+1L]]<-data.frame(celltype=meta$top[t],n=n,scenario=s,iteration=i,method=labels[k],error=e[k]);z<-items[[i]]$p;bad<-which(is.finite(z)&(z<0|z>1),arr.ind=TRUE);if(nrow(bad))for(j in seq_len(nrow(bad)))bounds[[length(bounds)+1L]]<-data.frame(celltype=meta$top[t],n=n,scenario=s,iteration=i,gene=rownames(z)[bad[j,1]],method=colnames(z)[bad[j,2]],p_exact=sprintf('%.17g',z[bad[j,1],bad[j,2]]))}
 for(kind in c('null','pow')){
  rows<-if(kind=='null')6:100 else 1:5;m<-do.call(rbind,lapply(items,function(x)x$p[rows,,drop=FALSE]));d<-file.path(out,paste0(kind,'_',s));dir.create(d,showWarnings=FALSE);dest<-file.path(d,paste(name,'n',n,'ct',t,'fc','0.5','csv',sep='.'));tmp<-tempfile(tmpdir=d);write.csv(m,tmp);if(!file.rename(tmp,dest))stop('CSV commit failed')
  for(k in 1:6){v<-m[,k];valid<-is.finite(v);stats[[length(stats)+1L]]<-data.frame(dataset=name,celltype=meta$top[t],n=n,scenario=s,kind=kind,method=labels[k],total=length(v),finite=sum(valid),nonfinite=sum(!valid),fraction_p_lt_0.05=if(any(valid))mean(v[valid]<0.05) else NA_real_)

  }
 }
}
write.csv(rbindlist(stats),file.path(out,'summary.csv'),row.names=FALSE);err<-if(length(errors))rbindlist(errors) else data.frame(celltype=character(),n=integer(),scenario=character(),iteration=integer(),method=character(),error=character());write.csv(err,file.path(out,'model-errors.csv'),row.names=FALSE);bound<-if(length(bounds))rbindlist(bounds) else data.frame(celltype=character(),n=integer(),scenario=character(),iteration=integer(),gene=character(),method=character(),p_exact=character());write.csv(bound,file.path(out,'out-of-range-p.csv'),row.names=FALSE)
cat('Completed 240 checkpoints and 48 result CSVs; model error records:',nrow(err),'p-value bound deviations:',nrow(bound),'\n')
