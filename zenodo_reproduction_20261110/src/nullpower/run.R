suppressPackageStartupMessages({library(Matrix);library(data.table);library(SingleCellExperiment);library(glmGamPoi);library(edgeR);library(Seurat);library(nebula);library(fixest);library(dplyr);library(tidyr);library(limma);library(Rcpp);library(rhdf5);library(filelock)})
args<-commandArgs(trailingOnly=TRUE);if(length(args)!=2L||!args[1]%in%c('yazar','reichart','bca')||!args[2]%in%c('smoke','smoke-cell','full'))stop('Usage: run.R yazar|reichart|bca smoke|smoke-cell|full')
name<-args[1];mode<-args[2];options(degpaper.nebula.ncore=1L,nullpower.fixef_none=name!='bca');source('src/nullpower/notebook-methods.R')
root<-file.path('data/nullpower-reproduction-20261104',name);prepared_path<-file.path(root,'prepared.rds');meta<-readRDS(prepared_path);fi<-file.info(meta$input);stopifnot(identical(as.numeric(fi$size),as.numeric(meta$input_info$size)),identical(as.numeric(fi$mtime),as.numeric(meta$input_info$mtime)),identical(unname(meta$prepare_hash),unname(tools::md5sum('src/nullpower/prepare.R'))))
out<-file.path(root,mode);dir.create(out,recursive=TRUE,showWarnings=FALSE);lock<-filelock::lock(file.path(out,'writer.lock'),timeout=0);if(is.null(lock))stop('Another writer holds this output lock.');on.exit(filelock::unlock(lock),add=TRUE)
scripts<-c('src/nullpower/prepare.R','src/nullpower/notebook-methods.R','src/nullpower/run.R');pkgs<-c('Matrix','data.table','SingleCellExperiment','glmGamPoi','edgeR','Seurat','nebula','fixest','dplyr','tidyr','limma','Rcpp','rhdf5','filelock')
manifest<-list(input=meta$input,input_info=meta$input_info,prepared_md5=unname(tools::md5sum(prepared_path)),script_md5=tools::md5sum(scripts),packages=sapply(pkgs,function(x)as.character(packageVersion(x))),R=R.version.string,seed=20261104L,ncore=1L,mode=mode,top=meta$top,eligible=meta$eligible)
manifest_path<-file.path(out,'manifest.rds');if(file.exists(manifest_path)){if(!identical(readRDS(manifest_path),manifest))stop('Incompatible output manifest; refusing to mix results.')}else{if(length(setdiff(list.files(out),'writer.lock')))stop('Output directory exists without a manifest.');saveRDS(manifest,manifest_path)}
getcounts<-function(cells,genes){
 cells<-sort(unique(cells));stopifnot(length(cells)>0);result<-list()
 for(block in unique((cells-1L)%/%10000L)){
  selected<-cells[(cells-1L)%/%10000L==block];start<-block*10000L+1L;end<-min(length(meta$donor),start+9999L);p<-meta$ptr[start:(end+1L)];lo<-p[1]+1;nn<-p[length(p)]-p[1]
  ix<-as.integer(h5read(meta$input,paste0(meta$counts_path,'/indices'),start=lo,count=nn))+1L;v<-as.numeric(h5read(meta$input,paste0(meta$counts_path,'/data'),start=lo,count=nn));col<-rep(seq_len(end-start+1L),diff(p));gi<-match(ix,genes);ci<-match(col,selected-start+1L);ok<-!is.na(gi)&!is.na(ci)
  result[[length(result)+1L]]<-sparseMatrix(i=gi[ok],j=ci[ok],x=v[ok],dims=c(length(genes),length(selected)))
 }
 mat<-do.call(cbind,result);colnames(mat)<-as.character(cells);rownames(mat)<-meta$feature[genes];mat
}
labels<-c('glmGamPoi (Pb)','edgeR (Pb)','limma (Pb)','glmGamPoi (cell)','robust GLM (cell)','NB GLMM')
for(t in if(mode!='full')1L else 1:6)for(n in if(mode!='full')10L else c(10L,50L))for(scenario in if(mode=='smoke')'subject' else if(mode=='smoke-cell')'cell' else c('cell','subject')){
 seed<-20261104L+match(name,c('yazar','reichart','bca'))*100000L+t*1000L+n*10L+match(scenario,c('cell','subject'));set.seed(seed)
 allowed<-which(meta$mean[,t]>0.1 & meta$mean[,7]>0.01);ngenes<-if(mode!='full')20L else 100L;stopifnot(length(allowed)>=ngenes);genes<-sample(allowed,ngenes);pool<-unique(meta$donor[meta$ct==meta$top[t] & meta$donor%in%meta$eligible]);stopifnot(length(pool)>0)
 for(iter in seq_len(if(mode!='full')1L else 10L)){
  dest<-file.path(out,sprintf('ct%d_n%d_%s_iter%02d.rds',t,n,scenario,iter));if(file.exists(dest)){old<-readRDS(dest);if(!identical(old$genes,genes)||old$seed!=seed)stop('Checkpoint gene or seed mismatch');assign('.Random.seed',old$rng_after,envir=.GlobalEnv);cat('Reusing',dest,'\n');next}
  donors<-sample(pool,n,replace=TRUE);selected<-lapply(donors,function(d){cells<-which(meta$ct==meta$top[t]&meta$donor==d);if(mode!='full')cells<-head(cells,30L);cells[as.logical(rbinom(length(cells),1,0.5))]});ids<-rep(seq_len(n),lengths(selected));cells<-unlist(selected,use.names=FALSE);stopifnot(length(cells)>0)
  raw<-getcounts(cells,genes);m<-raw[,match(cells,as.integer(colnames(raw))),drop=FALSE];df<-data.frame(id=ids,donor_id=donors[ids]);colnames(m)<-rownames(df)
  df$tx_cell<-rbinom(nrow(df),1,0.5);if(scenario=='subject'){urn<-c(rep(1L,as.integer(n/2)),rep(0L,n-as.integer(n/2)));tx<-sample(urn,n,replace=FALSE);df$tx_cell<-rep(tx,times=lengths(selected))}
  prior_null<-m[-(1:5),,drop=FALSE];prior_x<-m@x+0;downCells(m,1:5,which(df$tx_cell==1),0.5);if(!identical(as.numeric(prior_null@x),as.numeric(m[-(1:5),,drop=FALSE]@x)))stop('Null gene values changed during downsampling')
  if(mode!='full' && identical(prior_x,as.numeric(m@x))) stop('Downsampling had no detectable effect in smoke test.')
  ans<-matrix(NA_real_,ngenes,6,dimnames=list(meta$feature[genes],labels));errors<-character(6);for(k in seq_along(list.func)){
   cat(name,mode,'type',t,'n',n,scenario,'iteration',iter,'method',labels[k],'cells',nrow(df),'\n');flush.console()
   z<-tryCatch(list.func[[k]](if(k==5)if(scenario=='cell') y~tx_cell|id else y~tx_cell else ~tx_cell,if(k==5)t(m) else m,df),error=function(e){errors[k]<<-conditionMessage(e);NULL})
   if(!is.null(z)){vals<-as.numeric(z[,4]);if(length(vals)==ngenes)ans[,k]<-vals else {rn<-rownames(z);ix<-match(rn,rownames(ans));if(!is.null(rn)&&!anyNA(ix)&&!anyDuplicated(ix))ans[ix,k]<-vals else errors[k]<-paste('Returned',length(vals),'rows; gene alignment unverified')}}
  }
  item<-list(p=ans,genes=genes,donors=donors,cell_indices=cells,group=df$tx_cell,seed=seed,rng_after=.Random.seed,errors=setNames(errors,labels),scenario=scenario,cell_type=meta$top[t]);tmp<-tempfile(pattern='checkpoint-',tmpdir=out);saveRDS(item,tmp);if(!file.rename(tmp,dest))stop('Checkpoint commit failed');cat('Saved',dest,'errors',sum(nzchar(errors)),'\n');flush.console()
 }
}
