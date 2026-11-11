suppressPackageStartupMessages({library(rhdf5);library(Matrix)})
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1,args[1]%in%c('yazar','reichart','bca'));name<-args[1]
input<-switch(name,yazar='data/cellxgene/yazar.h5ad',reichart='data/cellxgene/reichart-all.h5ad',bca='data/cellxgene/reed-global.h5ad');counts<-if(name=='yazar') 'X' else 'raw/X'
root<-file.path('data/nullpower-reproduction-20261104',name);dir.create(root,recursive=TRUE,showWarnings=FALSE);dest<-file.path(root,'prepared.rds');if(file.exists(dest))stop('Prepared file already exists; refusing to overwrite.')
readcat<-function(key){v<-h5read(input,paste0('obs/',key,'/categories'));i<-h5read(input,paste0('obs/',key,'/codes'));out<-rep(NA_character_,length(i));out[i>=0]<-v[i[i>=0]+1L];out}
ct<-readcat('cell_type');donor<-readcat('donor_id');freq<-sort(table(ct),decreasing=TRUE);top<-names(freq)[freq>=freq[6]];stopifnot(length(top)==6L,!anyNA(ct),!anyNA(donor))
tab<-table(donor,factor(ct,levels=top));cutoff<-if(name=='bca')100L else 50L;eligible<-rownames(tab)[rowSums(tab>cutoff)==6];shape<-h5readAttributes(input,paste0('/',counts))$shape;ptr<-h5read(input,paste0(counts,'/indptr'));N<-length(ct);ng<-as.integer(shape[2]);total<-matrix(0,ng,7);nums<-numeric(7);block<-10000L
cat(name,'cells',N,'genes',ng,'donors',length(unique(donor)),'eligible',length(eligible),'\n');flush.console()
for(start in seq.int(1L,N,by=block)){
 end<-min(N,start+block-1L);p<-ptr[start:(end+1L)];lo<-p[1]+1;nn<-p[length(p)]-p[1];if(nn==0)next
 ix<-as.integer(h5read(input,paste0(counts,'/indices'),start=lo,count=nn))+1L;v<-as.numeric(h5read(input,paste0(counts,'/data'),start=lo,count=nn));if(any(!is.finite(v))||any(v<0)||any(v!=floor(v)))stop('Invalid count in block ',start)
 mat<-sparseMatrix(i=ix,j=rep(seq_len(end-start+1L),diff(p)),x=v,dims=c(ng,end-start+1L));selected<-donor[start:end]%in%eligible
 total[,7]<-total[,7]+Matrix::rowSums(mat);nums[7]<-nums[7]+ncol(mat)
 for(j in 1:6){z<-selected & ct[start:end]==top[j];if(any(z)){total[,j]<-total[,j]+Matrix::rowSums(mat[,z,drop=FALSE]);nums[j]<-nums[j]+sum(z)}}
 if(start%%100000L==1L){cat(name,'processed',end,'/',N,'\n');flush.console()}
}
stopifnot(all(nums>0));feature<-h5read(input,if(name=='yazar')'var/_index' else 'raw/var/_index')
prepared<-list(name=name,input=normalizePath(input),input_info=file.info(input)[,c('size','mtime')],counts_path=counts,ct=ct,donor=donor,top=top,eligible=eligible,ptr=ptr,shape=shape,feature=feature,mean=sweep(total,2,nums,'/'),n=nums,cutoff=cutoff,prepare_hash=tools::md5sum('src/nullpower/prepare.R'))
tmp<-tempfile(tmpdir=root);saveRDS(prepared,tmp);if(!file.rename(tmp,dest))stop('Could not commit preparation file');cat('Prepared:',dest,'\n')
