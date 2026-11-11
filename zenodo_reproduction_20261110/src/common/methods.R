# Model implementations from perturb/kang.ipynb.
selectCol <- function(mat, j.col){    
    x.col.dense <- rep(0,nrow(mat))
    p.begin <- mat@p[j.col]+1
    p.end <- mat@p[j.col+1]
    i.col <- mat@i[p.begin:p.end]+1 # i counts from 0
    x.col <- mat@x[p.begin:p.end]
    x.col.dense[i.col] <- x.col
    return(x.col.dense)
    }

selectCols <- function(mat, j.cols){    
    return(sapply(j.cols, selectCol, mat=mat))
    }

fixest.mult <- function(formula, count, df){
    df.result <- data.frame(matrix(nrow=0, ncol=4))
    colnames(df.result) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)')
    for (j in 1:ncol(count)){
        df$y <- selectCol(count, j)
        fit <- fixest::fepois(formula, vcov='hetero', data=df)
        df.result[colnames(count)[j],] <- fixest::coeftable(fit)['tx_cell',] # fixed effect o/x 에 따라 다르게 들어가야함
        } 
    return(df.result)
    }

nebula.mult <- function(formula, count, df){
    pred <- model.matrix(formula, data=df)
    sid <- df$id
    fit.nebula <- nebula::nebula(
        count,
        sid,
        pred=pred,
        cpc=0,
        mincp=0,
        ncore=getOption("degpaper.nebula.ncore", 2L)
        )
    fit.result <- fit.nebula$summary
    rownames(fit.result) <- fit.result$gene
    return(
            fit.result %>%
            mutate(
                Estimate=logFC_tx_cell,
                'Std. Error'=se_tx_cell,
                't value'=logFC_tx_cell/se_tx_cell,
                'Pr(>|t|)'=p_tx_cell
                ) %>%
            select(Estimate, 'Std. Error', 't value', 'Pr(>|t|)')
        )
    }

glmgp.mult <- function(formula, count, df){
    sce.obj <- SingleCellExperiment::SingleCellExperiment(list(counts=count), colData=df)
    sce.pb <- glmGamPoi::pseudobulk(
        sce.obj,
        group_by=vars(id, tx_cell),
        verbose=FALSE
        )

    fit <- glmGamPoi::glm_gp(sce.pb, design=~1+tx_cell)
    test <- glmGamPoi::test_de(fit, reduced_design=~1)
    
    beta <- fit$Beta[,'tx_cell']
    pval <- test$pval
    tval <- qnorm(1-pval/2) * sign(beta)
    se <- beta/tval
    result <- cbind(beta, se, tval, pval)
    colnames(result) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)')
    return(result)
    }

glmgp.cell.mult <- function(formula, count, df){
    sce.obj <- SingleCellExperiment::SingleCellExperiment(list(counts=count), colData=df)
    fit <- glmGamPoi::glm_gp(sce.obj, design=~1+tx_cell, on_disk=FALSE, size_factors=FALSE)
    test <- glmGamPoi::test_de(fit, reduced_design=~1)
    
    beta <- fit$Beta[,'tx_cell']
    pval <- test$pval
    tval <- qnorm(1-pval/2) * sign(beta)
    se <- beta/tval
    result <- cbind(beta, se, tval, pval)
    colnames(result) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)')
    return(result)
    }

edger.mult <- function(formula, count, df){
    sce.obj <- SingleCellExperiment::SingleCellExperiment(list(counts=count), colData=df)
    sce.pb <- glmGamPoi::pseudobulk(
        sce.obj,
        group_by=vars(id, tx_cell),
        verbose=FALSE
        )

    design <- model.matrix(~1+tx_cell, data=colData(sce.pb))
    edger.obj <- edgeR::DGEList(counts(sce.pb))
    edger.obj <- edgeR::estimateDisp(edger.obj, design)
    fit <- edgeR::glmQLFit(y=edger.obj, design=design)
    test <- edgeR::glmTreat(fit, coef=2)

    beta <- test$coefficients[,'tx_cell']
    pval <- test$table[,'PValue']
    tval <- qnorm(1-pval/2) * sign(beta)
    se <- beta/tval

    result <- cbind(beta, se, tval, pval)
    colnames(result) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)')
    
    return(result)
    }

limma.mult <- function(formula, count, df){
    sce.obj <- SingleCellExperiment::SingleCellExperiment(list(counts=count), colData=df)
    sce.pb <- glmGamPoi::pseudobulk(
        sce.obj,
        group_by=vars(id, tx_cell),
        verbose=FALSE
        )

    design <- model.matrix(~1+tx_cell, data=colData(sce.pb))
    edger.obj <- edgeR::DGEList(counts(sce.pb))
    v <- limma::voom(edger.obj, design)
    vfit <- limma::lmFit(v, design)
    efit <- limma::eBayes(vfit)
    
    beta <- efit$coefficients[,'tx_cell'] * log(2)
    pval <- efit$p.value[,'tx_cell']
    tval <- qnorm(1-pval/2) * sign(beta)
    se <- beta/tval

    result <- cbind(beta, se, tval, pval)
    colnames(result) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)')
    
    return(result)
    }
