sigMeasuresProj <- function(results, annotation, contrastlist, corMethod = c('pearson', 'kendall', 'spearman'), resampling = c('none', 'permutation', 'bootstrap'), adjustMethod = 'BH', verbose = T, useall = T, uniqueonly = F, wilcox = F) {
    if(verbose) message("Finding significantly correlated measurements")
    corMethod <- match.arg(corMethod)
    resampling <- match.arg(resampling)
    relpcs <- sort(as.numeric(unique(unlist(lapply(results$pcs, rownames)))))
    relpcs.rle <- rle(sort(as.numeric(unlist(lapply(results$pcs, rownames)))))
    rl <- lapply(vector("list", length(unique(unlist(lapply(contrastlist, function(e) return(e[1])))))), function(e) character(0))
    names(rl) <- colnames(annotation)[colnames(annotation) %in% unlist(lapply(contrastlist, function(e) return(e[1])))]
    if(length(relpcs) == 0) return(rl)
    rlcount <- 1
    for(e in which(colnames(annotation) %in% unique(unlist(lapply(contrastlist, function(e) return(e[1])))))) {
        if(class(annotation[,e]) != "factor") {
            rv <- list("high-low" = character(0))
            if(nrow(results$pcs[[e]]) > 0) {
                glmnet.coefs <- rep(0, ncol(results$pcasp) + 1)
                cpcs <- as.numeric(rownames(results$pcs[[e]])) + 1
                glmnet.coefs[cpcs] <- results$pcs[[e]]
                subsp.vals <- cbind(as.numeric(annotation[, e]), results$pcasp)
                zero.pt <- colMeans(subsp.vals)[cpcs]
                neg1.pt <- glmnet.coefs * c(1, rep(-1, ncol(results$pcasp)))
                proj.pts <- t(apply(subsp.vals, 1, function(v) projectPoint(glmnet.coefs, neg1.pt, v)))[, cpcs, drop = F]
                proj.vals <- ifelse(as.matrix(dist(rbind(subsp.vals[which.max(subsp.vals[,1]), cpcs], proj.pts)))[-1, 1] < as.matrix(dist(rbind(subsp.vals[which.min(subsp.vals[,1]), cpcs], proj.pts)))[-1, 1], 1, -1) * as.matrix(dist(rbind(zero.pt, proj.pts)))[-1, 1]
                cor.vals <- foreach(j=1:ncol(results$pcaobj$call$X), .combine = rbind, .packages = c('stats'), .export = 'mcor.test') %dopar% mcor.test(proj.vals, results$pcaobj$call$X[,j], method = corMethod, resampling = resampling, wilcox = wilcox)
                cor.vals <- cbind(cor.vals, p.adjust(cor.vals[,2], method = adjustMethod))
                colnames(cor.vals) <- c(corMethod, "p", adjustMethod)
                rownames(cor.vals) <- colnames(results$pcaobj$call$X)
                rv[[1]] <- cor.vals
            }
            rl[[rlcount]] <- rv
        } else {
            grp.pairs <- gtools::combinations(length(levels(annotation[,e])), 2)
            contrasts <- apply(grp.pairs, 1, function(i) paste(levels(annotation[,e])[i[1]], levels(annotation[,e])[i[2]], sep='-'))
            curcons <- which(unlist(lapply(contrastlist, function(e) e[1])) == colnames(annotation)[e])
            inccons <- unlist(lapply(strsplit(contrasts, "-"), function(e) any(unlist(lapply(contrastlist[curcons], function(clcc) return(any(all(unlist(clcc)[-1] == e), all(unlist(clcc)[-1] == rev(e)))))))))
            grp.pairs <- grp.pairs[inccons, , drop = F]
            contrasts <- contrasts[inccons]
            rv <- lapply(vector("list", length(contrasts)), function(e) character(0))
            names(rv) <- contrasts
            crpcs <- results$pcs[[e]]
            if(uniqueonly) {
                crpcs <- crpcs[relpcs.rle$lengths[relpcs.rle$values %in% as.numeric(rownames(crpcs))] == 1, , drop = F]
            }
            if(nrow(crpcs) > 0) {
                vs <- matrix(0, nrow = 2, ncol = ncol(results$pcasp))
                cpcs <- as.numeric(rownames(crpcs))
                for(p in 1:nrow(grp.pairs)) {
                    vs[1, cpcs] <- crpcs[,grp.pairs[p, 1]]
                    vs[2, cpcs] <- crpcs[,grp.pairs[p, 2]]
                    zero.pt <- colMeans(vs)[cpcs]
                    proj.pts <- t(apply(results$pcasp, 1, function(v) projectPoint(vs[1,], vs[2,], v)))[, cpcs, drop = F]
                    cx <- results$pcaobj$call$X
                    if(!useall) {
                        proj.pts <- proj.pts[as.integer(annotation[,e]) %in% grp.pairs[p,], , drop = F]
                        cx <- cx[as.integer(annotation[,e]) %in% grp.pairs[p,], , drop = F]
                    }
                    proj.vals <- ifelse(as.matrix(dist(rbind(vs[1, cpcs], proj.pts)))[-1, 1] < as.matrix(dist(rbind(vs[2, cpcs], proj.pts)))[-1, 1], 1, -1) * as.matrix(dist(rbind(zero.pt, proj.pts)))[-1, 1]
                    if(wilcox) {
                        cor.vals <- t(apply(cx, 2, function(j) { cv <- wilcox.test(j[as.integer(annotation[,e]) %in% grp.pairs[p,1]], j[as.integer(annotation[,e]) %in% grp.pairs[p,2]], conf.int = T); return(c(cv$estimate, cv$p.value)) }))
                    } else {
                        cor.vals <- foreach(j=1:ncol(results$pcaobj$call$X), .combine = rbind, .packages = c('stats'), .export = 'mcor.test') %dopar% mcor.test(proj.vals, results$pcaobj$call$X[,j], method = corMethod, resampling = resampling)
                    }
                    cor.vals <- cbind(cor.vals, p.adjust(cor.vals[,2], method = adjustMethod))
                    colnames(cor.vals) <- c(corMethod, "p", adjustMethod)
                    rownames(cor.vals) <- colnames(results$pcaobj$call$X)
                    rv[[p]] <- cor.vals
                }
            }
            rl[[rlcount]] <- rv
        }
        rlcount <- rlcount + 1
    }
    return(rl)
}

mcor.test <- function(x, y, method, resampling) {
    rv <- cor.test(x, y, method = method)
    if(resampling == "bootstrap") {
        rv$p.value <- 1-ecdf(abs(na.omit(replicate(10000, suppressWarnings(cor(sample(x, replace = T), sample(y, replace = T), method = method))))))(abs(rv$estimate))
    } else if(resampling == "permutation") {
        rv$p.value <- 1-ecdf(abs(replicate(10000, cor(sample(x), sample(y), method = method))))(abs(rv$estimate))
    }
    rv <- c(rv$estimate, rv$p.value)
    return(rv)
}

# deprecated
sigMeasures <- function(results, annotation, ntop = 5000, alpha = .05, verbose = T) {
    if(verbose) message("Finding significantly correlated measurements")
    relpcs <- sort(as.numeric(unique(unlist(lapply(results$pcs, rownames)))))
    rl <- lapply(vector("list", ncol(annotation)), function(e) character(0))
    if(length(relpcs) == 0) return(rl)
    rld.pca.dd <- FactoMineR::dimdesc(results$pcaobj, axes=relpcs, proba=alpha)
    pc.bbs <- lapply(results$pcs, function(a) {
        if(nrow(a) == 0) return(numeric(0))
        return(diff(apply(a, 1, range)))
    })
    cor.lb <- 1/sqrt((nrow(annotation) - 2)/qt(alpha/2, nrow(annotation) - 2, lower.tail = F)^2 + 1)
    for(e in 1:ncol(annotation)) {
        if(class(annotation[,e]) != "factor") {
            rv <- lapply(vector("list", 2), function(e) character(0))
            names(rv) <- c("high-low", "low-high")
            if(nrow(results$pcs[[e]]) > 0) {
                sig.coef <- summary(lm(annotation[,e] ~ results$pcasp[,as.integer(rownames(results$pcs[[e]]))]))$coefficients[-1,]
                sig.list <- sapply(1:nrow(results$pcs[[e]]), function(i) {
                    if(sig.coef[i,4] < alpha) {
                        sig.cor <- rld.pca.dd[[match(as.integer(rownames(results$pcs[[e]])[i]), relpcs)]]$quanti
                        if(sig.coef[i,1] > 0) return(list(first=rownames(sig.cor)[sig.cor[,1] >= cor.lb], second=rownames(sig.cor)[sig.cor[,1] <= -cor.lb]))
                        else return(list(first=rownames(sig.cor)[sig.cor[,1] <= -cor.lb], second=rownames(sig.cor)[sig.cor[,1] >= cor.lb]))
                    }
                    return(list(first=character(0), second=character(0)))
                })
                ui <- c()
                ui1 <- c()
                for(i in seq(1,nrow(sig.list),by=2)) {
                    ui <- union(ui, unlist(sig.list[i,]))
                    ui1 <- union(ui1, unlist(sig.list[i+1,]))
                }
                rv[[1]] <- setdiff(ui, ui1)
                rv[[2]] <- setdiff(ui1, ui)
            }
            rl[[e]] <- rv
        }
        else {
            grp.pairs <- gtools::combinations(length(levels(annotation[,e])), 2)
            contrasts <- c(apply(grp.pairs, 1, function(i) c(paste(levels(annotation[,e])[i[2]], levels(annotation[,e])[i[1]], sep='-'), paste(levels(annotation[,e])[i[1]], levels(annotation[,e])[i[2]], sep='-'))))
            rv <- lapply(vector("list", length(contrasts)), function(e) character(0))
            names(rv) <- contrasts
            if(nrow(results$pcs[[e]]) > 0) {
                sig.list <- sapply(1:nrow(results$pcs[[e]]), function(i) {
                    hsd <- TukeyHSD(aov(results$pcasp[,as.numeric(rownames(results$pcs[[e]])[i])] ~ annotation[,e]), deparse(substitute(annotation[,e])))[[1]][,4]
                    return(sapply(1:length(hsd), function(p) {
                        if(hsd[p] < alpha) {
                            rel.dist <- abs(diff(results$pcs[[e]][i,grp.pairs[p,]]))/pc.bbs[[e]][i]
                            pw.lb <- 1 - rel.dist*(1 - cor.lb)
                            cor.ord <- order(c(abs(results$pcs[[e]][i,grp.pairs[p,1]]), abs(results$pcs[[e]][i,grp.pairs[p,2]])), decreasing=T)
                            sig.cor <- rld.pca.dd[[match(as.integer(rownames(results$pcs[[e]])[i]), relpcs)]]$quanti
                            if(results$pcs[[e]][i,grp.pairs[p,2]] > results$pcs[[e]][i,grp.pairs[p,1]]) return(list(first=rownames(sig.cor)[sig.cor[,1] >= pw.lb], second=rownames(sig.cor)[sig.cor[,1] <= -pw.lb]))
                            else return(list(first=rownames(sig.cor)[sig.cor[,1] <= -pw.lb], second=rownames(sig.cor)[sig.cor[,1] >= pw.lb]))
                       }
                        return(list(first=character(0), second=character(0)))
                    }))
                })
                for(i in seq(1,nrow(sig.list),by=2)) {
                    ui <- unique(unlist(sig.list[i,]))
                    ui1 <- unique(unlist(sig.list[i+1,]))
                    rv[[i]] <- setdiff(ui, ui1)
                    rv[[i+1]] <- setdiff(ui1, ui)
                }
            }
            rl[[e]] <- rv
        }
    }
    return(rl)
}

# deprecated
plotVars2 <- function(results, annotation) {
    vars <- round(100*unlist(lapply(results$pcs, function(p) sum(results$pcavar[p]))))
    varpcs <- results$pcs
    vn <- paste(colnames(annotation), paste(vars, "%", sep=''), sep='\n')
    cl <- sapply(1:length(vn), function(i) gtools::combinations(length(vn), i))
    vdn <- unlist(lapply(cl, function(e) apply(e, 1, function(i) paste(vn[i], collapse='&'))))
    vdv <- rep(0, sum(unlist(lapply(cl, nrow))))
    curind <- 1
    for(e in length(cl):1) {
        for(i in nrow(cl[[e]]):1) {
            ov <- rle(sort(unlist(sapply(cl[[e]][i,], function(p) varpcs[[p]]))))
            jvex <- ov$values[ov$lengths == length(cl[[e]][i,])]
            vdv[curind] <- sum(results$pcavar[jvex])
            for(p in cl[[e]][i,]) varpcs[[p]] <- varpcs[[p]][!(varpcs[[p]] %in% jvex)]
            curind <- curind + 1
        }
    }
    vdv <- rev(vdv)
    names(vdv) <- vdn
    plot(venneuler::venneuler(sort(vdv[vdv > 0], decreasing=T)))
}

#' @import stats
htsPCA <- function(rd, l, ntop, verbose, prenorm = FALSE, autosel = TRUE, scaleVar = TRUE, minVar = 1) {
    rv <- list()
    if(prenorm) {
        tdds <- rd
    } else {
        if(all(as.integer(rd) == rd) && min(rd) >= 0) {
            sf <- DESeq2::estimateSizeFactorsForMatrix(rd)
            dds <- DESeq2::DESeqDataSetFromMatrix(countData = rd, colData = l, design = ~ 1)
            if(verbose) message("Normalizing data")
            if((max(sf)/min(sf) >= 4) && (nrow(l) < 200)) {
                tdds <- DESeq2::rlog(dds)
            } else {
                tdds <- DESeq2::varianceStabilizingTransformation(dds)
            }
            tdds <- SummarizedExperiment::assay(tdds)
        } else
            tdds <- vsn::justvsn(rd)
    }
    if(autosel) {
        message("Finding high variance features")
        tdds.means <- apply(tdds, 1, mean)
        tdds.vars <- apply(tdds, 1, var)
        th <- ifelse(length(tdds.means) > 1000, 'approximate', 'exact')
        tdds.lsd <- msir::loess.sd(tdds.means, tdds.vars, qnorm(.975), span = .75, family = 'symmetric', trace.hat = th)
        tdds <- tdds[order(tdds.means), ]
        sel.feat <- (tdds.vars[order(tdds.means)] > tdds.lsd$upper)
        rv$as <- (tdds.vars[order(tdds.means)] > tdds.lsd$upper)
    } else if(ntop < nrow(tdds)) {
        message("Finding high variance features")
        var.order <- order(apply(tdds, 1, var), decreasing=T)
        rv$varord <- var.order
        sel.feat <- sort(head(var.order, n = ntop))
    } else
        sel.feat <- 1:nrow(tdds)
    sel.rows <- tdds[sel.feat,]
    zrows <- apply(sel.rows, 1, function(i) return(sum(i == 0) == length(i)))
    if(sum(zrows) > 0) {
        warning("Some selected measurements are 0 for all samples. Removing these...")
        sel.rows <- sel.rows[!zrows,]
    }
    if(verbose) message("Performing PCA")
    pcs <- FactoMineR::PCA(t(sel.rows), graph=F, ncp=ncol(sel.rows), scale.unit = scaleVar)
    varex <- pcs$eig$percentage
    rv$dat <- tdds
    rv$varex <- varex[varex >= minVar]
    rv$pcasp <- pcs$ind$coord[, varex >= minVar]
    rv$pcaobj <- pcs
    return(rv)
}

pcaShell <- function(dat, pcs) {
    rv <- list()
    varex <- pcs$eig$percentage
    rv$dat <- dat
    rv$varex <- varex[varex >= 1]
    rv$pcasp <- pcs$ind$coord[, varex >= 1]
    rv$pcaobj <- pcs
    return(rv)
}

totVarEx <- function(pcsp, ann, csp, grp) {
    if(class(ann[,grp]) == "factor") return(discVarEx(pcsp, ann, csp, grp))
    return(contVarEx(ann[,grp,drop=F], pcsp[,csp]))
}

discVarEx <- function(pcsp, ann, csp, grp) {
    return(1 - sum(unlist(sapply(levels(ann[, grp]), function(l) {
        cmean <- mean(pcsp[ann[, grp] == l, csp]);
        sapply(pcsp[ann[, grp] == l, csp], function(i) sum((i - cmean)^2))
    })))/sum(scale(pcsp[, csp], scale = F)^2))
}

# deprecated
#newMean <- function(x) .Internal(mean(x))

#' @import cluster
newSil <- function(x, dist) {
    n <- length(x)
    k <- length(unique(x))
    .C(sildist, d = dist, n, x, k, diC = numeric(n * k), counts = integer(k), si = numeric(n), neighbor = integer(n), ismat = FALSE)$si
}

quickcor <- function(x, y) .Call(C_cor, x, y, 4, FALSE)

environment(quickcor) <- asNamespace('stats')

environment(newSil) <- asNamespace('cluster')

# deprecated
#environment(newMean) <- asNamespace('base')

# deprecated
#aicc <- function(fit, scale = 0, k = 2, ...) {
#    aic <- extractAIC(fit = fit, scale = scale, k = k, ...)
#    return(c(aic[1], aic[2] + 2*aic[1]*(aic[1]+1)/(length(fit$residuals) - aic[1] - 1)))
#}

contVarEx <- function(x, y) apply(x, 2, function(j) quickcor(j, y)^2)

# deprecated
#forSelAICc <- function(x, y) {
#    sel.feat <- c()
#    paicc <- Inf
#    while(length(sel.feat) < ncol(x)) {
#        cinds <- (1:ncol(x))[!(1:ncol(x) %in% sel.feat)]
#        caicc <- sapply(cinds, function(i) {
#            rv <- aicc(lm(y ~ x[,c(sel.feat, i)]))[2]
#            if(is.finite(rv))
#                return(rv)
#            return(Inf)
#        })
#        caicc.wm <- which.min(caicc)
#        if(caicc[caicc.wm] >= paicc)
#            break
#        sel.feat <- c(sel.feat, cinds[caicc.wm])
#        paicc <- caicc[caicc.wm]
#    }
#    return(sel.feat)
#}

#' @import foreach
forSelParRepDisc <- function(x, l, alpha, eff.size, disc.thresh) {
    pccol.range <- 1:ncol(x)
    rv <- list()
    cann <- NULL
    j <- NULL
    rv$pcs <- foreach(cann=1:ncol(l), .packages = c('cluster', 'stats'), .export = 'newSil') %dopar% {
        if(class(l[,cann]) != "factor") {
            cfolds <- floor(nrow(x)/3)
            cfolds <- ifelse(cfolds < 5, nrow(x), min(cfolds, 10))
            foo <- glmnet::glmnet(x, as.numeric(l[,cann]), lambda = glmnet::cv.glmnet(x, as.numeric(l[,cann]), nfolds = cfolds, family='gaussian')$lambda.1se, family='gaussian')
            return(which(coef(foo)[, 1] != 0) - 1)
        }
        int.classes <- as.integer(l[,cann])
        if(disc.thresh < 1) {
            ips <- foreach(j=1:ncol(x), .combine = c, .packages = c('cluster', 'stats'), .export = 'newSil') %dopar% ecdf(replicate(1000, mean(newSil(int.classes, dist(x[sample(nrow(x)), j])))))(mean(newSil(int.classes, dist(x[, j]))))
            disc.feat <- which((1 - ips) >= disc.thresh)
        } else {
            disc.feat <- c()
        }
        sel.feat <- c()
        while(length(c(disc.feat, sel.feat)) < ncol(x)) {
            cinds <- pccol.range[!(pccol.range %in% c(disc.feat, sel.feat))]
            pcsel <- x[,sel.feat]
            pvals <- matrix(foreach(j=1:length(cinds), .combine=cbind, .packages = c('cluster', 'stats'), .export = 'newSil') %dopar% {
                ccol <- cinds[j]
                csil <- mean(newSil(int.classes, dist(x[,c(sel.feat,ccol)])))
                perms <- replicate(10000, {
                    mean(newSil(int.classes, dist(cbind(pcsel, x[sample(nrow(x)),ccol]))))
                })
                c(1-ecdf(perms)(csil), csil-mean(perms))
            }, nrow=2)
            sig.inds <- ((pvals[1,] <= alpha) & (pvals[2,] >= eff.size))
            if(sum(sig.inds) >= 1) {
                sel.feat <- c(sel.feat, (cinds[order(pvals[2,], decreasing=T)])[which.max(sig.inds[order(pvals[2,], decreasing=T)])])
            } else {
                break
            }
            disc.feat <- c(disc.feat, cinds[which(pvals[1,] >= disc.thresh)])
        }
        return(sel.feat)
    }
    names(rv$pcs) <- colnames(l)
    rv$pcs <- lapply(seq_along(rv$pcs), function(e) {
        if(length(rv$pcs[[e]]) == 0) return(matrix(0, nrow = 0, ncol = 0))
        if(class(l[,e]) != "factor") {
            cfolds <- floor(nrow(x)/3)
            cfolds <- ifelse(cfolds < 5, nrow(x), min(cfolds, 10))
            return(matrix(coef(suppressWarnings(glmnet::glmnet(x, as.numeric(l[,e]), lambda = glmnet::cv.glmnet(x, as.numeric(l[,e]), nfolds = cfolds, family='gaussian')$lambda.1se, family='gaussian')))[rv$pcs[[e]] + 1], ncol = 1, dimnames = list(rv$pcs[[e]])))
        }
        return(matrix(sapply(levels(l[,e]), function(i) sapply(rv$pcs[[e]], function(j) mean(x[l[,e] == i, j]))), ncol = length(levels(l[,e])), dimnames = list(rv$pcs[[e]], levels(l[,e]))))
    })
    return(rv)
}

# deprecated
#findPCs <- function(x, l, n, ord) {
#    rv <- list()
#    pcsets <- gtools::combinations(ncol(x), n)
#    vals <- sapply(1:ncol(l), function(j) {
#        sil.inds <- apply(pcsets, 1, function(pcs) c(summary(cluster::silhouette(as.integer(l[,j]), dist(x[ord, pcs])))$avg.width))
#        return(c(which.max(sil.inds), max(sil.inds), summary(svm(x[ord, pcsets[which.max(sil.inds),]], l[,j], cross = nrow(x)))$tot.accuracy))
##        return(c(which.max(sil.inds), max(sil.inds), max(summary(svm(x[ord, pcsets[which.max(sil.inds),]], l[,j], cross = nrow(x)))$tot.accuracy)))
#    })
#    rv$pcs <- vals[1,]
#    rv$sil <- vals[2,]
#    rv$cv <- vals[3,]
#    return(rv)
#}

# deprecated
#iterateDims <- function(x, l, ord) {
#    dv <- sapply(1:min(3,ncol(x)), function(d) {
#        r <- findPCs(x, l, d, ord)
#        return(c(r$cv, r$sil, r$pcs))
#    })
#    pcsets <- vector("list", choose(ncol(x), min(3, ncol(x))))
#    for(i in 1:min(3,ncol(x))) {
#        pcsets[[i]] <- gtools::combinations(ncol(x), i)
#    }
#    rv <- list()
#    rv$cv <- rep(0, ncol(l))
#    rv$sil <- rep(0, ncol(l))
#    rv$pcs <- vector("list", ncol(l))
#    for(i in 1:ncol(l)) {
###        max.sil <- which.max(dv[ncol(l) + i,])
###        max.dims <- which.max(dv[i,max.sil:3]) + (max.sil - 1)
#        max.dims <- which.max(dv[i,])
#        rv$cv[i] <- dv[i, max.dims]
#        rv$sil[i] <- dv[ncol(l) + i, max.dims]
#        cpcs <- pcsets[[max.dims]][dv[2*ncol(l) + i, max.dims],]
#        rv$pcs[[i]] <- sapply(levels(l[,i]), function(g) sapply(cpcs, function(pc) mean(x[l[,i] == g, pc]))) #{
##            if(sum(l[,i] == g) == 1) return(mean(x[l[,i] == g, pc]))
##            pc.t <- t.test(x[l[,i] == g, pc])
##            if(signif(pc.t$p.value, 2) <= .1) return(pc.t$estimate)
##            return(0)
##        }))
#        if(!is.matrix(rv$pcs[[i]])) {
#            rv$pcs[[i]] <- matrix(rv$pcs[[i]], nrow=1)
#            colnames(rv$pcs[[i]]) <- levels(l[,i])
#        }
#        rownames(rv$pcs[[i]]) <- cpcs
#    }
#    return(rv)
#}

# deprecated
#permuteLabels <- function(x, l, sil, cv, b, verbose) {
#    if(verbose) message("Permuting labels")
#    rv <- list()
#    orderings <- t(replicate(b, sample(nrow(l))))
##    pcsets <- combinations(21, n)
#    rv$orderings <- orderings
#    if(verbose) message("Finding PCs for permuted data")
#   # vals <- apply(orderings, 1, function(ord) {
#    vals <- foreach(ord=iter(orderings, by='row'), .combine=cbind, .inorder=F) %dopar% {
#        r <- iterateDims(x, l, ord)
#        return(c(r$sil, r$cv))
#   # })
#    }
#    rv$sil.p <- rep(0, ncol(l))
#    rv$cv.p <- rep(0, ncol(l))
#    for(i in 1:ncol(l)) {
#        rv$sil.p[i] <- sum(vals[i, ] >= sil[i])/b
#        rv$cv.p[i] <- sum(vals[ncol(l) + i, ] >= cv[i])/b
#    }
##    vals <- t(apply(orderings, 1,
##        function(ord) sapply(1:ncol(l),
##            function(j) {
##                sil.inds <- apply(pcsets, 1,
##                    function(pcs) c(summary(silhouette(as.integer(l[,j]), dist(x[ord, c(pcs)])))$avg.width))
##                    return(c(max(sil.inds), max(summary(svm(x[ord, pcsets[which.max(sil.inds),]], l[,j], cross = nrow(x)))$tot.accuracy)))
##            }
##        )
##    ))
##    rv$sil <- t(apply(vals, 1, function(i) i[seq(1, length(i), by=2)]))
##    rv$cv <- t(apply(vals, 1, function(i) i[seq(2, length(i), by=2)]))
#    return(rv)
#}

projectPoint <- function(v1, v2, p) {
    return(v1 + sum((p - v1)*(v2 - v1))/sum((v2 - v1)*(v2 - v1))*(v2 - v1))
}

# Adapted from FactoMineR::reconst
reconstMatrix <- function(res, ncp = NULL) {
    if (is.null(ncp)) 
        ncp <- 1:ncol(res$ind$coord)
    coord.var = as.matrix(res$var$coord)[, ncp, drop = F]
    coord.ind = as.matrix(res$ind$coord)[, ncp, drop = F]
    hatX = coord.ind %*% t(sweep(coord.var, 2, sqrt(res$eig[ncp, 1]), FUN = "/"))
    hatX = sweep(hatX, 2, res$call$ecart.type, FUN = "*")
    hatX = sweep(hatX, 2, res$call$centre, FUN = "+")
    return(hatX)
}

# deprecated
#projectPoint2d <- function(v1, v2, p) {
#    m <- (v2[2] - v1[2])/(v2[1] - v1[1])
#    b <- v1[2] - v1[1]*m
#    mp <- -1/m
#    bp <- p[2] - p[1]*mp
#    px <- (bp - b)/(m - mp)
#    py <- m*px + b
#    return(c(px, py))
#}
