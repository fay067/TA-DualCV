    set.seed(100)
    library(hash)
    library(mice)
    library(GPfit)
    library(doParallel)
    library(foreach)
    library(abind)
    source('evalPtTensorImpImportTrTe.R')
    source('miceCrossSectionInit.R')
    source('selfNorm.R')
    source('selfDeNorm.R')
    source('ciSelfDeNorm.R')
    source(fncf)
    
    m=10 # number of iterations
    trte='tr'
    nimp=50 # multiple number of ncores
    ncores=5 # number of cores
    fncf=fncf
    ci.thr=0.2
    norm="self"
    sd=11242015
    maxit=10
    mincor=0.5
    cl=NULL
    param=T
    nug_thres=20
    w1_1 = 0.5 # weight of CFP
    w1_2 = 0.5 # weight of CTP
    n_lab = 13 # number of features
    n_time = 24 # should be median of number of time steps
    fncf='mimicConfig_release.R'
    source(fncf)
    
    dnroot <- "./data/mimic"
    
    pts.tr = read.csv(sprintf('%s/pts.tr.csv', dnroot), header=F)[,1]
    
    tnagt.tr = list()
    dn.tnagt = sprintf('%s/data_with_missing', dnroot)
    for (i in pts.tr) {
        fn = sprintf('%s/%d.csv', dn.tnagt, i)
        tnagt.tr[[i]] = t(read.csv(fn))
    }
    names(tnagt.tr) = as.character(1:length(tnagt.tr))
    
    naidx = read.csv(sprintf('%s/naidx.csv', dnroot), stringsAsFactors=F)
    print(length(tnagt.tr))
    h = hash()
    h[['tna']] = tnagt.tr
    h[['naidx']] = naidx
    
    print(timeidx)
    cl.init = F;
    if (is.null(cl)) {
        cl=makeCluster(ncores, type="FORK", outfile="") #
        registerDoParallel(cl)
        cl.init = T
        clusterSetRNGStream(cl, sd)
    }

    qtl = c(0, 0.25, 0.5, 0.75, 1)
    V = h[['tna']]; Y = V; naidx = h[['naidx']];
    ci = V; pts = names(V)
    cat(sprintf('na tensor len %d\n', length(V)))

    hmax = hash(); hmin = hash(); hstd = hash()
    Jl = dim(V)[2]; n = length(V)
    if (norm=='self') {
        V = selfNorm(V, hmax, hmin, fncf=fncf)
    }else if (norm=='std') {
        V = stdNorm(V, hstd, fncf=fncf)
    }
    r = miceCrossSectionInit(V, nimp=nimp, sd=sd, maxit=maxit, mincor=mincor, cl=cl, ncores=ncores, fncf=fncf)
    V = r$V; V.sd = r$std
    for (pt in pts) {
        ci[[pt]] = 1.96*V.sd[[pt]]
    }
    if (norm=="self") {
        V.raw = selfDeNorm(V, hmax, hmin, fncf=fncf)
        ci.raw = ciSelfDeNorm(ci, hmax, hmin, fncf=fncf)
    }else if (norm=="std") {
        V.raw = stdDeNorm(V, hstd, fncf=fncf)
        ci.raw = ciStdDeNorm(ci, hstd)
    }else {
        V.raw = V
        ci.raw = ci
    }
    
    for (c in 1:m) {# for m iters
        nleft = 0; ci.v = c()
        for (ipt in pts) {
            for (iv in tests) {
                selna = is.na(Y[[ipt]][iv,])
                sel = selna & (ci[[ipt]][iv,] > ci.thr)
                V[[ipt]][iv,sel] = NA
                nleft = nleft + sum(sel)
                ci.v = c(ci.v, ci[[ipt]][iv, selna])
            }
        }
        cat(sprintf('iter=%d, conf int quantile\n', c))
        print(quantile(ci.v, probs=qtl, na.rm=T), digits=4)
        ## cross sectional impute for all variables
        ## concatenate all time slices
        x = c() # columns are variables
        x_l = c() # used for saving rows number of each patient
        x_l_s = c(0) # used for saving acc rows number of each patient
        for (xsub in V) {
            x_l = c(x_l, ncol(xsub))
            x_l_s = c(x_l_s, nrow(x))
            x = rbind(x, t(xsub[tests,]))
        }

        imp <- foreach(no = 1:ncores, .combine=ibind, .packages="mice") %dopar% {
            mice(x, m=nimp/ncores, mincor=mincor, maxit=maxit, printFlag=F) # may want to increase the maxit and monitors the convergence, which statistics? may even consider parallelizing
        }
        xt = array(NA, dim=c(dim(x)[1], dim(x)[2], nimp))
        for (i in 1:nimp) {
            xt[,,i] = as.matrix(complete(imp, i))
        }
        
        # # cross sectional impute for all time slices
        # # Concatenate all variables
        x_2 = list() # total table: columns are time slices
        count = 1 # used for initialize x_2
        for (xsub_2 in V) {
            # make current xsub_2 to a 13tests*17time table
            curr_col = ncol(xsub_2)
            if (curr_col < n_time) {
                for (iix in (curr_col+1):n_time) {
                    xsub_2 = cbind(xsub_2, NA)
                }
            }
            else if (curr_col > n_time) {
                xsub_2 = xsub_2[, 1:n_time]
            }
            if (count == 1) { # initialize x_2
                for (varb in 1:length(tests)) {
                    x_2[[varb]] = xsub_2[tests[varb], 1:n_time]
                }
            }
            else {
                for (varb in 1:length(tests)) {
                    x_2[[varb]] = rbind(x_2[[varb]], xsub_2[tests[varb], 1:n_time])
                }
            }
            count = count + 1
        }
        
        xt_2 = array(NA, dim=c(1, n_time, nimp)) # xt_2 is the total table used for combining
        for (temp_i in 1:n_lab) { # we have (n_lab) tables totally, each table contains row as Xi and col as time slices
            imp_2 <- foreach(no = 1:ncores, .combine=ibind, .packages="mice") %dopar% {
                mice(x_2[[temp_i]], m=nimp/ncores, mincor=mincor, maxit=maxit, printFlag=F)
            }
            
            temp_xt_2 = array(NA, dim=c(length(V), n_time, nimp)) # temp_xt_2 is used to save imputed table for each temp_i
            for (i in 1:nimp) {
                temp_xt_2[,,i] = as.matrix(complete(imp_2, i))
            }
            if (temp_i == 1) {
                xt_2 = temp_xt_2
            }
            else {
                xt_2 = abind(xt_2, temp_xt_2, along = 1)
            }
        }

        # combine CFP and CTP with weights
        for (temp_v in 1:dim(xt_2)[1]) {
            for (temp_t in 1:dim(xt_2)[2]) {
                trans_v = ceiling(temp_v / length(V))
                curr_pts = temp_v %% length(V)
                if (curr_pts == 0) { curr_pts = length(V) }
                trans_t = x_l_s[curr_pts] + temp_t
                if (temp_t <= x_l[curr_pts]) {
                    xt[trans_t, trans_v,] = w1_1 * xt[trans_t, trans_v,] + w1_2 * xt_2[temp_v, temp_t,]
                }
            }
        }
        
        ## reconstruct V1
        V1 = V; V1.sd = V; cursor = 0
        if (param) {
            xm = apply(xt, c(1,2), mean)
            sdm = apply(xt, c(1,2), sd)
            for (pt in names(V)) {
                batch = dim(V1[[pt]])[2]
                ## rnames = rownames(V[[pt]])[-1]; cnames = colnames(V[[pt]])
                V1[[pt]][tests,] = t(xm[cursor+1:batch,]); # rownames(V1[[pt]]) = rnames; colnames(V1[[pt]]) = cnames
                V1.sd[[pt]][tests,] = t(sdm[cursor+1:batch,]); # rownames(V1.sd[[pt]]) = rnames; colnames(V1.sd[[pt]]) = cnames
                cursor = cursor + batch
            }
        }else {
            ## reconstruct V1m, with multiple imputations
            V1m = list()
            for (pt in names(V)) {
                vm = V[[pt]]
                batch = dim(vm)[2]
                V1m[[pt]] = array(NA, dim=c(dim(vm)[1], dim(vm)[2], nimp))
                dimnames(V1m[[pt]])[[1]] = rownames(vm)
                dimnames(V1m[[pt]])[[2]] = colnames(vm)
                for (i in 1:nimp) {
                    V1m[[pt]][tests,,i] = t(xt[cursor+1:batch,,i]);
                }
                V1[[pt]][tests,] = apply(V1m[[pt]][tests,,], c(1,2), mean)
                V1.sd[[pt]][tests,] = apply(V1m[[pt]][tests,,], c(1,2), sd)
                cursor = cursor + batch
            }
        }

        cat(sprintf('iter=%d, gp\n', c))
        ## parallel computation of V2 and V
        vci <- foreach (i=1:length(V), .packages='GPfit') %dopar% {
            pt = pts[i]
            ev = V[[pt]]
            ev2 = V[[pt]]
            x = ev[timeidx,]
            x = (x - min(x))/(max(x)-min(x))
            eci = ci[[pt]]
            esd2 = V1.sd[[pt]]
            ediff = list()
            for (test in tests) {
                ediff[[test]] = c()
                y = ev[test,]; itr = !is.na(y); ite = is.na(y)
                ytr = y[itr]; xtr = x[itr]
                xte = x[ite]
                if (length(xte)>0) {
                    gpmod = GP_fit(xtr, ytr, nug_thres=nug_thres)
                    gppred = predict(gpmod, xte)
                    mu2 = gppred$Y_hat; ev2[test, ite] = mu2
                    sigma2 = sqrt(gppred$MSE); esd2[test, ite] = sigma2

                    if (param) {
                        mu1 = V1[[pt]][test, ite]
                        sigma1 = V1.sd[[pt]][test, ite]
                        w1 = sigma2/ (sigma1 + sigma2)
                        w2 = 1-w1
                        mu = w1*mu1 + w2*mu2
                        ev[test, ite] = mu
                        ediff[[test]] = abs(mu1 - mu2)
                        eci[test, ite] = 1.96*sqrt(((mu1-mu)^2+(mu2-mu)^2)/2)
                    }else { # non-parametric
                        j = 0
                        for (i in which(ite)) {
                            j = j+1
                            V1.v = V1m[[pt]][test,i,]
                            if (sigma2[j]==0) {
                                nimpw = nimp
                                cat(sprintf('zero gp sd for %s, %s, %s\n', pt, test, i))
                            }else {
                                nimpw = ceiling(nimp*sd(V1.v)/sigma2[j])
                            }
                            V2.v = rnorm(nimpw, mean=mu2[j], sd=sigma2[j])
                            if (sigma2[j]==0) {
                                V12.v = V2.v
                            }else {
                                V12.v = c(V1.v, V2.v)
                            }
                            ev[test, i] = mean(V12.v)
                            eci[test, i] = 1.96*sd(V12.v)
                            if (eci[test,i]>1.2) {
                                cat(sprintf('big ci %.3f at %s, %s, %s\n', eci[test, i], pt, test, i))
                            }
                            ediff[[test]] = abs(mean(V1.v) - mean(V2.v))
                        }
                    }
                }
            }
            er = list(v=ev, ci=eci, diff=ediff, v2=ev2, sd2=esd2)
            er
        }
        V = lapply(vci, function(a) {a$v}); names(V) = pts
        ci = lapply(vci, function(a) {a$ci}); names(ci) = pts

        V2 = lapply(vci, function(a) {a$v2}); names(V2) = pts
        V2.sd = lapply(vci, function(a) {a$sd2}); names(V2.sd) = pts

        ptdiff = lapply(vci, function(a) {a$diff}); names(ptdiff) = pts
        vdiff = list()
        for (test in tests) {
            vdiff[[test]] = c()
            for (pt in pts) {
                vdiff[[test]] = c(vdiff[[test]], ptdiff[[pt]][[test]])
            }
            cat(sprintf('iter=%d, %s diff\n', c, test))
            print(quantile(vdiff[[test]], qtl, na.rm=T), digits=4)
        }

        if (norm=='self') {
            V.raw = selfDeNorm(V, hmax, hmin, fncf=fncf)
            V1.raw = selfDeNorm(V1, hmax, hmin, fncf=fncf)
            V2.raw = selfDeNorm(V2, hmax, hmin, fncf=fncf)
            ci.raw = ciSelfDeNorm(ci, hmax, hmin, fncf=fncf)
        }else if (norm=='std') {
            V.raw = stdDeNorm(V, hstd, fncf=fncf)
            ci.raw = ciStdDeNorm(ci, hstd)
        }else {
            V.raw = V
            ci.raw = ci
        }
        
        res = list(t.imp=V.raw, tn.imp=V, ci.imp=ci.raw, cin.imp=ci, V1.raw=V1.raw, V2.raw=V2.raw, V1.sd=V1.sd, V2.sd=V2.sd, h=h)
        fnres = sprintf(fnres.tmp, trte, c)
        save(res, file=fnres)
    }
    if (cl.init) {
        stopCluster(cl)
    }
    
    load("./data/data_with_missing/micegp_log/tr_res_iter2.RData")
    for (i in seq(1,length(V), by=1)) {write.csv(t(res$t.imp[[i]]), sprintf("%d.csv", i), row.names=FALSE)}
    

