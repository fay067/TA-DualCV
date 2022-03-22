evalPtTensorImpImportTrTe <- function(timp, tgt, ci, naidx, V1.raw=NULL, V2.raw=NULL, V1.sd=NULL, V2.sd=NULL, metric='RMSE range norm', err.out=T, trte='tr', iter=1, fncf='mghtsConfig.R') {
    ## metric can be RMSE std norm, RMSE range norm, MAPE
    ## naidx = as.data.frame(which(is.na(tna), arr.ind=T))
    library(hash)
    source(fncf)
    qtl = (1:10)/10
    
    if (metric=='RMSE std norm') {
        hstd = hash()
        for (test in tests) {
            hstd[[test]] = c()
            for (pt in names(tgt)) {
                hstd[[test]] = c(hstd[[test]], tgt[[pt]][test,])
            }
            hstd[[test]] = sd(hstd[[test]], na.rm=T)
        }
    }
    print(tests)
    nrmse = rep(0, length(tests))
    names(nrmse) = tests
    nrmse.v = list()
    output.v = list()
    for (test in tests) {
        nrmse.v[[test]] = c()
        output.v[[test]] = c()
    }
    n = dim(naidx)[1]
    cfun <- function(a, b){return(a+b)}
    # for(i in 1:n) {
    nrmse <- foreach(i = 1:n, .combine="cfun")%dopar%{
        ipt = naidx[i,1]; iv = naidx[i,2]; it = naidx[i,3]
        eimp = timp[[ipt]][iv, it]; eraw = tgt[[ipt]][iv, it]
        ev1 = V1.raw[[ipt]][iv, it]; ev2 = V2.raw[[ipt]][iv, it]
        esd1 = V1.sd[[ipt]][iv, it]; esd2 = V2.sd[[ipt]][iv, it]
        eci = ci[[ipt]][iv, it]
        if (metric=='RMSE std norm') {
            update = ((eimp - eraw) / hstd[[iv]])^2
        }else if (metric=='RMSE range norm') {
            update = ((eimp - eraw) / (max(tgt[[ipt]][iv,], na.rm=T) - min(tgt[[ipt]][iv,], na.rm=T)))^2
        }else if (metric=='MAPE') {
            update = abs(eimp - eraw) / eraw # elem-wise div
        }
        if (is.na(update) | is.infinite(update)) {
            cat(sprintf('%s, %s, %s\n', ipt, iv, it))                
        }else if (update > .6 & err.out) {
            if (is.null(V2.raw)) {
                cat(sprintf('pt %s, v %s, t %s, err %f, gt %f, imp %f\n', ipt, iv, it, update, eraw, eimp))
            }else {
                cat(sprintf('pt %s, v %s, t %s, err %f, gt %f, imp %f, v1 %f, v2 %f, sd1 %f, sd2 %f\n', ipt, iv, it, sqrt(update), eraw, eimp, ev1, ev2, esd1, esd2))
            }
        }
        
        update_one_hot <- rep(0, 13)
        update_one_hot[as.numeric(gsub("X", "", iv))] <- update
        
        return(update_one_hot)
        
    }
    names(nrmse) = tests
    
    for(test in tests){
        if (metric=='RMSE std norm' | metric=='RMSE range norm') {
            nrmse[test] = sqrt(nrmse[test] / sum(naidx$test==test))
        }else if (metric=='MAPE') {
            nrmse[test] = nrmse[test] / sum(naidx$test==test)
        }
        cat(sprintf('%s n%s\n', test, metric))
    }
    return (nrmse);
}
