ciSelfDeNorm <- function(ci, hmax, hmin, fncf='mghtsConfig.R') {
    source(fncf)
    ci.raw = ci; 
    cat(sprintf('CI self denormalizing\n'))
    for (ipt in names(ci)) {
        for (iv in tests) {
            k = sprintf('%s|%s', ipt, iv)
            kmax = hmax[[k]]
            kmin = hmin[[k]]
            sel = !is.na(ci[[ipt]][iv,])
            ci.raw[[ipt]][iv,sel] = ci[[ipt]][iv,sel] * (kmax - kmin)
        }
    }
    
    cat(sprintf('done\n'))
    return (ci.raw)
}
