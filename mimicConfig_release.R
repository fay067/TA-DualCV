
dnroot = "./data/test_with_missing"
tests = c('X1', 'X2', 'X3', 'X4', 'X5', 'X6', 'X7', 'X8', 'X9', 'X10', 'X11', 'X12', 'X13')

fnlog.tmp = sprintf('%s/micegp_log/%%s_output_iter%%s.RData', dnroot)
fnres.tmp = sprintf('%s/micegp_log/%%s_res_iter%%s.RData', dnroot)
fnkm.tmp = sprintf('%s/micegp_log/%%s_km_iter%%s.RData', dnroot)

timeidx = 'time'


