from planck import planck
from numpy import logspace, sum, log10, log
from matplotlib.pyplot import loglog, show, xlim

temp=100.                                           # temperature in K
mum=1e-4                                            # micrometer in cgs
log_lambda_0=log10(9.*mum)                          # lower limit
n=16                                                # interval split
lamb=logspace(log_lambda_0,log_lambda_0+2.0,n+1)    # log scale
loglog(lamb,planck(temp,lamb))                      # loglog plot
xlim(lamb[0]*0.9,lamb[n]*1.1)                       # tight margins
show()                                              # show plot
delta=log(lamb[1]/lamb[0])                          # delta ln(lambda)
print sum(planck(temp,lamb))*delta                  # approx integral