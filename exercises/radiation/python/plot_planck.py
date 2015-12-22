import planck
reload(planck)
from cgs import cgs
from planck import planck
from numpy import logspace, sum, log10, log, pi
from matplotlib.pyplot import loglog, show, xlim, draw, clf

temp=100.                                           # temperature in K
log_lambda_0=log10(7e-4)                            # lower limit
n=128                                               # interval split
lamb=logspace(log_lambda_0,log_lambda_0+3.0,n+1)    # log scale
clf()                                               # clear figure
loglog(lamb,planck(temp,lamb))                      # loglog plot
xlim(lamb[0]*0.9,lamb[n]*1.1)                       # tight margins
draw()                                              # redraw
show()                                              # show plot
delta=log(lamb[1]/lamb[0])                          # delta ln(lambda)
print sum(planck(temp,lamb))*delta                  # approx integral
stefan = 5.670367e-5                                # Stefan's constant in cgs
h=6.626e-27
k=1.38064e-16
c=2.9997e10
stefan = 2.0*pi**5*(k/h)**3*k/(15.*c**2)
print pi*sum(planck(temp,lamb))*delta/(stefan*temp**4)