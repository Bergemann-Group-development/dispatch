#!/usr/bin/env python

"""
make_plot.py
Created by Mikolaj on 2022-02-14.
"""

import os, sys
import pylab as plt


plt.style.use('dark_background')

# dispatch_ver = 'dispatch'
dispatch_ver = 'dispatch2'

dispatch_python_path = os.getenv('DISPATCH_DIR', '../../../../')
sys.path.append(os.path.join(dispatch_python_path,dispatch_ver + '/utilities/python/'))

import dispatch
import dispatch.select as ds
import dispatch.graphics as dg

d = dispatch.snapshot(1,data="data/ws1")

dr = ds.haver(d,iv='d')
dt = ds.haver(d,iv='T')
dp = ds.haver(d,iv='pg')

fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
plt.suptitle("Profiles after Δt = out_time",va="center",ha="center",size=15.)
fig.set_size_inches(14, 4.5)

ax1.plot(dr[0],dr[1],c='r',label='var(d)')
ax1.set_xlabel('z-axis')
ax1.set_ylabel('Density [code-units]')
ax1.legend(frameon=False)
ax1.set_title("Density")

ax2.plot(dt[0],dt[1],c='r',label='var(T)')
ax2.set_xlabel('z-axis')
ax2.set_ylabel('Temperature [code-units]')
ax2.legend(frameon=False)
ax2.set_title("Temperature")

ax3.plot(dp[0],dp[1],c='r',label='var(pg)')
ax3.set_xlabel('z-axis')
ax3.set_ylabel('Pressure [code-units]')
ax3.legend(frameon=False)
ax3.set_title("Pressure")

fig.set_size_inches(16,5)
plt.savefig("ws.1.png",dpi=150)
plt.show()