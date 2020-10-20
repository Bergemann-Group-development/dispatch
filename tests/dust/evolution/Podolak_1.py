# %%
import numpy as np
import astropy.units as u
from astropy import constants as const
import matplotlib.pyplot as plt
from pylab import rcParams
rcParams['figure.figsize'] = 15, 5

# %%
%%time
print("Initializing starting conditions")
T = 196 *u.K
cs = (np.sqrt(const.k_B*T/(2.37*1.008*u.u))).decompose()
r = 1 * u.au
msolar = 0.5 *u.solMass
Omega_K = np.sqrt(const.G*msolar/r**3)
Hp = (cs / Omega_K).decompose()
sigma_d = 0.18 * u.g/u.cm**2
#rho_d = (sigma_d/(np.sqrt(2*np.pi)*Hp)).decompose()
rho_s = 1.6 *u.g * u.cm**(-3)
a_s = 1*u.micron
m_s = (4/3*np.pi*a_s**3 * rho_s).decompose()
u_f = 1*u.m/u.s
trans_width = 0.001*u.m/u.s
alpha = 1e-3
Cdt = 0.1
xi = 1.83

dt = 0.1 * 60*60*24*365 * u.s
t_end = 1e3 * 60*60*24*365 * u.s
save_time = 1 *u.yr
eps_time = 1.10

# read from figure A.2 in Binsteil+ 2010. a(t) = 1e-3 at t = 6e4 yr, and then plugged in to their eq B.5

# rho_d = (2/5 * ((1e-3 * u.cm)**(5/2) - a_s**(5/2)) * rho_s * np.pi /
#          ((12. * const.k_B * T/rho_s)**(1/2) * 6e4 * u.yr)).decompose()
# rho_d
rho_d =( (1e4*u.cm)*2*rho_s*np.sqrt(np.pi)/sigma_d * 1/(8*1e3*u.yr)/cs * sigma_d /(np.sqrt(2*np.pi)) ).decompose()

#print("rho_d :", rho_d)

# %%
%%time
print("Calculating bins")

nbin = 200 # 40 for brown 200 for turb as per birnstiel
a = np.logspace(-4, 2, nbin) * u.cm
st = rho_s * a * np.pi/(2*sigma_d*100)
rho_d = (sigma_d/ (np.sqrt(2*np.pi)* Hp * min(1,np.sqrt(alpha/min(st[0],1/2)*(1+st[0]**2))) )).decompose()

m = 4/3 * a**3 * rho_s * np.pi
n = np.zeros(len(m)) * u.m**(-3)
n[0] = (rho_d / m[0]).decompose()
ms = np.add.outer(m,m)

# %%
%%time
print("Calculating velocities and sigma_col")
ubm = (np.sqrt(8 * const.k_B*T * (np.add.outer(m, m))/(np.pi*np.outer(m, m)))).decompose()
utm = cs * np.sqrt(2*alpha*st)
utm[st>1] = cs*np.sqrt(2*alpha/st[st>1])

utm = np.maximum.outer(utm,utm)

np.max(utm) /np.max(ubm)

du = ubm + utm
#du = ubm
sigma_col = np.pi * np.add.outer(a, a)**2

# %%
%%time
print("Setting up neighbour, probability matrices etc. ")

size = len(a)
mn = np.zeros((size,size))
mm = np.zeros((size,size))
pf = np.zeros((size,size))

crat_ind = (np.divide.outer(m,m) > 10) | (np.divide.outer(m,m) < 0.1)
SS = np.zeros((size,size,size))

for i in range(len(m)):
    for j in range(len(m)):
        kn = (np.where(m<ms[i][j])[0])[-1]
        #print(i,j,kn)
        #mm[i][j] = (m[m<ms[i][j]])[-1].value
        #if len(m[m>ms[i][j]]) > 0:

        if crat_ind[i][j]:
            #cratering
            SS[:min(i,j)+1,i,j] = a[:min(i,j)+1] ** (-3*xi)
            Stot = np.sum(m.value * SS[:,i,j] ) / (m[min(i,j)].value *2) 
            SS[:min(i,j)+1,i,j] /= Stot
            m_crat =  np.abs(m[i].value-m[j].value)
            l = np.where(m.value < m_crat )[0][-1]
            eps = (m_crat - m[l].value)/(m[l+1].value - m[l].value)
            SS[l,i,j] = 1-eps
            SS[l+1,i,j] = eps            

        else:
            #fragmentation
            SS[:max(i,j)+1,i,j] = a[:max(i,j)+1] ** (-3*xi)
            Stot = np.sum(m.value * SS[:,i,j] ) / (m[i].value + m[j].value) 
            SS[:max(i,j)+1,i,j] /= Stot

        if kn+1 < nbin:
            mm[i][j] = m[kn].value
            mn[i][j] = m[kn+1].value
            #mn[i][j] = (m[m>ms[i][j]])[0].value
        else:
            sigma_col[i][j] = 0
            #mn[i][j] = 0

        if du[i][j] < u_f-trans_width:
            pf[i][j] = 0
        elif du[i][j] > u_f:
            pf[i][j] = 1
        else:
            pf[i][j] = 1-(u_f-du[i][j])/trans_width

pc = 1-pf


# %%
%%time
print("Setting up coagulation matrix")

ep = (mn - ms.value) / (mn - mm)
ep[sigma_col==0] = 0

C = np.zeros((size,size,size))
for k in range(size):
    imm = m[k].value == mm
    imn = m[k].value == mn
    C[k][imm] = ep[imm]
    C[k][imn] = 1 - ep[imn]

K = (np.multiply( np.multiply(du, sigma_col), pc)).decompose()
L = (np.multiply( np.multiply(du, sigma_col), pf)).decompose()

# %%
# pf = (du > u_f) * 1
# pf[du < u_f-trans_width] = 0
# elseind = (du < u_f) | (du > u_f-trans_width)
# pf[elseind] = (1-(u_f-du)/trans_width)


# %%
%%time
print("Time evolving the number density")

n = n.value
K = K.value
L = L.value
dt = dt.value
m = m.to(u.kg).value
t_end = t_end.value

x = np.array([])
y = np.array([])
y2 = np.array([])

t = 0
n_save = n
i = 0
while t < t_end:
    if np.any(n<0) or np.abs(np.sum(n*m)/rho_d.value - 1.) > 1e-2:
        if np.any(n<0):
            print("Negative number particles, aborting", np.argmin(n), np.min(n))
        if np.abs(np.sum(n*m)/rho_d.value - 1.) > 1e-4:
            print("Mass not conserved, aborting", (np.sum(n*m)/rho_d.value - 1.))
        print(f"Time elapsed: {t/(60*60*24*365)} yr")
        break

    dn = np.zeros(size)
    dnk_frag = np.zeros(size)
    S = np.zeros(size)
    J = np.zeros((size,size))
    nn = np.outer(n,n)
    Q = np.multiply(K,nn)
    Q_frag = np.multiply(L,nn)
    for k in range(size):
        #dnk_frag[k] = 1/2 * np.sum( np.multiply( Q_frag, SS[k])) - np.sum(Q_frag[k])
        KC_k = np.multiply(K, C[k])
        LS_k = np.multiply(L, SS[k])
        #QC = np.multiply(KC_k, nn)

        Kn_k = np.multiply(K[k], n)
        #dn[k] = (0.5 * np.sum(QC) - np.sum(Kn_k*n[k]))   # Explicit method
        #S[k] = (0.5 * np.sum(QC) - np.sum(Kn_k*n[k]))
        S[k] = (0.5 * np.sum(np.multiply(Q, C[k])) - np.sum(Q[k])) + \
               (0.5 * np.sum(np.multiply(Q_frag, SS[k])) - np.sum(Q_frag[k]))

        J[k] = (np.dot(KC_k, n) - n[k]*K[k]) + \
               (np.dot(LS_k, n) - n[k]*L[k])
        J[k][k] -= np.sum(Kn_k)

    dn = np.linalg.solve((np.identity(size)/dt - J), S)
    #dn = (S + dnk_frag) * dt 
    n += dn  # (u.s * u.m**3)**(-1)
    t += dt

    if t > save_time.value:
        i += 1
        save_time *= eps_time
        n_save = np.append(n_save, n)
        y = np.append(y, np.argmax(n*m) )
        y2 = np.append(y2, np.argmax(n))
        x = np.append(x,t)
        if i % 10==0:
            print("Time left:", (t_end-t)/(60*60*24*365), n[0], "Stepping: ", dt/(60*60*24*365))


print("Lost mass to total ratio : ")
print(((sum(n*m)-rho_d.value)/rho_d.value))

fig, ax = plt.subplots(1,2)
ax[1].bar(range(0, len(n)), n*m)
ax[1].set_title('Density (n*m)')
ax[0].bar(range(0, len(n)), n)
ax[0].set_title('Number density (n)')

n_save = np.reshape(n_save, (i+1,len(n)))

# %%

x = x/(60*60*24*365)
y = a[y.astype(int)]
y2 = a[y2.astype(int)]
fig, ax = plt.subplots(1, 2)
ax[0].loglog(x,y2)
ax[1].loglog(x,y)
ax[0].grid(True)
ax[0].set_title("Number density (n)")
ax[1].set_title("Density (n*m)")
ax[1].grid(True)

fig.text(0.5, 0.04, r"t [yr]", fontsize="16", ha='center')
fig.text(0.04, 0.5, r"peak position [cm]", fontsize="16", va='center', rotation='vertical')


# %%

t_space = np.linspace(0, t_end, 10000) * u.s

a_Brown_anal = ((5/2 * rho_d/(rho_s * np.pi) * np.sqrt(12*const.k_B * T/rho_s) * t_space + a[0]**(5/2) ) **(2/5)).to(u.cm)
#a_Turb_anal = a[0] * np.exp(sigma_d / (sigma_d*100) * np.pi/np.sqrt(2) * t_space.value)
#a_Turb_anal = sigma_d * Omega_K/ ( 2*rho_s * np.sqrt(np.pi)) * t_space + a[0]


# %%

fig, ax = plt.subplots()
ax.loglog(t_space.to(u.yr),a_Brown_anal)
#ax.loglog(t_space.to(u.yr),a_Turb_anal)

ax.loglog(x,y)
ax.grid(True)

# %%
%%time
# Fragmenting and cratering
np.set_printoptions(threshold=sys.maxsize)
