pro plot_phase_space, i, w
default,i,30
if n_elements(w) gt 0 then wset,w
p=rp(i)
plot, (p.r[*,2,0]+500) mod 48, p.p[*,2,0], psym=3, xst=3, title=str(i)
END