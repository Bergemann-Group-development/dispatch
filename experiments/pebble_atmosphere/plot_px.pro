FUNCTION file, id, no
  return, 'data/'+str(id,form='(i4.4)')+'_'+str(no,form='(i4.4)')+'.dat'
END

PRO plot_px, no, id, m

default,id,230
default,no,1
default,m,36
n=m-4
v=2

open,a,file(id,no),dim=[1,1,1]*m,/quiet
x=findgen(m)-(m-1.)/2.
plot,x,a[*,m/2,m/2,v],xr=[-20,20+n],psy=-1

open,a,file(id+1,no),dim=[1,1,1]*m,/quiet
oplot,x+n,a[*,m/2,m/2,v],psy=-6

END