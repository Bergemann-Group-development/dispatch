FUNCTION file, id, no, data
  return, str(data)+'/'+str(id,form='(i4.4)')+'_'+str(no,form='(i4.4)')+'.dat'
END

PRO plot_d, no, id, m, data=data

default,id,230
default,no,1
default,m,36
default,data,'data'
n=m-4

open,a,file(id,no,data),dim=[1,1,1]*m,/quiet
x=findgen(m)-(m-1.)/2.
plot,x,a[*,m/2,m/2,0],/ylog,yr=[1e-12,1e-8],xr=[-20,20+n],psy=-1

open,a,file(id+1,no,data),dim=[1,1,1]*m,/quiet
oplot,x+n,a[*,m/2,m/2,0],psy=-1

END
