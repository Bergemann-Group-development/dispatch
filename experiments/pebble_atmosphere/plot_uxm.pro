FUNCTION file, data, id, no
  return, data+'/'+str(id,form='(i4.4)')+'_'+str(no,form='(i4.4)')+'.dat'
END

; For a staggered quantity the physical boundary lies at x[2] (IDL index) on
; the left hand side, and at x[m-2] on the right hand side

PRO plot_uxm

window,2,xsize=1500,ys=600
no=100
v=2
data='data/peb01'

x=fltarr(36,7)
x[*,6]=findgen(36)-17.5
x[*,5]=x[*,6]
x[*,5]=x[*,5]-x[34,5]+x[2,6]
for i=4,0,-1 do begin
  x[*,i]=x[*,i+1]*3
  x[*,i]=x[*,i]-x[34,i]+x[2,i+1]
end

plot,-[x[0,0],x[12,6]],[0,1]*1,/xlog,/nodat
j=0
id=94
for i=0,5 do begin
  ;print,file(data,id,no)
  open,a,file(data,id,no),dim=[1,1,1]*36,/quiet
  ddn=exp(xdn(alog(a[0])))
  oplot,-x[1:*,i],-a[1:*,17,17,v]/ddn[1:*,17,17],line=j
  oplot,-[1,1]*x[2,i],[-1,1],line=1
  id=id+27
end
id=id-26
open,a,file(data,id,no),dim=[1,1,1]*36,/quiet
ddn=exp(xdn(alog(a[0])))
oplot,-x[1:*,i],-a[1:*,17,17,v]/ddn[1:*,17,17],line=j
oplot,-[1,1]*x[2,i],[-1,1],line=1
j=j+2

END