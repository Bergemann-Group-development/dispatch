FUNCTION file, data, id, no
  return, data+'/'+str(id,form='(i4.4)')+'_'+str(no,form='(i4.4)')+'.dat'
END

; For a staggered quantity the physical boundary lies at x[2] (IDL index) on
; the left hand side, and at x[m-2] on the right hand side

PRO plot_uxm1
ws,2
no=10
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

i0=4

j=6
id=94
for i=i0,i0+1 do begin
  ;print,file(data,id,no)
  id=94+i*27
  open,a,file(data,id,no),dim=[1,1,1]*36,/quiet
  ddn=exp(xdn(alog(a[0])))
  ux=a[1:*,17,17,2]/ddn[1:*,17,17]
  print,ux
  if i eq i0 then begin
    plot,-x[1:*,i],ux,psym=-j,xr=-[x[29,i0],x[12,i0+1]]
  end else begin
    oplot,-x[1:*,i],ux,psym=-j
  end
  oplot,-[1,1]*x[2,i],[-1,1],line=1
  j=j-5
end

END