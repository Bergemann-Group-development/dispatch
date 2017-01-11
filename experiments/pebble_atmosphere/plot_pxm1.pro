FUNCTION file, data, id, no
  return, data+'/'+str(id,form='(i4.4)')+'_'+str(no,form='(i4.4)')+'.dat'
END

; For a staggered quantity the physical boundary lies at x[2] (IDL index) on
; the left hand side, and at x[m-2] on the right hand side

PRO plot_pxm1
ws,2
no=1
v=2
data='data'

x=fltarr(36,7)
x[*,6]=findgen(36)-17.5
x[*,5]=x[*,6]
x[*,5]=x[*,5]-x[34,5]+x[2,6]
for i=4,0,-1 do begin
  x[*,i]=x[*,i+1]*3
  x[*,i]=x[*,i]-x[34,i]+x[2,i+1]
end

i0=0
min=1e9
max=-min
for i=i0,i0+1 do begin
  if i eq i0 then id=67+i*27
  if id lt 220 then id=id+27 else id=id+1
  print,'id:',id
  open,a,file(data,id,no),dim=[1,1,1]*36,/quiet
  px=a[*,17,17,2]
  max=max>max(px)
  min=min<min(px)
end

j=6
for i=i0,i0+1 do begin
  if i eq i0 then id=67+i*27
  if id lt 220 then id=id+27 else id=id+1
  print,file(data,id,no)
  open,a,file(data,id,no),dim=[1,1,1]*36,/quiet
  px=a[*,17,17,2]
  if i eq i0 then begin
    plot,-x[*,i],px,psym=-j,xr=-[x[29,i0],x[12,i0+1]],yr=[min,max]
  end else begin
    oplot,-x[*,i],a[*,17,17,v],psym=-j
    oplot,-x[*,i],a[*,18,17,v],psym=-j
    oplot,-x[*,i],a[*,17,18,v],psym=-j
    oplot,-x[*,i],a[*,18,18,v],psym=-j
  end
  oplot,-[1,1]*x[2,i],[-1,1],line=1
  j=j-5
end

END