FUNCTION file, data, id, no
  return, data+'/'+str(id,form='(i4.4)')+'_'+str(no,form='(i4.4)')+'.dat'
END

; For a centered quantity the physical boundary lies at x[1.5] (IDL index) on
; the left hand side, and at x[m-2.5] on the right hand side

PRO plot_dm
ws,2
no=0
v=1
data='data/peb01'

x=fltarr(36,7)
x[*,6]=findgen(36)-17.5
x[*,5]=x[*,6]
x[*,5]=x[*,5]-0.5*(x[33,5]+x[34,5])+0.5*(x[1,6]+x[2,6])
for i=4,0,-1 do begin
  x[*,i]=x[*,i+1]*3
  x[*,i]=x[*,i]-0.5*(x[33,i]+x[34,i])+0.5*(x[1,i+1]+x[2,i+1])
end

plot,-[x[0,0],x[12,6]],[1e-14,1e-6],/xlog,/ylog,/nodat
j=0
for no=0,100,100 do begin
  id=94
  for i=0,5 do begin
    ;print,file(data,id,no)
    open,a,file(data,id,no),dim=[1,1,1]*36,/quiet
    oplot,-x[*,i],a[*,17,17,v],line=j
    oplot,-[1,1]*0.5*(x[1,i]+x[2,i]),[1e-20,1],line=1
    id=id+27
  end
  id=id-26
  open,a,file(data,id,no),dim=[1,1,1]*36,/quiet
  oplot,-x[*,i],a[*,17,17,v],line=j
  oplot,-[1,1]*0.5*(x[1,i]+x[2,i]),[1e-20,1],line=1
  j=j+2
end

END