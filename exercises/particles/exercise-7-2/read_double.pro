;openr, lun, 'data/variable9_rk4.dat', /get_lun, /f77
openr, lun, 'data/fixed10_rk4.dat', /get_lun, /f77

n = 0l
mrec=0l
ntot=0l
while (not eof(lun)) do begin
  readu, lun, n
  for i=1l,n do readu, lun
  mrec++
  ntot+=n
endwhile

rewind, lun
print, 'ntot, mrec :', ntot, mrec

r = dblarr(3,mrec,n)
v = dblarr(3,mrec,n)
t = dblarr(mrec,n)
m = dblarr(mrec,n)
s = dblarr(7,n)
tt = 0.0
ss = dblarr(7)
i=0l
while (not eof(lun)) do begin
  readu, lun, n
  for j=0,n-1 do begin
    readu, lun, tt, ss
    s[*,j]=ss
  endfor
  t[i,*]=tt & r[0:2,i,*]=s[0:2,*] & v[0:2,i,*]=s[3:5,*] & m[i,*]=s[6,*]
  i++
endwhile

free_lun, lun

;window, xsize=800, ysize=800
size=10.
;plot, [-1,1]*size, [-1,1]*size, /nodata
for j=0,n-1 do begin
  ;oplot, r[0,*,j], r[1,*,j], psym=3 ;, symsize=.2
end

END
