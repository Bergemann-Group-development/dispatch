openr, lun, 'data/single_rk3.dat', /get_lun, /f77

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

r = fltarr(3,mrec,n)
v = fltarr(3,mrec,n)
t = fltarr(mrec,n)
m = fltarr(mrec,n)
s= fltarr(8,n)
ss = fltarr(8)
i=0l
while (not eof(lun)) do begin
  readu, lun, n
  for j=0,n-1 do begin
    readu, lun, ss
    s[*,j]=ss
  endfor
  t[i,*]=s[0,*] & r[0:2,i,*]=s[1:3,*] & v[0:2,i,*]=s[4:6,*] & m[i,*]=s[7,*]
  i++
endwhile

free_lun, lun

window, xsize=800, ysize=800
size=10.
plot, [-1,1]*size, [-1,1]*size, /nodata
for j=0,n-1 do begin
  oplot, r[0,*,j], r[1,*,j], psym=3 ;, symsize=.2
end

END
