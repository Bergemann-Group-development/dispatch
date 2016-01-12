openr, lun, 'data/particles.dat', /get_lun, /f77

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

r = dblarr(3,mrec*n)
v = dblarr(3,mrec*n)
t = dblarr(mrec*n)  
m = dblarr(mrec*n)
s= dblarr(8,n)
ss = dblarr(8)
i=0l
while (not eof(lun)) do begin
  readu, lun, n
  for j=0,n-1 do begin
    readu, lun, ss
    s[*,j]=ss
  endfor
  t[i:i+n-1]=s[0,*] & r[*,i:i+n-1]=s[1:3,*] & v[*,i:i+n-1]=s[4:6,*] & m[i:i+n-1]=s[7,*]
  i++ 
endwhile

free_lun, lun

plot, r[0,*], r[1,*], ps=3

END
