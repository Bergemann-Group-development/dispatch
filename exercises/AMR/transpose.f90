real:: a(4,3,2)

do iz=1,2
do iy=1,3
  do ix=1,4
    a(ix,iy,iz) = ix+10*iy+100*iz
  end do
  print*,(a(ix,iy,iz), ix=1,4)
end do
print*,''
end do

open (1,file='transpose.dat',access='direct',recl=sizeof(a),status='unknown')
write (1,rec=1) a
close (1)

END
