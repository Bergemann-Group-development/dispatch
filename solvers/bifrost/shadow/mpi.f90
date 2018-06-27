MODULE parallel_mod

! **
! ** IMPORTANT NOTICE:
! ** This module cannot use the io-module, because that would cause
! ** circular dependencies in the compilation process.
! **

 USE trace_mod
 USE params_mod

 USE mpi
 !%
 logical :: io_master,io_mesh
!$omp threadprivate (io_master,io_mesh)
 logical :: exist_xup,exist_yup,exist_zup
!$omp threadprivate (exist_xup,exist_yup,exist_zup)
 logical :: exist_xdn,exist_ydn,exist_zdn
!$omp threadprivate (exist_xdn,exist_ydn,exist_zdn)
 logical :: exist_xupyup,exist_xupydn,exist_xdnyup,exist_xdnydn
!$omp threadprivate (exist_xupyup,exist_xupydn,exist_xdnyup,exist_xdnydn)
 logical :: exist_xupzup,exist_xupzdn,exist_xdnzup,exist_xdnzdn
!$omp threadprivate (exist_xupzup,exist_xupzdn,exist_xdnzup,exist_xdnzdn)
 logical :: exist_yupzup,exist_yupzdn,exist_ydnzup,exist_ydnzdn
!$omp threadprivate (exist_yupzup,exist_yupzdn,exist_ydnzup,exist_ydnzdn)
 logical :: exist_xupyupzup,exist_xupyupzdn
!$omp threadprivate (exist_xupyupzup,exist_xupyupzdn)
 logical :: exist_xupydnzup,exist_xupydnzdn
!$omp threadprivate (exist_xupydnzup,exist_xupydnzdn)
 logical :: exist_xdnyupzup,exist_xdnyupzdn
!$omp threadprivate (exist_xdnyupzup,exist_xdnyupzdn)
 logical :: exist_xdnydnzup,exist_xdnydnzdn
!$omp threadprivate (exist_xdnydnzup,exist_xdnydnzdn)
 logical,dimension(3) :: periodic
!$omp threadprivate (periodic)
 logical :: init_bcgrid=.false.
!$omp threadprivate (init_bcgrid)
 logical :: reorder=.true.
!$omp threadprivate (reorder)
 logical :: manual_cuts=.false. ! if .true. dim will be read from the mhd.in otherwise dims_create3 will be called
!$omp threadprivate (manual_cuts)
 integer(kind=4) :: ierr,size_a_buf
!$omp threadprivate (ierr,size_a_buf)
 integer(kind=4) :: ngridpoints,gridsize,nproc,ndim
!$omp threadprivate (ngridpoints,gridsize,nproc,ndim)
 integer(kind=4) :: mycolid,myhorid,ncol,nhor
!$omp threadprivate (mycolid,myhorid,ncol,nhor)
 integer(kind=4) :: xsb,xs,xe,xeb
!$omp threadprivate (xsb,xs,xe,xeb)
 integer(kind=4) :: ysb,ys,ye,yeb
!$omp threadprivate (ysb,ys,ye,yeb)
 integer(kind=4) :: zsb,zs,ze,zeb
!$omp threadprivate (zsb,zs,ze,zeb)
 integer(kind=4) :: zsbtmp,zstmp,zetmp,zebtmp
!$omp threadprivate (zsbtmp,zstmp,zetmp,zebtmp)
 integer(kind=4) :: zsbold,zsold,zeold,zebold
!$omp threadprivate (zsbold,zsold,zeold,zebold)
 integer(kind=4) :: xr,xrb,yr,yrb,zr,zrb
!$omp threadprivate (xr,xrb,yr,yrb,zr,zrb)
 integer(kind=4) :: R_CART,C_CART,R_MASTER
!$omp threadprivate (R_CART,C_CART,R_MASTER)
 integer(kind=4) :: MASTER_UB,MASTER_LB,COMM_LB,COMM_UB
!$omp threadprivate (MASTER_UB,MASTER_LB,COMM_LB,COMM_UB)

 integer(kind=4) :: R_XUP,R_XDN,R_YUP,R_YDN,R_ZUP,R_ZDN
!$omp threadprivate (R_XUP,R_XDN,R_YUP,R_YDN,R_ZUP,R_ZDN)
 integer(kind=4) :: R_XUPYDN,R_XUPYUP,R_XUPZUP,R_XUPZDN
!$omp threadprivate (R_XUPYDN,R_XUPYUP,R_XUPZUP,R_XUPZDN)
 integer(kind=4) :: R_XDNYDN,R_XDNYUP,R_XDNZUP,R_XDNZDN
!$omp threadprivate (R_XDNYDN,R_XDNYUP,R_XDNZUP,R_XDNZDN)
 integer(kind=4) :: R_YUPZUP,R_YUPZDN,R_YDNZUP,R_YDNZDN
!$omp threadprivate (R_YUPZUP,R_YUPZDN,R_YDNZUP,R_YDNZDN)
 integer(kind=4) :: R_XUPYUPZUP,R_XUPYUPZDN,R_XUPYDNZUP,R_XUPYDNZDN
!$omp threadprivate (R_XUPYUPZUP,R_XUPYUPZDN,R_XUPYDNZUP,R_XUPYDNZDN)
 integer(kind=4) :: R_XDNYUPZUP,R_XDNYUPZDN,R_XDNYDNZUP,R_XDNYDNZDN
!$omp threadprivate (R_XDNYUPZUP,R_XDNYUPZDN,R_XDNYDNZUP,R_XDNYDNZDN)

 integer(kind=4),dimension(3) :: dims,coord
!$omp threadprivate (dims,coord)
 integer(kind=4),dimension(:),allocatable :: x0,xmin,xmax,x1
!$omp threadprivate (x0,xmin,xmax,x1)
 integer(kind=4),dimension(:),allocatable :: y0,ymin,ymax,y1
!$omp threadprivate (y0,ymin,ymax,y1)
 integer(kind=4),dimension(:),allocatable :: z0,zmin,zmax,z1
!$omp threadprivate (z0,zmin,zmax,z1)
 integer(kind=4),dimension(:),allocatable :: xra,yra,zra
!$omp threadprivate (xra,yra,zra)
 integer(kind=4),dimension(:,:),allocatable :: griddat
!$omp threadprivate (griddat)

 integer(kind=4),dimension(MPI_STATUS_SIZE) :: status
!$omp threadprivate (status)
 integer(KIND=4) :: OSC_REAL, OSC_DOUBLE_PRECISION, OSC_INTEGER
!$omp threadprivate (OSC_REAL,OSC_DOUBLE_PRECISION,OSC_INTEGER)
 integer(KIND=4) :: OSC_2INTEGER,OSC_2REAL,OSC_2DOUBLE_PRECISION
!$omp threadprivate (OSC_2INTEGER,OSC_2REAL,OSC_2DOUBLE_PRECISION)

 integer(kind=4) :: T_FARRAY_REAL4,T_IOARRAY_REAL4
!$omp threadprivate (T_FARRAY_REAL4,T_IOARRAY_REAL4)
 integer(kind=4) :: T_LARRAY_REAL4,T_DARRAY_REAL4
!$omp threadprivate (T_LARRAY_REAL4,T_DARRAY_REAL4)
 integer(kind=4) :: T_SENDXUP_REAL4,T_RECVXUP_REAL4,T_SENDXDN_REAL4,T_RECVXDN_REAL4
!$omp threadprivate (T_SENDXUP_REAL4,T_RECVXUP_REAL4,T_SENDXDN_REAL4,T_RECVXDN_REAL4)
 integer(kind=4) :: T_SENDYUP_REAL4,T_RECVYUP_REAL4,T_SENDYDN_REAL4,T_RECVYDN_REAL4
!$omp threadprivate (T_SENDYUP_REAL4,T_RECVYUP_REAL4,T_SENDYDN_REAL4,T_RECVYDN_REAL4)
 integer(kind=4) :: T_SENDZUP_REAL4,T_RECVZUP_REAL4,T_SENDZDN_REAL4,T_RECVZDN_REAL4
!$omp threadprivate (T_SENDZUP_REAL4,T_RECVZUP_REAL4,T_SENDZDN_REAL4,T_RECVZDN_REAL4)
 integer(kind=4) :: T_DATA,T_IOLXYPLANE
!$omp threadprivate (T_DATA,T_IOLXYPLANE)
 integer,dimension(:),allocatable :: a_buf
!$omp threadprivate (a_buf)
 integer,dimension(2000) :: b_sending_request
!$omp threadprivate (b_sending_request)
 integer(KIND=4),dimension(2000,MPI_STATUS_SIZE) :: b_sending_status
!$omp threadprivate (b_sending_status)
 integer :: b_sending_counter
!$omp threadprivate (b_sending_counter)
 real,dimension(:),allocatable :: mpi_buf
!$omp threadprivate (mpi_buf)

 ! column communicator, rank, and communicator size
 integer(kind=4) :: C_COL,R_COL,S_COL
!$omp threadprivate (C_COL,R_COL,S_COL)
 integer(kind=4) :: C_HOR,R_HOR,S_HOR
!$omp threadprivate (C_HOR,R_HOR,S_HOR)

 character(len=72) :: id_mpi= &

     &"$Id$"
!$omp threadprivate (id_mpi)

 interface update
    module procedure update2,update4,update8,updatec !,update16
 end interface update

 interface update_x
    module procedure update2_x,update4_x,update8_x,updatec_x !,update16_x
 end interface update_x

 interface update_xy
    module procedure update2_xy,update4_xy,update8_xy,updatec_xy !,update16_xy
 end interface update_xy

interface update_y
    module procedure update2_y,update4_y,update8_y,updatec_y !,update16_y
 end interface update_y

interface update_z
    module procedure update2_z,update4_z,update8_z,updatec_z !,update16_z
 end interface update_z

interface update2d_xy
    module procedure update2d2_xy,update2d4_xy,update2d8_xy,updatec2d_xy !,update16_xy
 end interface update2d_xy

interface comm_x
   module procedure comm2_x,comm4_x2,comm8_x,commc_x
end interface comm_x

interface comm_y
   module procedure comm2_y,comm4_y2,comm8_y,commc_y
end interface comm_y

interface comm_z
   module procedure comm2_z,comm4_z2,comm8_z,commc_z
end interface comm_z

interface comm_xy
   module procedure comm2_xy,comm4_xy,comm8_xy,commc_xy
end interface comm_xy

interface comm_xz
   module procedure comm2_xz,comm4_xz,comm8_xz,commc_xz
end interface comm_xz

interface comm_yz
   module procedure comm2_yz,comm4_yz,comm8_yz,commc_yz
end interface comm_yz

interface comm_xyz
   module procedure comm2_xyz,comm4_xyz,comm8_xyz,commc_xyz
end interface comm_xyz

interface comm2d_xy
   module procedure comm2d2_xy,comm2d4_xy,comm2d8_xy,commc2d_xy
end interface comm2d_xy

interface adjust_array2
   module procedure adjust_array2_4,adjust_array2_8
end interface adjust_array2

Contains

 subroutine init_parallel

!  USE params_mod

  implicit none

  integer :: R_MASTER_CART,int_test,intkind,realkind,doublekind,i,j,xx,yy,zz
  real :: real_test
  double precision :: double_test
  integer,dimension(:),allocatable :: arr1,arr2
  logical,dimension(3) :: belongs
  integer,dimension(3) :: sizeb,sizel,sizes
  integer(Kind=MPI_ADDRESS_KIND) :: llb=0
  integer(kind=4) :: extent,extentsdn,extentsup,extentrup,extentrdn,m,k,nsize,nr,nn
! ** set default values here
  periodic(:) = (/.true.,.true.,.false./)
  ndim = 3

  call trace%begin ('mpi.f90::init_parallel')

  if (firstit) then

     call MPI_INIT(ierr)

     ! ** Getting compiler define size of reals, doubles etc

     intkind=kind(int_test)
     realkind=kind(real_test)
     doublekind=kind(double_test)

     if (intkind == 2) then
        OSC_INTEGER=MPI_INTEGER2
     elseif (intkind == 4) then
        OSC_INTEGER=MPI_INTEGER4
     elseif (intkind == 8) then
        OSC_INTEGER=MPI_INTEGER8
     else
        write (*,*) 'Warning : integer kind is not 2,4 or 8, will most likely crash MPI'
     endif
     nn=2
     call MPI_TYPE_CONTIGUOUS(nn,MPI_INTEGER,OSC_2INTEGER,ierr)
     call MPI_TYPE_COMMIT(OSC_2INTEGER,ierr)

     if (realkind == 2) then

        !    OSC_REAL=MPI_REAL2
     elseif (realkind == 4) then
        OSC_REAL=MPI_REAL4
     elseif (realkind == 8) then
        OSC_REAL=MPI_REAL8
     else
        write (*,*) 'Warning : real kind is not 2,4 or 8, will most likely crash MPI'
     endif
     nn=2
     call MPI_TYPE_CONTIGUOUS(nn,OSC_REAL,OSC_2REAL,ierr)
     call MPI_TYPE_COMMIT(OSC_2REAL,ierr)

     if (doublekind == 4) then
        OSC_DOUBLE_PRECISION=MPI_REAL4
     elseif (doublekind == 8) then
        OSC_DOUBLE_PRECISION=MPI_REAL8
     elseif (doublekind == 16) then
        OSC_DOUBLE_PRECISION=MPI_DOUBLE_COMPLEX
        write(*,*) 'The definition of data-type: double precision is 16 bits, that might produce problems'
     else
        write (*,*) 'Warning : double kind is not 4,8 or 16, will most likely crash MPI'
     endif
     nn=2
     call MPI_TYPE_CONTIGUOUS(nn,OSC_DOUBLE_PRECISION,OSC_2DOUBLE_PRECISION,ierr)
     call MPI_TYPE_COMMIT(OSC_2DOUBLE_PRECISION,ierr)

     ! ** Getting CPU id number (myid) and number of processors

     call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
     call MPI_COMM_SIZE(MPI_COMM_WORLD,gridsize,ierr)

! ** Setting processor 0 equal to master
     R_MASTER=0

     if (myid.eq.R_MASTER) then
        write (*,*)
        write (*,FMT=ilo) 'mpi     '
        write (*,*) id_mpi
     end if
  end if
  reread = 0
  call read_params("parallel")

  if (manual_cuts .and. (product(dims) .ne. gridsize)) then
    if (myid == R_MASTER) write(*,*) '[Error] manual grid cuts do not match gridsize'
    if (myid == R_MASTER) write(*,*) '[Error] manual_cuts will be set to .false.'
    manual_cuts = .false.
  endif

  if (firstit) then
     if((gridsize.gt.1)) then
        if (.not. manual_cuts) call dims_create3
     else
        dims=1
     endif
  endif

  if (myid == R_MASTER) then
     if ((firstit).or.(reread.eq.1)) then
        write (*,FMT=icol) 'periodic_x',periodic(1)
        write (*,FMT=icol) 'periodic_y',periodic(2)
        write (*,FMT=icol) 'periodic_z',periodic(3)
        write (*,FMT=icol) 'manual_cuts',manual_cuts
        write (*,FMT=icoi) 'ndim',ndim
        write (*,FMT=icoi) 'gridsize',gridsize
        write (*,FMT=icoi3x) 'dims',dims
        write (*,FMT=icol) 'reorder',reorder
     endif
  endif

  if (firstit) then
     call MPI_CART_CREATE(MPI_COMM_WORLD,3,dims,periodic,reorder, &
          C_CART,ierr)
     call MPI_COMM_RANK(C_CART,R_CART,ierr)

     ! ** Resetting everything to only use C_CART instead of MPI_COMM_WORLD

     call MPI_COMM_SIZE(C_CART,nproc,ierr)
     !  print*,'testing size',nproc

     if (myid.eq.R_MASTER) then
        R_MASTER_CART = R_CART
     endif
     call MPI_BCAST(R_MASTER_CART,1,MPI_INTEGER,R_MASTER,MPI_COMM_WORLD,ierr)
     R_MASTER = R_MASTER_CART
     myid = R_CART

     call MPI_CART_COORDS(C_CART,R_CART,3,coord,ierr)

     !  allocate(status_arr(MPI_STATUS_SIZE,0:gridsize))

     ! ** Send this info to master **

     allocate(griddat(9,0:gridsize-1))

     call MPI_ALLGATHER(coord,3,MPI_INTEGER, &
          griddat(1:3,:),3,MPI_INTEGER,C_CART,ierr)

     call find_neighb

     call setup_grid

     call MPI_ALLGATHER((/xs,xe,ys,ye,zs,ze/),6,MPI_INTEGER, &
          griddat(4:9,:),6,MPI_INTEGER,C_CART,ierr)


     ! initialize column communicator
     belongs(:) = (/.false.,.false.,.true./)
     call MPI_Cart_sub(C_CART, belongs, C_COL, ierr)
     call MPI_Comm_size(C_COL,S_COL,ierr)
     call MPI_Comm_rank(C_COL,R_COL,ierr)
     call MPI_comm_size(C_COL,ncol,ierr)
     call MPI_comm_rank(C_COL,mycolid,ierr)

     ! initialize horizontal communicator
     belongs(:) = (/.true.,.true.,.false./)
     call MPI_Cart_sub(C_CART, belongs, C_HOR, ierr)
     call MPI_Comm_size(C_HOR,S_HOR,ierr)
     call MPI_Comm_rank(C_HOR,R_HOR,ierr)
     call MPI_comm_size(C_HOR,nhor,ierr)
     call MPI_comm_rank(C_HOR,myhorid,ierr)


     !  deallocate(status)

     !***  do i=0,gridsize-1
     !***    if (myid.eq.i) then
     !***      print*,'X :',myid,xsb,xs,xe,xeb
     !***      print*,'Y :',myid,ysb,ys,ye,yeb
     !***      print*,'Z :',myid,zsb,zs,ze,zeb
     !***    end if
     !***    call MPI_BARRIER(C_CART,ierr)
     !***  end do
     !***
     !***  call MPI_BARRIER(C_CART,ierr)
     !***
     !***  if (myid == R_MASTER) then
     !***
     !***    do i=0,gridsize-1
     !***      print*,' '
     !***      print*,'Master X:',i,x0(i),xmin(i),xmax(i),x1(i)
     !***      print*,'Master Y:',i,y0(i),ymin(i),ymax(i),y1(i)
     !***      print*,'Master Z:',i,z0(i),zmin(i),zmax(i),z1(i)
     !***    end do
     !***  end if
     !***  call MPI_BARRIER(C_CART,ierr)

     !  call MPI_PACK_SIZE((xeb-xsb+1)*(yeb-ysb+1)*2*32,OSC_REAL,C_CART,size_a_buf,ierr)
     !  size_a_buf=size_a_buf+2*32*MPI_BSEND_OVERHEAD
     !  allocate(a_buf(size_a_buf))
     !  call MPI_BUFFER_ATTACH(a_buf,size_a_buf,ierr)

     !  print*,myid,'making types',xrb,xr,mb-1
     call barrier(C_CART,ierr)
     ! Full data array - All data for a single variable
     !   Simply a row with nx*ny*nz elements of type OSC_REAL
     nn=mx*my*mz
!     call MPI_TYPE_CONTIGUOUS(nn,OSC_REAL,T_FARRAY_REAL4,ierr)
!     call MPI_TYPE_COMMIT(T_FARRAY_REAL4,ierr)
     ! Creating a section of full data array belonging to each mpi-thread
!     call MPI_TYPE_SIZE(OSC_REAL,extent,ierr)
!     allocate(arr1(xr*yr*zr))
!     arr1(:)=OSC_REAL
!     allocate(arr2(xr*yr*zr))
!     m=0
!     do k=zs-1,ze-1
!        do j=ys-1,ye-1
!           do i=xs-1,xe-1
!              m=m+1
!              arr2(m)=(k*mx*my+j*mx+i)*extent
!           end do
!        end do
!     end do
     ! call MPI_TYPE_INDEXED(xr*yr*zr,arr1,arr2,OSC_REAL,T_DATA,ierr)
     ! call MPI_TYPE_COMMIT(T_DATA,ierr)
     !  call MPI_TYPE_CONTIGIGUOUS(mx*my,OSC_REAL,T_IOGXYPLANE,ierr)
     !  call MPI_TYPE_COMMIT(T_IOGXYPLANE,ierr)
     call MPI_TYPE_CREATE_SUBARRAY(2,(/mx,my/),(/xr,yr/),(/xs-1,ys-1/),MPI_ORDER_FORTRAN,OSC_REAL,T_IOLXYPLANE,ierr)
     call MPI_TYPE_COMMIT(T_IOLXYPLANE,ierr)

     if (myid == R_MASTER) print*,'*************',mx,my,mz,'***********'
     sizeb=(/mx,my,mz/)
     sizel=(/xr,yr,zr/)
     sizes=(/xs,ys,zs/)-1
     call MPI_TYPE_CREATE_SUBARRAY(3,sizeb,sizel,sizes,MPI_ORDER_FORTRAN,OSC_REAL,T_IOARRAY_REAL4,ierr)
     call MPI_TYPE_COMMIT(T_IOARRAY_REAL4,ierr)

     nn=xrb*yrb*zrb
     call MPI_TYPE_CONTIGUOUS(nn,OSC_REAL,T_LARRAY_REAL4,ierr)
     call MPI_TYPE_COMMIT(T_LARRAY_REAL4,ierr)
     call MPI_TYPE_SIZE(T_LARRAY_REAL4,extent,ierr)
     !  if (myid == R_MASTER) print*,'T_LARRAY_REAL4: ',xrb,yrb,zrb,xrb*yrb*zrb,extent/4
     sizeb=(/xrb,yrb,zrb/)
     sizel=(/xr,yr,zr/)
     sizes=(/mb,mb,mb/)
     if (xr == xrb) sizes(1)=0
     if (yr == yrb) sizes(2)=0
     if (zr == zrb) sizes(3)=0
     !  print*,'testing X',myid,ndim,sizeb(1),sizel(1),sizes(1)
     !  print*,'testing Y',myid,ndim,sizeb(2),sizel(2),sizes(2)
     !  print*,'testing Z',myid,ndim,sizeb(3),sizel(3),sizes(3)
     !  call MPI_TYPE_CREATE_SUBARRAY(ndim,sizeb,(/1,yr,1/),(/0,mb,0/),MPI_ORDER_FORTRAN,T_LARRAY_REAL4,T_YCOL,ierr)
     call MPI_TYPE_CREATE_SUBARRAY(ndim,sizeb,sizel,sizes,MPI_ORDER_FORTRAN,OSC_REAL,T_DARRAY_REAL4,ierr)
     llb=0
     call MPI_TYPE_COMMIT(T_DARRAY_REAL4,ierr)
     call MPI_TYPE_SIZE(T_DARRAY_REAL4,extent,ierr)
     !  if (myid == R_MASTER) print*,'data part of array:',xr,yr,zr,xr*yr*zr,extent/4
     ! ****** X UP
     if (exist_xup) then
        !    print*,myid,'making X up types'
        !    call barrier(C_CART,ierr)

        sizel=(/mb,yr,zr/)
        if ((.not.periodic(3)).and.(.not.exist_zdn)) sizel(3)=sizel(3)+mb
        if ((.not.periodic(3)).and.(.not.exist_zup)) sizel(3)=sizel(3)+mb

!    print*,'sending 0 - ',xrb-1,xrb-mb,xrb-2*mb,xrb-2*mb-1+mb,xr-1
        sizes=(/xrb-2*mb,mb,mb/)
        if (.not.exist_zdn) sizes(3)=0
        if (yr == yrb) sizes(2)=0
        if (zr == zrb) sizes(3)=0
        call MPI_TYPE_CREATE_SUBARRAY(ndim,sizeb,sizel,sizes,MPI_ORDER_FORTRAN,OSC_REAL,T_SENDXUP_REAL4,ierr)

        !    print*,'recieving',mb,xrb,-2*mb
        sizes=(/xrb-mb,mb,mb/)
        if ((.not.periodic(3)).and.(.not.exist_zdn)) sizes(3)=0
        if (yr == yrb) sizes(2)=0
        if (zr == zrb) sizes(3)=0
        call MPI_TYPE_CREATE_SUBARRAY(ndim,sizeb,sizel,sizes,MPI_ORDER_FORTRAN,OSC_REAL,T_RECVXUP_REAL4,ierr)
        call MPI_TYPE_COMMIT(T_SENDXUP_REAL4,ierr)
        call MPI_TYPE_COMMIT(T_RECVXUP_REAL4,ierr)
        call MPI_TYPE_SIZE(T_SENDXUP_REAL4,extentsup,ierr)
        call MPI_TYPE_SIZE(T_RECVXUP_REAL4,extentrup,ierr)
     end if
     ! ******* X DN
     if (exist_xdn) then
        !    print*,myid,'making X dn types'
        !    call barrier(C_CART,ierr)
        sizel=(/mb,yr,zr/)
        if ((.not.periodic(3)).and.(.not.exist_zdn)) sizel(3)=sizel(3)+mb
        if ((.not.periodic(3)).and.(.not.exist_zup)) sizel(3)=sizel(3)+mb

        sizes=(/mb,mb,mb/)
        if ((.not.periodic(3)).and.(.not.exist_zdn)) sizes(3)=0
        if (yr == yrb) sizes(2)=0
        if (zr == zrb) sizes(3)=0
        call MPI_TYPE_CREATE_SUBARRAY(ndim,sizeb,sizel,sizes,MPI_ORDER_FORTRAN,OSC_REAL,T_SENDXDN_REAL4,ierr)
        sizes=(/0,mb,mb/)
        if ((.not.periodic(3)).and.(.not.exist_zdn)) sizes(3)=0
        if (yr == yrb) sizes(2)=0
        if (zr == zrb) sizes(3)=0
        call MPI_TYPE_CREATE_SUBARRAY(ndim,sizeb,sizel,sizes,MPI_ORDER_FORTRAN,OSC_REAL,T_RECVXDN_REAL4,ierr)
        call MPI_TYPE_COMMIT(T_SENDXDN_REAL4,ierr)
        call MPI_TYPE_COMMIT(T_RECVXDN_REAL4,ierr)
        call MPI_TYPE_SIZE(T_SENDXDN_REAL4,extentsdn,ierr)
        call MPI_TYPE_SIZE(T_RECVXDN_REAL4,extentrdn,ierr)
     end if
     if (exist_xup.and.exist_xdn) then
        if (extentsdn /= extentrup) print*,myid,'error in sizes send xup recvxdn',extentsdn/4,extentrup/4
        if (extentsup /= extentrdn) print*,myid,'error in sizes send xdn recvxup',extentsup/4,extentrdn/4
     end if

     ! ******* Y UP
     if (exist_yup) then
        !    print*,myid,'making Y up types'
        !    call barrier(C_CART,ierr)

       sizel=(/xrb,mb,zr/)
       if ((.not.periodic(3)).and.(.not.exist_zdn)) sizel(3)=sizel(3)+mb
       if ((.not.periodic(3)).and.(.not.exist_zup)) sizel(3)=sizel(3)+mb
       sizes=(/0,yrb-2*mb,mb/)
        !    if (xr == xrb) sizes(1)=0
       if ((.not.periodic(3)).and.(.not.exist_zdn)) sizes(3)=0
       if (zr == zrb) sizes(3)=0
        call MPI_TYPE_CREATE_SUBARRAY(ndim,sizeb,sizel,sizes,MPI_ORDER_FORTRAN,OSC_REAL,T_SENDYUP_REAL4,ierr)

        sizes=(/0,yrb-mb,mb/)
        if ((.not.periodic(3)).and.(.not.exist_zdn)) sizes(3)=0
        !    if (xr == xrb) sizes(1)=0
        if (zr == zrb) sizes(3)=0
        call MPI_TYPE_CREATE_SUBARRAY(ndim,sizeb,sizel,sizes,MPI_ORDER_FORTRAN,OSC_REAL,T_RECVYUP_REAL4,ierr)
        call MPI_TYPE_COMMIT(T_SENDYUP_REAL4,ierr)
        call MPI_TYPE_COMMIT(T_RECVYUP_REAL4,ierr)
        call MPI_TYPE_SIZE(T_SENDYUP_REAL4,extentsup,ierr)
        call MPI_TYPE_SIZE(T_RECVYUP_REAL4,extentrup,ierr)
     end if
     ! ******* Y DN
     if (exist_ydn) then
        !    print*,myid,'making Y dn types'
        !    call barrier(C_CART,ierr)
       sizel=(/xrb,mb,zr/)
       if ((.not.periodic(3)).and.(.not.exist_zdn)) sizel(3)=sizel(3)+mb
       if ((.not.periodic(3)).and.(.not.exist_zup)) sizel(3)=sizel(3)+mb
       sizes=(/0,mb,mb/)
       if ((.not.periodic(3)).and.(.not.exist_zdn)) sizes(3)=0
        !    if (xr == xrb) sizes(1)=0
        if (zr == zrb) sizes(3)=0
        call MPI_TYPE_CREATE_SUBARRAY(ndim,sizeb,sizel,sizes,MPI_ORDER_FORTRAN,OSC_REAL,T_SENDYDN_REAL4,ierr)
        sizes=(/0,0,mb/)
        if ((.not.periodic(3)).and.(.not.exist_zdn)) sizes(3)=0
        if (zr == zrb) sizes(3)=0
        call MPI_TYPE_CREATE_SUBARRAY(ndim,sizeb,sizel,sizes,MPI_ORDER_FORTRAN,OSC_REAL,T_RECVYDN_REAL4,ierr)
        call MPI_TYPE_COMMIT(T_SENDYDN_REAL4,ierr)
        call MPI_TYPE_COMMIT(T_RECVYDN_REAL4,ierr)
        call MPI_TYPE_SIZE(T_SENDYDN_REAL4,extentsdn,ierr)
        call MPI_TYPE_SIZE(T_RECVYDN_REAL4,extentrdn,ierr)
     end if
     if (exist_yup.and.exist_ydn) then
        if (extentsdn /= extentrup) print*,myid,'error in sizes send yup recvydn',extentsdn/4,extentrup/4
        if (extentsup /= extentrdn) print*,myid,'error in sizes send ydn recvyup',extentsup/4,extentrdn/4
     end if
     ! ******* Z UP
     if (exist_zup) then
        !    print*,myid,'making Z up types'
        !    call barrier(C_CART,ierr)
        sizel=(/xrb,yrb,mb/)
        sizes=(/0,0,zrb-2*mb/)
        call MPI_TYPE_CREATE_SUBARRAY(ndim,sizeb,sizel,sizes,MPI_ORDER_FORTRAN,OSC_REAL,T_SENDZUP_REAL4,ierr)
        sizes=(/0,0,zrb-mb/)
        call MPI_TYPE_CREATE_SUBARRAY(ndim,sizeb,sizel,sizes,MPI_ORDER_FORTRAN,OSC_REAL,T_RECVZUP_REAL4,ierr)
        call MPI_TYPE_COMMIT(T_SENDZUP_REAL4,ierr)
        call MPI_TYPE_COMMIT(T_RECVZUP_REAL4,ierr)
        call MPI_TYPE_SIZE(T_SENDZUP_REAL4,extentsup,ierr)
        call MPI_TYPE_SIZE(T_RECVZUP_REAL4,extentrup,ierr)
     end if
     ! ******* Z DN
     if (exist_zdn) then
        !    print*,myid,'making Z dn types'
        !    call barrier(C_CART,ierr)
        sizel=(/xrb,yrb,mb/)
        sizes=(/0,0,mb/)
        call MPI_TYPE_CREATE_SUBARRAY(ndim,sizeb,sizel,sizes,MPI_ORDER_FORTRAN,OSC_REAL,T_SENDZDN_REAL4,ierr)
        sizes=(/0,0,0/)
        call MPI_TYPE_CREATE_SUBARRAY(ndim,sizeb,sizel,sizes,MPI_ORDER_FORTRAN,OSC_REAL,T_RECVZDN_REAL4,ierr)
        call MPI_TYPE_COMMIT(T_SENDZDN_REAL4,ierr)
        call MPI_TYPE_COMMIT(T_RECVZDN_REAL4,ierr)
        call MPI_TYPE_SIZE(T_SENDZDN_REAL4,extentsdn,ierr)
        call MPI_TYPE_SIZE(T_RECVZDN_REAL4,extentrdn,ierr)
     end if
     if (exist_zup.and.exist_zdn) then
        if (extentsdn /= extentrup) print*,myid,'error in sizes send zup recvzdn',extentsdn/4,extentrup/4
        if (extentsup /= extentrdn) print*,myid,'error in sizes send zdn recvzup',extentsup/4,extentrdn/4
     end if
     !  print*,myid,'done types'
     !  call barrier(C_CART,ierr)

     ! Allocating mpi-buffer (If not needed - delete this part)

     nsize=30*2*(mb*xrb*yrb+mb*xrb*zrb+mb*yrb*zrb)
     allocate(mpi_buf(nsize))
     call MPI_TYPE_SIZE(OSC_REAL,nr,ierr)
     call MPI_BUFFER_ATTACH(mpi_buf,nsize*nr,ierr)

  endif
  call trace%end
end subroutine init_parallel

! ***********************************************************************

 subroutine output_params_parallel (lun,ioerr)

  implicit none

  integer :: lun,ioerr,tmp
  logical :: op

  call trace%begin ('mpi.f90::output_params_parallel')

  inquire (unit=lun,opened=op)

  if (op) then
    write (lun,FMT=iloo) 'parallel'
    write (lun,FMT=idout) id_mpi
    if (periodic(1)) then
      tmp=1
    else
      tmp=0
    endif
    write (lun,FMT=icoi) 'periodic_x',tmp
    if (periodic(2)) then
      tmp=1
    else
      tmp=0
    endif
    write (lun,FMT=icoi) 'periodic_y',tmp
    if (periodic(3)) then
      tmp=1
    else
      tmp=0
    endif
    write (lun,FMT=icoi) 'periodic_z',tmp
    write (lun,FMT=icoi) 'ndim',ndim
    if (reorder) then
      tmp=1
    else
      tmp=0
    endif
    if (manual_cuts) then ! write to files only if relevant
        write (lun,FMT=icoi) 'manual_cuts',manual_cuts
        write (lun,FMT=icoi3x) 'dims',dims
    endif
    write (lun,FMT=icoi) 'reorder',tmp
    ioerr=0
    write (lun,*)
  else
    ioerr=1
  end if

  call trace%end
 end subroutine output_params_parallel

! ***********************************************************************

 subroutine find_neighb

  implicit none

  integer(KIND=4),dimension(3) :: dimstmp

  call trace%begin ('mpi.f90::find_neighb')

  R_XUP=MPI_PROC_NULL
  R_XDN=MPI_PROC_NULL
  R_YUP=MPI_PROC_NULL
  R_YDN=MPI_PROC_NULL
  R_ZUP=MPI_PROC_NULL
  R_ZDN=MPI_PROC_NULL
  R_XUPYUP=MPI_PROC_NULL
  R_XUPYDN=MPI_PROC_NULL
  R_XDNYUP=MPI_PROC_NULL
  R_XDNYDN=MPI_PROC_NULL
  R_XUPZUP=MPI_PROC_NULL
  R_XUPZDN=MPI_PROC_NULL
  R_XDNZUP=MPI_PROC_NULL
  R_XDNZDN=MPI_PROC_NULL
  R_YUPZUP=MPI_PROC_NULL
  R_YUPZDN=MPI_PROC_NULL
  R_YDNZUP=MPI_PROC_NULL
  R_YDNZDN=MPI_PROC_NULL
  R_XUPYUPZUP=MPI_PROC_NULL
  R_XUPYUPZDN=MPI_PROC_NULL
  R_XUPYDNZUP=MPI_PROC_NULL
  R_XUPYDNZDN=MPI_PROC_NULL
  R_XDNYUPZUP=MPI_PROC_NULL
  R_XDNYUPZDN=MPI_PROC_NULL
  R_XDNYDNZUP=MPI_PROC_NULL
  R_XDNYDNZDN=MPI_PROC_NULL

  exist_xup=.false.
  exist_xdn=.false.
  exist_yup=.false.
  exist_ydn=.false.
  exist_zup=.false.
  exist_zdn=.false.

  exist_xupyup=.false.
  exist_xupydn=.false.
  exist_xdnyup=.false.
  exist_xdnydn=.false.
  exist_xupzup=.false.
  exist_xupzdn=.false.
  exist_xdnzup=.false.
  exist_xdnzdn=.false.
  exist_yupzup=.false.
  exist_yupzdn=.false.
  exist_ydnzup=.false.
  exist_ydnzdn=.false.

  exist_xupyupzup=.false.
  exist_xupyupzdn=.false.
  exist_xupydnzup=.false.
  exist_xupydnzdn=.false.
  exist_xdnyupzup=.false.
  exist_xdnyupzdn=.false.
  exist_xdnydnzup=.false.
  exist_xdnydnzdn=.false.

! ** Ranks for xup,xdn

  call MPI_CART_SHIFT(C_CART,0, 1,R_XDN,R_XUP,ierr)

! ** Ranks for yup,ydn

  call MPI_CART_SHIFT(C_CART,1, 1,R_YDN,R_YUP,ierr)

! ** Ranks for zup,zdn

  call MPI_CART_SHIFT(C_CART,2, 1,R_ZDN,R_ZUP,ierr)

! ** Checking which neighbors exist

  if ((R_XUP /= MPI_PROC_NULL).and.(dims(1) /= 1)) then
    exist_xup=.true.
  end if
  if ((R_XDN /= MPI_PROC_NULL).and.(dims(1) /= 1)) then
    exist_xdn=.true.
  end if

  if ((R_YUP /= MPI_PROC_NULL).and.(dims(2) /= 1)) then
    exist_yup=.true.
  end if
  if ((R_YDN /= MPI_PROC_NULL).and.(dims(2) /= 1)) then
    exist_ydn=.true.
  end if

  if ((R_ZUP /= MPI_PROC_NULL).and.(dims(3) /= 1)) then
    exist_zup=.true.
  end if
  if ((R_ZDN /= MPI_PROC_NULL).and.(dims(3) /= 1)) then
    exist_zdn=.true.
  end if


! ** Cant do cart_shift so must find manually. Everything has been set
! ** to PROC_NULL per default

! ** Checking side-corners

  if (exist_xup.and.exist_yup) then
    dimstmp=coord+(/1,1,0/)
    dimstmp(1)=mod(dimstmp(1),dims(1))
    dimstmp(2)=mod(dimstmp(2),dims(2))
    call MPI_CART_RANK(C_CART,dimstmp,R_XUPYUP,ierr)
    exist_xupyup=.true.
  endif

  if (exist_xup.and.exist_ydn) then
    dimstmp=coord+(/1,-1,0/)
    dimstmp(1)=mod(dimstmp(1),dims(1))
    dimstmp(2)=mod(dimstmp(2)+dims(2),dims(2))
    call MPI_CART_RANK(C_CART,dimstmp,R_XUPYDN,ierr)
    exist_xupydn=.true.
  end if

  if (exist_xdn.and.exist_yup) then
    dimstmp=coord+(/-1,1,0/)
    dimstmp(1)=mod(dimstmp(1),dims(1))
    dimstmp(2)=mod(dimstmp(2),dims(2))
    call MPI_CART_RANK(C_CART,dimstmp,R_XDNYUP,ierr)
    exist_xdnyup=.true.
  end if

  if (exist_xdn.and.exist_ydn) then
    dimstmp=coord+(/-1,-1,0/)
    dimstmp(1)=mod(dimstmp(1),dims(1))
    dimstmp(2)=mod(dimstmp(2)+dims(2),dims(2))
    call MPI_CART_RANK(C_CART,dimstmp,R_XDNYDN,ierr)
    exist_xdnydn=.true.
  end if

  if (exist_xup.and.exist_zup) then
    dimstmp=coord+(/1,0,1/)
    dimstmp(1)=mod(dimstmp(1),dims(1))
    dimstmp(3)=mod(dimstmp(3),dims(3))
    call MPI_CART_RANK(C_CART,dimstmp,R_XUPZUP,ierr)
    exist_xupzup=.true.
  endif

  if (exist_xup.and.exist_zdn) then
    dimstmp=coord+(/1,0,-1/)
    dimstmp(1)=mod(dimstmp(1),dims(1))
    dimstmp(3)=mod(dimstmp(3)+dims(3),dims(3))
    call MPI_CART_RANK(C_CART,dimstmp,R_XUPZDN,ierr)
    exist_xupzdn=.true.
  end if

  if (exist_xdn.and.exist_zup) then
    dimstmp=coord+(/-1,0,1/)
    dimstmp(1)=mod(dimstmp(1),dims(1))
    dimstmp(3)=mod(dimstmp(3),dims(3))
    call MPI_CART_RANK(C_CART,dimstmp,R_XDNZUP,ierr)
    exist_xdnzup=.true.
  end if

  if (exist_xdn.and.exist_zdn) then
    dimstmp=coord+(/-1,0,-1/)
    dimstmp(1)=mod(dimstmp(1),dims(1))
    dimstmp(3)=mod(dimstmp(3)+dims(3),dims(3))
    call MPI_CART_RANK(C_CART,dimstmp,R_XDNZDN,ierr)
    exist_xdnzdn=.true.
  end if

  if (exist_yup.and.exist_zup) then
    dimstmp=coord+(/0,1,1/)
    dimstmp(2)=mod(dimstmp(2),dims(2))
    dimstmp(3)=mod(dimstmp(3),dims(3))
    call MPI_CART_RANK(C_CART,dimstmp,R_YUPZUP,ierr)
    exist_yupzup=.true.
  end if

  if (exist_yup.and.exist_zdn) then
    dimstmp=coord+(/0,1,-1/)
    dimstmp(2)=mod(dimstmp(2),dims(2))
    dimstmp(3)=mod(dimstmp(3)+dims(3),dims(3))
    call MPI_CART_RANK(C_CART,dimstmp,R_YUPZDN,ierr)
    exist_yupzdn=.true.
  end if

  if (exist_ydn.and.exist_zup) then
    dimstmp=coord+(/0,-1,1/)
    dimstmp(2)=mod(dimstmp(2)+dims(2),dims(2))
    dimstmp(3)=mod(dimstmp(3),dims(3))
    call MPI_CART_RANK(C_CART,dimstmp,R_YDNZUP,ierr)
    exist_ydnzup=.true.
  endif

  if (exist_ydn.and.exist_zdn) then
    dimstmp=coord+(/0,-1,-1/)
    dimstmp(2)=mod(dimstmp(2)+dims(2),dims(2))
    dimstmp(3)=mod(dimstmp(3)+dims(3),dims(3))
    call MPI_CART_RANK(C_CART,dimstmp,R_YDNZDN,ierr)
    exist_ydnzdn=.true.
  end if

! ** Checking corners

  if (exist_xup.and.exist_yup.and.exist_zup) then
    dimstmp=coord+(/1,1,1/)
    dimstmp(1)=mod(dimstmp(1),dims(1))
    dimstmp(2)=mod(dimstmp(2),dims(2))
    dimstmp(3)=mod(dimstmp(3),dims(3))
    call MPI_CART_RANK(C_CART,dimstmp,R_XUPYUPZUP,ierr)
    exist_xupyupzup=.true.
  end if

  if (exist_xup.and.exist_yup.and.exist_zdn) then
    dimstmp=coord+(/1,1,-1/)
    dimstmp(1)=mod(dimstmp(1),dims(1))
    dimstmp(2)=mod(dimstmp(2),dims(2))
    dimstmp(3)=mod(dimstmp(3)+dims(3),dims(3))
    call MPI_CART_RANK(C_CART,dimstmp,R_XUPYUPZDN,ierr)
    exist_xupyupzdn=.true.
  endif

  if (exist_xup.and.exist_ydn.and.exist_zup) then
    dimstmp=coord+(/1,-1,1/)
    dimstmp(1)=mod(dimstmp(1),dims(1))
    dimstmp(2)=mod(dimstmp(2)+dims(2),dims(2))
    dimstmp(3)=mod(dimstmp(3),dims(3))
    call MPI_CART_RANK(C_CART,dimstmp,R_XUPYDNZUP,ierr)
    exist_xupydnzup=.true.
  end if

  if (exist_xup.and.exist_ydn.and.exist_zdn) then
    dimstmp=coord+(/1,-1,-1/)
    dimstmp(1)=mod(dimstmp(1),dims(1))
    dimstmp(2)=mod(dimstmp(2)+dims(2),dims(2))
    dimstmp(3)=mod(dimstmp(3)+dims(3),dims(3))
    call MPI_CART_RANK(C_CART,dimstmp,R_XUPYDNZDN,ierr)
    exist_xupydnzdn=.true.
  end if

  if (exist_xdn.and.exist_yup.and.exist_zup) then
    dimstmp=coord+(/-1,1,1/)
    dimstmp(1)=mod(dimstmp(1),dims(1))
    dimstmp(2)=mod(dimstmp(2),dims(2))
    dimstmp(3)=mod(dimstmp(3),dims(3))
    call MPI_CART_RANK(C_CART,dimstmp,R_XDNYUPZUP,ierr)
    exist_xdnyupzup=.true.
  end if

  if (exist_xdn.and.exist_yup.and.exist_zdn) then
    dimstmp=coord+(/-1,1,-1/)
    dimstmp(1)=mod(dimstmp(1),dims(1))
    dimstmp(2)=mod(dimstmp(2),dims(2))
    dimstmp(3)=mod(dimstmp(3)+dims(3),dims(3))
    call MPI_CART_RANK(C_CART,dimstmp,R_XDNYUPZDN,ierr)
    exist_xdnyupzdn=.true.
  end if

  if (exist_xdn.and.exist_ydn.and.exist_zup) then
    dimstmp=coord+(/-1,-1,1/)
    dimstmp(1)=mod(dimstmp(1),dims(1))
    dimstmp(2)=mod(dimstmp(2)+dims(2),dims(2))
    dimstmp(3)=mod(dimstmp(3),dims(3))
    call MPI_CART_RANK(C_CART,dimstmp,R_XDNYDNZUP,ierr)
    exist_xdnydnzup=.true.
  end if

  if (exist_xdn.and.exist_ydn.and.exist_zdn) then
    dimstmp=coord+(/-1,-1,-1/)
    dimstmp(1)=mod(dimstmp(1),dims(1))
    dimstmp(2)=mod(dimstmp(2)+dims(2),dims(2))
    dimstmp(3)=mod(dimstmp(3)+dims(3),dims(3))
    call MPI_CART_RANK(C_CART,dimstmp,R_XDNYDNZDN,ierr)
    exist_xdnydnzdn=.true.
  end if

  call trace%end
 end subroutine find_neighb

! ***********************************************************************

 subroutine setup_grid

  implicit none

  integer :: i,ii,j,k,kk,n,xra_tmp,yra_tmp,zra_tmp
  character(len=20) :: fmt,fmti

  call trace%begin ('mpi.f90::setup_grid')

  if (myid.eq.R_MASTER) then

! ** Find subregions for all threads
    allocate(x0(0:gridsize-1))
    allocate(xmin(0:gridsize-1))
    allocate(xmax(0:gridsize-1))
    allocate(x1(0:gridsize-1))

    allocate(y0(0:gridsize-1))
    allocate(ymin(0:gridsize-1))
    allocate(ymax(0:gridsize-1))
    allocate(y1(0:gridsize-1))

    allocate(z0(0:gridsize-1))
    allocate(zmin(0:gridsize-1))
    allocate(zmax(0:gridsize-1))
    allocate(z1(0:gridsize-1))

    allocate(xra(0:gridsize-1))
    allocate(yra(0:gridsize-1))
    allocate(zra(0:gridsize-1))

    xra_tmp=int(mx/(1.0*dims(1))+0.5)
    yra_tmp=int(my/(1.0*dims(2))+0.5)
    zra_tmp=int(mz/(1.0*dims(3))+0.5)

    do i=0,gridsize-1

      xmin(i)=max0(1,xra_tmp*griddat(1,i)+1)
      xmax(i)=min0(mx,xra_tmp*(griddat(1,i)+1))
      if (griddat(1,i).eq.(dims(1)-1)) then
        xmax(i)=mx
      endif

      ymin(i)=max0(1,yra_tmp*griddat(2,i)+1)
      ymax(i)=min0(my,yra_tmp*(griddat(2,i)+1))
      if (griddat(2,i).eq.(dims(2)-1)) then
        ymax(i)=my
      endif

      zmin(i)=max0(1,zra_tmp*griddat(3,i)+1)
      zmax(i)=min0(mz,zra_tmp*(griddat(3,i)+1))
      if (griddat(3,i).eq.(dims(3)-1)) then
        zmax(i)=mz
      endif

      xra(i)=xmax(i)-xmin(i)+1
      yra(i)=ymax(i)-ymin(i)+1
      zra(i)=zmax(i)-zmin(i)+1

! ** If only one cpu in periodic direction, no boundaries are needed

      if ((dims(1) == 1).and.(periodic(1))) then
        x0(i)=xmin(i)
        x1(i)=xmax(i)
      else
        x0(i)=xmin(i)-mb
        x1(i)=xmax(i)+mb
      end if

      if ((dims(2) == 1).and.(periodic(2))) then
        y0(i)=ymin(i)
        y1(i)=ymax(i)
      else
        y0(i)=ymin(i)-mb
        y1(i)=ymax(i)+mb
      end if

      if ((dims(3) == 1).and.(periodic(3))) then
        z0(i)=zmin(i)
        z1(i)=zmax(i)
      else
        z0(i)=zmin(i)-mb
        z1(i)=zmax(i)+mb
      end if

    end do

    xsb=x0(myid)
    xs=xmin(myid)
    xe=xmax(myid)
    xeb=x1(myid)
    xr=xe-xs+1
    xrb=xeb-xsb+1

    ysb=y0(myid)
    ys=ymin(myid)
    ye=ymax(myid)
    yeb=y1(myid)
    yr=ye-ys+1
    yrb=yeb-ysb+1

    zsb=z0(myid)
    zs=zmin(myid)
    ze=zmax(myid)
    zeb=z1(myid)
    zr=ze-zs+1
    zrb=zeb-zsb+1

    if(mx.lt.1000 .and. my.lt.1000 .and. mz.lt.1000) then
      fmt='(/A12,I3," - ",I3)'
    else
      fmt='(/A12,I4," - ",I4)'
    endif
    do i=0,dims(1)-1
      ii=0
      do while (ii.le.(gridsize-1))
        if (griddat(1,ii).eq.i) then
          if (i.eq.0) then
            write(*,fmt,ADVANCE='NO') 'X-Splits : ',xmin(ii),xmax(ii)
          else
            write(*,fmt,ADVANCE='NO') '           ',xmin(ii),xmax(ii)
          end if
          ii=gridsize
        else
          ii=ii+1
        end if
      end do
    end do
    write (*,*) ' '

    do i=0,dims(2)-1
      ii=0
      do while (ii.le.(gridsize-1))
        if (griddat(2,ii).eq.i) then
          if (i.eq.0) then
            write(*,fmt,ADVANCE='NO') 'Y-Splits : ',ymin(ii),ymax(ii)
          else
            write(*,fmt,ADVANCE='NO') '           ',ymin(ii),ymax(ii)
          end if
          ii=gridsize
        else
          ii=ii+1
        end if
      end do
    end do
    write (*,*) ' '

    do i=0,dims(3)-1
      ii=0
      do while (ii.le.(gridsize-1))
        if (griddat(3,ii).eq.i) then
          if (i.eq.0) then
            write(*,fmt,ADVANCE='NO') 'Z-Splits : ',zmin(ii),zmax(ii)
          else
            write(*,fmt,ADVANCE='NO') '           ',zmin(ii),zmax(ii)
          end if
          ii=gridsize
        else
          ii=ii+1
        end if
      end do
    end do

    write (*,*) ' '
    write (*,*) ' '

    if(gridsize < 1000) then
      fmti="(' ',I3)"
    else
      fmti="(' ',I4)"
    endif
    do k=dims(3)-1,0,-1
      do j=dims(2)-1,0,-1
        fmt="(A1)"
        do kk=0,(j-1)*2+5
          write (*,fmt,ADVANCE='NO') " "
        end do
        fmt=fmti
        do i=0,dims(1)-1
          do n=0,gridsize-1
            if (i.eq.(griddat(1,n)).and.(j.eq.griddat(2,n)).and.(k.eq.griddat(3,n))) write (*,fmt,ADVANCE='NO') n
          end do
        end do
        fmt="(A1,/)"
        write (*,fmt) " "
      end do
      write (*,*) " "
    end do

    write (*,*) ' '
    write (*,*) 'Z |    / Y'
    write (*,*) '  |   /'
    write (*,*) '  |  /'
    write (*,*) '  | /'
    write (*,*) '  |_______ X'
    write (*,*) ' '

  else

! ** Find subregion for individual thread

    xr=int(mx/(1.0*dims(1))+0.5)
    yr=int(my/(1.0*dims(2))+0.5)
    zr=int(mz/(1.0*dims(3))+0.5)

    xs=max0(1,xr*coord(1)+1)
    xe=min0(mx,xr*(coord(1)+1))

    if (coord(1) == (dims(1)-1)) xe=mx

    if ((dims(1) == 1).and.(periodic(1))) then
      xsb=xs
      xeb=xe
    else
      xsb=xs-mb
      xeb=xe+mb
    end if

    ys=max0(1,yr*coord(2)+1)
    ye=min0(my,yr*(coord(2)+1))
    if (coord(2) == (dims(2)-1)) ye=my

    if ((dims(2) == 1).and.(periodic(2))) then
      ysb=ys
      yeb=ye
    else
      ysb=ys-mb
      yeb=ye+mb
    end if

    zs=max0(1,zr*coord(3)+1)
    ze=min0(mz,zr*(coord(3)+1))

    if (coord(3) == (dims(3)-1)) ze=mz

    if ((dims(3) == 1).and.(periodic(3))) then
      zsb=zs
      zeb=ze
    else
      zsb=zs-mb
      zeb=ze+mb
    end if

  endif

!*** resetting xr,yr,zr to ranges of internal box

  xr=xe-xs+1
  yr=ye-ys+1
  zr=ze-zs+1

  xrb=xeb-xsb+1
  yrb=yeb-ysb+1
  zrb=zeb-zsb+1

  ngridpoints=xr*yr*zr

  call trace%end
 end subroutine setup_grid

! ********************************************************

subroutine setup_bcgrid

  implicit none
  integer :: i,j,k
  integer, dimension(:), allocatable :: ub_id,lb_id
  integer :: G_CART,G_UB,G_LB
  logical :: debug_boundary=.false.

  if(init_bcgrid)return

  call trace%begin ('mpi.f90::setup_bcgrid')

  allocate(ub_id(dims(1)*dims(2)))
  allocate(lb_id(dims(1)*dims(2)))
  j=0
  k=0
  do i=0,gridsize-1
    if (griddat(3,i).eq.0) then
      j=j+1
      lb_id(j)=i
    endif
    if (griddat(3,i).eq.(dims(3)-1)) then
      k=k+1
      ub_id(k)=i
    endif
  enddo
  if (myid.eq.R_MASTER.and.debug_boundary) then
    print*,'ub',ub_id
    print*,'lb',lb_id
  endif

  call MPI_COMM_GROUP(C_CART,G_CART,ierr)
  call MPI_GROUP_INCL(G_CART,j,lb_id,G_LB,ierr)
  call MPI_GROUP_INCL(G_CART,k,ub_id,G_UB,ierr)
  call MPI_COMM_CREATE(C_CART,G_LB,COMM_LB,ierr)
  call MPI_COMM_CREATE(C_CART,G_UB,COMM_UB,ierr)
  MASTER_UB=minval(ub_id)
  MASTER_LB=minval(lb_id)
  call MPI_GROUP_FREE(G_CART,ierr)
  call MPI_GROUP_FREE(G_LB,ierr)
  call MPI_GROUP_FREE(G_UB,ierr)

  deallocate(ub_id,lb_id)

  init_bcgrid=.true.

  call trace%end
end subroutine setup_bcgrid

! ***********************************************************************

 subroutine barrier (who,id)

  implicit none

  integer(kind=4) :: who,id

  call MPI_BARRIER(who,ierr)

  call trace%end
 end subroutine barrier

! ***********************************************************************

!!$ subroutine dims_create
!!$
!!$  implicit none
!!$
!!$  integer :: counter,i,j,k,ncpus
!!$  integer,dimension(1) :: tmploc
!!$  real,dimension(3) :: weight
!!$  real :: maxw
!!$
!!$! **
!!$! ** Had to redo this built in MPI routine, since it automatically
!!$! ** sets most processors on the first coordinate. Very
!!$! ** unfortunately for the stagger, now puts least cpu's in first
!!$! ** non-periodic direction
!!$! **
!!$
!!$!  dims(:)=0
!!$!  call MPI_DIMS_CREATE(gridsize,3,dims,ierr)
!!$
!!$!
!!$! Have seen mpi_dims_create being very bad at this
!!$!
!!$
!!$  dims(:)=1
!!$  weight(:)=(1.0*(/mx,my,mz/))/minval((/mx,my,mz/))/dims(:)
!!$
!!$  counter=1
!!$  k=0
!!$
!!$  do
!!$! ** Finding location of maximum weight
!!$    if (k == 0) then
!!$      tmploc=maxloc(weight(:))
!!$
!!$      k=tmploc(1)
!!$    end if
!!$
!!$! ** adding one cpu to this direction
!!$
!!$    dims(k)=dims(k)+1
!!$
!!$! ** Counting number of cpus
!!$    ncpus=1
!!$    do j=1,ndim
!!$      ncpus=ncpus*dims(j)
!!$    end do
!!$
!!$! ** Testing if we reached correct result
!!$
!!$    if (ncpus == gridsize) then
!!$!      if (myid == R_MASTER) print*,'Finishing of with ',dims
!!$      exit
!!$! ** If were above value when adding one more cpu in the direction with fewest
!!$! ** then remove extra cpu in direction k and change scaling
!!$    elseif (( gridsize < ((1.0*ncpus)/maxval(dims)*(maxval(dims)+1))).and.(counter < 4)) then
!!$!      if (myid == R_MASTER) print*,'tried',dims,'must correct'
!!$      dims(k)=dims(k)-1
!!$      maxw=0.0
!!$      do i=1,ndim
!!$        if (i == k) cycle
!!$        maxw=max(maxw,weight(i))
!!$      end do
!!$      weight(k)=maxw-1e-10
!!$      counter=counter+1
!!$      k=mod(k,ndim)+1
!!$    elseif ((gridsize == ((1.0*ncpus)/maxval(dims)*(maxval(dims)+1))).and.(counter < 4)) then
!!$      tmploc=maxloc(dims)
!!$      k=tmploc(1)
!!$!      if (myid == R_MASTER) print*,'making sure direction',k,' is chosen'
!!$    elseif (counter < 4) then
!!$!      if (myid == R_MASTER) print*,'Found no problem with ',dims
!!$      weight(k)=weight(k)*(dims(k)-1)/dims(k)-1.0e-10
!!$!      if (myid == R_MASTER) print*,'New weights ',weight
!!$      k=0
!!$      counter=1
!!$    elseif (counter >= 4) then
!!$      if (myid == R_MASTER) print*,'Could not solve problem',dims
!!$      dims(:)=0
!!$      call MPI_DIMS_CREATE(gridsize,ndim,dims,ierr)
!!$      exit
!!$    end if
!!$  end do
!!$
!!$ end subroutine dims_create
!!$
!!$! ***********************************************************************
!!$
!!$ subroutine dims_create2
!!$
!!$  USE params_mod
!!$
!!$  implicit none
!!$
!!$  integer,dimension(3,gridsize) :: cdims
!!$  real,dimension(3,gridsize) :: diff,grprcpu
!!$  real,dimension(gridsize,gridsize,gridsize) :: errval,grprcpudiff
!!$  real :: grprcpuref
!!$  integer :: i,j,k
!!$
!!$! ** checking which combinations produces ok results
!!$! ** grid cells per cpu : grprcpu
!!$! ** how many grid cells per cpu : cdims
!!$! ** number describing the unevenness of point
!!$! **    distribution to each cpu : diff
!!$
!!$  do i=1,gridsize
!!$
!!$    grprcpu(1,i)=(1.0*mx)/(1.0*i)
!!$    cdims(1,i)=int(grprcpu(1,i))
!!$    diff(1,i)=abs((1.0*mx)/(1.0*cdims(1,i))-grprcpu(1,i))
!!$
!!$    grprcpu(2,i)=(1.0*my)/(1.0*i)
!!$    cdims(2,i)=max(int(grprcpu(2,i)),1)
!!$    diff(2,i)=abs((1.0*my)/(1.0*cdims(2,i))-grprcpu(2,i))
!!$
!!$    grprcpu(3,i)=(1.0*mz)/(1.0*i)
!!$    cdims(3,i)=max(int(grprcpu(3,i)),1)
!!$    diff(3,i)=abs((1.0*mz)/(1.0*cdims(3,i))-grprcpu(3,i))
!!$
!!$  end do
!!$
!!$  do i=1,gridsize
!!$    do j=1,gridsize
!!$      do k=1,gridsize
!!$        if ((i+j+k) == 3) then
!!$          errval(i,j,k)=1e30
!!$        else
!!$          errval(i,j,k)=diff(1,i)/grprcpu(1,i)+ &
!!$                        diff(2,j)/grprcpu(2,j)+ &
!!$                        diff(3,k)/grprcpu(3,k)
!!$        end if
!!$        grprcpudiff(i,j,k)=abs(grprcpu(1,i)-grprcpu(2,j))+ &
!!$                           abs(grprcpu(1,i)-grprcpu(3,k))+ &
!!$                           abs(grprcpu(2,j)-grprcpu(3,k))
!!$      end do
!!$    end do
!!$  end do
!!$
!!$
!!$  grprcpuref=maxval(errval)
!!$  errval(1,1,1)=grprcpuref
!!$  grprcpuref=1e6
!!$
!!$  do k=1,gridsize
!!$    do j=1,gridsize
!!$      do i=1,gridsize
!!$        if (i*j*k == gridsize) then
!!$
!!$          if (grprcpudiff(i,j,k)*errval(i,j,k) < grprcpuref) then
!!$            dims(:)=(/i,j,k/)
!!$            grprcpuref=grprcpudiff(i,j,k)*errval(i,j,k)
!!$          end if
!!$        end if
!!$      end do
!!$    end do
!!$  end do
!!$
!!$ end subroutine dims_create2
!!$

! ***********************************************************************

 subroutine dims_create3

  USE params_mod

  implicit none

  integer,dimension(3,gridsize) :: cdims
  real,dimension(3,gridsize) :: diff,grprcpu
!  real,dimension(gridsize,gridsize,gridsize) :: errval,grprcpudiff
  real :: grprcpudiff
  real :: grprcpuref,tmp
  integer :: i,j,k

! ** checking which combinations produces ok results
! ** grid cells per cpu : grprcpu
! ** how many grid cells per cpu : cdims
! ** number describing the unevenness of point
! **    distribution to each cpu : diff

  call trace%begin ('mpi.f90::dims_create3')

  do i=1,gridsize

    grprcpu(1,i)=(1.0*mx)/(1.0*i)
    cdims(1,i)=int(grprcpu(1,i))
    if (cdims(1,i) == 0) then
      diff(1,i) = 1e20
    else
      diff(1,i)=abs((1.0*mx)/(1.0*cdims(1,i))-grprcpu(1,i))
    endif

    grprcpu(2,i)=(1.0*my)/(1.0*i)
    cdims(2,i)=max(int(grprcpu(2,i)),1)
    if (cdims(2,i) == 0) then
      diff(2,i) = 1e20
    else
      diff(2,i)=abs((1.0*my)/(1.0*cdims(2,i))-grprcpu(2,i))
    endif

    grprcpu(3,i)=(1.0*mz)/(1.0*i)
    cdims(3,i)=max(int(grprcpu(3,i)),1)
    if (cdims(3,i) == 0) then
      diff(3,i) = 1e20
    else
      diff(3,i)=abs((1.0*mz)/(1.0*cdims(3,i))-grprcpu(3,i))
    endif

  end do

  grprcpuref=1e6
  tmp=1e30

  do i=1,gridsize
    do j=1,gridsize/i
      do k=1,gridsize/(i*j)
        if (i*j*k == gridsize) then
          tmp= diff(1,i)/grprcpu(1,i)+ &
               diff(2,j)/grprcpu(2,j)+ &
               diff(3,k)/grprcpu(3,k)
          grprcpudiff=abs(grprcpu(1,i)-grprcpu(2,j))+ &
                      abs(grprcpu(1,i)-grprcpu(3,k))+ &
                      abs(grprcpu(2,j)-grprcpu(3,k))
          if (grprcpudiff*tmp < grprcpuref) then
            dims(:) = (/i,j,k/)
            grprcpuref=grprcpudiff*tmp
          endif
        end if
      end do
    end do
  end do

  call trace%end
 end subroutine dims_create3


! ***********************************************************************

 subroutine resolve_params_parallel (name,value,reread)

  implicit none

  character(len=*) :: name, value
  integer ::strl,strl2,reread

  call trace%begin ('mpi.f90::resolve_params_parallel', 1)

!  write(*,*) "DEBUG: handle: ", id, value
  strl=len_trim(value)
  strl2=len_trim(name)

  select case (name)
  case ('periodic_x')
    call read_nochg_lgk(periodic(1),value(1:strl),reread)

  case ('periodic_y')
     call read_nochg_lgk(periodic(2),value(1:strl),reread)

  case ('periodic_z')
     call read_nochg_lgk(periodic(3),value(1:strl),reread)

  case ('manual_cuts')
     call read_nochg_lgk(manual_cuts,value(1:strl),reread)

  case ('ndim')
    call read_nochg_int(ndim,value(1:strl),reread)

  case ('dims')
    call read_nochg_int_arr(dims,value(1:strl),reread,3)

  case ('reorder')
    call read_nochg_lgk(reorder,value(1:strl),reread)

  end select

  call trace%end
 end subroutine resolve_params_parallel

! ***********************************************************************
 subroutine update4 (f)

  implicit none

  real(kind=4),dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::update4')

  call comm_x(f)
  call comm_y(f)
  call comm_z(f)

!  call comm_xy(f)
!  call comm_xz(f)
!  call comm_yz(f)

!  call comm_xyz(f)

!  call barrier(C_CART,ierr)
  call trace%end
 end subroutine update4

! ***********************************************************************
subroutine update8 (f)

  implicit none

  real(kind=8),dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::update8')

  call comm_x(f)
  call comm_y(f)
  call comm_z(f)

!  call comm_xy(f)
!  call comm_xz(f)
!  call comm_yz(f)

!  call comm_xyz(f)

  call trace%end
 end subroutine update8

! ***********************************************************************
subroutine update2 (f)

  implicit none

  integer,dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::update2')


  call comm_x(f)
  call comm_y(f)
  call comm_z(f)

!  call comm_xy(f)
!  call comm_xz(f)
!  call comm_yz(f)

!  call comm_xyz(f)

  call trace%end
 end subroutine update2
! ***********************************************************************
 subroutine updatec (f)

  implicit none

  complex,dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::updatec')


  call comm_x(f)
  call comm_y(f)
  call comm_z(f)

!  call comm_xy(f)
!  call comm_xz(f)
!  call comm_yz(f)

!  call comm_xyz(f)


!  call barrier(C_CART,ierr)
  call trace%end
 end subroutine updatec

! ************************************************************************
 subroutine update4_x (f)

  implicit none

  real(kind=4),dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::update4_x')

  call comm_x(f)

  call comm_xy(f)
  call comm_xz(f)

  call comm_xyz(f)

  call trace%end
 end subroutine update4_x

! ************************************************************************
 subroutine update8_x (f)

  implicit none

  real(kind=8),dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::update8_x')

  call comm_x(f)

  call comm_xy(f)
  call comm_xz(f)

  call comm_xyz(f)

  call trace%end
 end subroutine update8_x

! ************************************************************************
 subroutine update2_x (f)

  implicit none

  integer,dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::update2_x')

  call comm_x(f)

  call comm_xy(f)
  call comm_xz(f)

  call comm_xyz(f)

  call trace%end
 end subroutine update2_x

! ************************************************************************
 subroutine updatec_x (f)

  implicit none

  complex,dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::updatec_x')

  call comm_x(f)

  call comm_xy(f)
  call comm_xz(f)

  call comm_xyz(f)

  call trace%end
 end subroutine updatec_x

! ************************************************************************
 subroutine update4_xy (f)

  implicit none

  real(kind=4),dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::update4_xy')

  call comm_xy(f)

  call comm_xyz(f)

  call trace%end
 end subroutine update4_xy

! ************************************************************************
 subroutine update8_xy (f)

  implicit none

  real(kind=8),dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::update8_xy')

  call comm_xy(f)

  call comm_xyz(f)

  call trace%end
 end subroutine update8_xy

! ************************************************************************
 subroutine update2_xy (f)

  implicit none

  integer,dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::update2_xy')

  call comm_xy(f)

  call comm_xyz(f)

  call trace%end
 end subroutine update2_xy

! ************************************************************************
 subroutine updatec_xy (f)

  implicit none

  complex,dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::updatec_xy')

  call comm_xy(f)

  call comm_xyz(f)

  call trace%end
 end subroutine updatec_xy

! ***********************************************************************
 subroutine update4_y (f)

  implicit none

  real(kind=4),dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::update4_y')

  call comm_y(f)

  call comm_xy(f)
  call comm_yz(f)

  call comm_xyz(f)


  call trace%end
 end subroutine update4_y

! ***********************************************************************
 subroutine update8_y (f)

  implicit none

  real(kind=8),dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::update8_y')

  call comm_y(f)

  call comm_xy(f)
  call comm_yz(f)

  call comm_xyz(f)


  call trace%end
 end subroutine update8_y

! ***********************************************************************
 subroutine update2_y (f)

  implicit none

  integer,dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::update2_y')

  call comm_y(f)

  call comm_xy(f)
  call comm_yz(f)

  call comm_xyz(f)

  call trace%end
 end subroutine update2_y

! ***********************************************************************
 subroutine updatec_y (f)

  implicit none

  complex,dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::updatec_y')

  call comm_y(f)

  call comm_xy(f)
  call comm_yz(f)

  call comm_xyz(f)


  call trace%end
 end subroutine updatec_y

! ***********************************************************************
 subroutine update4_z (f)

  implicit none

  real(kind=4),dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::update4_z')

  call comm_z(f)

  call comm_xz(f)
  call comm_yz(f)

  call comm_xyz(f)


  call trace%end
 end subroutine update4_z

! ***********************************************************************
 subroutine update8_z (f)

  implicit none

  real(kind=8),dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::update8_z')

  call comm_z(f)

  call comm_xz(f)
  call comm_yz(f)

  call comm_xyz(f)

  call trace%end
 end subroutine update8_z

! ***********************************************************************
 subroutine update2_z (f)

  implicit none

  integer,dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::update2_z')

  call comm_z(f)

  call comm_xz(f)
  call comm_yz(f)

  call comm_xyz(f)

  call trace%end
 end subroutine update2_z

! ***********************************************************************
 subroutine updatec_z (f)

  implicit none

  complex,dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::updatec_z')

  call comm_z(f)

  call comm_xz(f)
  call comm_yz(f)

  call comm_xyz(f)


  call trace%end
 end subroutine updatec_z


! ************************************************************************
 subroutine update2d4_xy (f)

  implicit none

  real(kind=4),dimension(xsb:xeb,ysb:yeb) :: f

  call trace%begin ('mpi.f90::update2d4_xy')

  call comm2d_xy(f)

  call trace%end
 end subroutine update2d4_xy

! ************************************************************************
 subroutine update2d8_xy (f)

  implicit none

  real(kind=8),dimension(xsb:xeb,ysb:yeb) :: f

  call trace%begin ('mpi.f90::update2d8_xy')

  call comm2d_xy(f)

  call trace%end
 end subroutine update2d8_xy

! ************************************************************************
 subroutine update2d2_xy (f)

  implicit none

  integer,dimension(xsb:xeb,ysb:yeb) :: f

  call trace%begin ('mpi.f90::update2d2_xy')

  call comm2d_xy(f)

  call trace%end
 end subroutine update2d2_xy

! ************************************************************************
 subroutine updatec2d_xy (f)

  implicit none

  complex,dimension(xsb:xeb,ysb:yeb) :: f

  call trace%begin ('mpi.f90::updatec2d_xy')

  call comm2d_xy(f)

  call trace%end
 end subroutine updatec2d_xy

! ***********************************************************************

 subroutine check_update (f)

  implicit none

  real,dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::check_update')

  call print_check(f,"X  ")

  call print_check(f,"Y  ")

  call print_check(f,"Z  ")

  call print_check(f,"XY ")

  call print_check(f,"XZ ")

  call print_check(f,"YZ ")

  call print_check(f,"XYZ")


  call trace%end
 end subroutine check_update

! ***********************************************************************
 subroutine comm4_x (f)

  implicit none

  real(kind=4),dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::comm4_x')

  if (exist_xup.or.exist_xdn) then

    call MPI_SENDRECV(f(xe-(mb-1) :xe         , &
         ys        :ye         , &
         zs        :ze        ), &
         mb*yr*zr,MPI_REAL,R_XUP,1  , &
         f(xsb     :xsb+(mb-1) , &
         ys        :ye         , &
         zs        :ze        ), &
         mb*yr*zr,MPI_REAL,R_XDN,1  , &
         C_CART,status,ierr)

    call MPI_SENDRECV(f(xs        :xs+(mb-1)  , &
         ys        :ye         , &
         zs        :ze        ), &
         mb*yr*zr,MPI_REAL,R_XDN,2  , &
         f(xeb-(mb-1):xeb        , &
         ys        :ye         , &
         zs        :ze        ), &
         mb*yr*zr,MPI_REAL,R_XUP,2  , &
         C_CART,status,ierr)

  end if

  call trace%end
 end subroutine comm4_x

! ***********************************************************************
 subroutine comm4_x2 (f)

  implicit none

  real(kind=4),dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f
!  integer(KIND=MPI_ADDRESS_KIND) :: address
  integer,dimension(MPI_STATUS_SIZE,2) :: status_xup,status_xdn
!  integer(Kind=4) :: req_recvxup,req_recvxdn,req_sendxdn,req_sendxup
  integer(kind=4),dimension(2) :: req_xdn,req_xup
  integer(kind=4) :: req_sxup,req_sxdn,req_rxup,req_rxdn

  call trace%begin ('mpi.f90::comm4_x2')

!  if (exist_xup) call MPI_IRECV(f,1,T_RECVXUP_REAL4,R_XUP,1,C_CART,req_xup(1),ierr)
!  if (exist_xdn) call MPI_IRECV(f,1,T_RECVXDN_REAL4,R_XDN,2,C_CART,req_xdn(1),ierr)
  if (exist_xup) call MPI_IRECV(f,1,T_RECVXUP_REAL4,R_XUP,1,C_CART,req_rxup,ierr)
  if (exist_xdn) call MPI_IRECV(f,1,T_RECVXDN_REAL4,R_XDN,2,C_CART,req_rxdn,ierr)

  if (exist_xup) then
!    call MPI_ISEND(f,1,T_SENDXUP_REAL4,R_XUP,2,C_CART,req_xup(2),ierr)
    call MPI_ISEND(f,1,T_SENDXUP_REAL4,R_XUP,2,C_CART,req_sxup,ierr)
!    call MPI_WAITALL(2,req_xup,status_xup,ierr)
  end if

  if (exist_xdn) then
    call MPI_ISEND(f,1,T_SENDXDN_REAL4,R_XDN,1,C_CART,req_sxdn,ierr)
!    call MPI_ISEND(f,1,T_SENDXDN_REAL4,R_XDN,1,C_CART,req_xdn(2),ierr)
!    call MPI_WAITALL(2,req_xdn,status_xdn,ierr)
  end if

!  if (exist_xup) call MPI_WAITALL(2,req_xup,status_xup,ierr)
  if (exist_xup) then
    call MPI_WAIT(req_rxup,status,ierr)
    call MPI_WAIT(req_sxup,status,ierr)
  end if
!  if (exist_xdn) call MPI_WAITALL(2,req_xdn,status_xdn,ierr)
  if (exist_xdn) then
    call MPI_WAIT(req_sxdn,status,ierr)
    call MPI_WAIT(req_rxdn,status,ierr)
  end if

!  call MPI_BARRIER(C_CART,ierr)
!  print*,myid,'out of comm_x'

  call trace%end
 end subroutine comm4_x2

! ***********************************************************************
 subroutine comm4_x3 (f)

  implicit none

  real(kind=4),dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f
!  integer(KIND=MPI_ADDRESS_KIND) :: address
  integer,dimension(MPI_STATUS_SIZE,2) :: status_xup,status_xdn
!  integer(Kind=4) :: req_recvxup,req_recvxdn,req_sendxdn,req_sendxup
  integer(kind=4),dimension(2) :: req_xdn,req_xup
  integer(kind=4) :: req_sxup,req_sxdn,req_rxup,req_rxdn

  call trace%begin ('mpi.f90::comm4_x3')

!  if (exist_xup) call MPI_IRECV(f,1,T_RECVXUP_REAL4,R_XUP,1,C_CART,req_xup(1),ierr)
!  if (exist_xdn) call MPI_IRECV(f,1,T_RECVXDN_REAL4,R_XDN,2,C_CART,req_xdn(1),ierr)
  if (exist_xup) call MPI_IRECV(f,1,T_RECVXUP_REAL4,R_XUP,1,C_CART,req_rxup,ierr)
  if (exist_xdn) call MPI_IRECV(f,1,T_RECVXDN_REAL4,R_XDN,2,C_CART,req_rxdn,ierr)

  if (exist_xup) then
!    call MPI_ISEND(f,1,T_SENDXUP_REAL4,R_XUP,2,C_CART,req_xup(2),ierr)
    call MPI_IBSEND(f,1,T_SENDXUP_REAL4,R_XUP,2,C_CART,req_sxup,ierr)
!    call MPI_WAITALL(2,req_xup,status_xup,ierr)
  end if

  if (exist_xdn) then
    call MPI_IBSEND(f,1,T_SENDXDN_REAL4,R_XDN,1,C_CART,req_sxdn,ierr)
!    call MPI_ISEND(f,1,T_SENDXDN_REAL4,R_XDN,1,C_CART,req_xdn(2),ierr)
!    call MPI_WAITALL(2,req_xdn,status_xdn,ierr)
  end if

!  if (exist_xup) call MPI_WAITALL(2,req_xup,status_xup,ierr)
  if (exist_xup) then
    call MPI_WAIT(req_rxup,status,ierr)
    call MPI_WAIT(req_sxup,status,ierr)
  end if
!  if (exist_xdn) call MPI_WAITALL(2,req_xdn,status_xdn,ierr)
  if (exist_xdn) then
    call MPI_WAIT(req_sxdn,status,ierr)
    call MPI_WAIT(req_rxdn,status,ierr)
  end if

!  call MPI_BARRIER(C_CART,ierr)
!  print*,myid,'out of comm_x'

  call trace%end
 end subroutine comm4_x3

! ***********************************************************************
 subroutine comm8_x (f)

  implicit none

  real(kind=8),dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::comm8_x')


  if (exist_xup.or.exist_xdn) then

    call MPI_SENDRECV(f(xe-(mb-1) :xe         , &
         ys        :ye         , &
         zs        :ze        ), &
         mb*yr*zr,MPI_DOUBLE_PRECISION,R_XUP,1  , &
         f(xsb     :xsb+(mb-1) , &
         ys        :ye         , &
         zs        :ze        ), &
         mb*yr*zr,MPI_DOUBLE_PRECISION,R_XDN,1  , &
         C_CART,status,ierr)

    call MPI_SENDRECV(f(xs        :xs+(mb-1)  , &
         ys        :ye         , &
         zs        :ze        ), &
         mb*yr*zr,MPI_DOUBLE_PRECISION,R_XDN,2  , &
         f(xeb-(mb-1):xeb        , &
         ys        :ye         , &
         zs        :ze        ), &
         mb*yr*zr,MPI_DOUBLE_PRECISION,R_XUP,2  , &
         C_CART,status,ierr)

  end if

  call trace%end
 end subroutine comm8_x

! ***********************************************************************
 subroutine comm2_x (f)

  implicit none

  integer,dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::comm2_x')


  if (exist_xup.or.exist_xdn) then

    call MPI_SENDRECV(f(xe-(mb-1) :xe         , &
         ys        :ye         , &
         zs        :ze        ), &
         mb*yr*zr,OSC_INTEGER,R_XUP,1  , &
         f(xsb     :xsb+(mb-1) , &
         ys        :ye         , &
         zs        :ze        ), &
         mb*yr*zr,OSC_INTEGER,R_XDN,1  , &
         C_CART,status,ierr)

    call MPI_SENDRECV(f(xs        :xs+(mb-1)  , &
         ys        :ye         , &
         zs        :ze        ), &
         mb*yr*zr,OSC_INTEGER,R_XDN,2  , &
         f(xeb-(mb-1):xeb        , &
         ys        :ye         , &
         zs        :ze        ), &
         mb*yr*zr,OSC_INTEGER,R_XUP,2  , &
         C_CART,status,ierr)

  end if

  call trace%end
 end subroutine comm2_x

! ***********************************************************************
 subroutine commc_x (f)

  implicit none

  complex,dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::commc_x')


  if (exist_xup.or.exist_xdn) then

    call MPI_SENDRECV(f(xe-(mb-1) :xe         , &
         ys        :ye         , &
         zs        :ze        ), &
         mb*yr*zr,MPI_COMPLEX,R_XUP,1  , &
         f(xsb     :xsb+(mb-1) , &
         ys        :ye         , &
         zs        :ze        ), &
         mb*yr*zr,MPI_COMPLEX,R_XDN,1  , &
         C_CART,status,ierr)

    call MPI_SENDRECV(f(xs        :xs+(mb-1)  , &
         ys        :ye         , &
         zs        :ze        ), &
         mb*yr*zr,MPI_COMPLEX,R_XDN,2  , &
         f(xeb-(mb-1):xeb        , &
         ys        :ye         , &
         zs        :ze        ), &
         mb*yr*zr,MPI_COMPLEX,R_XUP,2  , &
         C_CART,status,ierr)

  end if

  call trace%end
 end subroutine commc_x

! ***********************************************************************
 subroutine comm4_y2 (f)

  implicit none

  real(kind=4),dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f
!  integer(kind=4) :: req_sendyup,req_recvyup,req_sendydn,req_recvydn
  integer(kind=4),dimension(2) :: req_ydn,req_yup
!  integer,dimension(MPI_STATUS_SIZE,2) :: status_ydn,status_yup
!  integer(kind=MPI_ADDRESS_KIND) :: address
  integer(kind=4) :: req_syup,req_sydn,req_ryup,req_rydn

  call trace%begin ('mpi.f90::comm4_y2')

!  if (exist_yup) call MPI_IRECV(f,1,T_RECVYUP_REAL4,R_YUP,3,C_CART,req_yup(1),ierr)
!  if (exist_ydn) call MPI_IRECV(f,1,T_RECVYDN_REAL4,R_YDN,4,C_CART,req_ydn(1),ierr)
  if (exist_yup) call MPI_IRECV(f,1,T_RECVYUP_REAL4,R_YUP,3,C_CART,req_ryup,ierr)
  if (exist_ydn) call MPI_IRECV(f,1,T_RECVYDN_REAL4,R_YDN,4,C_CART,req_rydn,ierr)

  if (exist_yup) then
!    call MPI_ISEND(f,1,T_SENDYUP_REAL4,R_YUP,4,C_CART,req_yup(2),ierr)
    call MPI_ISEND(f,1,T_SENDYUP_REAL4,R_YUP,4,C_CART,req_syup,ierr)
!    call MPI_WAITALL(2,req_yup,status_yup,ierr)
  end if
  if (exist_ydn) then
!    call MPI_ISEND(f,1,T_SENDYDN_REAL4,R_YDN,3,C_CART,req_ydn(2),ierr)
    call MPI_ISEND(f,1,T_SENDYDN_REAL4,R_YDN,3,C_CART,req_sydn,ierr)
!    call MPI_WAITALL(2,req_ydn,status_ydn,ierr)
  end if

  if (exist_yup) then
    call MPI_WAIT(req_ryup,status,ierr)
    call MPI_WAIT(req_syup,status,ierr)
  endif
  if (exist_ydn) then
    call MPI_WAIT(req_rydn,status,ierr)
    call MPI_WAIT(req_sydn,status,ierr)
  end if
!  call MPI_GET_ADDRESS(f,address,ierr)

  call trace%end
 end subroutine comm4_y2

! ***********************************************************************
 subroutine comm4_y3 (f)

  implicit none

  real(kind=4),dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f
!  integer(kind=4) :: req_sendyup,req_recvyup,req_sendydn,req_recvydn
  integer(kind=4),dimension(2) :: req_ydn,req_yup
!  integer,dimension(MPI_STATUS_SIZE,2) :: status_ydn,status_yup
!  integer(kind=MPI_ADDRESS_KIND) :: address
  integer(kind=4) :: req_syup,req_sydn,req_ryup,req_rydn

  call trace%begin ('mpi.f90::comm4_y3')

!  if (exist_yup) call MPI_IRECV(f,1,T_RECVYUP_REAL4,R_YUP,3,C_CART,req_yup(1),ierr)
!  if (exist_ydn) call MPI_IRECV(f,1,T_RECVYDN_REAL4,R_YDN,4,C_CART,req_ydn(1),ierr)
  if (exist_yup) call MPI_IRECV(f,1,T_RECVYUP_REAL4,R_YUP,3,C_CART,req_ryup,ierr)
  if (exist_ydn) call MPI_IRECV(f,1,T_RECVYDN_REAL4,R_YDN,4,C_CART,req_rydn,ierr)

  if (exist_yup) then
!    call MPI_ISEND(f,1,T_SENDYUP_REAL4,R_YUP,4,C_CART,req_yup(2),ierr)
    call MPI_IBSEND(f,1,T_SENDYUP_REAL4,R_YUP,4,C_CART,req_syup,ierr)
!    call MPI_WAITALL(2,req_yup,status_yup,ierr)
  end if
  if (exist_ydn) then
!    call MPI_ISEND(f,1,T_SENDYDN_REAL4,R_YDN,3,C_CART,req_ydn(2),ierr)
    call MPI_IBSEND(f,1,T_SENDYDN_REAL4,R_YDN,3,C_CART,req_sydn,ierr)
!    call MPI_WAITALL(2,req_ydn,status_ydn,ierr)
  end if

  if (exist_yup) then
    call MPI_WAIT(req_ryup,status,ierr)
    call MPI_WAIT(req_syup,status,ierr)
  endif
  if (exist_ydn) then
    call MPI_WAIT(req_rydn,status,ierr)
    call MPI_WAIT(req_sydn,status,ierr)
  end if
!  call MPI_GET_ADDRESS(f,address,ierr)

  call trace%end
 end subroutine comm4_y3

! ***********************************************************************

 subroutine comm4_y (f)

  implicit none

  real(kind=4),dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::comm4_y')

  if (exist_yup.or.exist_ydn) then

      call MPI_SENDRECV(f(xsb       :xeb              , &
           ye-(mb-1):ye              , &
           zs       :ze       )      , &
           xrb*mb*zr,MPI_REAL,R_YUP,3, &
           f(xsb       :xeb              , &
           ysb         :ysb+(mb-1)         , &
           zs       :ze       )      , &
           xrb*mb*zr,MPI_REAL,R_YDN,3, &
           C_CART,status,ierr)


      call MPI_SENDRECV(f(xsb       :xeb              , &
           ys       :ys+(mb-1)       , &
           zs       :ze       )      , &
           xrb*mb*zr,MPI_REAL,R_YDN,4          , &
           f(xsb       :xeb              , &
           yeb-(mb-1)  :yeb                , &
           zs       :ze       )      , &
           xrb*mb*zr,MPI_REAL,R_YUP,4           , &
           C_CART,status,ierr)

    end if

  call trace%end
   end subroutine comm4_y

! ***********************************************************************
 subroutine comm8_y (f)

  implicit none

  real(kind=8),dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::comm8_y')


  if (exist_yup.or.exist_ydn) then

      call MPI_SENDRECV(f(xs       :xe              , &
           ye-(mb-1):ye              , &
           zs       :ze       )      , &
           xr*mb*zr,MPI_DOUBLE_PRECISION,R_YUP,3, &
           f(xs       :xe              , &
           ysb         :ysb+(mb-1)         , &
           zs       :ze       )      , &
           xr*mb*zr,MPI_DOUBLE_PRECISION,R_YDN,3, &
           C_CART,status,ierr)


      call MPI_SENDRECV(f(xs       :xe              , &
           ys       :ys+(mb-1)       , &
           zs       :ze       )      , &
           xr*mb*zr,MPI_DOUBLE_PRECISION,R_YDN,4          , &
           f(xs       :xe              , &
           yeb-(mb-1)  :yeb                , &
           zs       :ze       )      , &
           xr*mb*zr,MPI_DOUBLE_PRECISION,R_YUP,4           , &
           C_CART,status,ierr)

    end if

  call trace%end
   end subroutine comm8_y

! ***********************************************************************
 subroutine comm2_y (f)

  implicit none

  integer,dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::comm2_y')


  if (exist_yup.or.exist_ydn) then

      call MPI_SENDRECV(f(xsb       :xeb              , &
           ye-(mb-1):ye              , &
           zs       :ze       )      , &
           xrb*mb*zr,OSC_INTEGER,R_YUP,3, &
           f(xsb       :xeb              , &
           ysb         :ysb+(mb-1)         , &
           zs       :ze       )      , &
           xrb*mb*zr,OSC_INTEGER,R_YDN,3, &
           C_CART,status,ierr)


      call MPI_SENDRECV(f(xsb       :xeb              , &
           ys       :ys+(mb-1)       , &
           zs       :ze       )      , &
           xrb*mb*zr,OSC_INTEGER,R_YDN,4          , &
           f(xsb       :xeb              , &
           yeb-(mb-1)  :yeb                , &
           zs       :ze       )      , &
           xrb*mb*zr,OSC_INTEGER,R_YUP,4           , &
           C_CART,status,ierr)

    end if

  call trace%end
   end subroutine comm2_y

! ***********************************************************************
 subroutine commc_y (f)

  implicit none

  complex,dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::commc_y')


  if (exist_yup.or.exist_ydn) then

      call MPI_SENDRECV(f(xsb       :xeb              , &
           ye-(mb-1):ye              , &
           zs       :ze       )      , &
           xrb*mb*zr,MPI_COMPLEX,R_YUP,3, &
           f(xsb       :xeb              , &
           ysb         :ysb+(mb-1)         , &
           zs       :ze       )      , &
           xrb*mb*zr,MPI_COMPLEX,R_YDN,3, &
           C_CART,status,ierr)


      call MPI_SENDRECV(f(xsb       :xeb              , &
           ys       :ys+(mb-1)       , &
           zs       :ze       )      , &
           xrb*mb*zr,MPI_COMPLEX,R_YDN,4          , &
           f(xsb       :xeb              , &
           yeb-(mb-1)  :yeb                , &
           zs       :ze       )      , &
           xrb*mb*zr,MPI_COMPLEX,R_YUP,4           , &
           C_CART,status,ierr)

    end if

  call trace%end
   end subroutine commc_y

! ***********************************************************************

 subroutine comm4_z2 (f)

  implicit none

  real(kind=4),dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f
  integer(kind=4) :: req_recvzdn,req_recvzup,req_sendzup,req_sendzdn
!  integer(kind=4),dimension(2) :: req_zdn,req_zup
!  integer,dimension(MPI_STATUS_SIZE,2) :: status_zdn,status_zup
!  integer(kind=MPI_ADDRESS_KIND) :: address
  integer(kind=4) :: req_szup,req_szdn,req_rzup,req_rzdn

  call trace%begin ('mpi.f90::comm4_z2')

!  if (exist_zup) call MPI_IRECV(f,1,T_RECVZUP_REAL4,R_ZUP,5,C_CART,req_zup(1),ierr)
!  if (exist_zdn) call MPI_IRECV(f,1,T_RECVZDN_REAL4,R_ZDN,6,C_CART,req_zdn(1),ierr)
  if (exist_zup) call MPI_IRECV(f,1,T_RECVZUP_REAL4,R_ZUP,5,C_CART,req_rzup,ierr)
  if (exist_zdn) call MPI_IRECV(f,1,T_RECVZDN_REAL4,R_ZDN,6,C_CART,req_rzdn,ierr)

  if (exist_zup) then
!    call MPI_ISEND(f,1,T_SENDZUP_REAL4,R_ZUP,6,C_CART,req_zup(2),ierr)
    call MPI_ISEND(f,1,T_SENDZUP_REAL4,R_ZUP,6,C_CART,req_szup,ierr)
!    call MPI_WAITALL(2,req_zup,status_zup,ierr)
  end if

  if (exist_zdn) then
!    call MPI_ISEND(f,1,T_SENDZDN_REAL4,R_ZDN,5,C_CART,req_zdn(2),ierr)
    call MPI_ISEND(f,1,T_SENDZDN_REAL4,R_ZDN,5,C_CART,req_szdn,ierr)
!    call MPI_WAITALL(2,req_zdn,status_zdn,ierr)
  end if

  if (exist_zup) then
    call MPI_WAIT(req_rzup,status,ierr)
    call MPI_WAIT(req_szup,status,ierr)
  end if
  if (exist_zdn) then
    call MPI_WAIT(req_rzdn,status,ierr)
    call MPI_WAIT(req_szdn,status,ierr)
  end if
!  call MPI_GET_ADDRESS(f,address,ierr)

  call trace%end
 end subroutine comm4_z2

! ***********************************************************************
 subroutine comm4_z3 (f)

  implicit none

  real(kind=4),dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f
  integer(kind=4) :: req_recvzdn,req_recvzup,req_sendzup,req_sendzdn
!  integer(kind=4),dimension(2) :: req_zdn,req_zup
!  integer,dimension(MPI_STATUS_SIZE,2) :: status_zdn,status_zup
!  integer(kind=MPI_ADDRESS_KIND) :: address
  integer(kind=4) :: req_szup,req_szdn,req_rzup,req_rzdn

  call trace%begin ('mpi.f90::comm4_z3')

!  if (exist_zup) call MPI_IRECV(f,1,T_RECVZUP_REAL4,R_ZUP,5,C_CART,req_zup(1),ierr)
!  if (exist_zdn) call MPI_IRECV(f,1,T_RECVZDN_REAL4,R_ZDN,6,C_CART,req_zdn(1),ierr)
  if (exist_zup) call MPI_IRECV(f,1,T_RECVZUP_REAL4,R_ZUP,5,C_CART,req_rzup,ierr)
  if (exist_zdn) call MPI_IRECV(f,1,T_RECVZDN_REAL4,R_ZDN,6,C_CART,req_rzdn,ierr)

  if (exist_zup) then
!    call MPI_ISEND(f,1,T_SENDZUP_REAL4,R_ZUP,6,C_CART,req_zup(2),ierr)
    call MPI_IBSEND(f,1,T_SENDZUP_REAL4,R_ZUP,6,C_CART,req_szup,ierr)
!    call MPI_WAITALL(2,req_zup,status_zup,ierr)
  end if

  if (exist_zdn) then
!    call MPI_ISEND(f,1,T_SENDZDN_REAL4,R_ZDN,5,C_CART,req_zdn(2),ierr)
    call MPI_IBSEND(f,1,T_SENDZDN_REAL4,R_ZDN,5,C_CART,req_szdn,ierr)
!    call MPI_WAITALL(2,req_zdn,status_zdn,ierr)
  end if

  if (exist_zup) then
    call MPI_WAIT(req_rzup,status,ierr)
    call MPI_WAIT(req_szup,status,ierr)
  end if
  if (exist_zdn) then
    call MPI_WAIT(req_rzdn,status,ierr)
    call MPI_WAIT(req_szdn,status,ierr)
  end if
!  call MPI_GET_ADDRESS(f,address,ierr)

  call trace%end
 end subroutine comm4_z3

! ***********************************************************************

 subroutine comm4_z (f)

  implicit none

  real(kind=4),dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::comm4_z')

    if (exist_zdn.or.exist_zup) then

      call MPI_SENDRECV( &
           f(xsb        :xeb              , &
             ysb        :yeb              , &
             ze-(mb-1) :ze          )   , &
           xrb*yrb*mb,MPI_REAL,R_ZUP,5, &
           f(xsb        :xeb              , &
             ysb        :yeb              , &
             zsb       :zsb+(mb-1)  )   , &
           xrb*yrb*mb,MPI_REAL,R_ZDN,5, &
           C_CART,status,ierr)


      call MPI_SENDRECV( &
           f(xsb        :xeb              , &
             ysb        :yeb              , &
             zs        :zs+(mb-1))      , &
           xrb*yrb*mb,MPI_REAL,R_ZDN,6, &
           f(xsb        :xeb              , &
             ysb        :yeb              , &
             zeb-(mb-1):zeb         )   , &
           xrb*yrb*mb,MPI_REAL,R_ZUP,6, &
           C_CART,status,ierr)

    end if

  call trace%end
   end subroutine comm4_z

! ***********************************************************************

 subroutine comm8_z (f)

  implicit none

  real(kind=8),dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::comm8_z')

    if (exist_zdn.or.exist_zup) then

      call MPI_SENDRECV(f(xsb       :xeb              , &
           ysb      :yeb             , &
           ze-(mb-1):ze       )      , &
           xrb*yrb*mb,MPI_DOUBLE_PRECISION,R_ZUP,5, &
           f(xsb      :xeb             , &
           ysb      :yeb             , &
           zsb         :zsb+(mb-1)  )      , &
           xrb*yrb*mb,MPI_DOUBLE_PRECISION,R_ZDN,5, &
           C_CART,status,ierr)


      call MPI_SENDRECV(f(xsb      :xeb             , &
           ysb      :yeb             , &
           zs       :zs+(mb-1))      , &
           xrb*yrb*mb,MPI_DOUBLE_PRECISION,R_ZDN,6, &
           f(xsb      :xeb             , &
           ysb      :yeb             , &
           zeb-(mb-1)  :zeb         )      , &
           xrb*yrb*mb,MPI_DOUBLE_PRECISION,R_ZUP,6, &
           C_CART,status,ierr)

    end if

  call trace%end
   end subroutine comm8_z

! ***********************************************************************

 subroutine comm2_z (f)

  implicit none

  integer,dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::comm2_z')

    if (exist_zdn.or.exist_zup) then

      call MPI_SENDRECV(f(xsb      :xeb             , &
           ysb      :yeb             , &
           ze-(mb-1):ze       )      , &
           xrb*yrb*mb,OSC_INTEGER,R_ZUP,5, &
           f(xsb      :xeb             , &
           ysb      :yeb             , &
           zsb         :zsb+(mb-1)  )      , &
           xrb*yrb*mb,OSC_INTEGER,R_ZDN,5, &
           C_CART,status,ierr)


      call MPI_SENDRECV(f(xsb      :xeb             , &
           ysb      :yeb             , &
           zs       :zs+(mb-1))      , &
           xrb*yrb*mb,OSC_INTEGER,R_ZDN,6, &
           f(xsb      :xeb             , &
           ysb      :yeb             , &
           zeb-(mb-1)  :zeb         )      , &
           xrb*yrb*mb,OSC_INTEGER,R_ZUP,6, &
           C_CART,status,ierr)

    end if

  call trace%end
   end subroutine comm2_z

! ***********************************************************************

 subroutine commc_z (f)

  implicit none

  complex,dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::commc_z')

    if (exist_zdn.or.exist_zup) then

      call MPI_SENDRECV( &
           f(xsb       :xeb             , &
             ysb       :yeb             , &
             ze-(mb-1) :ze          )   , &
           xrb*yrb*mb,MPI_COMPLEX,R_ZUP,5, &
           f(xsb       :xeb             , &
             ysb       :yeb             , &
             zsb       :zsb+(mb-1)  )   , &
           xrb*yrb*mb,MPI_COMPLEX,R_ZDN,5, &
           C_CART,status,ierr)


      call MPI_SENDRECV( &
           f(xsb       :xeb             , &
             ysb       :yeb             , &
             zs        :zs+(mb-1))      , &
           xrb*yrb*mb,MPI_COMPLEX,R_ZDN,6, &
           f(xsb       :xeb             , &
             ysb       :yeb             , &
             zeb-(mb-1):zeb         )   , &
           xrb*yrb*mb,MPI_COMPLEX,R_ZUP,6, &
           C_CART,status,ierr)

    end if

  call trace%end
   end subroutine commc_z

!************************************************************************

   subroutine comm4_xy (f)

    implicit none

    real(kind=4),dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::comm4_xy')

    if ((mx /= 1).and.(my /= 1)) then
      if (exist_xupyup.or.exist_xdnydn) then

        call MPI_SENDRECV(f(xe-(mb-1):xe               , &
             ye-(mb-1):ye               , &
             zs       :ze       )       , &
             mb*mb*zr,MPI_REAL,R_XUPYUP,7, &
             f(xsb         :xsb+(mb-1)          , &
             ysb         :ysb+(mb-1)          , &
             zs       :ze       )       , &
             mb*mb*zr,MPI_REAL,R_XDNYDN,7, &
             C_CART,status,ierr)

        call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
             ys       :ys+(mb-1)        , &
             zs       :ze       )       , &
             mb*mb*zr,MPI_REAL,R_XDNYDN,8, &
             f(xeb-(mb-1)  :xeb                 , &
             yeb-(mb-1)  :yeb                 , &
             zs       :ze       )       , &
             mb*mb*zr,MPI_REAL,R_XUPYUP,8, &
             C_CART,status,ierr)
      end if

      if ((exist_xup.and.exist_ydn).or.(exist_xdn.and.exist_yup)) then

        call MPI_SENDRECV(f(xe-(mb-1):xe               , &
             ys       :ys+(mb-1)        , &
             zs       :ze       )       , &
             mb*mb*zr,MPI_REAL,R_XUPYDN,9, &
             f(xsb         :xsb+(mb-1)          , &
             yeb-(mb-1)  :yeb                 , &
             zs       :ze       )       , &
             mb*mb*zr,MPI_REAL,R_XDNYUP,9, &
             C_CART,status,ierr)

        call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
             ye-(mb-1):ye               , &
             zs       :ze       )       , &
             mb*mb*zr,MPI_REAL,R_XDNYUP,10, &
             f(xeb-(mb-1)  :xeb                 , &
             ysb         :ysb+(mb-1)          , &
             zs       :ze       )       , &
             mb*mb*zr,MPI_REAL,R_XUPYDN,10, &
             C_CART,status,ierr)

      end if
    end if
  call trace%end
   end subroutine comm4_xy

! ***********************************************************************

 subroutine comm8_xy (f)

  implicit none

  real(kind=8),dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::comm8_xy')

    if ((exist_xup.and.exist_yup).or.(exist_xdn.and.exist_ydn)) then

      call MPI_SENDRECV(f(xe-(mb-1):xe               , &
           ye-(mb-1):ye               , &
           zs       :ze       )       , &
           mb*mb*zr,MPI_DOUBLE_PRECISION,R_XUPYUP,7, &
           f(xsb         :xsb+(mb-1)          , &
           ysb         :ysb+(mb-1)          , &
           zs       :ze       )       , &
           mb*mb*zr,MPI_DOUBLE_PRECISION,R_XDNYDN,7, &
           C_CART,status,ierr)

      call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
           ys       :ys+(mb-1)        , &
           zs       :ze       )       , &
           mb*mb*zr,MPI_DOUBLE_PRECISION,R_XDNYDN,8, &
           f(xeb-(mb-1)  :xeb                 , &
           yeb-(mb-1)  :yeb                 , &
           zs       :ze       )       , &
           mb*mb*zr,MPI_DOUBLE_PRECISION,R_XUPYUP,8, &
           C_CART,status,ierr)

    endif

    if ((exist_xup.and.exist_ydn).or.(exist_xdn.and.exist_yup)) then

      call MPI_SENDRECV(f(xe-(mb-1):xe               , &
           ys       :ys+(mb-1)        , &
           zs       :ze       )       , &
           mb*mb*zr,MPI_DOUBLE_PRECISION,R_XUPYDN,9, &
           f(xsb         :xsb+(mb-1)          , &
           yeb-(mb-1)  :yeb                 , &
           zs       :ze       )       , &
           mb*mb*zr,MPI_DOUBLE_PRECISION,R_XDNYUP,9, &
           C_CART,status,ierr)

      call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
           ye-(mb-1):ye               , &
           zs       :ze       )       , &
           mb*mb*zr,MPI_DOUBLE_PRECISION,R_XDNYUP,10, &
           f(xeb-(mb-1)  :xeb                 , &
           ysb         :ysb+(mb-1)          , &
           zs       :ze       )       , &
           mb*mb*zr,MPI_DOUBLE_PRECISION,R_XUPYDN,10, &
           C_CART,status,ierr)

    end if

  call trace%end
   end subroutine comm8_xy

! ***********************************************************************

 subroutine comm2_xy (f)

  implicit none

  integer,dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::comm2_xy')

    if ((exist_xup.and.exist_yup).or.(exist_xdn.and.exist_ydn)) then

      call MPI_SENDRECV(f(xe-(mb-1):xe               , &
           ye-(mb-1):ye               , &
           zs       :ze       )       , &
           mb*mb*zr,OSC_INTEGER,R_XUPYUP,7, &
           f(xsb         :xsb+(mb-1)          , &
           ysb         :ysb+(mb-1)          , &
           zs       :ze       )       , &
           mb*mb*zr,OSC_INTEGER,R_XDNYDN,7, &
           C_CART,status,ierr)

      call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
           ys       :ys+(mb-1)        , &
           zs       :ze       )       , &
           mb*mb*zr,OSC_INTEGER,R_XDNYDN,8, &
           f(xeb-(mb-1)  :xeb                 , &
           yeb-(mb-1)  :yeb                 , &
           zs       :ze       )       , &
           mb*mb*zr,OSC_INTEGER,R_XUPYUP,8, &
           C_CART,status,ierr)

    endif

    if ((exist_xup.and.exist_ydn).or.(exist_xdn.and.exist_yup)) then

      call MPI_SENDRECV(f(xe-(mb-1):xe               , &
           ys       :ys+(mb-1)        , &
           zs       :ze       )       , &
           mb*mb*zr,OSC_INTEGER,R_XUPYDN,9, &
           f(xsb         :xsb+(mb-1)          , &
           yeb-(mb-1)  :yeb                 , &
           zs       :ze       )       , &
           mb*mb*zr,OSC_INTEGER,R_XDNYUP,9, &
           C_CART,status,ierr)

      call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
           ye-(mb-1):ye               , &
           zs       :ze       )       , &
           mb*mb*zr,OSC_INTEGER,R_XDNYUP,10, &
           f(xeb-(mb-1)  :xeb                 , &
           ysb         :ysb+(mb-1)          , &
           zs       :ze       )       , &
           mb*mb*zr,OSC_INTEGER,R_XUPYDN,10, &
           C_CART,status,ierr)

    end if

  call trace%end
   end subroutine comm2_xy

!************************************************************************

 subroutine commc_xy (f)

  implicit none

  complex,dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::commc_xy')

    if ((exist_xup.and.exist_yup).or.(exist_xdn.and.exist_ydn)) then

      call MPI_SENDRECV(f(xe-(mb-1):xe               , &
           ye-(mb-1):ye               , &
           zs       :ze       )       , &
           mb*mb*zr,MPI_COMPLEX,R_XUPYUP,7, &
           f(xsb         :xsb+(mb-1)          , &
           ysb         :ysb+(mb-1)          , &
           zs       :ze       )       , &
           mb*mb*zr,MPI_COMPLEX,R_XDNYDN,7, &
           C_CART,status,ierr)

      call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
           ys       :ys+(mb-1)        , &
           zs       :ze       )       , &
           mb*mb*zr,MPI_COMPLEX,R_XDNYDN,8, &
           f(xeb-(mb-1)  :xeb                 , &
           yeb-(mb-1)  :yeb                 , &
           zs       :ze       )       , &
           mb*mb*zr,MPI_COMPLEX,R_XUPYUP,8, &
           C_CART,status,ierr)

    endif

    if ((exist_xup.and.exist_ydn).or.(exist_xdn.and.exist_yup)) then

      call MPI_SENDRECV(f(xe-(mb-1):xe               , &
           ys       :ys+(mb-1)        , &
           zs       :ze       )       , &
           mb*mb*zr,MPI_COMPLEX,R_XUPYDN,9, &
           f(xsb         :xsb+(mb-1)          , &
           yeb-(mb-1)  :yeb                 , &
           zs       :ze       )       , &
           mb*mb*zr,MPI_COMPLEX,R_XDNYUP,9, &
           C_CART,status,ierr)

      call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
           ye-(mb-1):ye               , &
           zs       :ze       )       , &
           mb*mb*zr,MPI_COMPLEX,R_XDNYUP,10, &
           f(xeb-(mb-1)  :xeb                 , &
           ysb         :ysb+(mb-1)          , &
           zs       :ze       )       , &
           mb*mb*zr,MPI_COMPLEX,R_XUPYDN,10, &
           C_CART,status,ierr)

    end if

  call trace%end
   end subroutine commc_xy

!************************************************************************

   subroutine comm2d4_xy (f)

    implicit none

    real(kind=4),dimension(xsb:xeb,ysb:yeb) :: f

  call trace%begin ('mpi.f90::comm2d4_xy')

    if ((mx /= 1).and.(my /= 1)) then
      if (exist_xupyup.or.exist_xdnydn) then

        call MPI_SENDRECV(f(xe-(mb-1):xe               , &
             ye-(mb-1):ye)       , &
             mb*mb,MPI_REAL,R_XUPYUP,7, &
             f(xsb         :xsb+(mb-1)          , &
             ysb         :ysb+(mb-1))       , &
             mb*mb,MPI_REAL,R_XDNYDN,7, &
             C_CART,status,ierr)

        call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
             ys       :ys+(mb-1))       , &
             mb*mb,MPI_REAL,R_XDNYDN,8, &
             f(xeb-(mb-1)  :xeb                 , &
             yeb-(mb-1)  :yeb)       , &
             mb*mb,MPI_REAL,R_XUPYUP,8, &
             C_CART,status,ierr)
      end if

      if ((exist_xup.and.exist_ydn).or.(exist_xdn.and.exist_yup)) then

        call MPI_SENDRECV(f(xe-(mb-1):xe               , &
             ys       :ys+(mb-1))       , &
             mb*mb,MPI_REAL,R_XUPYDN,9, &
             f(xsb         :xsb+(mb-1)          , &
             yeb-(mb-1)  :yeb)       , &
             mb*mb,MPI_REAL,R_XDNYUP,9, &
             C_CART,status,ierr)

        call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
             ye-(mb-1):ye)       , &
             mb*mb,MPI_REAL,R_XDNYUP,10, &
             f(xeb-(mb-1)  :xeb                 , &
             ysb         :ysb+(mb-1))       , &
             mb*mb,MPI_REAL,R_XUPYDN,10, &
             C_CART,status,ierr)

      end if
    end if
  call trace%end
   end subroutine comm2d4_xy

! ***********************************************************************

 subroutine comm2d8_xy (f)

  implicit none

  real(kind=8),dimension(xsb:xeb,ysb:yeb) :: f

  call trace%begin ('mpi.f90::comm2d8_xy')

    if ((exist_xup.and.exist_yup).or.(exist_xdn.and.exist_ydn)) then

      call MPI_SENDRECV(f(xe-(mb-1):xe               , &
           ye-(mb-1):ye)       , &
           mb*mb,MPI_DOUBLE_PRECISION,R_XUPYUP,7, &
           f(xsb         :xsb+(mb-1)          , &
           ysb         :ysb+(mb-1))       , &
           mb*mb,MPI_DOUBLE_PRECISION,R_XDNYDN,7, &
           C_CART,status,ierr)

      call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
           ys       :ys+(mb-1))       , &
           mb*mb,MPI_DOUBLE_PRECISION,R_XDNYDN,8, &
           f(xeb-(mb-1)  :xeb                 , &
           yeb-(mb-1)  :yeb)       , &
           mb*mb,MPI_DOUBLE_PRECISION,R_XUPYUP,8, &
           C_CART,status,ierr)

    endif

    if ((exist_xup.and.exist_ydn).or.(exist_xdn.and.exist_yup)) then

      call MPI_SENDRECV(f(xe-(mb-1):xe               , &
           ys       :ys+(mb-1))       , &
           mb*mb,MPI_DOUBLE_PRECISION,R_XUPYDN,9, &
           f(xsb         :xsb+(mb-1)          , &
           yeb-(mb-1)  :yeb)       , &
           mb*mb,MPI_DOUBLE_PRECISION,R_XDNYUP,9, &
           C_CART,status,ierr)

      call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
           ye-(mb-1):ye )       , &
           mb*mb,MPI_DOUBLE_PRECISION,R_XDNYUP,10, &
           f(xeb-(mb-1)  :xeb                 , &
           ysb         :ysb+(mb-1))       , &
           mb*mb,MPI_DOUBLE_PRECISION,R_XUPYDN,10, &
           C_CART,status,ierr)

    end if

  call trace%end
   end subroutine comm2d8_xy

! ***********************************************************************

 subroutine comm2d2_xy (f)

  implicit none

  integer,dimension(xsb:xeb,ysb:yeb) :: f

  call trace%begin ('mpi.f90::comm2d2_xy')

    if ((exist_xup.and.exist_yup).or.(exist_xdn.and.exist_ydn)) then

      call MPI_SENDRECV(f(xe-(mb-1):xe               , &
           ye-(mb-1):ye)       , &
           mb*mb,OSC_INTEGER,R_XUPYUP,7, &
           f(xsb         :xsb+(mb-1)          , &
           ysb         :ysb+(mb-1))       , &
           mb*mb,OSC_INTEGER,R_XDNYDN,7, &
           C_CART,status,ierr)

      call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
           ys       :ys+(mb-1) )       , &
           mb*mb,OSC_INTEGER,R_XDNYDN,8, &
           f(xeb-(mb-1)  :xeb                 , &
           yeb-(mb-1)  :yeb)       , &
           mb*mb,OSC_INTEGER,R_XUPYUP,8, &
           C_CART,status,ierr)

    endif

    if ((exist_xup.and.exist_ydn).or.(exist_xdn.and.exist_yup)) then

      call MPI_SENDRECV(f(xe-(mb-1):xe               , &
           ys       :ys+(mb-1))       , &
           mb*mb,OSC_INTEGER,R_XUPYDN,9, &
           f(xsb         :xsb+(mb-1)          , &
           yeb-(mb-1)  :yeb )       , &
           mb*mb,OSC_INTEGER,R_XDNYUP,9, &
           C_CART,status,ierr)

      call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
           ye-(mb-1):ye)       , &
           mb*mb,OSC_INTEGER,R_XDNYUP,10, &
           f(xeb-(mb-1)  :xeb                 , &
           ysb         :ysb+(mb-1))       , &
           mb*mb,OSC_INTEGER,R_XUPYDN,10, &
           C_CART,status,ierr)

    end if

  call trace%end
   end subroutine comm2d2_xy

!************************************************************************

 subroutine commc2d_xy (f)

  implicit none

  complex,dimension(xsb:xeb,ysb:yeb) :: f

  call trace%begin ('mpi.f90::commc2d_xy')

    if ((exist_xup.and.exist_yup).or.(exist_xdn.and.exist_ydn)) then

      call MPI_SENDRECV(f(xe-(mb-1):xe               , &
           ye-(mb-1):ye )       , &
           mb*mb,MPI_COMPLEX,R_XUPYUP,7, &
           f(xsb         :xsb+(mb-1)          , &
           ysb         :ysb+(mb-1))       , &
           mb*mb,MPI_COMPLEX,R_XDNYDN,7, &
           C_CART,status,ierr)

      call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
           ys       :ys+(mb-1))       , &
           mb*mb,MPI_COMPLEX,R_XDNYDN,8, &
           f(xeb-(mb-1)  :xeb                 , &
           yeb-(mb-1)  :yeb)       , &
           mb*mb,MPI_COMPLEX,R_XUPYUP,8, &
           C_CART,status,ierr)

    endif

    if ((exist_xup.and.exist_ydn).or.(exist_xdn.and.exist_yup)) then

      call MPI_SENDRECV(f(xe-(mb-1):xe               , &
           ys       :ys+(mb-1) )       , &
           mb*mb,MPI_COMPLEX,R_XUPYDN,9, &
           f(xsb         :xsb+(mb-1)          , &
           yeb-(mb-1)  :yeb )       , &
           mb*mb,MPI_COMPLEX,R_XDNYUP,9, &
           C_CART,status,ierr)

      call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
           ye-(mb-1):ye )       , &
           mb*mb,MPI_COMPLEX,R_XDNYUP,10, &
           f(xeb-(mb-1)  :xeb                 , &
           ysb         :ysb+(mb-1))       , &
           mb*mb,MPI_COMPLEX,R_XUPYDN,10, &
           C_CART,status,ierr)

    end if

  call trace%end
   end subroutine commc2d_xy

! ***********************************************************************


 subroutine comm4_xz (f)

  implicit none

  real(kind=4),dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::comm4_xz')

  if ((mx /= 1).and.(mz /= 1)) then
    if ((exist_xup.and.exist_zup).or.(exist_xdn.and.exist_zdn)) then

    call MPI_SENDRECV(f(xe-(mb-1):xe               , &
         ys       :ye               , &
         ze-(mb-1):ze       )       , &
         mb*yr*mb,MPI_REAL,R_XUPZUP,11, &
         f(xsb         :xsb+(mb-1)          , &
         ys       :ye               , &
         zsb         :zsb+(mb-1)  )       , &
         mb*yr*mb,MPI_REAL,R_XDNZDN,11, &
         C_CART,status,ierr)

    call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
         ys       :ye               , &
         zs       :zs+(mb-1))       , &
         mb*yr*mb,MPI_REAL,R_XDNZDN,12, &
         f(xeb-(mb-1)  :xeb                 , &
         ys       :ye               , &
         zeb-(mb-1)  :zeb         )       , &
         mb*yr*mb,MPI_REAL,R_XUPZUP,12, &
         C_CART,status,ierr)

  end if


  if ((exist_xup.and.exist_zdn).or.(exist_xdn.and.exist_zup)) then

    call MPI_SENDRECV(f(xe-(mb-1):xe               , &
         ys       :ye               , &
         zs       :zs+(mb-1))       , &
         mb*yr*mb,MPI_REAL,R_XUPZDN,13, &
         f(xsb         :xsb+(mb-1)          , &
         ys       :ye               , &
         zeb-(mb-1)  :zeb         )       , &
         mb*yr*mb,MPI_REAL,R_XDNZUP,13, &
         C_CART,status,ierr)

    call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
         ys       :ye               , &
         ze-(mb-1):ze       )       , &
         mb*yr*mb,MPI_REAL,R_XDNZUP,14, &
         f(xeb-(mb-1)  :xeb                 , &
         ys       :ye               , &
         zsb         :zsb+(mb-1)  )       , &
         mb*yr*mb,MPI_REAL,R_XUPZDN,14, &
         C_CART,status,ierr)

  end if
endif
  call trace%end
 end subroutine comm4_xz

! ***********************************************************************

 subroutine comm8_xz (f)

  implicit none

  real(kind=8),dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::comm8_xz')

  if ((exist_xup.and.exist_zup).or.(exist_xdn.and.exist_zdn)) then

    call MPI_SENDRECV(f(xe-(mb-1):xe               , &
         ys       :ye               , &
         ze-(mb-1):ze       )       , &
         mb*yr*mb,MPI_DOUBLE_PRECISION,R_XUPZUP,11, &
         f(xsb         :xsb+(mb-1)          , &
         ys       :ye               , &
         zsb         :zsb+(mb-1)  )       , &
         mb*yr*mb,MPI_DOUBLE_PRECISION,R_XDNZDN,11, &
         C_CART,status,ierr)

    call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
         ys       :ye               , &
         zs       :zs+(mb-1))       , &
         mb*yr*mb,MPI_DOUBLE_PRECISION,R_XDNZDN,12, &
         f(xeb-(mb-1)  :xeb                 , &
         ys       :ye               , &
         zeb-(mb-1)  :zeb         )       , &
         mb*yr*mb,MPI_DOUBLE_PRECISION,R_XUPZUP,12, &
         C_CART,status,ierr)

  end if


  if ((exist_xup.and.exist_zdn).or.(exist_xdn.and.exist_zup)) then

    call MPI_SENDRECV(f(xe-(mb-1):xe               , &
         ys       :ye               , &
         zs       :zs+(mb-1))       , &
         mb*yr*mb,MPI_DOUBLE_PRECISION,R_XUPZDN,13, &
         f(xsb         :xsb+(mb-1)          , &
         ys       :ye               , &
         zeb-(mb-1)  :zeb         )       , &
         mb*yr*mb,MPI_DOUBLE_PRECISION,R_XDNZUP,13, &
         C_CART,status,ierr)

    call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
         ys       :ye               , &
         ze-(mb-1):ze       )       , &
         mb*yr*mb,MPI_DOUBLE_PRECISION,R_XDNZUP,14, &
         f(xeb-(mb-1)  :xeb                 , &
         ys       :ye               , &
         zsb         :zsb+(mb-1)  )       , &
         mb*yr*mb,MPI_DOUBLE_PRECISION,R_XUPZDN,14, &
         C_CART,status,ierr)

  end if

  call trace%end
 end subroutine comm8_xz

! ***********************************************************************

 subroutine comm2_xz (f)

  implicit none

  integer,dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::comm2_xz')

  if ((exist_xup.and.exist_zup).or.(exist_xdn.and.exist_zdn)) then

    call MPI_SENDRECV(f(xe-(mb-1):xe               , &
         ys       :ye               , &
         ze-(mb-1):ze       )       , &
         mb*yr*mb,OSC_INTEGER,R_XUPZUP,11, &
         f(xsb         :xsb+(mb-1)          , &
         ys       :ye               , &
         zsb         :zsb+(mb-1)  )       , &
         mb*yr*mb,OSC_INTEGER,R_XDNZDN,11, &
         C_CART,status,ierr)

    call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
         ys       :ye               , &
         zs       :zs+(mb-1))       , &
         mb*yr*mb,OSC_INTEGER,R_XDNZDN,12, &
         f(xeb-(mb-1)  :xeb                 , &
         ys       :ye               , &
         zeb-(mb-1)  :zeb         )       , &
         mb*yr*mb,OSC_INTEGER,R_XUPZUP,12, &
         C_CART,status,ierr)

  end if


  if ((exist_xup.and.exist_zdn).or.(exist_xdn.and.exist_zup)) then

    call MPI_SENDRECV(f(xe-(mb-1):xe               , &
         ys       :ye               , &
         zs       :zs+(mb-1))       , &
         mb*yr*mb,OSC_INTEGER,R_XUPZDN,13, &
         f(xsb         :xsb+(mb-1)          , &
         ys       :ye               , &
         zeb-(mb-1)  :zeb         )       , &
         mb*yr*mb,OSC_INTEGER,R_XDNZUP,13, &
         C_CART,status,ierr)

    call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
         ys       :ye               , &
         ze-(mb-1):ze       )       , &
         mb*yr*mb,OSC_INTEGER,R_XDNZUP,14, &
         f(xeb-(mb-1)  :xeb                 , &
         ys       :ye               , &
         zsb         :zsb+(mb-1)  )       , &
         mb*yr*mb,OSC_INTEGER,R_XUPZDN,14, &
         C_CART,status,ierr)

  end if

  call trace%end
 end subroutine comm2_xz

! ***********************************************************************

 subroutine commc_xz (f)

  implicit none

  complex,dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::commc_xz')

  if ((exist_xup.and.exist_zup).or.(exist_xdn.and.exist_zdn)) then

    call MPI_SENDRECV(f(xe-(mb-1):xe               , &
         ys       :ye               , &
         ze-(mb-1):ze       )       , &
         mb*yr*mb,MPI_COMPLEX,R_XUPZUP,11, &
         f(xsb         :xsb+(mb-1)          , &
         ys       :ye               , &
         zsb         :zsb+(mb-1)  )       , &
         mb*yr*mb,MPI_COMPLEX,R_XDNZDN,11, &
         C_CART,status,ierr)

    call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
         ys       :ye               , &
         zs       :zs+(mb-1))       , &
         mb*yr*mb,MPI_COMPLEX,R_XDNZDN,12, &
         f(xeb-(mb-1)  :xeb                 , &
         ys       :ye               , &
         zeb-(mb-1)  :zeb         )       , &
         mb*yr*mb,MPI_COMPLEX,R_XUPZUP,12, &
         C_CART,status,ierr)

  end if


  if ((exist_xup.and.exist_zdn).or.(exist_xdn.and.exist_zup)) then

    call MPI_SENDRECV(f(xe-(mb-1):xe               , &
         ys       :ye               , &
         zs       :zs+(mb-1))       , &
         mb*yr*mb,MPI_COMPLEX,R_XUPZDN,13, &
         f(xsb         :xsb+(mb-1)          , &
         ys       :ye               , &
         zeb-(mb-1)  :zeb         )       , &
         mb*yr*mb,MPI_COMPLEX,R_XDNZUP,13, &
         C_CART,status,ierr)

    call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
         ys       :ye               , &
         ze-(mb-1):ze       )       , &
         mb*yr*mb,MPI_COMPLEX,R_XDNZUP,14, &
         f(xeb-(mb-1)  :xeb                 , &
         ys       :ye               , &
         zsb         :zsb+(mb-1)  )       , &
         mb*yr*mb,MPI_COMPLEX,R_XUPZDN,14, &
         C_CART,status,ierr)

  end if

  call trace%end
 end subroutine commc_xz

! ***********************************************************************

 subroutine comm4_yz (f)

  implicit none

  real(kind=4),dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::comm4_yz')

  if ((my /= 1) .and.(mz /= 1)) then
    if ((exist_yup.and.exist_zup).or.(exist_ydn.and.exist_zdn)) then

    call MPI_SENDRECV(f(xs       :xe               , &
         ye-(mb-1):ye               , &
         ze-(mb-1):ze       )       , &
         xr*mb*mb,MPI_REAL,R_YUPZUP,15, &
         f(xs       :xe               , &
         ysb         :ysb+(mb-1)          , &
         zsb         :zsb+(mb-1)  )       , &
         xr*mb*mb,MPI_REAL,R_YDNZDN,15, &
         C_CART,status,ierr)

    call MPI_SENDRECV(f(xs       :xe               , &
         ys       :ys+(mb-1)        , &
         zs       :zs+(mb-1))       , &
         xr*mb*mb,MPI_REAL,R_YDNZDN,16, &
         f(xs       :xe               , &
         yeb-(mb-1)  :yeb                 , &
         zeb-(mb-1)  :zeb         )       , &
         xr*mb*mb,MPI_REAL,R_YUPZUP,16 , &
         C_CART,status,ierr)

  end if

  if ((exist_yup.and.exist_zdn).or.(exist_ydn.and.exist_zup)) then

    call MPI_SENDRECV(f(xs       :xe               , &
         ye-(mb-1):ye               , &
         zs       :zs+(mb-1))       , &
         xr*mb*mb,MPI_REAL,R_YUPZDN,17, &
         f(xs       :xe               , &
         ysb         :ysb+(mb-1)          , &
         zeb-(mb-1)  :zeb         )       , &
         xr*mb*mb,MPI_REAL,R_YDNZUP,17, &
         C_CART,status,ierr)

    call MPI_SENDRECV(f(xs       :xe               , &
         ys       :ys+(mb-1)        , &
         ze-(mb-1):ze       )       , &
         xr*mb*mb,MPI_REAL,R_YDNZUP,18                   , &
         f(xs       :xe               , &
         yeb-(mb-1)  :yeb                 , &
         zsb         :zsb+(mb-1)  )       , &
         xr*mb*mb,MPI_REAL,R_YUPZDN,18                   , &
         C_CART,status,ierr)

  end if
end if
  call trace%end
 end subroutine comm4_yz

! ***********************************************************************

 subroutine comm8_yz (f)

  implicit none

  real(kind=8),dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::comm8_yz')

  if ((exist_yup.and.exist_zup).or.(exist_ydn.and.exist_zdn)) then

    call MPI_SENDRECV(f(xs       :xe               , &
         ye-(mb-1):ye               , &
         ze-(mb-1):ze       )       , &
         xr*mb*mb,MPI_DOUBLE_PRECISION,R_YUPZUP,15, &
         f(xs       :xe               , &
         ysb         :ysb+(mb-1)          , &
         zsb         :zsb+(mb-1)  )       , &
         xr*mb*mb,MPI_DOUBLE_PRECISION,R_YDNZDN,15, &
         C_CART,status,ierr)

    call MPI_SENDRECV(f(xs       :xe               , &
         ys       :ys+(mb-1)        , &
         zs       :zs+(mb-1))       , &
         xr*mb*mb,MPI_DOUBLE_PRECISION,R_YDNZDN,16, &
         f(xs       :xe               , &
         yeb-(mb-1)  :yeb                 , &
         zeb-(mb-1)  :zeb         )       , &
         xr*mb*mb,MPI_DOUBLE_PRECISION,R_YUPZUP,16 , &
         C_CART,status,ierr)

  end if

  if ((exist_yup.and.exist_zdn).or.(exist_ydn.and.exist_zup)) then

    call MPI_SENDRECV(f(xs       :xe               , &
         ye-(mb-1):ye               , &
         zs       :zs+(mb-1))       , &
         xr*mb*mb,MPI_DOUBLE_PRECISION,R_YUPZDN,17, &
         f(xs       :xe               , &
         ysb         :ysb+(mb-1)          , &
         zeb-(mb-1)  :zeb         )       , &
         xr*mb*mb,MPI_DOUBLE_PRECISION,R_YDNZUP,17, &
         C_CART,status,ierr)

    call MPI_SENDRECV(f(xs       :xe               , &
         ys       :ys+(mb-1)        , &
         ze-(mb-1):ze       )       , &
         xr*mb*mb,MPI_DOUBLE_PRECISION,R_YDNZUP,18                   , &
         f(xs       :xe               , &
         yeb-(mb-1)  :yeb                 , &
         zsb         :zsb+(mb-1)  )       , &
         xr*mb*mb,MPI_DOUBLE_PRECISION,R_YUPZDN,18                   , &
         C_CART,status,ierr)

  end if

  call trace%end
 end subroutine comm8_yz

! ***********************************************************************

 subroutine comm2_yz (f)

  implicit none

  integer,dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::comm2_yz')

  if ((exist_yup.and.exist_zup).or.(exist_ydn.and.exist_zdn)) then

    call MPI_SENDRECV(f(xs       :xe               , &
         ye-(mb-1):ye               , &
         ze-(mb-1):ze       )       , &
         xr*mb*mb,OSC_INTEGER,R_YUPZUP,15, &
         f(xs       :xe               , &
         ysb         :ysb+(mb-1)          , &
         zsb         :zsb+(mb-1)  )       , &
         xr*mb*mb,OSC_INTEGER,R_YDNZDN,15, &
         C_CART,status,ierr)

    call MPI_SENDRECV(f(xs       :xe               , &
         ys       :ys+(mb-1)        , &
         zs       :zs+(mb-1))       , &
         xr*mb*mb,OSC_INTEGER,R_YDNZDN,16, &
         f(xs       :xe               , &
         yeb-(mb-1)  :yeb                 , &
         zeb-(mb-1)  :zeb         )       , &
         xr*mb*mb,OSC_INTEGER,R_YUPZUP,16 , &
         C_CART,status,ierr)

  end if

  if ((exist_yup.and.exist_zdn).or.(exist_ydn.and.exist_zup)) then

    call MPI_SENDRECV(f(xs       :xe               , &
         ye-(mb-1):ye               , &
         zs       :zs+(mb-1))       , &
         xr*mb*mb,OSC_INTEGER,R_YUPZDN,17, &
         f(xs       :xe               , &
         ysb         :ysb+(mb-1)          , &
         zeb-(mb-1)  :zeb         )       , &
         xr*mb*mb,OSC_INTEGER,R_YDNZUP,17, &
         C_CART,status,ierr)

    call MPI_SENDRECV(f(xs       :xe               , &
         ys       :ys+(mb-1)        , &
         ze-(mb-1):ze       )       , &
         xr*mb*mb,OSC_INTEGER,R_YDNZUP,18                   , &
         f(xs       :xe               , &
         yeb-(mb-1)  :yeb                 , &
         zsb         :zsb+(mb-1)  )       , &
         xr*mb*mb,OSC_INTEGER,R_YUPZDN,18                   , &
         C_CART,status,ierr)

  end if

  call trace%end
 end subroutine comm2_yz

! ***********************************************************************

 subroutine commc_yz (f)

  implicit none

  complex,dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::commc_yz')

  if ((exist_yup.and.exist_zup).or.(exist_ydn.and.exist_zdn)) then

    call MPI_SENDRECV(f(xs       :xe               , &
         ye-(mb-1):ye               , &
         ze-(mb-1):ze       )       , &
         xr*mb*mb,MPI_COMPLEX,R_YUPZUP,15, &
         f(xs       :xe               , &
         ysb         :ysb+(mb-1)          , &
         zsb         :zsb+(mb-1)  )       , &
         xr*mb*mb,MPI_COMPLEX,R_YDNZDN,15, &
         C_CART,status,ierr)

    call MPI_SENDRECV(f(xs       :xe               , &
         ys       :ys+(mb-1)        , &
         zs       :zs+(mb-1))       , &
         xr*mb*mb,MPI_COMPLEX,R_YDNZDN,16, &
         f(xs       :xe               , &
         yeb-(mb-1)  :yeb                 , &
         zeb-(mb-1)  :zeb         )       , &
         xr*mb*mb,MPI_COMPLEX,R_YUPZUP,16 , &
         C_CART,status,ierr)

  end if

  if ((exist_yup.and.exist_zdn).or.(exist_ydn.and.exist_zup)) then

    call MPI_SENDRECV(f(xs       :xe               , &
         ye-(mb-1):ye               , &
         zs       :zs+(mb-1))       , &
         xr*mb*mb,MPI_COMPLEX,R_YUPZDN,17, &
         f(xs       :xe               , &
         ysb         :ysb+(mb-1)          , &
         zeb-(mb-1)  :zeb         )       , &
         xr*mb*mb,MPI_COMPLEX,R_YDNZUP,17, &
         C_CART,status,ierr)

    call MPI_SENDRECV(f(xs       :xe               , &
         ys       :ys+(mb-1)        , &
         ze-(mb-1):ze       )       , &
         xr*mb*mb,MPI_COMPLEX,R_YDNZUP,18                   , &
         f(xs       :xe               , &
         yeb-(mb-1)  :yeb                 , &
         zsb         :zsb+(mb-1)  )       , &
         xr*mb*mb,MPI_COMPLEX,R_YUPZDN,18                   , &
         C_CART,status,ierr)

  end if

  call trace%end
 end subroutine commc_yz

! ***********************************************************************

 subroutine comm4_xyz (f)

  implicit none

  real(kind=4),dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::comm4_xyz')

  if (exist_xupyupzup.or.exist_xdnydnzdn) then

      call MPI_SENDRECV(f(xe-(mb-1):xe               , &
           ye-(mb-1):ye               , &
           ze-(mb-1):ze       )       , &
           mb*mb*mb,MPI_REAL,R_XUPYUPZUP,19, &
           f(xsb         :xsb+(mb-1)          , &
           ysb         :ysb+(mb-1)          , &
           zsb         :zsb+(mb-1)  )       , &
           mb*mb*mb,MPI_REAL,R_XDNYDNZDN,19, &
           C_CART,status,ierr)

      call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
           ys       :ys+(mb-1)        , &
           zs       :zs+(mb-1))       , &
           mb*mb*mb,MPI_REAL,R_XDNYDNZDN,20, &
           f(xeb-(mb-1)  :xeb                 , &
           yeb-(mb-1)  :yeb                 , &
           zeb-(mb-1)  :zeb         )       , &
           mb*mb*mb,MPI_REAL,R_XUPYUPZUP,20 , &
           C_CART,status,ierr)

    end if

    if ((exist_xup.and.exist_yup.and.exist_zdn).or.(exist_xdn.and.exist_ydn.and.exist_zup)) then

      call MPI_SENDRECV(f(xe-(mb-1):xe               , &
           ye-(mb-1):ye               , &
           zs       :zs+(mb-1))       , &
           mb*mb*mb,MPI_REAL,R_XUPYUPZDN,21, &
           f(xsb         :xsb+(mb-1)           , &
           ysb         :ysb+(mb-1)           , &
           zeb-(mb-1)  :zeb         )        , &
           mb*mb*mb,MPI_REAL,R_XDNYDNZUP,21, &
           C_CART,status,ierr)

      call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
           ys       :ys+(mb-1)        , &
           ze-(mb-1):ze       )       , &
           mb*mb*mb,MPI_REAL,R_XDNYDNZUP,22, &
           f(xeb-(mb-1)  :xeb                 , &
           yeb-(mb-1)  :yeb                 , &
           zsb         :zsb+(mb-1)  )       , &
           mb*mb*mb,MPI_REAL,R_XUPYUPZDN,22, &
           C_CART,status,ierr)

    end if

    if ((exist_xup.and.exist_ydn.and.exist_zup).or.(exist_xdn.and.exist_yup.and.exist_zdn)) then

      call MPI_SENDRECV(f(xe-(mb-1):xe               , &
           ys       :ys+(mb-1)        , &
           ze-(mb-1):ze       )       , &
           mb*mb*mb,MPI_REAL,R_XUPYDNZUP,23, &
           f(xsb         :xsb+(mb-1)          , &
           yeb-(mb-1)  :yeb                 , &
           zsb         :zsb+(mb-1)  )       , &
           mb*mb*mb,MPI_REAL,R_XDNYUPZDN,23, &
           C_CART,status,ierr)

      call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
           ye-(mb-1):ye               , &
           zs       :zs+(mb-1))       , &
           mb*mb*mb,MPI_REAL,R_XDNYUPZDN,24, &
           f(xeb-(mb-1)  :xeb                 , &
           ysb         :ysb+(mb-1)          , &
           zeb-(mb-1)  :zeb         )       , &
           mb*mb*mb,MPI_REAL,R_XUPYDNZUP,24, &
           C_CART,status,ierr)

    end if

    if ((exist_xup.and.exist_ydn.and.exist_zdn).or.(exist_xdn.and.exist_yup.and.exist_zup)) then

      call MPI_SENDRECV(f(xe-(mb-1):xe               , &
           ys       :ys+(mb-1)        , &
           zs       :zs+(mb-1))       , &
           mb*mb*mb,MPI_REAL,R_XUPYDNZDN,25, &
           f(xsb         :xsb+(mb-1)          , &
           yeb-(mb-1)  :yeb                 , &
           zeb-(mb-1)  :zeb         )       , &
           mb*mb*mb,MPI_REAL,R_XDNYUPZUP,25, &
           C_CART,status,ierr)

      call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
           ye-(mb-1):ye               , &
           ze-(mb-1):ze       )       , &
           mb*mb*mb,MPI_REAL,R_XDNYUPZUP,26, &
           f(xeb-(mb-1)  :xeb                 , &
           ysb         :ysb+(mb-1)          , &
           zsb         :zsb+(mb-1)  )       , &
           mb*mb*mb,MPI_REAL,R_XUPYDNZDN,26, &
           C_CART,status,ierr)

    end if

  call trace%end
   end subroutine comm4_xyz

! ***********************************************************************

 subroutine comm8_xyz (f)

  implicit none

  real(kind=8),dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::comm8_xyz')

  if ((exist_xup.and.exist_yup.and.exist_zup).or.(exist_xdn.and.exist_ydn.and.exist_zdn)) then

    call MPI_SENDRECV(f(xe-(mb-1):xe               , &
         ye-(mb-1):ye               , &
         ze-(mb-1):ze       )       , &
         mb*mb*mb,MPI_DOUBLE_PRECISION,R_XUPYUPZUP,19, &
         f(xsb         :xsb+(mb-1)          , &
         ysb         :ysb+(mb-1)          , &
         zsb         :zsb+(mb-1)  )       , &
         mb*mb*mb,MPI_DOUBLE_PRECISION,R_XDNYDNZDN,19, &
         C_CART,status,ierr)

    call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
         ys       :ys+(mb-1)        , &
         zs       :zs+(mb-1))       , &
         mb*mb*mb,MPI_DOUBLE_PRECISION,R_XDNYDNZDN,20, &
         f(xeb-(mb-1)  :xeb                 , &
         yeb-(mb-1)  :yeb                 , &
         zeb-(mb-1)  :zeb         )       , &
         mb*mb*mb,MPI_DOUBLE_PRECISION,R_XUPYUPZUP,20 , &
         C_CART,status,ierr)

  end if

  if ((exist_xup.and.exist_yup.and.exist_zdn).or.(exist_xdn.and.exist_ydn.and.exist_zup)) then

    call MPI_SENDRECV(f(xe-(mb-1):xe               , &
         ye-(mb-1):ye               , &
         zs       :zs+(mb-1))       , &
         mb*mb*mb,MPI_DOUBLE_PRECISION,R_XUPYUPZDN,21, &
         f(xsb         :xsb+(mb-1)           , &
         ysb         :ysb+(mb-1)           , &
         zeb-(mb-1)  :zeb         )        , &
         mb*mb*mb,MPI_DOUBLE_PRECISION,R_XDNYDNZUP,21, &
         C_CART,status,ierr)

    call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
         ys       :ys+(mb-1)        , &
         ze-(mb-1):ze       )       , &
         mb*mb*mb,MPI_DOUBLE_PRECISION,R_XDNYDNZUP,22, &
         f(xeb-(mb-1)  :xeb                 , &
         yeb-(mb-1)  :yeb                 , &
         zsb         :zsb+(mb-1)  )       , &
         mb*mb*mb,MPI_DOUBLE_PRECISION,R_XUPYUPZDN,22, &
         C_CART,status,ierr)

  end if

  if ((exist_xup.and.exist_ydn.and.exist_zup).or.(exist_xdn.and.exist_yup.and.exist_zdn)) then

    call MPI_SENDRECV(f(xe-(mb-1):xe               , &
         ys       :ys+(mb-1)        , &
         ze-(mb-1):ze       )       , &
         mb*mb*mb,MPI_DOUBLE_PRECISION,R_XUPYDNZUP,23, &
         f(xsb         :xsb+(mb-1)          , &
         yeb-(mb-1)  :yeb                 , &
         zsb         :zsb+(mb-1)  )       , &
         mb*mb*mb,MPI_DOUBLE_PRECISION,R_XDNYUPZDN,23, &
         C_CART,status,ierr)

    call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
         ye-(mb-1):ye               , &
         zs       :zs+(mb-1))       , &
         mb*mb*mb,MPI_DOUBLE_PRECISION,R_XDNYUPZDN,24, &
         f(xeb-(mb-1)  :xeb                 , &
         ysb         :ysb+(mb-1)          , &
         zeb-(mb-1)  :zeb         )       , &
         mb*mb*mb,MPI_DOUBLE_PRECISION,R_XUPYDNZUP,24, &
         C_CART,status,ierr)

  end if

  if ((exist_xup.and.exist_ydn.and.exist_zdn).or.(exist_xdn.and.exist_yup.and.exist_zup)) then

    call MPI_SENDRECV(f(xe-(mb-1):xe               , &
         ys       :ys+(mb-1)        , &
         zs       :zs+(mb-1))       , &
         mb*mb*mb,MPI_DOUBLE_PRECISION,R_XUPYDNZDN,25, &
         f(xsb         :xsb+(mb-1)          , &
         yeb-(mb-1)  :yeb                 , &
         zeb-(mb-1)  :zeb         )       , &
         mb*mb*mb,MPI_DOUBLE_PRECISION,R_XDNYUPZUP,25, &
         C_CART,status,ierr)

    call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
         ye-(mb-1):ye               , &
         ze-(mb-1):ze       )       , &
         mb*mb*mb,MPI_DOUBLE_PRECISION,R_XDNYUPZUP,26, &
         f(xeb-(mb-1)  :xeb                 , &
         ysb         :ysb+(mb-1)          , &
         zsb         :zsb+(mb-1)  )       , &
         mb*mb*mb,MPI_DOUBLE_PRECISION,R_XUPYDNZDN,26, &
         C_CART,status,ierr)

  end if

  call trace%end
 end subroutine comm8_xyz

! ***********************************************************************

 subroutine comm2_xyz (f)

  implicit none

  integer,dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::comm2_xyz')

  if ((exist_xup.and.exist_yup.and.exist_zup).or.(exist_xdn.and.exist_ydn.and.exist_zdn)) then

    call MPI_SENDRECV(f(xe-(mb-1):xe               , &
         ye-(mb-1):ye               , &
         ze-(mb-1):ze       )       , &
         mb*mb*mb,OSC_INTEGER,R_XUPYUPZUP,19, &
         f(xsb         :xsb+(mb-1)          , &
         ysb         :ysb+(mb-1)          , &
         zsb         :zsb+(mb-1)  )       , &
         mb*mb*mb,OSC_INTEGER,R_XDNYDNZDN,19, &
         C_CART,status,ierr)

    call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
         ys       :ys+(mb-1)        , &
         zs       :zs+(mb-1))       , &
         mb*mb*mb,OSC_INTEGER,R_XDNYDNZDN,20, &
         f(xeb-(mb-1)  :xeb                 , &
         yeb-(mb-1)  :yeb                 , &
         zeb-(mb-1)  :zeb         )       , &
         mb*mb*mb,OSC_INTEGER,R_XUPYUPZUP,20 , &
         C_CART,status,ierr)

  end if

  if ((exist_xup.and.exist_yup.and.exist_zdn).or.(exist_xdn.and.exist_ydn.and.exist_zup)) then

    call MPI_SENDRECV(f(xe-(mb-1):xe               , &
         ye-(mb-1):ye               , &
         zs       :zs+(mb-1))       , &
         mb*mb*mb,OSC_INTEGER,R_XUPYUPZDN,21, &
         f(xsb         :xsb+(mb-1)           , &
         ysb         :ysb+(mb-1)           , &
         zeb-(mb-1)  :zeb         )        , &
         mb*mb*mb,OSC_INTEGER,R_XDNYDNZUP,21, &
         C_CART,status,ierr)

    call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
         ys       :ys+(mb-1)        , &
         ze-(mb-1):ze       )       , &
         mb*mb*mb,OSC_INTEGER,R_XDNYDNZUP,22, &
         f(xeb-(mb-1)  :xeb                 , &
         yeb-(mb-1)  :yeb                 , &
         zsb         :zsb+(mb-1)  )       , &
         mb*mb*mb,OSC_INTEGER,R_XUPYUPZDN,22, &
         C_CART,status,ierr)

  end if

  if ((exist_xup.and.exist_ydn.and.exist_zup).or.(exist_xdn.and.exist_yup.and.exist_zdn)) then

    call MPI_SENDRECV(f(xe-(mb-1):xe               , &
         ys       :ys+(mb-1)        , &
         ze-(mb-1):ze       )       , &
         mb*mb*mb,OSC_INTEGER,R_XUPYDNZUP,23, &
         f(xsb         :xsb+(mb-1)          , &
         yeb-(mb-1)  :yeb                 , &
         zsb         :zsb+(mb-1)  )       , &
         mb*mb*mb,OSC_INTEGER,R_XDNYUPZDN,23, &
         C_CART,status,ierr)

    call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
         ye-(mb-1):ye               , &
         zs       :zs+(mb-1))       , &
         mb*mb*mb,OSC_INTEGER,R_XDNYUPZDN,24, &
         f(xeb-(mb-1)  :xeb                 , &
         ysb         :ysb+(mb-1)          , &
         zeb-(mb-1)  :zeb         )       , &
         mb*mb*mb,OSC_INTEGER,R_XUPYDNZUP,24, &
         C_CART,status,ierr)

  end if

  if ((exist_xup.and.exist_ydn.and.exist_zdn).or.(exist_xdn.and.exist_yup.and.exist_zup)) then

    call MPI_SENDRECV(f(xe-(mb-1):xe               , &
         ys       :ys+(mb-1)        , &
         zs       :zs+(mb-1))       , &
         mb*mb*mb,OSC_INTEGER,R_XUPYDNZDN,25, &
         f(xsb         :xsb+(mb-1)          , &
         yeb-(mb-1)  :yeb                 , &
         zeb-(mb-1)  :zeb         )       , &
         mb*mb*mb,OSC_INTEGER,R_XDNYUPZUP,25, &
         C_CART,status,ierr)

    call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
         ye-(mb-1):ye               , &
         ze-(mb-1):ze       )       , &
         mb*mb*mb,OSC_INTEGER,R_XDNYUPZUP,26, &
         f(xeb-(mb-1)  :xeb                 , &
         ysb         :ysb+(mb-1)          , &
         zsb         :zsb+(mb-1)  )       , &
         mb*mb*mb,OSC_INTEGER,R_XUPYDNZDN,26, &
         C_CART,status,ierr)

  end if

  call trace%end
 end subroutine comm2_xyz

! ***********************************************************************

 subroutine commc_xyz (f)

  implicit none

  complex,dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::commc_xyz')

  if ((exist_xup.and.exist_yup.and.exist_zup).or.(exist_xdn.and.exist_ydn.and.exist_zdn)) then

    call MPI_SENDRECV(f(xe-(mb-1):xe               , &
         ye-(mb-1):ye               , &
         ze-(mb-1):ze       )       , &
         mb*mb*mb,MPI_COMPLEX,R_XUPYUPZUP,19, &
         f(xsb         :xsb+(mb-1)          , &
         ysb         :ysb+(mb-1)          , &
         zsb         :zsb+(mb-1)  )       , &
         mb*mb*mb,MPI_COMPLEX,R_XDNYDNZDN,19, &
         C_CART,status,ierr)

    call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
         ys       :ys+(mb-1)        , &
         zs       :zs+(mb-1))       , &
         mb*mb*mb,MPI_COMPLEX,R_XDNYDNZDN,20, &
         f(xeb-(mb-1)  :xeb                 , &
         yeb-(mb-1)  :yeb                 , &
         zeb-(mb-1)  :zeb         )       , &
         mb*mb*mb,MPI_COMPLEX,R_XUPYUPZUP,20 , &
         C_CART,status,ierr)

  end if

  if ((exist_xup.and.exist_yup.and.exist_zdn).or.(exist_xdn.and.exist_ydn.and.exist_zup)) then

    call MPI_SENDRECV(f(xe-(mb-1):xe               , &
         ye-(mb-1):ye               , &
         zs       :zs+(mb-1))       , &
         mb*mb*mb,MPI_COMPLEX,R_XUPYUPZDN,21, &
         f(xsb         :xsb+(mb-1)           , &
         ysb         :ysb+(mb-1)           , &
         zeb-(mb-1)  :zeb         )        , &
         mb*mb*mb,MPI_COMPLEX,R_XDNYDNZUP,21, &
         C_CART,status,ierr)

    call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
         ys       :ys+(mb-1)        , &
         ze-(mb-1):ze       )       , &
         mb*mb*mb,MPI_COMPLEX,R_XDNYDNZUP,22, &
         f(xeb-(mb-1)  :xeb                 , &
         yeb-(mb-1)  :yeb                 , &
         zsb         :zsb+(mb-1)  )       , &
         mb*mb*mb,MPI_COMPLEX,R_XUPYUPZDN,22, &
         C_CART,status,ierr)

  end if

  if ((exist_xup.and.exist_ydn.and.exist_zup).or.(exist_xdn.and.exist_yup.and.exist_zdn)) then

    call MPI_SENDRECV(f(xe-(mb-1):xe               , &
         ys       :ys+(mb-1)        , &
         ze-(mb-1):ze       )       , &
         mb*mb*mb,MPI_COMPLEX,R_XUPYDNZUP,23, &
         f(xsb         :xsb+(mb-1)          , &
         yeb-(mb-1)  :yeb                 , &
         zsb         :zsb+(mb-1)  )       , &
         mb*mb*mb,MPI_COMPLEX,R_XDNYUPZDN,23, &
         C_CART,status,ierr)

    call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
         ye-(mb-1):ye               , &
         zs       :zs+(mb-1))       , &
         mb*mb*mb,MPI_COMPLEX,R_XDNYUPZDN,24, &
         f(xeb-(mb-1)  :xeb                 , &
         ysb         :ysb+(mb-1)          , &
         zeb-(mb-1)  :zeb         )       , &
         mb*mb*mb,MPI_COMPLEX,R_XUPYDNZUP,24, &
         C_CART,status,ierr)

  end if

  if ((exist_xup.and.exist_ydn.and.exist_zdn).or.(exist_xdn.and.exist_yup.and.exist_zup)) then

    call MPI_SENDRECV(f(xe-(mb-1):xe               , &
         ys       :ys+(mb-1)        , &
         zs       :zs+(mb-1))       , &
         mb*mb*mb,MPI_COMPLEX,R_XUPYDNZDN,25, &
         f(xsb         :xsb+(mb-1)          , &
         yeb-(mb-1)  :yeb                 , &
         zeb-(mb-1)  :zeb         )       , &
         mb*mb*mb,MPI_COMPLEX,R_XDNYUPZUP,25, &
         C_CART,status,ierr)

    call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
         ye-(mb-1):ye               , &
         ze-(mb-1):ze       )       , &
         mb*mb*mb,MPI_COMPLEX,R_XDNYUPZUP,26, &
         f(xeb-(mb-1)  :xeb                 , &
         ysb         :ysb+(mb-1)          , &
         zsb         :zsb+(mb-1)  )       , &
         mb*mb*mb,MPI_COMPLEX,R_XUPYDNZDN,26, &
         C_CART,status,ierr)

  end if

  call trace%end
 end subroutine commc_xyz

! ***********************************************************************
 subroutine comm_bdry_x (f)

  implicit none

  real(kind=4),dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::comm_bdry_x')


  if (exist_xup.or.exist_xdn) then

    if ((.not.exist_zdn).and.(zsb /= zs)) then

      call MPI_SENDRECV(f(xe-(mb-1) :xe         , &
           ys        :ye         , &
           zs-mb     :zs-1      ), &
           mb*yr*mb,MPI_REAL,R_XUP,1  , &
           f(xsb     :xsb+(mb-1) , &
           ys        :ye         , &
           zs-mb     :zs-1      ), &
           mb*yr*mb,MPI_REAL,R_XDN,1  , &
           C_CART,status,ierr)

      call MPI_SENDRECV(f(xs        :xs+(mb-1)  , &
           ys        :ye         , &
           zs-mb     :zs-1      ), &
           mb*yr*mb,MPI_REAL,R_XDN,2  , &
           f(xeb-(mb-1):xeb      , &
           ys        :ye         , &
           zs-mb     :zs-1      ), &
           mb*yr*mb,MPI_REAL,R_XUP,2  , &
           C_CART,status,ierr)

    endif

    if ((.not.exist_zup).and.(zeb /= ze)) then

      call MPI_SENDRECV(f(xe-(mb-1) :xe         , &
           ys        :ye         , &
           ze+1      :ze+mb     ), &
           mb*yr*mb,MPI_REAL,R_XUP,3  , &
           f(xsb     :xsb+(mb-1) , &
           ys        :ye         , &
           ze+1      :ze+mb     ), &
           mb*yr*mb,MPI_REAL,R_XDN,3  , &
           C_CART,status,ierr)

      call MPI_SENDRECV(f(xs        :xs+(mb-1)  , &
           ys        :ye         , &
           ze+1      :ze+mb     ), &
           mb*yr*mb,MPI_REAL,R_XDN,4  , &
           f(xeb-(mb-1):xeb      , &
           ys        :ye         , &
           ze+1      :ze+mb     ), &
           mb*yr*mb,MPI_REAL,R_XUP,4  , &
           C_CART,status,ierr)

    end if

  end if

  call trace%end
 end subroutine comm_bdry_x
! ***********************************************************************
 subroutine comm_bdry_y (f)

  implicit none

  real(kind=4),dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::comm_bdry_y')

  if (exist_yup.or.exist_ydn) then

      if ((.not.exist_zdn).and.(zsb /= zs)) then

        call MPI_SENDRECV(f(xs       :xe              , &
                         ye-(mb-1)   :ye              , &
                         zs-mb       :zs-1     )      , &
               xr*mb*mb,MPI_REAL,R_YUP,5, &
                          f(xs       :xe              , &
                          ysb        :ysb+(mb-1)      , &
                          zs-mb      :zs-1       )    , &
               xr*mb*mb,MPI_REAL,R_YDN,5, &
               C_CART,status,ierr)


        call MPI_SENDRECV(f(xs       :xe              , &
                            ys       :ys+(mb-1)       , &
                            zs-mb    :zs-1     )      , &
               xr*mb*mb,MPI_REAL,R_YDN,6          , &
                          f(xs       :xe              , &
                         yeb-(mb-1)  :yeb             , &
                         zs-mb       :zs-1     )      , &
               xr*mb*mb,MPI_REAL,R_YUP,6           , &
               C_CART,status,ierr)
      endif

      if ((.not.exist_zup).and.(zeb /= ze)) then

        call MPI_SENDRECV(f(xs       :xe              , &
                         ye-(mb-1)   :ye              , &
                         ze+1        :ze+mb    )      , &
               xr*mb*mb,MPI_REAL,R_YUP,7, &
                          f(xs       :xe              , &
                          ysb        :ysb+(mb-1)      , &
                          ze+1       :ze+mb      )    , &
               xr*mb*mb,MPI_REAL,R_YDN,7, &
               C_CART,status,ierr)


        call MPI_SENDRECV(f(xs       :xe              , &
                            ys       :ys+(mb-1)       , &
                            ze+1     :ze+mb    )      , &
               xr*mb*mb,MPI_REAL,R_YDN,8          , &
                          f(xs       :xe              , &
                         yeb-(mb-1)  :yeb             , &
                         ze+1        :ze+mb    )      , &
               xr*mb*mb,MPI_REAL,R_YUP,8           , &
               C_CART,status,ierr)
      endif

    end if

  call trace%end
   end subroutine comm_bdry_y

!************************************************************************

subroutine comm_bdry_xy (f)

  implicit none

  real(kind=4),dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f

  call trace%begin ('mpi.f90::comm_bdry_xy')

  if(.not.exist_zdn) then

    if (((exist_xup.and.exist_yup).or.(exist_xdn.and.exist_ydn)).and.(zsb /= zs)) then

      call MPI_SENDRECV(f(xe-(mb-1):xe               , &
           ye-(mb-1):ye               , &
           zs-mb    :zs-1     )       , &
           mb*mb*mb,MPI_REAL,R_XUPYUP,9, &
           f(xsb         :xsb+(mb-1)          , &
           ysb         :ysb+(mb-1)          , &
           zs-mb       :zs-1     )       , &
           mb*mb*mb,MPI_REAL,R_XDNYDN,9, &
           C_CART,status,ierr)

      call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
           ys       :ys+(mb-1)        , &
           zs-mb    :zs-1       )       , &
           mb*mb*mb,MPI_REAL,R_XDNYDN,10, &
           f(xeb-(mb-1)  :xeb                 , &
           yeb-(mb-1)  :yeb                 , &
           zs-mb       :zs-1     )       , &
           mb*mb*mb,MPI_REAL,R_XUPYUP,10, &
           C_CART,status,ierr)

    endif

    if (((exist_xup.and.exist_ydn).or.(exist_xdn.and.exist_yup)).and.(zsb /= zs)) then

      call MPI_SENDRECV(f(xe-(mb-1):xe               , &
           ys       :ys+(mb-1)        , &
           zs-mb     :zs-1       )       , &
           mb*mb*mb,MPI_REAL,R_XUPYDN,11, &
           f(xsb         :xsb+(mb-1)          , &
           yeb-(mb-1)  :yeb                 , &
           zs-mb       :zs-1     )       , &
           mb*mb*mb,MPI_REAL,R_XDNYUP,11, &
           C_CART,status,ierr)

      call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
           ye-(mb-1):ye               , &
           zs-mb    :zs-1       )       , &
           mb*mb*mb,MPI_REAL,R_XDNYUP,12, &
           f(xeb-(mb-1)  :xeb                 , &
           ysb         :ysb+(mb-1)          , &
           zs-mb       :zs-1     )       , &
           mb*mb*mb,MPI_REAL,R_XUPYDN,12, &
           C_CART,status,ierr)

    end if

  endif

  if(.not.exist_zup) then

    if (((exist_xup.and.exist_yup).or.(exist_xdn.and.exist_ydn)).and.(zeb /= ze)) then

      call MPI_SENDRECV(f(xe-(mb-1):xe               , &
           ye-(mb-1):ye               , &
           ze+1    :ze+mb     )       , &
           mb*mb*mb,MPI_REAL,R_XUPYUP,13, &
           f(xsb         :xsb+(mb-1)          , &
           ysb         :ysb+(mb-1)          , &
           ze+1       :ze+mb     )       , &
           mb*mb*mb,MPI_REAL,R_XDNYDN,13, &
           C_CART,status,ierr)

      call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
           ys       :ys+(mb-1)        , &
           ze+1       :ze+mb       )       , &
           mb*mb*mb,MPI_REAL,R_XDNYDN,14, &
           f(xeb-(mb-1)  :xeb                 , &
           yeb-(mb-1)  :yeb                 , &
           ze+1       :ze+mb     )       , &
           mb*mb*mb,MPI_REAL,R_XUPYUP,14, &
           C_CART,status,ierr)

    endif

    if (((exist_xup.and.exist_ydn).or.(exist_xdn.and.exist_yup)).and.(zeb /= ze)) then

      call MPI_SENDRECV(f(xe-(mb-1):xe               , &
           ys          :ys+(mb-1)        , &
           ze+1        :ze+mb       )       , &
           mb*mb*mb,MPI_REAL,R_XUPYDN,15, &
           f(xsb       :xsb+(mb-1)          , &
           yeb-(mb-1)  :yeb                 , &
           ze+1        :ze+mb     )       , &
           mb*mb*mb,MPI_REAL,R_XDNYUP,15, &
           C_CART,status,ierr)

      call MPI_SENDRECV(f(xs       :xs+(mb-1)        , &
           ye-(mb-1) :ye               , &
           ze+1      :ze+mb       )       , &
           mb*mb*mb,MPI_REAL,R_XDNYUP,16, &
           f(xeb-(mb-1)  :xeb                 , &
           ysb         :ysb+(mb-1)          , &
           ze+1       :ze+mb     )       , &
           mb*mb*mb,MPI_REAL,R_XUPYDN,16, &
           C_CART,status,ierr)

    end if

  endif




  call trace%end
end subroutine comm_bdry_xy

! ***********************************************************************

subroutine update_mg (f,sx,ex,sy,ey,sz,ez,sxgp,exgp,sygp,eygp,szgp,ezgp)

 implicit none

 integer :: sx,ex,sy,ey,sz,ez,sxgp,exgp,sygp,eygp,szgp,ezgp
 real,dimension(sx-sxgp:ex+exgp,sy-sygp:ey+eygp,sz-szgp:ez+ezgp) :: f

  call trace%begin ('mpi.f90::update_mg')

 call comm_x_mg(f,sx,ex,sy,ey,sz,ez,sxgp,exgp,sygp,eygp,szgp,ezgp)
 call comm_y_mg(f,sx,ex,sy,ey,sz,ez,sxgp,exgp,sygp,eygp,szgp,ezgp)
 call comm_z_mg(f,sx,ex,sy,ey,sz,ez,sxgp,exgp,sygp,eygp,szgp,ezgp)

 call comm_xy_mg(f,sx,ex,sy,ey,sz,ez,sxgp,exgp,sygp,eygp,szgp,ezgp)
 call comm_xz_mg(f,sx,ex,sy,ey,sz,ez,sxgp,exgp,sygp,eygp,szgp,ezgp)
 call comm_yz_mg(f,sx,ex,sy,ey,sz,ez,sxgp,exgp,sygp,eygp,szgp,ezgp)

 call comm_xyz_mg(f,sx,ex,sy,ey,sz,ez,sxgp,exgp,sygp,eygp,szgp,ezgp)

 call barrier(c_cart,ierr)

  call trace%end
end subroutine update_mg

! *************************************************

subroutine comm_x_mg (f,sx,ex,sy,ey,sz,ez,sxgp,exgp,sygp,eygp,szgp,ezgp)

implicit none

integer :: sx,ex,sy,ey,sz,ez,sxgp,exgp,sygp,eygp,szgp,ezgp
real, dimension(sx-sxgp:ex+exgp,sy-sygp:ey+eygp,sz-szgp:ez+ezgp) :: f

  call trace%begin ('mpi.f90::comm_x_mg')


if (exist_xup.or.exist_xdn) then

  call MPI_SENDRECV(f(ex-exgp+1        :ex                        , &
                      sy               :ey                        , &
                      sz               :ez                       ), &
       exgp*(ey-sy+1)*(ez-sz+1),OSC_REAL,R_XUP,31 , &
                    f(sx-sxgp          :sx-1                      , &
                      sy               :ey                        , &
                      sz               :ez                       ), &
       sxgp*(ey-sy+1)*(ez-sz+1),OSC_REAL,R_XDN,31 , &
       C_CART,status,ierr)

  call MPI_SENDRECV(f(sx               :sx+sxgp-1                 , &
                      sy               :ey                        , &
                      sz               :ez                       ), &
       exgp*(ey-sy+1)*(ez-sz+1),OSC_REAL,R_XDN,32 , &
                    f(ex+1             :ex+exgp                   , &
                      sy               :ey                        , &
                      sz               :ez                       ), &
       sxgp*(ey-sy+1)*(ez-sz+1),OSC_REAL,R_XUP,32 , &
       C_CART,status,ierr)

end if

  call trace%end
end subroutine comm_x_mg

! ***********************************************************************

subroutine comm_y_mg (f,sx,ex,sy,ey,sz,ez,sxgp,exgp,sygp,eygp,szgp,ezgp)

implicit none

integer :: sx,ex,sy,ey,sz,ez,sxgp,exgp,sygp,eygp,szgp,ezgp
real, dimension(sx-sxgp:ex+exgp,sy-sygp:ey+eygp,sz-szgp:ez+ezgp) :: f


if (exist_yup.or.exist_ydn) then

  call MPI_SENDRECV(f(sx               :ex                        , &
                      ey-eygp+1        :ey                        , &
                      sz               :ez                       ), &
       eygp*(ex-sx+1)*(ez-sz+1),OSC_REAL,R_YUP,33 , &
                    f(sx               :ex                        , &
                      sy-sygp          :sy-1                      , &
                      sz               :ez                       ), &
       sygp*(ex-sx+1)*(ez-sz+1),OSC_REAL,R_YDN,33 , &
       C_CART,status,ierr)

   call MPI_SENDRECV(f(sx               :ex                        , &
                      sy               :sy+sygp-1                 , &
                      sz               :ez                       ), &
       eygp*(ex-sx+1)*(ez-sz+1),OSC_REAL,R_YDN,34 , &
                    f(sx               :ex                        , &
                      ey+1             :ey+eygp                   , &
                      sz               :ez                       ), &
       sygp*(ex-sx+1)*(ez-sz+1),OSC_REAL,R_YUP,34 , &
       C_CART,status,ierr)

end if

  call trace%end
end subroutine comm_y_mg

! ***********************************************************************

subroutine comm_z_mg (f,sx,ex,sy,ey,sz,ez,sxgp,exgp,sygp,eygp,szgp,ezgp)

implicit none

integer :: sx,ex,sy,ey,sz,ez,sxgp,exgp,sygp,eygp,szgp,ezgp
real, dimension(sx-sxgp:ex+exgp,sy-sygp:ey+eygp,sz-szgp:ez+ezgp) :: f

  call trace%begin ('mpi.f90::comm_z_mg')

if (exist_zup.or.exist_zdn) then

  call MPI_SENDRECV(f(sx               :ex                        , &
                      sy               :ey                        , &
                      ez-ezgp+1        :ez                       ), &
       ezgp*(ex-sx+1)*(ey-sy+1),OSC_REAL,R_ZUP,35 , &
                    f(sx               :ex                        , &
                      sy               :ey                        , &
                      sz-szgp          :sz-1                     ), &
       szgp*(ex-sx+1)*(ey-sy+1),OSC_REAL,R_ZDN,35 , &
       C_CART,status,ierr)

  call MPI_SENDRECV(f(sx               :ex                        , &
                      sy               :ey                        , &
                      sz               :sz+szgp-1                ), &
       ezgp*(ex-sx+1)*(ey-sy+1),OSC_REAL,R_ZDN,36 , &
                    f(sx               :ex                        , &
                      sy               :ey                        , &
                      ez+1             :ez+ezgp                  ), &
       szgp*(ex-sx+1)*(ey-sy+1),OSC_REAL,R_ZUP,36 , &
       C_CART,status,ierr)

end if

  call trace%end
end subroutine comm_z_mg

!************************************************************************

subroutine comm_xy_mg (f,sx,ex,sy,ey,sz,ez,sxgp,exgp,sygp,eygp,szgp,ezgp)

implicit none

integer :: sx,ex,sy,ey,sz,ez,sxgp,exgp,sygp,eygp,szgp,ezgp
real, dimension(sx-sxgp:ex+exgp,sy-sygp:ey+eygp,sz-szgp:ez+ezgp) :: f

if ((exist_xup.and.exist_yup).or.(exist_xdn.and.exist_ydn)) then

  call MPI_SENDRECV(f(ex-exgp+1:ex                               , &
                      ey-eygp+1:ey                               , &
                      sz       :ez                       )       , &
       exgp*eygp*(ez-sz+1),OSC_REAL,R_XUPYUP,37, &
                    f(sx-sxgp          :sx-1                     , &
                      sy-sygp          :sy-1                     , &
                      sz               :ez               )       , &
       sxgp*sygp*(ez-sz+1),OSC_REAL,R_XDNYDN,37, &
       C_CART,status,ierr)

  call MPI_SENDRECV(f(sx       :sx+sxgp-1                        , &
                      sy       :sy+sygp-1                        , &
                      sz       :ez       )       , &
       sxgp*sygp*(ez-sz+1),OSC_REAL,R_XDNYDN,38, &
                    f(ex+1     :ex+exgp                          , &
                      ey+1     :ey+eygp                          , &
                      sz       :ez       )                       , &
       exgp*eygp*(ez-sz+1),OSC_REAL,R_XUPYUP,38, &
       C_CART,status,ierr)

endif

if ((exist_xup.and.exist_ydn).or.(exist_xdn.and.exist_yup)) then

  call MPI_SENDRECV(f(ex-exgp+1:ex                               , &
                      sy       :sy+sygp-1                        , &
                      sz       :ez                       )       , &
       exgp*sygp*(ez-sz+1),OSC_REAL,R_XUPYDN,39, &
                    f(sx-sxgp  :sx-1                             , &
                      ey+1     :ey+eygp                          , &
                      sz       :ez                       )       , &
       sxgp*eygp*(ez-sz+1),OSC_REAL,R_XDNYUP,39, &
       C_CART,status,ierr)

  call MPI_SENDRECV(f(sx       :sx+sxgp-1                        , &
                      ey-eygp+1:ey                               , &
                      sz       :ez                       )       , &
       sxgp*eygp*(ez-sz+1),OSC_REAL,R_XDNYUP,40, &
                    f(ex+1     :ex+exgp                          , &
                      sy-sygp  :sy-1                             , &
                      sz       :ez                       )       , &
       exgp*sygp*(ez-sz+1),OSC_REAL,R_XUPYDN,40, &
       C_CART,status,ierr)

end if

  call trace%end
end subroutine comm_xy_mg

! ***********************************************************************

subroutine comm_xz_mg (f,sx,ex,sy,ey,sz,ez,sxgp,exgp,sygp,eygp,szgp,ezgp)

implicit none

integer :: sx,ex,sy,ey,sz,ez,sxgp,exgp,sygp,eygp,szgp,ezgp
real, dimension(sx-sxgp:ex+exgp,sy-sygp:ey+eygp,sz-szgp:ez+ezgp) :: f

  call trace%begin ('mpi.f90::comm_xz_mg')

if ((exist_xup.and.exist_zup).or.(exist_xdn.and.exist_zdn)) then

  call MPI_SENDRECV(f(ex-exgp+1:ex                               , &
                      sy       :ey                               , &
                      ez-ezgp+1:ez                       )       , &
       exgp*(ey-sy+1)*ezgp,OSC_REAL,R_XUPZUP,41, &
                    f(sx-sxgp  :sx-1                             , &
                      sy       :ey                               , &
                      sz-szgp  :sz-1  )                          , &
       sxgp*(ey-sy+1)*szgp,OSC_REAL,R_XDNZDN,41, &
       C_CART,status,ierr)

  call MPI_SENDRECV(f(sx       :sx+sxgp-1                        , &
                      sy       :ey                               , &
                      sz       :sz+szgp-1                )       , &
       sxgp*(ey-sy+1)*szgp,OSC_REAL,R_XDNZDN,42, &
                    f(ex+1     :ex+exgp                          , &
                      sy       :ey                               , &
                      ez+1     :ez+ezgp                  )       , &
       exgp*(ey-sy+1)*ezgp,OSC_REAL,R_XUPZUP,42, &
       C_CART,status,ierr)

end if

if ((exist_xup.and.exist_zdn).or.(exist_xdn.and.exist_zup)) then

  call MPI_SENDRECV(f(ex-exgp+1:ex                               , &
                      sy       :ey                               , &
                      sz       :sz+szgp-1                )       , &
       exgp*(ey-sy+1)*szgp,OSC_REAL,R_XUPZDN,43, &
                    f(sx-sxgp  :sx-1                             , &
                      sy       :ey                               , &
                      ez+1     :ez+ezgp                  )       , &
       sxgp*(ey-sy+1)*ezgp,OSC_REAL,R_XDNZUP,43, &
       C_CART,status,ierr)

  call MPI_SENDRECV(f(sx       :sx+sxgp-1                        , &
                      sy       :ey                               , &
                      ez-ezgp+1:ez                       )       , &
       sxgp*(ey-sy+1)*ezgp,OSC_REAL,R_XDNZUP,44, &
                    f(ex+1     :ex+exgp                          , &
                      sy       :ey                               , &
                      sz-szgp  :sz-1                     )       , &
       exgp*(ey-sy+1)*szgp,OSC_REAL,R_XUPZDN,44, &
       C_CART,status,ierr)

end if

  call trace%end
end subroutine comm_xz_mg

! ***********************************************************************

subroutine comm_yz_mg (f,sx,ex,sy,ey,sz,ez,sxgp,exgp,sygp,eygp,szgp,ezgp)

implicit none

integer :: sx,ex,sy,ey,sz,ez,sxgp,exgp,sygp,eygp,szgp,ezgp
real, dimension(sx-sxgp:ex+exgp,sy-sygp:ey+eygp,sz-szgp:ez+ezgp) :: f

  call trace%begin ('mpi.f90::comm_yz_mg')

if ((exist_yup.and.exist_zup).or.(exist_ydn.and.exist_zdn)) then

  call MPI_SENDRECV(f(sx       :ex                               , &
                      ey-eygp+1:ey                               , &
                      ez-ezgp+1:ez                       )       , &
       (ex-sx+1)*eygp*ezgp,OSC_REAL,R_YUPZUP,45, &
                    f(sx       :ex                               , &
                      sy-sygp  :sy-1                             , &
                      sz-szgp  :sz-1  )                          , &
       (ex-sx+1)*sygp*szgp,OSC_REAL,R_YDNZDN,45, &
       C_CART,status,ierr)

  call MPI_SENDRECV(f(sx       :ex                               , &
                      sy       :sy+sygp-1                        , &
                      sz       :sz+szgp-1                )       , &
       (ex-sx+1)*sygp*szgp,OSC_REAL,R_YDNZDN,46, &
                    f(sx       :ex                               , &
                      ey+1     :ey+eygp                          , &
                      ez+1     :ez+ezgp                  )       , &
       (ex-sx+1)*eygp*ezgp,OSC_REAL,R_YUPZUP,46, &
       C_CART,status,ierr)

end if

if ((exist_yup.and.exist_zdn).or.(exist_ydn.and.exist_zup)) then

  call MPI_SENDRECV(f(sx       :ex                               , &
                      ey-eygp+1:ey                               , &
                      sz       :sz+szgp-1                )       , &
       (ex-sx+1)*eygp*szgp,OSC_REAL,R_YUPZDN,47, &
                    f(sx       :ex                               , &
                      sy-sygp  :sy-1                             , &
                      ez+1     :ez+ezgp                  )       , &
       (ex-sx+1)*sygp*ezgp,OSC_REAL,R_YDNZUP,47, &
       C_CART,status,ierr)

  call MPI_SENDRECV(f(sx       :ex                               , &
                      sy       :sy+sygp-1                        , &
                      ez-ezgp+1:ez                       )       , &
       (ex-sx+1)*sygp*ezgp,OSC_REAL,R_YDNZUP,48, &
                    f(sx       :ex                               , &
                      ey+1     :ey+eygp                          , &
                      sz-szgp  :sz-1                     )       , &
       (ex-sx+1)*eygp*szgp,OSC_REAL,R_YUPZDN,48, &
       C_CART,status,ierr)

end if

  call trace%end
end subroutine comm_yz_mg


! ***********************************************************************

subroutine comm_xyz_mg (f,sx,ex,sy,ey,sz,ez,sxgp,exgp,sygp,eygp,szgp,ezgp)

implicit none

integer :: sx,ex,sy,ey,sz,ez,sxgp,exgp,sygp,eygp,szgp,ezgp
real, dimension(sx-sxgp:ex+exgp,sy-sygp:ey+eygp,sz-szgp:ez+ezgp) :: f

  call trace%begin ('mpi.f90::comm_xyz_mg')

if ((exist_xup.and.exist_yup.and.exist_zup).or.(exist_xdn.and.exist_ydn.and.exist_zdn)) then

  call MPI_SENDRECV(f(ex-exgp+1:ex               , &
                      ey-eygp+1:ey               , &
                      ez-ezgp+1:ez       )       , &
       exgp*eygp*ezgp,OSC_REAL,R_XUPYUPZUP,49, &
                    f(sx-sxgp:sx-1          , &
                      sy-sygp:sy-1          , &
                      sz-szgp:sz-1  )       , &
       sxgp*sygp*szgp,OSC_REAL,R_XDNYDNZDN,49, &
       C_CART,status,ierr)

  call MPI_SENDRECV(f(sx:sx+sxgp-1        , &
                      sy:sy+sygp-1        , &
                      sz:sz+szgp-1)       , &
       sxgp*sygp*szgp,OSC_REAL,R_XDNYDNZDN,50, &
                    f(ex+1:ex+exgp                 , &
                      ey+1:ey+eygp                 , &
                      ez+1:ez+ezgp         )       , &
       exgp*eygp*ezgp,OSC_REAL,R_XUPYUPZUP,50 , &
       C_CART,status,ierr)

end if

if ((exist_xup.and.exist_yup.and.exist_zdn).or.(exist_xdn.and.exist_ydn.and.exist_zup)) then

  call MPI_SENDRECV(f(ex-exgp+1:ex               , &
                      ey-eygp+1:ey               , &
                      sz:sz+szgp-1)       , &
       exgp*eygp*szgp,OSC_REAL,R_XUPYUPZDN,51, &
                    f(sx-sxgp:sx-1           , &
                      sy-sygp:sy-1           , &
                      ez+1:ez+ezgp         )        , &
       sxgp*sygp*ezgp,OSC_REAL,R_XDNYDNZUP,51, &
       C_CART,status,ierr)

  call MPI_SENDRECV(f(sx:sx+sxgp-1        , &
                      sy:sy+sygp-1        , &
                      ez-ezgp+1:ez       )       , &
       sxgp*sygp*ezgp,OSC_REAL,R_XDNYDNZUP,52, &
                    f(ex+1:ex+exgp                 , &
                      ey+1:ey+eygp                 , &
                      sz-szgp:sz-1  )       , &
       exgp*eygp*szgp,OSC_REAL,R_XUPYUPZDN,52, &
       C_CART,status,ierr)

end if

if ((exist_xup.and.exist_ydn.and.exist_zup).or.(exist_xdn.and.exist_yup.and.exist_zdn)) then

  call MPI_SENDRECV(f(ex-exgp+1:ex               , &
                      sy:sy+sygp-1        , &
                      ez-ezgp+1:ez       )       , &
       exgp*sygp*ezgp,OSC_REAL,R_XUPYDNZUP,53, &
                    f(sx-sxgp:sx-1          , &
                      ey+1:ey+eygp                 , &
                      sz-szgp:sz-1  )       , &
       szgp*eygp*szgp,OSC_REAL,R_XDNYUPZDN,53, &
       C_CART,status,ierr)

  call MPI_SENDRECV(f(sx:sx+sxgp-1        , &
                      ey-eygp+1:ey               , &
                      sz:sz+szgp-1)       , &
       sxgp*eygp*szgp,OSC_REAL,R_XDNYUPZDN,54, &
                    f(ex+1:ex+exgp                 , &
                      sy-sygp:sy-1          , &
                      ez+1:ez+ezgp         )       , &
       exgp*sygp*ezgp,OSC_REAL,R_XUPYDNZUP,54, &
       C_CART,status,ierr)

end if

if ((exist_xup.and.exist_ydn.and.exist_zdn).or.(exist_xdn.and.exist_yup.and.exist_zup)) then

  call MPI_SENDRECV(f(ex-exgp+1:ex               , &
                      sy:sy+sygp-1        , &
                      sz:sz+szgp-1)       , &
       exgp*sygp*szgp,OSC_REAL,R_XUPYDNZDN,55, &
                    f(sx-sxgp:sx-1          , &
                      ey+1:ey+eygp                 , &
                      ez+1:ez+ezgp         )       , &
       szgp*eygp*ezgp,OSC_REAL,R_XDNYUPZUP,55, &
       C_CART,status,ierr)

  call MPI_SENDRECV(f(sx:sx+sxgp-1        , &
                      ey-eygp+1:ey               , &
                      ez-ezgp+1:ez       )       , &
       sxgp*eygp*ezgp,OSC_REAL,R_XDNYUPZUP,56, &
                    f(ex+1:ex+exgp                 , &
                      sy-sygp:sy-1          , &
                      sz-szgp:sz-1  )       , &
       exgp*sygp*szgp,OSC_REAL,R_XUPYDNZDN,56, &
       C_CART,status,ierr)

end if

  call trace%end
end subroutine comm_xyz_mg

! ***********************************************************************

 subroutine print_check (f,loc)

  implicit none

  real,dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f
  character(len=3) :: loc
  integer :: i,nr,nr2,j1,j2,j
  character(len=10),dimension(8) :: pos
  real,dimension(8) :: sum_res
  logical :: go_through

  call trace%begin ('mpi.f90::print_check')

  nr=0
  if (index(loc,'X').ne.0) nr=nr+1
  if (index(loc,'Y').ne.0) nr=nr+1
  if (index(loc,'Z').ne.0) nr=nr+1

  if (nr.eq.1) then
    nr=2
    if (scan(loc,"X").ne.0) then
      pos(1)="XUP"
      pos(2)="XDN"
    elseif (scan(loc,"Y").ne.0) then
      pos(1)="YUP"
      pos(2)="YDN"
    elseif (scan(loc,"Z").ne.0) then
      pos(1)="ZUP"
      pos(2)="ZDN"
    end if
  else if (nr.eq.2) then
    nr=4
    if (index(loc,"X").ne.0) then
      if (index(loc,"Y").ne.0) then
        pos(1)="XUPYUP"
        pos(2)="XUPYDN"
        pos(3)="XDNYUP"
        pos(4)="XDNYDN"
      endif
      if (index(loc,"Z").ne.0) then
        pos(1)="XUPZUP"
        pos(2)="XUPZDN"
        pos(3)="XDNZUP"
        pos(4)="XDNZDN"
      end if
    else
      pos(1)="YUPZUP"
      pos(2)="YUPZDN"
      pos(3)="YDNZUP"
      pos(4)="YDNZDN"
    end if
  else if (nr.eq.3) then
    nr=8
    pos(1)="XUPYUPZUP"
    pos(2)="XUPYUPZDN"
    pos(3)="XUPYDNZUP"
    pos(4)="XUPYDNZDN"
    pos(5)="XDNYUPZUP"
    pos(6)="XDNYUPZDN"
    pos(7)="XDNYDNZUP"
    pos(8)="XDNYDNZDN"
  end if


  nr2=min(nr,4)
  go_through=.true.
  j1=1
  j2=nr2

  do while (go_through)

    if (myid.eq.R_MASTER) then
      print*,repeat('-',4*16+10)
      write(*,'("          ",4A16)') pos(j1:j2)
      print*
    end if

    do i=1,gridsize
      if (myid.eq.R_MASTER) then
        do j=j1,j2
          sum_res(j)=esa_sum(f,pos(j),i)
        end do
        write(*,'("Master",I3," ",4E16.8)') i,sum_res(j1:j2)
      end if
      call MPI_BARRIER(C_CART,ierr)
      if (myid.eq.i) then
        do j=j1,j2
          sum_res(j)=esa_sum(f,pos(j),i)
        end do
        write(*,'("Slave ",I3," ",4E16.8)') i,sum_res(j1:j2)
      endif
      call MPI_BARRIER(C_CART,ierr)
      if (myid.eq.R_MASTER) then
        print*
      endif
    end do

    if (nr.gt.j2) then
      go_through=.true.
      j1=j2+1
      j2=nr
      nr2=(j2-j1+1)
    else
      go_through=.false.
    end if
  end do

  call trace%end
 end subroutine print_check

! ***********************************************************************

 function esa_sum(f,loc,i)

! ** (Extract-Sub-Array)

  implicit none

  real,dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: f
  character(len=9),intent(IN) :: loc
  real :: esa_sum
  integer :: xx0,xx1,yy0,yy1,zz0,zz1,xstart,ystart,zstart
  integer,intent(in) :: i
  logical :: invx,invy,invz

  call trace%begin ('mpi.f90::print_check')

  xstart=index(loc,"X")
  ystart=index(loc,"Y")
  zstart=index(loc,"Z")

  invx=.false.
  invy=.false.
  invz=.false.

  if (xstart.ne.0) invx=.true.
  if (ystart.ne.0) invy=.true.
  if (zstart.ne.0) invz=.true.


  if (invx) then
    if (loc(xstart:xstart+2).eq.'XUP') then
      if (myid.eq.R_MASTER) then
        xx0=mod(x1(i)-(mb-1)-1+mx,mx)+1
        xx1=mod(x1(i)       -1+mx,mx)+1
      else
        xx0=xeb-(mb-1)
        xx1=xeb
      end if
    else if (loc(xstart:xstart+2).eq.'XDN') then
      if (myid.eq.R_MASTER) then
        xx0=mod(x0(i)       -1+mx,mx)+1
        xx1=mod(x0(i)+(mb-1)-1+mx,mx)+1
      else
        xx0=xsb
        xx1=xsb+(mb-1)
      end if
    end if
  else
    if (myid.eq.R_MASTER) then
      xx0=xmin(i)
      xx1=xmax(i)
    else
      xx0=xs
      xx1=xe
    end if
  end if

  if (invy) then
    if (loc(ystart:ystart+2).eq.'YUP') then
      if (myid.eq.R_MASTER) then
        yy0=mod(y1(i)-(mb-1)-1+my,my)+1
        yy1=mod(y1(i)       -1+my,my)+1
      else
        yy0=yeb-(mb-1)
        yy1=yeb
      end if
    else if (loc(ystart:ystart+2).eq.'YDN') then
      if (myid.eq.R_MASTER) then
        yy0=mod(y0(i)       -1+my,my)+1
        yy1=mod(y0(i)+(mb-1)-1+my,my)+1
      else
        yy0=ysb
        yy1=ysb+(mb-1)
      end if
    end if
  else
    if (myid.eq.R_MASTER) then
      yy0=ymin(i)
      yy1=ymax(i)
    else
      yy0=ys
      yy1=ye
    end if
  end if

  if (invz) then
    if (loc(zstart:zstart+2).eq.'ZUP') then
      if (myid.eq.R_MASTER) then
        zz0=mod(z1(i)-(mb-1)-1+mz,mz)+1
        zz1=mod(z1(i)       -1+mz,mz)+1
      else
        zz0=zeb-(mb-1)
        zz1=zeb
      end if
    else if (loc(zstart:zstart+2).eq.'ZDN') then
      if (myid.eq.R_MASTER) then
        zz0=mod(z0(i)       -1+mz,mz)+1
        zz1=mod(z0(i)+(mb-1)-1+mz,mz)+1
      else
        zz0=zsb
        zz1=zsb+(mb-1)
      end if
    end if
  else
    if (myid.eq.R_MASTER) then
      zz0=zmin(i)
      zz1=zmax(i)
    else
      zz0=zs
      zz1=ze
    end if
  end if

  esa_sum=sum(f(xx0:xx1,yy0:yy1,zz0:zz1))

 end function esa_sum

! ***********************************************************************

 subroutine adjust_array (arr,direction,arrnew,tag)

  implicit none

  real,dimension(:,:,:) :: arr
  character(len=6) :: direction
  logical :: sendup,senddn,recvup,recvdn
  real,dimension(xsb:xeb,ysb:yeb,zsbtmp:zebtmp) :: arrnew
  integer(kind=4) :: tag

  call trace%begin ('mpi.f90::adjust_array')

  sendup = .false.
  senddn = .false.
  recvup = .false.
  recvdn = .false.

  select case (direction)
  case ('sendup')
    sendup = .true.
  case ('senddn')
    senddn = .true.
  case ('recvup')
    recvup = .true.
  case ('recvdn')
    recvdn = .true.
  end select


! ** Sending recieving arr

  if (sendup) then
!    print*,'in sendup',myid
!    print*,'in sendup,checking bounds',lbound(arr,3),ubound(arr,3),1,zebtmp-(zsbtmp-1)
    arrnew(:,:,:)=arr(:,:,zsbtmp-(zsbtmp-1):zebtmp-(zsbtmp-1))
    call MPI_IBSEND(arr(:,:,zetmp+1-(zsbtmp-1)),xrb*yrb,OSC_REAL, &
                   R_ZUP,1013,C_CART,b_sending_request(b_sending_counter),ierr)
    b_sending_counter=b_sending_counter+1
  end if

  if (recvdn) then
!    print*,'in recvdn',myid
    arrnew(:,:,zsbtmp+1:zebtmp)=arr(:,:,:)
    call MPI_RECV(arrnew(:,:,zstmp),xrb*yrb,OSC_REAL, &
                  R_ZDN,1013,C_CART,status,ierr)
  end if

  if (senddn) then
    arrnew(:,:,:)=arr(:,:,zsbtmp+1-(zsbtmp-1):zebtmp-(zsbtmp-1))
    call MPI_IBSEND(arr(:,:,zstmp-1-(zsbtmp-1)),xrb*yrb,OSC_REAL,R_ZDN,1014,C_CART,b_sending_request(b_sending_counter),ierr)
    b_sending_counter=b_sending_counter+1
  end if

  if (recvup) then
    arrnew(:,:,zsbtmp:zebtmp-1)=arr(:,:,:)
    call MPI_RECV(arrnew(:,:,zetmp),xrb*yrb,OSC_REAL,R_ZUP,1014,C_CART,status,ierr)
  end if

  call trace%end
 end subroutine adjust_array

! ***********************************************************************

! ***********************************************************************

 subroutine adjust_array2_4 (arr,sendup,senddn,recvup,recvdn,arrnew)

  implicit none

  real,dimension(xsb:xeb,ysb:yeb,zsbold:zebold) :: arr
  logical :: sendup,senddn,recvup,recvdn
  real,dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: arrnew
  integer(kind=4) :: request_up,request_dn,lzs,lze
  real,dimension(xsb:xeb,ysb:yeb) :: buff_recvup,buff_recvdn

  call trace%begin ('mpi.f90::adjust_array2_4')

! ** Sending recieving arr

  if (recvup) then
    call MPI_IRECV(buff_recvup,xrb*yrb,OSC_REAL,R_ZUP,1013,C_CART,request_up,ierr)
  end if

  if (recvdn) then
    call MPI_IRECV(buff_recvdn,xrb*yrb,OSC_REAL,R_ZDN,1014,C_CART,request_dn,ierr)
  end if

  if (senddn) then
    call MPI_SEND(arr(:,:,zsold),xrb*yrb,OSC_REAL,R_ZDN,1013,C_CART,ierr)
  end if

  if (sendup) then
    call MPI_SEND(arr(:,:,zeold),xrb*yrb,OSC_REAL,r_ZUP,1014,C_CART,ierr)
  end if

  lzs=zs
  lze=ze
  if ((.not.sendup).and.(.not.recvup)) lze=zeb
  if ((.not.senddn).and.(.not.recvdn)) lzs=zsb
  arrnew(:,:,lzs:lze)=arr(:,:,lzs:lze)

  if (recvup) then
    call MPI_WAIT(request_up,status,ierr)
    arrnew(:,:,ze)=buff_recvup
  end if

  if (recvdn) then
    call MPI_WAIT(request_dn,status,ierr)
    arrnew(:,:,zs)=buff_recvdn
  end if


!call update(arrnew)
!***  if (sendup) then
!***!    print*,'in sendup',myid
!***    print*,'in sendup,checking bounds',lbound(arr,3),ubound(arr,3),1,zebtmp-(zsbtmp-1)
!***    arrnew(:,:,:)=arr(:,:,zsbtmp-(zsbtmp-1):zebtmp-(zsbtmp-1))
!***    call MPI_IBSEND(arr(:,:,zetmp+1-(zsbtmp-1)),xrb*yrb,OSC_REAL, &
!***                   R_ZUP,1013,C_CART,b_sending_request(b_sending_counter),ierr)
!***    b_sending_counter=b_sending_counter+1
!***  end if
!***
!***  if (recvdn) then
!***!    print*,'in recvdn',myid
!***    arrnew(:,:,zsbtmp+1:zebtmp)=arr(:,:,:)
!***    call MPI_RECV(arrnew(:,:,zstmp),xrb*yrb,OSC_REAL, &
!***                  R_ZDN,1013,C_CART,status,ierr)
!***  end if
!***
!***  if (senddn) then
!***    arrnew(:,:,:)=arr(:,:,zsbtmp+1-(zsbtmp-1):zebtmp-(zsbtmp-1))
!***    call MPI_IBSEND(arr(:,:,zstmp-1-(zsbtmp-1)),xrb*yrb,OSC_REAL,R_ZDN,1014,C_CART,b_sending_request(b_sending_counter),ierr)
!***    b_sending_counter=b_sending_counter+1
!***  end if
!***
!***  if (recvup) then
!***    arrnew(:,:,zsbtmp:zebtmp-1)=arr(:,:,:)
!***    call MPI_RECV(arrnew(:,:,zetmp),xrb*yrb,OSC_REAL,R_ZUP,1014,C_CART,status,ierr)
!***  end if

  call trace%end
 end subroutine adjust_array2_4

! ***********************************************************************

 subroutine adjust_array2_8 (arr,sendup,senddn,recvup,recvdn,arrnew)

  implicit none

  real(kind=8),dimension(xsb:xeb,ysb:yeb,zsbold:zebold) :: arr
  logical :: sendup,senddn,recvup,recvdn
  real(kind=8),dimension(xsb:xeb,ysb:yeb,zsb:zeb) :: arrnew
  integer(kind=4) :: request_up,request_dn
  real(kind=8),dimension(xsb:xeb,ysb:yeb) :: buff_recvup,buff_recvdn

  call trace%begin ('mpi.f90::adjust_array2_8')

! ** Sending recieving arr

  if (recvup) then
    call MPI_IRECV(buff_recvup,xrb*yrb,OSC_DOUBLE_PRECISION,R_ZUP,1013,C_CART,request_up,ierr)
  end if

  if (recvdn) then
    call MPI_IRECV(buff_recvdn,xrb*yrb,OSC_DOUBLE_PRECISION,R_ZDN,1014,C_CART,request_dn,ierr)
  end if

  if (senddn) then
    call MPI_SEND(arr(:,:,zsold),xrb*yrb,OSC_DOUBLE_PRECISION,R_ZDN,1013,C_CART,ierr)
  end if

  if (sendup) then
    call MPI_SEND(arr(:,:,zeold),xrb*yrb,OSC_DOUBLE_PRECISION,r_ZUP,1014,C_CART,ierr)
  end if


  arrnew(:,:,zs:ze)=arr(:,:,zs:ze)

  if (recvup) then
    call MPI_WAIT(request_up,status,ierr)
    arrnew(:,:,ze)=buff_recvup
  end if

  if (recvdn) then
    call MPI_WAIT(request_dn,status,ierr)
    arrnew(:,:,zs)=buff_recvdn
  end if


  call trace%end
 end subroutine adjust_array2_8
 ! *********************************************************************



 subroutine allocate_arrays_parallel (reallocation,sendup,senddn,recvup,recvdn)

  implicit none

  logical,optional,intent(in) :: reallocation,sendup,senddn,recvup,recvdn

  call trace%begin ('mpi.f90::allocate_arrays_parallel')

  call trace%end
 end subroutine allocate_arrays_parallel

! ***********************************************************************

 subroutine finalize_parallel

  implicit none

  integer :: ierr

  call MPI_FINALIZE(ierr)

 end subroutine finalize_parallel

end MODULE parallel_mod


!  1D
!  none needed
!
!  2D
!  if (exist_xupyup) then
!    standard
!  elseif (exist_xup).and.(.not.exist_yup) then
!    communicate xupyup to xup
!  elseif (.not.exist_xup).and.(exist_yup) then
!    communicate xupyup to yup
!  endif
!
!  3D
!  if (exist_xupyupzup) then
!    standard
!  elsif (exist_xupyup).and.(.not.exist_zup) then
!    communicat xupyupzup to xupyup
!  elseif (exist_xupzup).and.(.not.exist_yup) then
!    communicate xupyupzup to xupzup
!  elseif (exist_yupzup).and.(.not.exist_xup) then
!    communicate xupyupzup to yupzup
!  endif
!
