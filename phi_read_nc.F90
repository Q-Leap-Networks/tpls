!Standalone tool to read in the phi_channel in netcdf format
!Converts to standard phi output

program phi_read_nc
  use netcdf
  implicit none
  
  ! This is the name of the data file we will read.
  character (len = 128) :: filename, cmd
  integer :: ncid
  integer :: argc
  integer :: iargc
  integer, parameter :: NDIMS = 3, NRECS = 1

  !number of entries for each dimension
  integer :: n_x, n_y, n_z

  ! The start and count arrays will tell the netCDF library where to
  ! read our data.
  integer :: start(NDIMS), count(NDIMS)

  !Create variables to hold data
  double precision, allocatable, dimension(:,:,:)  :: phi2
  double precision, allocatable, dimension(:,:,:)  :: phi
 ! real :: lats(NLATS), lons(NLONS)
 ! integer :: lon_varid, lat_varid

  !Variable names
  character (len = *), parameter :: x="x"
  character (len = *), parameter :: y="y"
  character (len = *), parameter :: z="z"

  integer :: varid
  integer :: dimid(NDIMS)

  !Get filename
  ! take filename from command-line argument if there is any
  call getarg(0, cmd)
  argc = IARGC()
  if (argc .NE. 1) then
     print*,'Usage: ',trim(cmd),' filename'
  endif
  call getarg(1, filename)

  !Open a file to print to
  open(unit=7, file='phi_postnetcdf.dat')

  ! Open the file. 
  call check( nf90_open(filename, nf90_nowrite, ncid) )
  write(*,*)"open"
  ! Get the varids of the 3dims.
  call check( nf90_inq_dimid(ncid, x, dimid(1) ) )
  write(*,*)"dimidz"
  call check( nf90_inq_dimid(ncid, y, dimid(2) ) )
  write(*,*)"dimidx"
  call check( nf90_inq_dimid(ncid, z, dimid(3) ) )
  write(*,*)"dimidy"

  !Get size of variables
  call check( nf90_inquire_dimension(ncid,dimid(1), len=n_x) )
  write(*,*)"inq_dim 1"
  call check( nf90_inquire_dimension(ncid,dimid(2), len=n_y) )
  write(*,*)"inq_dim 2"
  call check( nf90_inquire_dimension(ncid,dimid(3), len=n_z) )
  write(*,*)"inq_dim 3"

  !allocate array to receive data
  allocate(phi2(0:n_x+1,0:n_y+1,0:n_z+1))
  allocate(phi(n_x,n_y,n_z))
  phi2 = 0

  ! Obtain variable id
  call check( nf90_inq_varid(ncid, "phi_out", varid) )
  write(*,*)"inq varid"

  ! Set start and count
  start = (/ 1, 1, 1 /)
  count = (/ n_x, n_y, n_z /)
  write(*,*) "dim x,y,z",n_x, n_y, n_z

  !Get data
  call check( nf90_get_var(ncid, varid, phi) )
  write(*,*)"get var"
  phi2(1:n_x,1:n_y,1:n_z)= phi(1:n_x,1:n_y,1:n_z)
  ! Close the file. 
  call check( nf90_close(ncid) )
  write(*,*)"close "

!!$  write(*,*) "phi2 check: ",phi2(1,1,1), phi2(n_x,n_y,n_z)
!!$  write(*,*) "phi(256,144,153):",phi2(256,144,153)
!!$  write(*,*) "phi(73,73,1)",phi2(73,73,1)
!!$  write(*,*) "phi(128,144,153):",phi2(128,144,153)
!!$  write(*,*) "phi(73,73,1):",phi2(73,73,1)
!!$  write(*,*) "phi(256,72,153):",phi2(256,72,153)
!!$  write(*,*) "phi(1,1,1):",phi2(1,1,1)
!!$  write(*,*) "phi(128,72,153):",phi2(128,72,153)
!!$  write(*,*) "phi(1,1,1):",phi2(1,1,1)

 write(*,*) "phi(128,144,151): ",phi(128,144,151)
 write(*,*) "phi(73,73,1): ",phi(73,73,1)
 write(*,*) "phi( 256,144,151): ",phi( 256,144,151)
 write(*,*) "phi( 73,73,1): ",phi( 73,73,1)
 write(*,*) "phi( 256, 72,151): ",phi( 256, 72,151)
 write(*,*) "phi( 1,1,1): ",phi( 1,1,1)
 write(*,*) "phi( 128,72,151): ",phi( 128,72,151)
 write(*,*) "phi(1,1,1): ",phi(1,1,1)

  ! If we got this far, everything worked as expected. Yipee! 
  print *,"*** SUCCESS reading ", filename

contains
  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check  
end program phi_read_nc

