!> MPI subroutines.
!!
!! @author Prashant Valluri, Lennon O Naraigh, Iain Bethune,
!! Toni Collis, David Scott, Peter Spelt.
!! @copyright (c) 2013-2015, Prashant Valluri, Lennon O Naraigh, 
!! Iain Bethune, Toni Collis, David Scott, Peter Spelt.
!! This program is distributed under the BSD License See LICENSE.txt
!! for details.

module tpls_mpi

  implicit none

  integer :: master_id
  parameter(master_id=0)

contains

  !> Is the given process ID the master?
  logical function is_master(my_id)
    integer,intent(in) :: my_id

    is_master = (master_id == my_id)
  end function is_master

  subroutine mpi_decomp_2d(sx,ex,sy,ey,n_local_x,n_local_y,maxl,maxm,coords,dims,Ndim)
    implicit none

    integer,intent(inout) :: sx,ex,sy,ey,n_local_x,n_local_y
    integer,intent(in) :: maxl,maxm,Ndim
    integer,intent(in) :: coords(Ndim),dims(Ndim)

    integer :: n_procs_x,n_procs_y,rank

    n_procs_x=dims(1)        
    rank=coords(1)
    n_local_x = (maxl-1)/n_procs_x
    sx = rank*n_local_x + 1
    ex = sx + n_local_x - 1

    n_procs_y=dims(2)        
    rank=coords(2)
    n_local_y = (maxm-1)/n_procs_y
    sy = rank*n_local_y + 1
    ey = sy + n_local_y - 1

    return
  end subroutine mpi_decomp_2d


  subroutine get_mpi_neighbours(neighbours,comm3d)     
    use mpi
    implicit none

    integer :: neighbours(6), N,E,S,W,Fr,Bk,comm3d,ierr
    parameter (N=1,E=2,S=3,W=4,Fr=5,Bk=6)

    neighbours(1:6)  =  MPI_PROC_NULL
    CALL  MPI_CART_SHIFT(comm3d,  0, 1, neighbours(W),  neighbours(E),  ierr)
    CALL  MPI_CART_SHIFT(comm3d,  1, 1, neighbours(S),  neighbours(N),  ierr)
    CALL  MPI_CART_SHIFT(comm3d,  2, 1, neighbours(Fr), neighbours(Bk), ierr)

    return
  end subroutine get_mpi_neighbours


  !> For passing 1st-order halos of non-augmented arrays WITH CORNERS
  !! (slightly modified code) 
  subroutine get_stride_p(stride_p_yz,stride_p_xz,sx,ex,sy,ey,sz,ez)
    use mpi
    implicit none

    integer :: stride_p_xz,stride_p_yz,type_y
    integer :: ex,ey,ez,sx,sy,sz,ierr,size_real

    ! For passing xz planes in the y direction (North-South).  This one is changed from before to allow
    ! the passing of vertex points.
    CALL  MPI_TYPE_VECTOR(    &
         ez-sz+1,              & ! nombre de blocs
         ex-sx+3,              & ! longueur d'un bloc (change this one for bigger blocks)
         (ex-sx+3)*(ey-sy+3),  & ! pas entre le debut de deux blocs consecutifs
         mpi_double_precision, &  
         stride_p_xz, ierr  )
    CALL  MPI_TYPE_COMMIT(  stride_p_xz,  ierr  )

    ! For passing yz planes in the x direction (East-West).  This one is ALSO changed from before.
    CALL MPI_TYPE_SIZE (mpi_double_precision,  size_real, ierr )
    CALL MPI_TYPE_HVECTOR ( &
         ey-sy+3, &
         1, &
         (ex-sx+3)*size_real,   &
         mpi_double_precision, & 
         type_y, ierr  )
    CALL MPI_TYPE_COMMIT (  type_y, ierr )

    CALL  MPI_TYPE_HVECTOR ( &
         ez-sz+1, &            
         1, &                  
         (ey-sy+3)*(ex-sx+3)*size_real, & 
         type_y,    &           
         stride_p_yz, ierr  )
    CALL MPI_TYPE_COMMIT (  stride_p_yz,  ierr  )

    return
  end subroutine get_stride_p


  !> For passing 1st-order halos of non-augmented arrays WITH CORNERS
  !! (slightly modified code) 
  Subroutine exchange2d(uu,stride_p_xz,stride_p_yz,neighbours,ex,ey,ez,sx,sy,sz,comm_topology)
    use mpi
    implicit none

    integer, intent(in) :: ex,ey,ez,sx,sy,sz
    integer, intent(in) :: stride_p_xz,stride_p_yz
    integer, intent(in) :: neighbours(6),comm_topology
    double precision, dimension(sx-1:ex+1,sy-1:ey+1,sz:ez), intent(inout) :: uu

    integer :: requests(4)
    integer :: statuses(mpi_status_size,4)
    integer :: tag1=100,tag2=200,ierr
    integer :: N,E,S,W
    parameter (N=1,E=2,S=3,W=4)

    ! Send neighbours "N" and receive neighbours "S" using the stride type stride_p_xz.  
    ! This part is changed to allow for sending extra strips of data.
    ! Send neighbours "S" and receive neighbours "N" using the stride type stride_p_xz

    CALL  MPI_IRECV( &
         uu(sx-1,sy-1,sz),  1,  stride_p_xz,  neighbours(S),  tag1, &
         comm_topology,  requests(1),  ierr )
    CALL  MPI_IRECV( &
         uu(sx-1,ey+1,sz),  1,  stride_p_xz,  neighbours(N),  tag2, &
         comm_topology,  requests(2),  ierr )

    CALL  MPI_ISEND( &
         uu(sx-1,ey,  sz),  1,  stride_p_xz,  neighbours(N),  tag1, &
         comm_topology,  requests(3),  ierr )
    CALL  MPI_ISEND( &
         uu(sx-1,sy,  sz),  1,  stride_p_xz,  neighbours(S),  tag2, &
         comm_topology,  requests(4),  ierr )

    CALL MPI_WAITALL(4, requests, statuses, ierr)

    ! This part is also changed to allow for sending extra strips of data.  
    ! Corners need to be exchanged twice becuase in the first (N-S) swap, 
    ! wrong information is received into corners.  
    ! The second (E-W) swap fixes this problem.

    ! Send neighbours "W" and receive neighbours "E" using the stride type stride_p_yz
    ! Send neighbours "E" and receive neighbours "W" using the stride type stride_p_yz

    CALL  MPI_IRECV( &
         uu(ex+1,sy-1,sz),  1,  stride_p_yz,  neighbours(E),  tag1, &
         comm_topology,  requests(1),  ierr )
    CALL  MPI_IRECV( &
         uu(sx-1,sy-1,sz),  1,  stride_p_yz,  neighbours(W),  tag2, &
         comm_topology,  requests(2),  ierr )

    CALL  MPI_ISEND( &
         uu(sx,  sy-1,sz),  1,  stride_p_yz,  neighbours(W),  tag1, &
         comm_topology,  requests(3),  ierr )
    CALL  MPI_ISEND( &
         uu(ex,  sy-1,sz),  1,  stride_p_yz,  neighbours(E),  tag2, &
         comm_topology,  requests(4),  ierr )

    CALL MPI_WAITALL(4, requests, statuses, ierr)

    return
  end subroutine exchange2d


  !> For passing 1st-order halos of AUGMENTED arrays WITH CORNERS
  !! (modified code) 
  subroutine get_stride_p_aug1(stride_p_aug1_yz,stride_p_aug1_xz,sx,ex,sy,ey,sz,ez)
    use mpi
    implicit none

    integer :: stride_p_aug1_xz,stride_p_aug1_yz,type_y
    integer :: ex,ey,ez,sx,sy,sz,ierr,size_real

    ! for passing xz planes in the y direction (North-South)
    CALL  MPI_TYPE_VECTOR(    &
         ez-sz+1,              & ! nombre de blocs
         ex-sx+3,              & ! longueur d'un bloc
         (ex-sx+5)*(ey-sy+5),  & ! pas entre le debut de deux blocs consecutifs
         mpi_double_precision, &  
         stride_p_aug1_xz, ierr  )
    CALL  MPI_TYPE_COMMIT(  stride_p_aug1_xz,  ierr  )        

    ! for passing yz planes in the x direction (East-West)
    CALL MPI_TYPE_SIZE (mpi_double_precision,  size_real, ierr )
    CALL MPI_TYPE_HVECTOR ( &
         ey-sy+3, &
         1, &
         (ex-sx+5)*size_real,   &
         mpi_double_precision, & 
         type_y, ierr  )
    CALL MPI_TYPE_COMMIT (  type_y, ierr )

    CALL  MPI_TYPE_HVECTOR ( &
         ez-sz+1, &            
         1, &                  
         (ey-sy+5)*(ex-sx+5)*size_real, & 
         type_y,    &           
         stride_p_aug1_yz, ierr  )
    CALL MPI_TYPE_COMMIT (  stride_p_aug1_yz,  ierr  )

    return
  end subroutine get_stride_p_aug1


  !> For passing 1st-order halos of augmented arrays WITH CORNERS
  !! (modified code) 
  Subroutine exchange2d_aug1(uu,stride_p_aug1_xz,stride_p_aug1_yz,neighbours,ex,ey,ez,sx,sy,sz,comm_topology)
    use mpi
    implicit none

    integer, intent(in) :: ex,ey,ez,sx,sy,sz
    integer, intent(in) :: stride_p_aug1_xz,stride_p_aug1_yz
    integer, intent(in) :: neighbours(6),comm_topology
    double precision, dimension(sx-2:ex+2,sy-2:ey+2,sz:ez), intent(inout) :: uu

    integer :: requests(4)
    integer :: statuses(mpi_status_size,4)
    integer :: tag1=100,tag2=200,ierr
    integer :: N,E,S,W
    parameter (N=1,E=2,S=3,W=4)

    ! Send neighbours "N" and receive neighbours "S" using the stride type stride_p_xz.  
    ! This one is changed to allow for sending extra strips of data.
    ! Send neighbours "S" and receive neighbours "N" using the stride type stride_p_aug1_xz

    CALL  MPI_IRECV( &
         uu(sx-1,sy-1,sz),  1,  stride_p_aug1_xz,  neighbours(S),  tag1, &
         comm_topology,  requests(1),  ierr )
    CALL  MPI_IRECV( &
         uu(sx-1,ey+1,sz),  1,  stride_p_aug1_xz,  neighbours(N),  tag2, &
         comm_topology,  requests(2),  ierr )

    CALL  MPI_ISEND( &
         uu(sx-1,ey,  sz),  1,  stride_p_aug1_xz,  neighbours(N),  tag1, &
         comm_topology,  requests(3),  ierr )
    CALL  MPI_ISEND( &
         uu(sx-1,sy,  sz),  1,  stride_p_aug1_xz,  neighbours(S),  tag2, &
         comm_topology,  requests(4),  ierr )

    CALL MPI_WAITALL(4, requests, statuses, ierr)

    ! Send neighbours "W" and receive neighbours "E" using the stride type stride_p_yz
    ! Send neighbours "E" and receive neighbours "W" using the stride type stride_p_aug1_yz

    CALL  MPI_IRECV( &
         uu(ex+1,sy-1,sz),  1,  stride_p_aug1_yz,  neighbours(E),  tag1, &
         comm_topology,  requests(1),  ierr )
    CALL  MPI_IRECV( &
         uu(sx-1,sy-1,sz),  1,  stride_p_aug1_yz,  neighbours(W),  tag2, &
         comm_topology,  requests(2),  ierr )

    CALL  MPI_ISEND( &
         uu(sx,  sy-1,sz),  1,  stride_p_aug1_yz,  neighbours(W),  tag1, &
         comm_topology,  requests(3),  ierr )
    CALL  MPI_ISEND( &
         uu(ex,  sy-1,sz),  1,  stride_p_aug1_yz,  neighbours(E),  tag2, &
         comm_topology,  requests(4),  ierr )

    CALL MPI_WAITALL(4, requests, statuses, ierr)

    return
  end subroutine exchange2d_aug1


  !> For passing 2nd-order halos of augmented arrays WITHOUT CORNERS
  !! (modified code) 
  subroutine get_stride_p_aug2(stride_p_aug2_yz,stride_p_aug2_xz,sx,ex,sy,ey,sz,ez)
    use mpi
    implicit none

    integer :: stride_p_aug2_xz,stride_p_aug2_yz,type_y
    integer :: ex,ey,ez,sx,sy,sz,ierr,size_real

    ! for passing xz planes in the y direction (North-South)
    CALL  MPI_TYPE_VECTOR(    &
         ez-sz+1,              & ! nombre de blocs
         ex-sx+3,              & ! longueur d'un bloc
         (ex-sx+5)*(ey-sy+5),  & ! pas entre le debut de deux blocs consecutifs
         mpi_double_precision, &  
         stride_p_aug2_xz, ierr  )
    CALL  MPI_TYPE_COMMIT(  stride_p_aug2_xz,  ierr  )        

    ! for passing yz planes in the x direction (East-West)
    CALL MPI_TYPE_SIZE (mpi_double_precision,  size_real, ierr )
    CALL MPI_TYPE_HVECTOR ( &
         ey-sy+3, &
         1, &
         (ex-sx+5)*size_real,   &
         mpi_double_precision, & 
         type_y, ierr  )
    CALL MPI_TYPE_COMMIT (  type_y, ierr )

    CALL  MPI_TYPE_HVECTOR ( &
         ez-sz+1, &            
         1, &                  
         (ey-sy+5)*(ex-sx+5)*size_real, & 
         type_y,    &           
         stride_p_aug2_yz, ierr  )
    CALL MPI_TYPE_COMMIT (  stride_p_aug2_yz,  ierr  )

    return
  end subroutine get_stride_p_aug2


  !> For passing 2nd-order halos of augmented arrays WITHOUT CORNERS
  !! (modified code) 
  Subroutine exchange2d_aug2(uu,stride_p_aug2_xz,stride_p_aug2_yz,neighbours,ex,ey,ez,sx,sy,sz,comm_topology)
    use mpi
    implicit none

    integer, intent(in) :: ex,ey,ez,sx,sy,sz
    integer, intent(in) :: stride_p_aug2_xz,stride_p_aug2_yz
    integer, intent(in) :: neighbours(6),comm_topology
    double precision, dimension(sx-2:ex+2,sy-2:ey+2,sz:ez), intent(inout) :: uu

    integer :: requests(4)
    integer :: statuses(mpi_status_size,4)
    integer :: tag1=100,tag2=200,ierr
    integer :: N,E,S,W
    parameter (N=1,E=2,S=3,W=4)

    ! Send neighbours "N" and receive neighbours "S" using the stride type stride_p_aug2_xz
    ! Send neighbours "S" and receive neighbours "N" using the stride type stride_p_aug2_xz

    CALL  MPI_IRECV( &
         uu(sx-1, sy-2,sz),  1,  stride_p_aug2_xz,  neighbours(S),  tag1, &
         comm_topology,  requests(1),  ierr )
    CALL  MPI_IRECV( &
         uu(sx-1, ey+2,sz),  1,  stride_p_aug2_xz,  neighbours(N),  tag2, &
         comm_topology,  requests(2),  ierr )

    CALL  MPI_ISEND( &
         uu(sx-1, ey-1,  sz),  1,  stride_p_aug2_xz,  neighbours(N),  tag1, &
         comm_topology,  requests(3),  ierr )
    CALL  MPI_ISEND( &
         uu(sx-1, sy+1,  sz),  1,  stride_p_aug2_xz,  neighbours(S),  tag2, &
         comm_topology,  requests(4),  ierr )

    CALL MPI_WAITALL(4, requests, statuses, ierr)

    ! Send neighbours "W" and receive neighbours "E" using the stride type stride_p_aug2_yz
    ! Send neighbours "E" and receive neighbours "W" using the stride type stride_p_aug2_yz

    CALL  MPI_IRECV( &
         uu(ex+2, sy-1,sz),  1,  stride_p_aug2_yz,  neighbours(E),  tag1, &
         comm_topology,  requests(1),  ierr )
    CALL  MPI_IRECV( &
         uu(sx-2, sy-1,sz),  1,  stride_p_aug2_yz,  neighbours(W),  tag2, &
         comm_topology,  requests(2),  ierr )

    CALL  MPI_ISEND( &
         uu(sx+1, sy-1,sz),  1,  stride_p_aug2_yz,  neighbours(W),  tag1, &
         comm_topology,  requests(3),  ierr )
    CALL  MPI_ISEND( &
         uu(ex-1, sy-1,sz),  1,  stride_p_aug2_yz,  neighbours(E),  tag2, &
         comm_topology,  requests(4),  ierr )

    CALL MPI_WAITALL(4, requests, statuses, ierr)

    return
  end subroutine exchange2d_aug2


  ! Doubly-Augmented arrays


  !> For passing 1st-order halos of DOUBLY-augmented arrays WITH CORNERS
  subroutine get_stride_p_augaug1(stride_p_augaug1_yz,stride_p_augaug1_xz,sx,ex,sy,ey,sz,ez)
    use mpi
    implicit none

    integer :: stride_p_augaug1_xz,stride_p_augaug1_yz,type_y
    integer :: ex,ey,ez,sx,sy,sz,ierr,size_real

    ! for passing xz planes in the y direction (North-South)
    CALL  MPI_TYPE_VECTOR(    &
         ez-sz+1,              & ! nombre de blocs
         ex-sx+3,              & ! longueur d'un bloc
         (ex-sx+7)*(ey-sy+7),  & ! pas entre le debut de deux blocs consecutifs
         mpi_double_precision, &  
         stride_p_augaug1_xz, ierr  )
    CALL  MPI_TYPE_COMMIT(  stride_p_augaug1_xz,  ierr  )        

    ! for passing yz planes in the x direction (East-West)
    CALL MPI_TYPE_SIZE (mpi_double_precision,  size_real, ierr )
    CALL MPI_TYPE_HVECTOR ( &
         ey-sy+3, &
         1, &
         (ex-sx+7)*size_real,   &
         mpi_double_precision, & 
         type_y, ierr  )
    CALL MPI_TYPE_COMMIT (  type_y, ierr )

    CALL  MPI_TYPE_HVECTOR ( &
         ez-sz+1, &            
         1, &                  
         (ey-sy+7)*(ex-sx+7)*size_real, & 
         type_y,    &           
         stride_p_augaug1_yz, ierr  )
    CALL MPI_TYPE_COMMIT (  stride_p_augaug1_yz,  ierr  )

    return
  end subroutine get_stride_p_augaug1


  !> For passing 1st-order halos of DOUBLY-augmented arrays WITH CORNERS
  Subroutine exchange2d_augaug1(uu,stride_p_augaug1_xz,stride_p_augaug1_yz,neighbours,ex,ey,ez,sx,sy,sz,comm_topology)
    use mpi
    implicit none

    integer, intent(in) :: ex,ey,ez,sx,sy,sz
    integer, intent(in) :: stride_p_augaug1_xz,stride_p_augaug1_yz
    integer, intent(in) :: neighbours(6),comm_topology
    double precision, dimension(sx-3:ex+3,sy-3:ey+3,sz:ez), intent(inout) :: uu

    integer :: requests(4)
    integer :: statuses(mpi_status_size,4)
    integer :: tag1=100,tag2=200,ierr
    integer :: N,E,S,W
    parameter (N=1,E=2,S=3,W=4)

    ! Send neighbours "N" and receive neighbours "S" using the stride type stride_p_xz.  
    ! This one is changed to allow for sending extra strips of data.
    ! Send neighbours "S" and receive neighbours "N" using the stride type stride_p_augaug1_xz

    CALL  MPI_IRECV( &
         uu(sx-1,sy-1,sz),  1,  stride_p_augaug1_xz,  neighbours(S),  tag1, &
         comm_topology,  requests(1),  ierr )
    CALL  MPI_IRECV( &
         uu(sx-1,ey+1,sz),  1,  stride_p_augaug1_xz,  neighbours(N),  tag2, &
         comm_topology,  requests(2),  ierr )

    CALL  MPI_ISEND( &
         uu(sx-1,ey,  sz),  1,  stride_p_augaug1_xz,  neighbours(N),  tag1, &
         comm_topology,  requests(3),  ierr )
    CALL  MPI_ISEND( &
         uu(sx-1,sy,  sz),  1,  stride_p_augaug1_xz,  neighbours(S),  tag2, &
         comm_topology,  requests(4),  ierr )

    CALL MPI_WAITALL(4, requests, statuses, ierr)

    ! Send neighbours "W" and receive neighbours "E" using the stride type stride_p_yz
    CALL  MPI_IRECV( &
         uu(ex+1,sy-1,sz),  1,  stride_p_augaug1_yz,  neighbours(E),  tag1, &
         comm_topology,  requests(1),  ierr )
    CALL  MPI_IRECV( &
         uu(sx-1,sy-1,sz),  1,  stride_p_augaug1_yz,  neighbours(W),  tag2, &
         comm_topology,  requests(2),  ierr )

    ! Send neighbours "E" and receive neighbours "W" using the stride type stride_p_augaug1_yz
    CALL  MPI_ISEND( &
         uu(sx,  sy-1,sz),  1,  stride_p_augaug1_yz,  neighbours(W),  tag1, &
         comm_topology,  requests(3),  ierr )
    CALL  MPI_ISEND( &
         uu(ex,  sy-1,sz),  1,  stride_p_augaug1_yz,  neighbours(E),  tag2, &
         comm_topology,  requests(4),  ierr )

    CALL MPI_WAITALL(4, requests, statuses, ierr)

    return
  end subroutine exchange2d_augaug1


  !> For passing 2nd-order halos of DOUBLY-augmented arrays WITH CORNERS
  subroutine get_stride_p_augaug2(stride_p_augaug2_yz,stride_p_augaug2_xz,sx,ex,sy,ey,sz,ez)
    use mpi
    implicit none

    integer :: stride_p_augaug2_xz,stride_p_augaug2_yz,type_y
    integer :: ex,ey,ez,sx,sy,sz,ierr,size_real

    ! for passing xz planes in the y direction (North-South)
    CALL  MPI_TYPE_VECTOR(    &
         ez-sz+1,              & ! nombre de blocs
         ex-sx+5,              & ! longueur d'un bloc
         (ex-sx+7)*(ey-sy+7),  & ! pas entre le debut de deux blocs consecutifs
         mpi_double_precision, &  
         stride_p_augaug2_xz, ierr  )
    CALL  MPI_TYPE_COMMIT(  stride_p_augaug2_xz,  ierr  )        

    ! for passing yz planes in the x direction (East-West)
    CALL MPI_TYPE_SIZE (mpi_double_precision,  size_real, ierr )
    CALL MPI_TYPE_HVECTOR ( &
         ey-sy+5, &
         1, &
         (ex-sx+7)*size_real,   &
         mpi_double_precision, & 
         type_y, ierr  )
    CALL MPI_TYPE_COMMIT (  type_y, ierr )

    CALL  MPI_TYPE_HVECTOR ( &
         ez-sz+1, &            
         1, &                  
         (ey-sy+7)*(ex-sx+7)*size_real, & 
         type_y,    &           
         stride_p_augaug2_yz, ierr  )
    CALL MPI_TYPE_COMMIT (  stride_p_augaug2_yz,  ierr  )

    return
  end subroutine get_stride_p_augaug2


  !> For passing 2nd-order halos of DOUBLY-augmented arrays WITH CORNERS
  Subroutine exchange2d_augaug2(uu,stride_p_augaug2_xz,stride_p_augaug2_yz,neighbours,ex,ey,ez,sx,sy,sz,comm_topology)
    use mpi
    implicit none

    integer, intent(in) :: ex,ey,ez,sx,sy,sz
    integer, intent(in) :: stride_p_augaug2_xz,stride_p_augaug2_yz
    integer, intent(in) :: neighbours(6),comm_topology
    double precision, dimension(sx-3:ex+3,sy-3:ey+3,sz:ez), intent(inout) :: uu

    integer :: requests(4)
    integer :: statuses(mpi_status_size,4)
    integer :: tag1=100,tag2=200,ierr
    integer :: N,E,S,W
    parameter (N=1,E=2,S=3,W=4)

    ! Send neighbours "N" and receive neighbours "S" using the stride type stride_p_augaug2_xz
    ! Send neighbours "S" and receive neighbours "N" using the stride type stride_p_augaug2_xz

    CALL  MPI_IRECV( &
         uu(sx-2, sy-2,sz),  1,  stride_p_augaug2_xz,  neighbours(S),  tag1, &
         comm_topology,  requests(1),  ierr )
    CALL  MPI_IRECV( &
         uu(sx-2, ey+2,sz),  1,  stride_p_augaug2_xz,  neighbours(N),  tag2, &
         comm_topology,  requests(2),  ierr )

    CALL  MPI_ISEND( &
         uu(sx-2, ey-1,  sz),  1,  stride_p_augaug2_xz,  neighbours(N),  tag1, &
         comm_topology,  requests(3),  ierr )
    CALL  MPI_ISEND( &
         uu(sx-2, sy+1,  sz),  1,  stride_p_augaug2_xz,  neighbours(S),  tag2, &
         comm_topology,  requests(4),  ierr )

    CALL MPI_WAITALL(4, requests, statuses, ierr)

    ! Send neighbours "W" and receive neighbours "E" using the stride type stride_p_augaug2_yz
    CALL  MPI_IRECV( &
         uu(ex+2, sy-2,sz),  1,  stride_p_augaug2_yz,  neighbours(E),  tag1, &
         comm_topology,  requests(1),  ierr )
    CALL  MPI_IRECV( &
         uu(sx-2, sy-2,sz),  1,  stride_p_augaug2_yz,  neighbours(W),  tag2, &
         comm_topology,  requests(2),  ierr )

    ! Send neighbours "E" and receive neighbours "W" using the stride type stride_p_augaug2_yz
    CALL  MPI_ISEND( &
         uu(sx+1, sy-2,sz),  1,  stride_p_augaug2_yz,  neighbours(W),  tag1, &
         comm_topology,  requests(3),  ierr )
    CALL  MPI_ISEND( &
         uu(ex-1, sy-2,sz),  1,  stride_p_augaug2_yz,  neighbours(E),  tag2, &
         comm_topology,  requests(4),  ierr )

    CALL MPI_WAITALL(4, requests, statuses, ierr)

    return
  end subroutine exchange2d_augaug2


  !> For passing 3rd-order halos of DOUBLY-augmented arrays WITHOUT CORNERS
  subroutine get_stride_p_augaug3(stride_p_augaug2_yz,stride_p_augaug2_xz,sx,ex,sy,ey,sz,ez)
    use mpi
    implicit none

    integer :: stride_p_augaug2_xz,stride_p_augaug2_yz,type_y
    integer :: ex,ey,ez,sx,sy,sz,ierr,size_real

    ! for passing xz planes in the y direction (North-South)
    CALL  MPI_TYPE_VECTOR(    &
         ez-sz+1,              & ! nombre de blocs
         ex-sx+5,              & ! longueur d'un bloc
         (ex-sx+7)*(ey-sy+7),  & ! pas entre le debut de deux blocs consecutifs
         mpi_double_precision, &  
         stride_p_augaug2_xz, ierr  )
    CALL  MPI_TYPE_COMMIT(  stride_p_augaug2_xz,  ierr  )        

    ! for passing yz planes in the x direction (East-West)
    CALL MPI_TYPE_SIZE (mpi_double_precision,  size_real, ierr )
    CALL MPI_TYPE_HVECTOR ( &
         ey-sy+5, &
         1, &
         (ex-sx+7)*size_real,   &
         mpi_double_precision, & 
         type_y, ierr  )
    CALL MPI_TYPE_COMMIT (  type_y, ierr )

    CALL  MPI_TYPE_HVECTOR ( &
         ez-sz+1, &            
         1, &                  
         (ey-sy+7)*(ex-sx+7)*size_real, & 
         type_y,    &           
         stride_p_augaug2_yz, ierr  )
    CALL MPI_TYPE_COMMIT (  stride_p_augaug2_yz,  ierr  )

    return
  end subroutine get_stride_p_augaug3


  !> For passing 2nd-order halos of DOUBLY-augmented arrays WITHOUT CORNERS
  Subroutine exchange2d_augaug3(uu,stride_p_augaug3_xz,stride_p_augaug3_yz,neighbours,ex,ey,ez,sx,sy,sz,comm_topology)
    use mpi
    implicit none

    integer, intent(in) :: ex,ey,ez,sx,sy,sz
    integer, intent(in) :: stride_p_augaug3_xz,stride_p_augaug3_yz
    integer, intent(in) :: neighbours(6),comm_topology
    double precision, dimension(sx-3:ex+3,sy-3:ey+3,sz:ez), intent(inout) :: uu

    integer :: requests(4)
    integer :: statuses(mpi_status_size,4)
    integer :: tag1=100,tag2=200,ierr
    integer :: N,E,S,W
    parameter (N=1,E=2,S=3,W=4)

    ! Send neighbours "N" and receive neighbours "S" using the stride type stride_p_augaug2_xz
    ! Send neighbours "S" and receive neighbours "N" using the stride type stride_p_augaug2_xz

    CALL  MPI_IRECV( &
         uu(sx-2, sy-3,sz),  1,  stride_p_augaug3_xz,  neighbours(S),  tag1, &
         comm_topology,  requests(1),  ierr )
    CALL  MPI_IRECV( &
         uu(sx-2, ey+3,sz),  1,  stride_p_augaug3_xz,  neighbours(N),  tag2, &
         comm_topology,  requests(2),  ierr )

    CALL  MPI_ISEND( &
         uu(sx-2, ey-2,  sz),  1,  stride_p_augaug3_xz,  neighbours(N),  tag1, &
         comm_topology,  requests(3),  ierr )
    CALL  MPI_ISEND( &
         uu(sx-2, sy+2,  sz),  1,  stride_p_augaug3_xz,  neighbours(S),  tag2, &
         comm_topology,  requests(4),  ierr )

    CALL MPI_WAITALL(4, requests, statuses, ierr)

    ! Send neighbours "W" and receive neighbours "E" using the stride type stride_p_augaug2_yz
    ! Send neighbours "E" and receive neighbours "W" using the stride type stride_p_augaug2_yz

    CALL  MPI_IRECV( &
         uu(ex+3, sy-2,sz),  1,  stride_p_augaug3_yz,  neighbours(E),  tag1, &
         comm_topology,  requests(1),  ierr )
    CALL  MPI_IRECV( &
         uu(sx-3, sy-2,sz),  1,  stride_p_augaug3_yz,  neighbours(W),  tag2, &
         comm_topology,  requests(2),  ierr )

    CALL  MPI_ISEND( &
         uu(sx+2, sy-2,sz),  1,  stride_p_augaug3_yz,  neighbours(W),  tag1, &
         comm_topology,  requests(3),  ierr )
    CALL  MPI_ISEND( &
         uu(ex-2, sy-2,sz),  1,  stride_p_augaug3_yz,  neighbours(E),  tag2, &
         comm_topology,  requests(4),  ierr )

    CALL MPI_WAITALL(4, requests, statuses, ierr)

    return
  end subroutine exchange2d_augaug3

end module tpls_mpi
