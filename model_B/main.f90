!      
!   Whole array Finite Difference Phase Field Code of Cahn-Hilliard Eq.
!
!
!   Author  :
!               Shahid Maqbool
!
!   Modified   :
!                    29 March 2023
!
!   To compile and run :
!                            check ReadMe file
!
!------------------------------------------------------------------------------


program fd_ch_test
  use Dislin
  implicit none


  ! ===========================================================================
  !                                parameters
  ! ===========================================================================


  integer ( kind = 4 ), parameter :: Nx = 128 , Ny = 128
  integer ( kind = 4 ) :: nsteps = 10000 , nprint = 50 , istep
  real ( kind = 8 )    :: dt = 0.01 , start , finish
  real ( kind = 8 )    :: c0 = 0.4 , mobility = 1.0 , grad_coef = 0.5
  real ( kind = 8 )    :: noise = 0.02 , A = 1.0
  real ( kind = 8 ), dimension ( Nx, Ny ) :: r, con, dfdcon, dummy_con


  call cpu_time ( start )


  ! ===========================================================================
  !                            initial microstucture
  ! ===========================================================================


  call random_number ( r )

  con = c0 + noise*( 0.5 - r )


  ! ===========================================================================
  !                  Setting initial dislin routins for animation 
  ! ===========================================================================


  call Metafl ( 'cons' )
  call Disini ( )


  ! ===========================================================================
  !                         evolution of microstructure 
  ! ===========================================================================


  time_loop: do istep = 1, nsteps


     dfdcon = A*( 2.0*con*( 1.0 - con )**2 &
          -2.0*con**2*( 1.0 - con ) )

     dummy_con = dfdcon - grad_coef*Laplacian( con )

     con = con + dt*mobility*Laplacian( dummy_con )


     ! adjust concentration in range

     where ( con >= 0.99999 )  con = 0.99999
     where ( con < 0.00001 )   con = 0.00001

     
     call Dislin_color_animation ( )

     
  end do time_loop


  call Disfin ( )
  call cpu_time ( finish )


  ! ===========================================================================
  !                     print computed time on the screen 
  ! ===========================================================================


  print*,'---------------------------------'
  print '("  Time       = ", f10.3," seconds." )', finish - start


contains


  ! ===========================================================================
  !                               Sub-programs
  ! ===========================================================================


  function Laplacian ( con_ )
    implicit none

    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( in )  :: con_ 
    real ( kind = 8 ), dimension ( Nx, Ny ) :: Laplacian
    integer ( kind = 4 ) :: i , j, jp, jm, ip, im, dx = 1, dy = 1


    do i = 1, Nx
       do j = 1, Ny

       jp = j + 1
       jm = j - 1

       ip = i + 1
       im = i - 1

       if ( im == 0 ) im = Nx
       if ( ip == ( Nx + 1 ) ) ip = 1
       if ( jm == 0 ) jm = Ny
       if ( jp == ( Ny + 1 ) ) jp = 1

       Laplacian(i,j) = ( phi_(ip,j) + phi_(im,j) + phi_(i,jm) + &
            phi_(i,jp) - 4.0*phi_(i,j))  / ( dx*dy )              

       end do
    end do


  end function Laplacian


  !----------------------------------------------------------------------------


  subroutine Dislin_color_animation ( )

    call autres ( Nx, Ny )
    if ( mod( istep, nprint ) .eq. 0 ) then 
       call erase ( )  
       call graf3 ( 0.d0, 128.d0, 0.d0, 32.d0, 0.d0, 128.d0,&
            & 0.d0, 32.d0, 0.0d0, 1.0d0, 0.0d0, 0.1d0 )
       call crvmat ( con, Nx, Ny, 1, 1 )   
       call endgrf
       call sendbf ( )
    end if

  end subroutine Dislin_color_animation


end program fd_ch_test
