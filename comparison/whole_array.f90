
!      
!   Whole array finite Difference Phase Field Code of Cahn-Hilliard Eq.
!
!        comparison code
!
!
!   Author  :
!               Shahid Maqbool
!
!   Modified   :
!                    29 March 2022
!
!   To compile and run :
!                            check ReadMe file
!
!------------------------------------------------------------------------------


program fd_ch_test
  implicit none


  ! ===========================================================================
  !                                parameters
  ! ===========================================================================


  integer ( kind = 4 ), parameter :: Nx = 256 , Ny = 256
  integer ( kind = 4 ) :: nsteps = 30000 , istep
  real ( kind = 8 )    :: dt = 0.01 , start , finish
  real ( kind = 8 )    :: c0 = 0.4 , mobility = 1.0 , grad_coef = 0.5
  real ( kind = 8 )    :: noise = 0.02 , A = 1.0
  real ( kind = 8 ), dimension ( Nx, Ny ) :: r, con, dfdcon, lap_con, dummy_con


  call cpu_time ( start )


  ! ===========================================================================
  !                            initial microstucture
  ! ===========================================================================


  call random_number ( r )

  con = c0 + noise*( 0.5 - r )


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

          
  end do time_loop

  call cpu_time ( finish )


  ! ===========================================================================
  !                     print computed time on the screen 
  ! ===========================================================================


  print*,'------------------------------------------------------'
  print '("  Whole array computed time = ", f10.3," seconds." )', finish - start


contains


  ! ===========================================================================
  !                               Sub-program
  ! ===========================================================================


  function Laplacian ( con_ )
    implicit none

    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( in )  :: con_ 
    real ( kind = 8 ), dimension ( Nx, Ny ) :: Laplacian
    integer ( kind = 4 ) :: i , j, jp, jm, ip, im, dx = 1, dy = 1


    do concurrent ( i = 1:Nx, j = 1:Ny )

          jp = j + 1
          jm = j - 1

          ip = i + 1
          im = i - 1

          if ( im == 0 ) im = Nx
          if ( ip == ( Nx + 1 ) ) ip = 1
          if ( jm == 0 ) jm = Ny
          if ( jp == ( Ny + 1 ) ) jp = 1

          
          Laplacian(i,j) = ( con_(ip,j) + con_(im,j) + con_(i,jm) + &
               & con_(i,jp) - 4.0*con_(i,j) ) / ( dx*dy )


    end do

  end function Laplacian


  
end program fd_ch_test
