
!      
!   Indexed array finite Difference Phase Field Code of Cahn-Hilliard Equation.
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
!-------------------------------------------------------------------------------


program fd_ch_test
  implicit none


  ! ===========================================================================
  !                                parameters
  ! ===========================================================================


  integer ( kind = 4 ), parameter :: Nx = 256 , Ny = 256, dx = 1, dy = 1
  integer ( kind = 4 ) :: nsteps = 30000 , istep
  real ( kind = 8 )    :: dt = 0.01 , start , finish
  real ( kind = 8 )    :: c0 = 0.4 , mobility = 1.0 , grad_coef = 0.5
  real ( kind = 8 )    :: noise = 0.02 , A = 1.0
  real ( kind = 8 ), dimension ( Nx, Ny ) :: r, con, dfdcon,  dummy_con
  real ( kind = 8 ), dimension ( Nx, Ny ) :: lap_con, lap_dummy
  integer ( kind = 4 ) :: i , j, jp, jm, ip, im


  call cpu_time ( start )


  ! ===========================================================================
  !                            initial microstucture
  ! ===========================================================================


  do i = 1 , Nx
     do j = 1, Ny

        call random_number ( r (i,j) )

        con(i,j) = c0 + noise*( 0.5 - r(i,j) )

     end do
  end do



  ! ===========================================================================
  !                         evolution of microstructure 
  ! ===========================================================================


  time_loop: do istep = 1, nsteps

     do concurrent ( i = 1:Nx, j = 1:Ny )

        dfdcon(i,j) = A*( 2.0*con(i,j)*( 1.0 - con(i,j) )**2 &
             - 2.0*con(i,j)**2*( 1.0 - con(i,j) ) )

        jp = j + 1
        jm = j - 1

        ip = i + 1
        im = i - 1

        if ( im == 0 ) im = Nx
        if ( ip == ( Nx + 1) ) ip = 1
        if ( jm == 0 ) jm = Ny
        if ( jp == ( Ny + 1) ) jp = 1

        lap_con(i,j)   = ( con(ip,j) + con(im,j) + con(i,jm) + con(i,jp) - &
             4.0*con(i,j) ) /( dx*dy )

        dummy_con(i,j) = dfdcon(i,j) - grad_coef*lap_con(i,j)

        lap_dummy(i,j) = ( dummy_con(ip,j) + dummy_con(im,j) + dummy_con(i,jm) &
             + dummy_con(i,jp) - 4.0*dummy_con(i,j) ) / ( dx*dy )

        con(i,j) =  con(i,j) + dt*mobility*lap_dummy(i,j)


     end do

     ! adjust concentration in range

     if ( con(i,j) >= 0.99999 )  con(i,j) = 0.99999
     if ( con(i,j) < 0.00001 )   con(i,j) = 0.00001


  end do time_loop

  call cpu_time ( finish )


  ! ===========================================================================
  !                     print computed time on the screen 
  ! ===========================================================================


  print*,'------------------------------------------------------'
  print '("  Indexed array compute time = ", f10.3," seconds." )', finish - start



end program fd_ch_test
