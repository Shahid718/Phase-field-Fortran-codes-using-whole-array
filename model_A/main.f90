!      
!   Whole array Finite Difference Phase Field Code of Allen-Cahn Equation.
!
!
!   Author :
!               Shahid Maqbool
!
!   Modified   :
!                    29 March 2023
!
!   To compile and run :
!                            check ReadMe file
!
!------------------------------------------------------------------------------


program fd_ac_test
  use Dislin
  implicit none


  ! ===========================================================================
  !                                parameters
  ! ===========================================================================


  integer ( kind = 4 ), parameter :: Nx = 128 , Ny = 128
  integer ( kind = 4 ) :: nsteps = 1500 , nprint = 10 , istep
  real ( kind = 8 )    :: dt = 0.01 , start , finish
  real ( kind = 8 )    :: phi_0 = 0.5 , mobility = 1.0 , grad_coef = 1.0
  real ( kind = 8 )    :: noise = 0.02 , A = 1.0
  real ( kind = 8 ) , dimension ( Nx, Ny ) :: r , phi, dfdphi  


  call cpu_time ( start )


  ! ===========================================================================
  !                            initial microstucture
  ! ===========================================================================


  call random_number ( r )

  phi = phi_0 + noise*( 0.5 - r )


  ! ===========================================================================
  !                       Setting dislin routines for animation 
  ! ===========================================================================

  
  call Metafl ( 'cons' )
  call Disini ( )

  
  ! ===========================================================================
  !                         evolution of microstructure 
  ! ===========================================================================


  time_loop: do istep = 1, nsteps


     dfdphi = A*( 2.0*phi*( 1.0 - phi )**2*( 1.0 - 2*phi ) )

     phi = phi - dt*mobility*( dfdphi - grad_coef*Laplacian( phi ) )


     ! adjust order parameter in range

     where ( phi >= 0.99999 )  phi = 0.99999
     where ( phi < 0.00001 )   phi = 0.00001
  

     call Dislin_color_animation ( )

     
  end do time_loop

  
  call Disfin ( )
  call cpu_time ( finish )


  ! ===========================================================================
  !              print computed time on the screen and dislin plot
  ! ===========================================================================


  print *,'---------------------------------'
  print '("  Time       = ", f10.3," seconds." )', finish - start


contains


  ! ===========================================================================
  !                               Sub-programs
  ! ===========================================================================


  function Laplacian ( phi_) 
    implicit none

    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( in ) :: phi_ 
    real ( kind = 8 ), dimension ( Nx, Ny ) :: Laplacian 
    integer ( kind = 4 ) :: i , j, jp, jm, ip, im , dx = 2, dy = 2


    do concurrent ( i = 1 : Nx, j= 1 : Ny )

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


  end function Laplacian


  !----------------------------------------------------------------------------


  subroutine Dislin_color_animation ()

    call autres ( Nx, Ny )
    if ( mod( istep, nprint ) .eq. 0 ) then 
       call erase ( )  
       call graf3 ( 0.d0, 128.d0, 0.d0, 32.d0, 0.d0, 128.d0,&
            & 0.d0, 32.d0, 0.05d0, 1.0d0, 0.05d0, 0.1d0 )
       call crvmat ( phi, Nx, Ny, 1, 1 )   
       call endgrf
       call sendbf ( )
    end if


  end subroutine Dislin_color_animation


end program fd_ac_test
