!      
!   Whole array Finite Difference PhaseField Code of CahnHilliard-AllenCahn Eqs.
!                                                                             
!
!   Author  :
!               Shahid Maqbool
!
!   Modified   :
!                    17 May 2023
!
!   To compile and run :
!                            check ReadMe file
!
!------------------------------------------------------------------------------


program model_ch_ac_test
  use Dislin
  implicit none


  ! ===========================================================================
  !                                parameters
  ! ===========================================================================

  integer ( kind = 4 ), parameter :: Nx = 128, Ny = 128 
  integer ( kind = 4 ) :: nsteps = 10000, nprint = 1000, istep, i, j
  real ( kind = 8 ) :: dt = 0.03, start, finish
  real ( kind = 8 ) :: A = 1.0, B = 1.0, D = 1.0, radius = 10.0
  real ( kind = 8 ) :: mobility_con = 0.5,  mobility_phi  = 0.5
  real ( kind = 8 ) :: grad_coef_con = 1.5, grad_coef_phi = 1.5 
  real ( kind = 8 ), dimension ( Nx, Ny ) :: con, phi, dfdcon,dfdphi, dummy_con


  call cpu_time ( start )


  ! ===========================================================================
  !                            initial microstucture
  ! ===========================================================================


  con = 0.02
  phi = 0.0

  do i = 1, Nx
     do j = 1, Ny


        if ( (i - Nx/2)*(i - Nx/2) + (j - Ny/2)*(j - Ny/2) < radius**2 ) then
           con(i,j) = 1.0
           phi(i,j) = 1.0
        endif


     end do
  end do


  ! ===========================================================================
  !                       Setting dislin routines for animation 
  ! ===========================================================================


  call Metafl ( 'cons' )
  call Disini ( ) 



  ! ===========================================================================
  !                         evolution of microstructure 
  ! ===========================================================================


  time_loop: do istep = 1, nsteps


     ! derivatives of free energy


     dfdcon = 2*A*con*(1 - ( phi**3*( 10 - 15*phi + 6*phi**2 ) ) ) & 
          - 2*B*(1 - con)*( phi**3*( 10 - 15*phi + 6*phi**2 ) )

     dfdphi = -A*con*con*( 3*phi**2*( 10 - 15*phi + 6*phi**2 ) + phi**3* &
          ( 12*phi - 15 )) + 2*B*(1 - con)*(1 - con)*( 3*phi**2*( 10 - 15*phi + &
          6*phi**2 ) + phi**3*( 12*phi - 15 )) + 2*D*phi*(1 - phi)*(1 - 2*phi )


     ! dummy array for concentration

     dummy_con = dfdcon - grad_coef_con*Laplacian( con )


     ! time integration

     con = con + dt*mobility_con*Laplacian( dummy_con )
     phi = phi - dt*mobility_phi*( dfdphi - grad_coef_phi*Laplacian( phi ) )



     ! adjust order parameter in range

     where ( phi >= 0.99999 )  phi = 0.99999
     where ( phi < 0.00001 )   phi = 0.00001


     call Dislin_color_animation ( )


  end do time_loop


  call Disfin ( )
  call cpu_time ( finish )


  ! ===========================================================================
  !              print computed time on the screen and dislin animation
  ! ===========================================================================


  print *,'---------------------------------'
  print '("  Time       = ", f10.3," seconds." )', finish - start


contains


  ! ===========================================================================
  !                               Sub-programs
  ! ===========================================================================


  function Laplacian ( order_parameter) 
    implicit none

    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( in ) :: order_parameter 
    real ( kind = 8 ), dimension ( Nx, Ny ) :: Laplacian 
    integer ( kind = 4 ) :: i , j, jp, jm, ip, im , dx = 1, dy = 1


    do concurrent ( i = 1 : Nx, j= 1 : Ny )

       jp = j + 1
       jm = j - 1

       ip = i + 1
       im = i - 1

       if ( im == 0 ) im = Nx
       if ( ip == ( Nx + 1 ) ) ip = 1
       if ( jm == 0 ) jm = Ny
       if ( jp == ( Ny + 1 ) ) jp = 1

       Laplacian(i,j) = ( order_parameter(ip,j) + order_parameter(im,j) + & 
            order_parameter(i,jm) + order_parameter(i,jp) - &
            4.0*order_parameter(i,j))  / ( dx*dy )              


    end do


  end function Laplacian


  !----------------------------------------------------------------------------


  subroutine Dislin_color_animation ()

    call autres ( Nx, Ny )
    if ( mod( istep, nprint ) .eq. 0 ) then 
       call erase ( )  
       call graf3 ( 0.d0, 128.d0, 0.d0, 32.d0, 0.d0, 128.d0,&
            & 0.d0, 32.d0, 0.0d0, 1.0d0, 0.0d0, 0.2d0 )
       call crvmat ( phi, Nx, Ny, 1, 1 )   
       call endgrf
       call sendbf ( )
    end if


  end subroutine Dislin_color_animation



end program model_ch_ac_test
