
!      
!   Whole array Finite Difference Phase Field Code of Dendrite Solidification.
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


program fd_Kobayashi_model_test
  use Dislin
  implicit none


  ! ===========================================================================
  !                                parameters
  ! ===========================================================================


  integer ( kind = 4 ), parameter :: Nx = 300 , Ny = 300 
  integer (kind = 4 ) :: nsteps = 2500 , nprint = 100 , istep
  real ( kind = 8 )   :: dtime  = 1.0e-4, start , finish, dx = 0.03 , dy = 0.03
  real ( kind = 8 )   :: tau   = 0.0003 , epsilonb = 0.01 , mu = 1.0
  real ( kind = 8 )   :: kappa = 1.8 , delta = 0.02 , aniso = 6.0 , alpha = 0.9
  real ( kind = 8 )   :: gama  = 10.0 , teq   = 1.0 , theta0= 0.2
  real ( kind = 8 )   :: pix   = 4.0*atan(1.0), phi_old, term1, term2, theta, m
  integer ( kind = 4 ):: i, j, ip, im, jp, jm
  real ( kind = 8 ), dimension( Nx, Ny ):: phi, tempr , lap_phi, lap_tempr
  real ( kind = 8 ), dimension( Nx, Ny ):: phidx, phidy , epsil, epsilon_deriv


  call cpu_time ( start )


  ! ===========================================================================
  !                            introduce initial nuclei
  ! ===========================================================================


  phi = 0.0
  tempr = 0.0

  phi ( 148:153, 148:153 ) = 1.0


  ! ===========================================================================
  !                 Setting initial dislin routines for animation 
  ! ===========================================================================


  call Metafl ( 'cons' )
  call Disini ( )


  ! ===========================================================================
  !                         evolution of microstructure 
  ! ===========================================================================


  time_loop: do istep = 1, nsteps

     do concurrent ( i = 1:Nx, j = 1:Ny )

        jp = j + 1
        jm = j - 1

        ip = i + 1
        im = i - 1

        if ( im == 0 ) im = Nx
        if ( ip == ( Nx + 1) ) ip = 1
        if ( jm == 0 ) jm = Ny
        if ( jp == ( Ny + 1) ) jp = 1

        ! laplacian

        lap_phi(i,j) = ( phi(ip,j) + phi(im,j) + phi(i,jm) + phi(i,jp)&
             & - 4.0*phi(i,j)) / ( dx*dy )
        lap_tempr(i,j) = ( tempr(ip,j) + tempr(im,j) + tempr(i,jm) + &
             & tempr(i,jp) - 4.0*tempr(i,j)) / ( dx*dy )

        ! gradients

        phidx(i,j) = ( phi(ip,j) - phi(im,j) ) / dx
        phidy(i,j) = ( phi(i,jp) - phi(i,jm) ) / dy

        ! angle

        theta  = atan2( phidy(i,j),phidx(i,j) )

        ! epsilon and its derivative

        epsil(i,j) = epsilonb*( 1.0 + delta*cos(aniso* &
             & ( theta - theta0 ) ) )
        epsilon_deriv(i,j) = -epsilonb*aniso*delta*sin &
             & ( aniso*( theta - theta0 ) )

     end do


     do concurrent ( i = 1:Nx, j = 1:Ny )

        jp = j + 1
        jm = j - 1

        ip = i + 1
        im = i - 1

        if ( im == 0 ) im = Nx
        if ( ip == ( Nx + 1) ) ip = 1
        if ( jm == 0 ) jm = Ny
        if ( jp == ( Ny + 1) ) jp = 1

        phi_old = phi(i,j)

        ! term1 and term2

        term1 = ( epsil(i,jp)*epsilon_deriv(i,jp)*phidx(i,jp) &
             & - epsil(i,jm)*epsilon_deriv(i,jm)*phidx(i,jm) ) / dy
        term2 = -( epsil(ip,j)*epsilon_deriv(ip,j)*phidy(ip,j) &
             & - epsil(im,j)*epsilon_deriv(im,j)*phidy(im,j) ) / dx

        ! factor m

        m = alpha/pix*atan( gama*( teq - tempr(i,j) ) )

        ! time integration

        phi(i,j) = phi(i,j) + ( dtime/tau )*( term1 + term2 +&
             & epsil(i,j)**2*lap_phi(i,j) ) + &
             & phi_old*( 1.0 - phi_old )*( phi_old -0.5 + m )
        tempr(i,j) = tempr(i,j) + dtime*lap_tempr(i,j) &
             & + kappa*( phi(i,j) - phi_old )

     end do

     
     call Dislin_color_animation ( )


  end do time_loop


  call Disfin ( )
  call cpu_time ( finish )


  ! ===========================================================================
  !                    prints computed time on the screen 
  ! ===========================================================================


  print*,'---------------------------------'
  print '("  Time       = ", f10.3," seconds." )', finish - start


contains


  ! ===========================================================================
  !                            Sub-program
  ! ===========================================================================


  subroutine Dislin_color_animation ( )

    call autres ( Nx, Ny )
    if ( mod( istep, nprint ) .eq. 0 ) then 
       call erase ( )  
       call graf3 ( 0.d0, 300.d0, 0.d0, 50.d0, 0.d0, 300.d0,&
            & 0.d0, 50.d0, 0.0d0, 1.2d0, 0.0d0, 0.2d0 )
       call crvmat ( phi, Nx, Ny, 1, 1 )   
       call endgrf
       call sendbf ( )
    end if

  end subroutine Dislin_color_animation



end program fd_Kobayashi_model_test
