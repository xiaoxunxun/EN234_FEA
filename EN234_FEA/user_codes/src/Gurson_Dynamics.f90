
!==========================SUBROUTINE Gurson_dynamic ==============================
subroutine Gurson_dynamic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords, length_coord_array, &                               ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                       ! Input variables
    n_state_variables, initial_state_variables, &                                                        ! Input variables
    updated_state_variables,element_residual,element_deleted)                                            ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Globals, only: TIME,DTIME
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only:  dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,(x2),(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
               
    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)
          
    logical, intent( inout )      :: element_deleted                                        ! Set to .true. to delete element

    ! Local Variables
    integer      :: n_points,kint,kk

    real (prec)  ::  strain(6), dstrain(3,3),dspin(3,3)              ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6),stress1(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]                       ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array), Bbard(6,length_dof_array),B_Bar(6,length_dof_array)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  J,noneed
    real (prec)  ::  el_vol,yta,dL_kk ,dyta
    real (prec)  ::  du_increment (3,length_coord_array/3),  du_total (3,length_coord_array/3)
    real (prec)  ::  dF (3,3), F_mid(3,3),dL_ij(3,3), invert_F_mid(3,3)
    real (prec)  ::  dNdy(length_coord_array/3,3), dNbardy(length_coord_array/3,3)
    real (prec)  ::  dL_bard(3,3),dR1(3,3),dR2(3,3),invert_dR1(3,3),dR(3,3)



    x = reshape(element_coords,(/3,length_coord_array/3/))
    du_increment = reshape(dof_increment,(/3,length_coord_array/3/))
    du_total =  reshape(dof_total,(/3,length_coord_array/3/))


    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27


    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0

    el_vol = 0.d0
    yta = 0.d0
    dyta = 0.d0
    dNbardy = 0.d0

    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

        dF (1:3,1:3) = matmul(du_increment(1:3,1:n_nodes), dNdx(1:n_nodes,1:3))

        F_mid (1:3,1:3) = eye3_D + matmul ((du_total(1:3,1:n_nodes)+0.5d0*du_increment(1:3,1:n_nodes)), dNdx(1:n_nodes,1:3))

        call invert_small(F_mid,invert_F_mid,J)

        dL_ij(1:3,1:3) = matmul(dF(1:3,1:3), invert_F_mid(1:3,1:3))
        dL_kk = dL_ij(1,1) +dL_ij(2,2)+dL_ij(3,3)
        yta = yta + J * w(kint)* determinant
        dyta = dyta + J * dL_kk * w(kint)* determinant

        dNdy(1:n_nodes,1:3) = matmul (dNdx(1:n_nodes,1:3), invert_F_mid(1:3,1:3))
        dNbardy(1:n_nodes,1:3) = dNbardy(1:n_nodes,1:3) + J*dNdy(1:n_nodes,1:3)*w(kint)*determinant

        el_vol = el_vol + w(kint)*determinant
    end do

    yta = 1.d0/el_vol* yta
    dyta = 1.d0/(el_vol*yta) * dyta
    dNbardy(1:n_nodes,1:3) =  1.d0/(el_vol*yta) * dNbardy(1:n_nodes,1:3)


    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

        dF (1:3,1:3) = matmul(du_increment(1:3,1:n_nodes), dNdx(1:n_nodes,1:3))

        F_mid (1:3,1:3) = eye3_D + matmul ((du_total(1:3,1:n_nodes)+0.5d0*du_increment(1:3,1:n_nodes)), dNdx(1:n_nodes,1:3))

        call invert_small(F_mid,invert_F_mid,J)

        dNdy(1:n_nodes,1:3) = matmul (dNdx(1:n_nodes,1:3), invert_F_mid(1:3,1:3))

        B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)

        Bbard = 0.d0

        do kk=1,n_nodes
            Bbard(1,3*kk-2:3*kk) = dNbardy(kk,1:3)-dNdy(kk,1:3)
            Bbard(2,3*kk-2:3*kk) = dNbardy(kk,1:3)-dNdy(kk,1:3)
            Bbard(3,3*kk-2:3*kk) = dNbardy(kk,1:3)-dNdy(kk,1:3)
        end do

        B_Bar = B + (1.d0/3.d0)* Bbard


        dL_ij(1:3,1:3) = matmul(dF(1:3,1:3), invert_F_mid(1:3,1:3))
        dL_kk = dL_ij(1,1) +dL_ij(2,2)+dL_ij(3,3)

        dL_bard(1:3,1:3) = matmul(dF(1:3,1:3), invert_F_mid(1:3,1:3)) + 1.d0/3.d0 * (dyta - dL_kk)*eye3_D

        dstrain(1:3,1:3) = 0.5d0 * (dL_bard(1:3,1:3) + transpose(dL_bard(1:3,1:3)))
        dspin(1:3,1:3) = 0.5d0 * (dL_bard(1:3,1:3) - transpose(dL_bard(1:3,1:3)))



        dR1(1:3,1:3)= eye3_D - 0.5*dspin(1:3,1:3)
        dR2(1:3,1:3)= eye3_D + 0.5*dspin(1:3,1:3)

        call invert_small(dR1,invert_dR1,noneed)
        dR = matmul (invert_dR1,dR2)


        call gurson(element_properties,n_properties,8, initial_state_variables(8*(kint-1)+1:8*kint), &
            updated_state_variables(8*(kint-1)+1:8*kint),dstrain,dR,stress1)


        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B_Bar),stress1)*w(kint)*determinant

    end do
  
    return
end subroutine Gurson_dynamic




subroutine gurson(element_properties,n_properties,n_state_variables, initial_state_variables,&
    updated_state_variables,dstrain,dR,stress1)
    use Types
    use ParamIO
    use Element_Utilities
    use Globals, only : TIME, DTIME

    implicit none

    integer, intent (in)  :: n_properties
    integer, intent (in)  :: n_state_variables

    real (prec), intent(in) :: element_properties(n_properties)
    real (prec), intent(in) :: initial_state_variables(n_state_variables)
    real (prec), intent(in) :: dstrain(3,3)
    real (prec), intent(in) :: dR (3,3)

    real (prec), intent(out) :: stress1(6)
    real (prec), intent(out) :: updated_state_variables(n_state_variables)

    integer      :: maxit,nit,k,n_points

    real (prec)  :: stress0(6),stressdR(6)
    real (prec)  :: ematrix, vf, ff_bar,P_star,f_star
    real (prec)  :: effect_dstrain(3,3), stress0_reshape(3,3), effect_stress(3,3)
    real (prec)  :: stressdR_reshape(3,3), S_star(3,3),EE(3,3)
    real (prec)  :: EEE(9)
    real (prec)  :: stress_star,Fi
    real (prec)  :: ev,nu,yy,z0,m,q1,q2,q3,fn,xn,sn,fc,ff
    real (prec)  :: stress1_1(3,3)
    real (prec)  :: dde,ddv
    real (prec)  :: dfdx,dfdy,F1,F2
    real (prec)  :: df1de,df1dv,df2de,df2dv
    real (prec)  :: dematrix, ematrix1,vf1
    real (prec)  :: tol,wnorm,err1,ddde,dddv


    ev = element_properties(1)
    nu = element_properties(2)
    yy = element_properties(3)
    z0 = element_properties(4)
    m  = element_properties(5)
    q1 = element_properties(6)
    q2 = element_properties(7)
    q3 = element_properties(8)
    fn = element_properties(9)
    xn = element_properties(10)
    sn = element_properties(11)
    fc = element_properties(12)
    ff = element_properties(13)


    stress0 = initial_state_variables(1:6)
    ematrix = initial_state_variables(7)
    vf = initial_state_variables(8)

    effect_dstrain(1:3,1:3) = dstrain(1:3,1:3) - 1.d0/3.d0 * eye3_D*(dstrain(1,1)+dstrain(2,2)+dstrain(3,3))
    stress0 = rotatesymvec(stress0,dR)
    stress0_reshape = 0.d0
    stress0_reshape(1,1) = stress0(1)
    stress0_reshape(2,2) = stress0(2)
    stress0_reshape(3,3) = stress0(3)
    stress0_reshape(1,2) = stress0(4)
    stress0_reshape(2,1) = stress0(4)
    stress0_reshape(1,3) = stress0(5)
    stress0_reshape(3,1) = stress0(5)
    stress0_reshape(2,3) = stress0(6)
    stress0_reshape(3,2) = stress0(6)


    effect_stress(1:3,1:3) = stress0_reshape(1:3,1:3) -1.d0/3.d0 * eye3_D*(stress0(1)+stress0(2)+stress0(3))

    S_star(1:3,1:3) = ev/(1.d0+nu)*effect_dstrain(1:3,1:3) + effect_stress(1:3,1:3)
    P_star =  1.d0/3.d0 * (stress0(1)+stress0(2)+stress0(3)) &
        + ev/(3.d0*(1.d0-2.d0*nu))*(dstrain(1,1)+dstrain(2,2)+dstrain(3,3))

    EE(1:3,1:3) = S_star(1:3,1:3)*S_star(1:3,1:3)
    EEE(1:9) = reshape(EE(1:3,1:3),(/9/))

    stress_star = sqrt (1.5d0*sum(EEE))

    ff_bar = (q1+sqrt(q1**2.d0-q3))/q3

    if (vf < fc) then
        f_star = vf

    else if( fc < vf .and. vf < ff) then
        f_star = fc + ((ff_bar-fc)/(ff-fc)*(vf-fc))

    else
        f_star = ff_bar

    end if


    Fi = (stress_star**2.d0/yy**2.d0) +2.d0*q1*f_star*cosh(1.5d0*q2*(P_star/yy))-(1.d0+q3*(f_star)**2.d0)

    if (Fi < 10**(-8)) then

        stress1_1(1:3,1:3) = S_star(1:3,1:3) +p_star*eye3_d

    else

        Fi = sqrt(Fi)


tol =0.0001d0

maxit =6

nit = 0
err1 =1.d0


dde = 0.d0
ddv = 0.d0


do while (err1>tol .and. nit<maxit)

nit = nit + 1


!h = (stress_star -1.5*ev/(1+nu)*dde)

!p = (p_star - ev/(3*(1-2*nu))*ddv)

!fi = ((h + 1.5*ev/(1+nu)*dde)**2/yy**2 + 2*q1*f_star*cosh(1.5*q2*(p + ev/(3*(1-2*nu))*ddv)/yy) - (1+q3*f_star**2))**0.5
!dfdx = ((2.d0*(stress_star -1.5*ev/(1+nu)*dde) + (3.d0*dde*ev)/(nu + 1.d0))/(2*(yy**2)*(((stress_star -1.5*ev/(1+nu)*dde) &
!        + (3.d0*dde*ev)/(2.d0*nu + 2.d0))**2/yy**2 &
!        - (f_star**2.d0)*q3 + 2.d0*f_star*q1*cosh((3.d0*q2*((p_star - ev/(3*(1-2*nu))*ddv) &
!        - (ddv*ev)/(6.d0*nu - 3.d0)))/(2.d0*yy)) - 1.d0)**(0.5d0)))
!
!
!dfdy = (3.d0*f_star*q1*q2*sinh((3.d0*q2*((p_star - ev/(3.d0*(1.d0-2.d0*nu))*ddv)&
!        - (ddv*ev)/(6.d0*nu - 3.d0)))/(2.d0*yy)))/(2.d0*yy*(((stress_star -1.5d0*ev/(1.d0+nu)*dde)&
!         + (3.d0*dde*ev)/(2.d0*nu + 2.d0))**2.d0/(yy**2.d0) &
!          - (f_star**2.d0)*q3 + 2.d0*f_star*q1*cosh((3.d0*q2*((p_star - ev/(3.d0*(1.d0-2.d0*nu))*ddv)&
!           - (ddv*ev)/(6.d0*nu - 3.d0)))/(2.d0*yy)) - 1.d0)**(0.5d0))




!fi = ((stress_star**2/yy**2) + 2*q1*f_star*cosh(1.5*q2*p_star /yy) - (1+q3*f_star**2))**0.5



dfdx = stress_star/((yy**2.d0)*((stress_star**2.d0)/(yy**2.d0) - (f_star**2.d0)*q3 &
+ 2.d0*f_star*q1*cosh((3.d0*p_star*q2)/(2.d0*yy)) - 1.d0)**(0.5d0))

dfdy = (3.d0*f_star*q1*q2*sinh((3.d0*p_star*q2)/(2.d0*yy)))/(2.d0*yy*((stress_star**2.d0)/(yy**2.d0) - (f_star**2.d0)*q3 &
+ 2.d0*f_star*q1*cosh((3.d0*p_star*q2)/(2.d0*yy)) - 1.d0)**(0.5d0))



df1de = (z0*(stress_star**2/(yy**4*(stress_star**2/yy**2 - f_star**2*q3 + 2*f_star*q1*cosh((3*p_star*q2)/(2*yy)) - 1))&
 + (f_star**2*q1**2*q2**2*sinh((3*p_star*q2)/(2*yy))**2)/(2*yy**2*(stress_star**2/yy**2 - f_star**2*q3 &
 + 2*f_star*q1*cosh((3*p_star*q2)/(2*yy)) - 1)))**(1/2))/DTIME

df1dv = 0

df2de = 0

df2dv = (z0*(stress_star**2/(yy**4*(stress_star**2/yy**2 - f_star**2*q3 + 2*f_star*q1*cosh((3*p_star*q2)/(2*yy)) - 1)) &
+ (f_star**2*q1**2*q2**2*sinh((3*p_star*q2)/(2*yy))**2)/(2*yy**2*(stress_star**2/yy**2 - f_star**2*q3 &
+ 2*f_star*q1*cosh((3*p_star*q2)/(2*yy)) - 1)))**(1/2))/DTIME


F1 = -(stress_star*((stress_star**2/yy**2 - f_star**2*q3 &
+ 2*f_star*q1*cosh((3*p_star*q2)/(2*yy)) - 1)**(1/2))**m)/(yy**2*(stress_star**2/yy**2 - f_star**2*q3&
 + 2*f_star*q1*cosh((3*p_star*q2)/(2*yy)) - 1)**(1/2))


F2 = -(3*f_star*q1*q2*sinh((3*p_star*q2)/(2*yy))*((stress_star**2/yy**2 &
- f_star**2*q3 + 2*f_star*q1*cosh((3*p_star*q2)/(2*yy)) - 1)**(1/2))**m)/(2*yy*(stress_star**2/yy**2 &
- f_star**2*q3 + 2*f_star*q1*cosh((3*p_star*q2)/(2*yy)) - 1)**(1/2))


! Fi = (((stress_star -1.5*ev/(1+nu)*dde) + 1.5*ev/(1+nu)*dde)**2/yy**2 &
!   + 2*q1*f_star*cosh(1.5*q2*((p_star - ev/(3*(1-2*nu))*ddv) + ev/(3*(1-2*nu))*ddv)/yy)&
!    - (1+q3*f_star**2))**0.5

!F1 = ((dfdx**2 +2/9*dfdy**2)**0.5)*(dde/DTIME*z0) - (dfdx)*(Fi**m)
!F2 = ((dfdx**2 +2/9*dfdy**2)**0.5)*(ddv/DTIME*z0) - (dfdy)*(Fi**m)
!F1 = (dde*z0*(stress_star**2.d0/((yy**4.d0)*((stress_star - (3.d0*dde*ev)/(2.d0*(nu + 1.d0))  &
!   +(3.d0*dde*ev)/(2.d0*nu + 2.d0))**2.d0/yy**2.d0 - (f_star**2.d0)*q3 + 2.d0*f_star*q1*cosh((3*p_star*q2)/(2*yy)) - 1)) &
!   + (f_star**2*q1**2*q2**2*sinh((3*p_star*q2)/(2*yy))**2)/(2*yy**2*((stress_star - (3*dde*ev)/(2*(nu + 1)) &
!   + (3*dde*ev)/(2*nu + 2))**2/yy**2 - f_star**2*q3 + 2*f_star*q1*cosh((3*p_star*q2)/(2*yy)) - 1)))**(1/2))/DTIME &
!   - (stress_star*((stress_star**2/yy**2 - f_star**2*q3 + 2*f_star*q1*cosh((3*p_star*q2)/(2*yy)) &
!   - 1)**(1/2))**m)/(yy**2*(stress_star**2/yy**2 - f_star**2*q3 &
!   + 2*f_star*q1*cosh((3*p_star*q2)/(2*yy)) - 1)**(1/2))
!
!F2 = (ddv*z0*(stress_star**2.d0/((yy**4.d0)*((stress_star - (3.d0*dde*ev)/(2.d0*(nu + 1.d0)) &
!     + (3.d0*dde*ev)/(2.d0*nu + 2.d0))**2.d0/yy**2.d0 - (f_star**2.d0)*q3 + 2.d0*f_star*q1*cosh((3*p_star*q2)/(2*yy)) - 1))&
!     + (f_star**2*q1**2*q2**2*sinh((3*p_star*q2)/(2*yy))**2)/(2*yy**2*((stress_star - (3*dde*ev)/(2*(nu + 1)) &
!     + (3*dde*ev)/(2*nu + 2))**2/yy**2 - f_star**2*q3 + 2*f_star*q1*cosh((3*p_star*q2)/(2*yy)) - 1)))**(1/2))/DTIME &
!     - (3*f_star*q1*q2*sinh((3*p_star*q2)/(2*yy))*((stress_star**2/yy**2 - f_star**2*q3 &
!     + 2*f_star*q1*cosh((3*p_star*q2)/(2*yy)) - 1)**(1/2))**m)/(2*yy*(stress_star**2/yy**2 - f_star**2*q3 &
!     + 2*f_star*q1*cosh((3*p_star*q2)/(2*yy)) - 1)**(1/2))
!
!
!df1de = (z0*(stress_star**2/(yy**4*((stress_star - (3*dde*ev)/(2*(nu + 1)) + (3*dde*ev)/(2*nu + 2))**2/yy**2 &
! - f_star**2*q3 + 2*f_star*q1*cosh((3*p_star*q2)/(2*yy)) - 1)) &
! + (f_star**2*q1**2*q2**2*sinh((3*p_star*q2)/(2*yy))**2)/(2*yy**2*((stress_star - (3*dde*ev)/(2*(nu + 1)) &
! + (3*dde*ev)/(2*nu + 2))**2/yy**2 - f_star**2*q3 + 2*f_star*q1*cosh((3*p_star*q2)/(2*yy)) - 1)))**(1/2))/DTIME &
!  - (dde*z0*((2*stress_star**2*((3*ev)/(2*nu + 2) - (3*ev)/(2*(nu + 1)))*(stress_star - (3*dde*ev)/(2*(nu + 1))&
!   + (3*dde*ev)/(2*nu + 2)))/(yy**6*(stress_star**2/yy**2 - f_star**2*q3 + 2*f_star*q1*cosh((3*p_star*q2)/(2*yy)) - 1)**2) &
!   + (f_star**2*q1**2*q2**2*sinh((3*p_star*q2)/(2*yy))**2*((3*ev)/(2*nu + 2) - (3*ev)/(2*(nu + 1)))*(stress_star&
!    - (3*dde*ev)/(2*(nu + 1)) + (3*dde*ev)/(2*nu + 2)))/(yy**4*(stress_star**2/yy**2 - f_star**2*q3 &
!    + 2*f_star*q1*cosh((3*p_star*q2)/(2*yy)) - 1)**2)))/(2*DTIME*(stress_star**2/(yy**4*(stress_star**2/yy**2 - f_star**2*q3 &
!    + 2*f_star*q1*cosh((3*p_star*q2)/(2*yy)) - 1))&
!     + (f_star**2*q1**2*q2**2*sinh((3*p_star*q2)/(2*yy))**2)/(2*yy**2*(stress_star**2/yy**2 &
!    - f_star**2*q3 + 2*f_star*q1*cosh((3*p_star*q2)/(2*yy)) - 1)))**(1/2)) + (stress_star*((stress_star**2/yy**2 &
!    - f_star**2*q3 + 2*f_star*q1*cosh((3*p_star*q2)/(2*yy)) - 1)**(1/2))**m*((3*ev)/(2*nu + 2)&
!     - (3*ev)/(2*(nu + 1)))*(stress_star - (3*dde*ev)/(2*(nu + 1)) + (3*dde*ev)/(2*nu + 2)))/(yy**4*(stress_star**2/yy**2 &
!     - f_star**2*q3 + 2*f_star*q1*cosh((3*p_star*q2)/(2*yy)) - 1)**(3/2))
!
!df1dv = 0.d0
!
!df2de = (3*f_star*q1*q2*sinh((3*p_star*q2)/(2*yy))*((stress_star**2/yy**2 &
!- f_star**2*q3 + 2*f_star*q1*cosh((3*p_star*q2)/(2*yy)) - 1)**(1/2))**m*((3*ev)/(2*nu + 2) &
!- (3*ev)/(2*(nu + 1)))*(stress_star - (3*dde*ev)/(2*(nu + 1)) + (3*dde*ev)/(2*nu + 2)))/(2*yy**3*(stress_star**2/yy**2&
! - f_star**2*q3 + 2*f_star*q1*cosh((3*p_star*q2)/(2*yy)) - 1)**(3/2)) &
! - (ddv*z0*((2*stress_star**2*((3*ev)/(2*nu + 2) - (3*ev)/(2*(nu + 1)))*(stress_star &
! - (3*dde*ev)/(2*(nu + 1)) + (3*dde*ev)/(2*nu + 2)))/(yy**6*(stress_star**2/yy**2 - f_star**2*q3 &
! + 2*f_star*q1*cosh((3*p_star*q2)/(2*yy)) - 1)**2) + (f_star**2*q1**2*q2**2*sinh((3*p_star*q2)/(2*yy))**2*((3*ev)/(2*nu + 2) &
! - (3*ev)/(2*(nu + 1)))*(stress_star - (3*dde*ev)/(2*(nu + 1)) + (3*dde*ev)/(2*nu + 2)))/(yy**4*(stress_star**2/yy**2 &
!  - f_star**2*q3 + 2*f_star*q1*cosh((3*p_star*q2)/(2*yy)) - 1)**2)))/(2*DTIME*(stress_star**2/(yy**4*(stress_star**2/yy**2 &
!  - f_star**2*q3 + 2*f_star*q1*cosh((3*p_star*q2)/(2*yy)) - 1)) + &
!  (f_star**2*q1**2*q2**2*sinh((3*p_star*q2)/(2*yy))**2)/(2*yy**2*(stress_star**2/yy**2 - f_star**2*q3 &
!  + 2*f_star*q1*cosh((3*p_star*q2)/(2*yy)) - 1)))**(1/2))
!
!
!df2dv = (z0*(stress_star**2/(yy**4*((stress_star - (3*dde*ev)/(2*(nu + 1)) + (3*dde*ev)/(2*nu + 2))**2/yy**2&
! - f_star**2*q3 + 2*f_star*q1*cosh((3*p_star*q2)/(2*yy)) - 1)) &
! + (f_star**2*q1**2*q2**2*sinh((3*p_star*q2)/(2*yy))**2)/(2*yy**2*((stress_star - (3*dde*ev)/(2*(nu + 1)) &
! + (3*dde*ev)/(2*nu + 2))**2/yy**2 - f_star**2*q3 + 2*f_star*q1*cosh((3*p_star*q2)/(2*yy)) - 1)))**(1/2))/DTIME

ddde = (df1dv*F2 - df2dv*F1)/(df1de*df2dv-df2de*df1dv)
dddv = -(df1de*F2 - df2de*F1)/(df1de*df2dv-df2de*df1dv)


!solve dde ddv
dde = dde + ddde
ddv = ddv + dddv

wnorm = dde*dde + ddv*ddv
err1 = ddde*ddde + dddv*dddv
err1 = sqrt (err1/wnorm)


 end do

stress1_1(1:3,1:3) = S_star(1:3,1:3) - dde/stress_star*(ev/(1.d0+nu))*1.5d0*S_star(1:3,1:3) &
  + (p_star-(ev/(3.d0*(1.d0-2.d0*nu)))*ddv)*eye3_d

end if


stress1(1) = stress1_1(1,1)
stress1(2) = stress1_1(2,2)
stress1(3) = stress1_1(3,3)
stress1(4) = stress1_1(1,2)
stress1(5) = stress1_1(1,3)
stress1(6) = stress1_1(2,3)

if (Fi< 10**(-8)) then
dematrix = 0.d0

else

dfdx = stress_star/(yy**2*(stress_star**2/yy**2 - f_star**2*q3 &
+ 2*f_star*q1*cosh((3*p_star*q2)/(2*yy)) - 1)**(1/2))

dfdy = (3*f_star*q1*q2*sinh((3*p_star*q2)/(2*yy)))/(2*yy*(stress_star**2/yy**2 - f_star**2*q3 &
+ 2*f_star*q1*cosh((3*p_star*q2)/(2*yy)) - 1)**(1/2))


dematrix = z0*DTIME/(1.d0-vf)*(Fi**m)*(((dfdx**2+2.d0/9.d0*dfdy**2))**(-0.5d0)) &
* (dfdx*((stress_star -1.5*ev/(1.d0+nu)*dde) +1.d0/3.d0*dfdy*(p_star - ev/(3.d0*(1.d0-2.d0*nu))*ddv)))

 end if

  ematrix1 = ematrix + dematrix

  vf1 = 1.d0 + (vf-1.d0)*exp(-ddv) + fn*dematrix/(sn*sqrt(2.d0*pi))*exp(-0.5d0*((ematrix-xn)/sn)**2)


   updated_state_variables(1:6) = stress1
   updated_state_variables(7)  = ematrix1
  updated_state_variables(8) = vf1


   end subroutine gurson




!==========================SUBROUTINE _Gurson_Dynamics ==============================
 subroutine fieldvars_Gurson_Dynamics(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords,length_coord_array, &                                ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                      ! Input variables
    n_state_variables, initial_state_variables,updated_state_variables, &                               ! Input variables
    n_field_variables,field_variable_names, &                                                           ! Field variable definition
    nodal_fieldvariables)      ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only: dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only: dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element
    integer, intent( in )         :: n_field_variables                                      ! # field variables

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    character (len=100), intent(in) :: field_variable_names(n_field_variables)

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,x2,(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
    real( prec ), intent( in )    :: updated_state_variables(n_state_variables)             ! State variables at end of time step

    real( prec ), intent( out )   :: nodal_fieldvariables(n_field_variables,n_nodes)        ! Nodal field variables

    ! Local Variables
    logical      :: strcmp

    integer      :: n_points,kint,k

    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real( prec ) ::  updated_state_variables_current(8)
    real( prec ) :: stress(6)
    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27



    call initialize_integration_points(n_points, n_nodes, xi, w)

   nodal_fieldvariables=0.d0

    do kint = 1, n_points

        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)

        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))

        call invert_small(dxdxi,dxidx,determinant)

  updated_state_variables_current(1:8)=updated_state_variables(((kint-1)*8+1):kint*8)
  stress(1:6)=updated_state_variables_current(1:6)


        ! In the code below the strcmp( string1, string2, nchar) function returns true if the first nchar characters in strings match
        do k = 1,n_field_variables
            if (strcmp(field_variable_names(k),'S11',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(1)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S22',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(2)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S33',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(3)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S12',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(4)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S13',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(5)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S23',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(6)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'Vf',2) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes)&
                 + updated_state_variables_current(8)*N(1:n_nodes)*determinant*w(kint)
            endif
        end do
end do

    return
end subroutine fieldvars_Gurson_Dynamics
