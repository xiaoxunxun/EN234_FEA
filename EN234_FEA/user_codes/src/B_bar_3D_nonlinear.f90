!     Subroutines for basic 3D nonlinear



!==========================SUBROUTINE el_nonlinelast ==============================
subroutine el_nonlinelast (lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties, element_coords, length_coord_array, &                      ! Input variables
    dof_increment, dof_total, length_dof_array, &                                                ! Input variables
    n_state_variables, initial_state_variables, &                                                ! Input variables
    updated_state_variables,element_stiffness,element_residual, fail)                          ! Output variables
    use Types
    use ParamIO
    !  use Globals, only: TIME,DTIME  For a time dependent problem uncomment this line to access the time increment and total time
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only : dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    use Element_Utilities, only : dNbardx => vol_avg_shape_function_derivatives_3D

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
  
    logical, intent( out )        :: fail                                                   ! Set to .true. to force a timestep cutback
    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_stiffness(length_dof_array,length_dof_array)   ! Element stiffness (ROW,COLUMN)
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)
          

    ! Local Variables
    integer      :: n_points,kint,k,kk

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                          ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  el_vol                            ! Element Volume
    real (prec)  ::  Bbard (6,length_dof_array)
    real (prec)  ::  B_Bar (6,length_dof_array)


    fail = .false.
    
    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

        el_vol = 0.d0
        dNbardx = 0.d0

        do kint = 1, n_points

        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        dNbardx(1:n_nodes,1) = dNbardx(1:n_nodes,1) + dNdx(1:n_nodes,1)*w(kint)*determinant
        dNbardx(1:n_nodes,2) = dNbardx(1:n_nodes,2) + dNdx(1:n_nodes,2)*w(kint)*determinant
        dNbardx(1:n_nodes,3) = dNbardx(1:n_nodes,3) + dNdx(1:n_nodes,3)*w(kint)*determinant
        el_vol = el_vol + w(kint)*determinant
        end do

        dNbardx(1:n_nodes,1:3) = (1.d0/el_vol)*dNbardx(1:n_nodes,1:3)


    element_residual = 0.d0
    element_stiffness = 0.d0
  
    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)


        Bbard = 0.d0

  !      if (element_identifier == 1002) then

        do kk=1,n_nodes
        Bbard(1,3*kk-2:3*kk) = dNbardx(kk,1:3)-dNdx(kk,1:3)
        Bbard(2,3*kk-2:3*kk) = dNbardx(kk,1:3)-dNdx(kk,1:3)
        Bbard(3,3*kk-2:3*kk) = dNbardx(kk,1:3)-dNdx(kk,1:3)
        end do
  !      endif

        B_Bar = B + (1.d0/3.d0)* Bbard

        strain = matmul(B_Bar,dof_total)
        dstrain = matmul(B_Bar,dof_increment)


call hypoelastic_material(strain, dstrain, &
 element_properties,n_properties,stress,D)


element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B_Bar),stress)*w(kint)*determinant

element_stiffness(1:3*n_nodes,1:3*n_nodes) = element_stiffness(1:3*n_nodes,1:3*n_nodes) &
            + matmul(transpose(B_Bar(1:6,1:3*n_nodes)),matmul(D,B_Bar(1:6,1:3*n_nodes)))*w(kint)*determinant

    end do
  
    return
end subroutine el_nonlinelast


!==========================SUBROUTINE el_linelast_3dbasic_dynamic ==============================
!subroutine nonlinear_3dbar_dynamic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
!    n_properties, element_properties,element_coords, length_coord_array, &                               ! Input variables
!    dof_increment, dof_total, length_dof_array,  &                                                       ! Input variables
!    n_state_variables, initial_state_variables, &                                                        ! Input variables
!    updated_state_variables,element_residual,element_deleted)                                            ! Output variables
!    use Types
!    use ParamIO
!    use Mesh, only : node
!    use Element_Utilities, only : N => shape_functions_3D
!    use Element_Utilities, only:  dNdxi => shape_function_derivatives_3D
!    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
!    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
!    use Element_Utilities, only : dxdxi => jacobian_3D
!    use Element_Utilities, only : initialize_integration_points
!    use Element_Utilities, only : calculate_shapefunctions
!    use Element_Utilities, only : invert_small
!    implicit none
!
!    integer, intent( in )         :: lmn                                                    ! Element number
!    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
!    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
!    integer, intent( in )         :: n_properties                                           ! # properties for the element
!    integer, intent( in )         :: length_coord_array                                     ! Total # coords
!    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
!    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element
!
!    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
!    !  type node
!    !      sequence
!    !      integer :: flag                          ! Integer identifier
!    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
!    !      integer :: n_coords                      ! Total no. coordinates for the node
!    !      integer :: dof_index                     ! Index of first DOF in dof array
!    !      integer :: n_dof                         ! Total no. of DOF for node
!    !   end type node
!    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element
!
!    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,(x2),(x3) for each node in turn
!    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
!    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment
!
!    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
!    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
!
!    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
!    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)
!
!    logical, intent( inout )      :: element_deleted                                        ! Set to .true. to delete element
!
!    ! Local Variables
!    integer      :: n_points,kint
!
!    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
!    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
!    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
!    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
!    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
!    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
!    real (prec)  :: E, xnu, D44, D11, D12              ! Material properties
!    !
!    !     Subroutine to compute element force vector for a linear elastodynamic problem
!    !     El props are:
!
!    !     element_properties(1)         Young's modulus
!    !     element_properties(2)         Poisson's ratio
!
!    x = reshape(element_coords,(/3,length_coord_array/3/))
!
!    if (n_nodes == 4) n_points = 1
!    if (n_nodes == 10) n_points = 4
!    if (n_nodes == 8) n_points = 8
!    if (n_nodes == 20) n_points = 27
!
!    call initialize_integration_points(n_points, n_nodes, xi, w)
!
!    element_residual = 0.d0
!
!    D = 0.d0
!    E = element_properties(1)
!    xnu = element_properties(2)
!    d44 = 0.5D0*E/(1+xnu)
!    d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
!    d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
!    D(1:3,1:3) = d12
!    D(1,1) = d11
!    D(2,2) = d11
!    D(3,3) = d11
!    D(4,4) = d44
!    D(5,5) = d44
!    D(6,6) = d44
!
!    !     --  Loop over integration points
!    do kint = 1, n_points
!        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
!        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
!        call invert_small(dxdxi,dxidx,determinant)
!        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
!        B = 0.d0
!        B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
!        B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
!        B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
!        B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
!        B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
!        B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
!        B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
!        B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
!        B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)
!
!        strain = matmul(B,dof_total)
!        dstrain = matmul(B,dof_increment)
!
!        stress = matmul(D,strain+dstrain)
!        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant
!
!    end do
!
!    return
!end subroutine nonlinear_3dbar_dynamic


!==========================SUBROUTINE fieldvars_linelast_3dbasic ==============================
subroutine fieldvars_nonlinear(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
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
    use Element_Utilities, only : dNbardx => vol_avg_shape_function_derivatives_3D

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
  
    integer      :: n_points,kint,k,kk

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                          ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  sdev(6)                           ! Deviatoric stress
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  p, smises                          ! Pressure and Mises stress
    real (prec)  ::  el_vol                            ! Element Volume
    real (prec)  ::  Bbard (6,length_dof_array)
    real (prec)  ::  B_Bar (6,length_dof_array)


    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    nodal_fieldvariables = 0.d0

        el_vol = 0.d0
        dNbardx = 0.d0

        do kint = 1, n_points

        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)

        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))

        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

        dNbardx(1:n_nodes,1) = dNbardx(1:n_nodes,1) + dNdx(1:n_nodes,1)*w(kint)*determinant
        dNbardx(1:n_nodes,2) = dNbardx(1:n_nodes,2) + dNdx(1:n_nodes,2)*w(kint)*determinant
        dNbardx(1:n_nodes,3) = dNbardx(1:n_nodes,3) + dNdx(1:n_nodes,3)*w(kint)*determinant
        el_vol = el_vol + w(kint)*determinant
        end do

        dNbardx(1:n_nodes,1:3) = (1.d0/el_vol)*dNbardx(1:n_nodes,1:3)


  
    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)

        Bbard = 0.d0

        do kk=1,n_nodes
        Bbard(1,3*kk-2:3*kk) = dNbardx(kk,1:3)-dNdx(kk,1:3)
        Bbard(2,3*kk-2:3*kk) = dNbardx(kk,1:3)-dNdx(kk,1:3)
        Bbard(3,3*kk-2:3*kk) = dNbardx(kk,1:3)-dNdx(kk,1:3)
        end do


        B_Bar = B + (1.d0/3.d0)* Bbard

        strain = matmul(B_Bar,dof_total)
        dstrain = matmul(B_Bar,dof_increment)


call hypoelastic_material(strain, dstrain, &
 element_properties,n_properties,stress,D)

        p = sum(stress(1:3))/3.d0
        sdev = stress
        sdev(1:3) = sdev(1:3)-p
        smises = dsqrt( dot_product(sdev(1:3),sdev(1:3)) + 2.d0*dot_product(sdev(4:6),sdev(4:6)) )*dsqrt(1.5d0)
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
            else if (strcmp(field_variable_names(k),'SMISES',6) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + smises*N(1:n_nodes)*determinant*w(kint)
            endif
        end do
 
    end do
  
    return
end subroutine fieldvars_nonlinear


subroutine hypoelastic_material(strain, dstrain, &
 element_properties,n_properties,stress,D)

   use Types
   use ParamIO

   implicit none

   integer, intent( in ) :: n_properties

   real (prec), intent( in )  :: strain(6),dstrain(6)
   real (prec), intent( in )  :: element_properties(n_properties)
   real (prec), intent( out ) :: stress(6)
   real (prec), intent( out ) :: D(6,6)


   real (prec) :: Et, Es
   real (prec) :: Nn1,Nn2,Nn3,Nn4,Nn5,Nn6
   real (prec) :: stresszero,strainzero,nzero,kzero
   real (prec) :: effect_stress
   real (prec) :: total_strain(6)
   real (prec) :: deviatoric_strain(6)
   real (prec) :: volumetric_strain
   real (prec) :: effect_strain
   real (prec) :: e_dyadic_e(6,6)
   real (prec) :: Aa(6,6)
   real (prec) :: Bb(6,6)

      stresszero = element_properties(1)
      strainzero = element_properties(2)
      nzero = element_properties(3)
      kzero = element_properties(4)

      total_strain = strain + dstrain


   volumetric_strain = total_strain(1)+total_strain(2)+total_strain(3)
   deviatoric_strain(1) = total_strain(1)-1.d0/3.d0*volumetric_strain
   deviatoric_strain(2) = total_strain(2)-1.d0/3.d0*volumetric_strain
   deviatoric_strain(3) = total_strain(3)-1.d0/3.d0*volumetric_strain
   deviatoric_strain(4) =  0.5*total_strain(4)
   deviatoric_strain(5) =  0.5*total_strain(5)
   deviatoric_strain(6) =  0.5*total_strain(6)


   effect_strain = dot_product(deviatoric_strain(1:3),deviatoric_strain(1:3)) &
                 + 2.d0*dot_product(deviatoric_strain(4:6),deviatoric_strain(4:6))

   effect_strain = dsqrt(2.d0*effect_strain/3.d0)

if (effect_strain == 0.d0) then
  Nn1 = nzero/(nzero-1.d0)-effect_strain/strainzero
  Nn2 = dsqrt((nzero**2+1.d0)/(nzero-1.d0)**2.d0-(nzero/(nzero-1.d0)-effect_strain/strainzero)**2)
  Et = Nn1/(strainzero*Nn2)
  Es = 0.d0

else if (effect_strain < strainzero ) then
   Nn3= (1.d0+nzero**2)/((nzero-1.d0)**2)
   Nn4= (nzero/(nzero-1.d0)-effect_strain/strainzero)**2
   effect_stress = stresszero * (dsqrt(Nn3-Nn4) -1.d0/(nzero-1.d0))

   Nn5 = nzero/(nzero-1.d0)-effect_strain/strainzero
   Nn6 = dsqrt((nzero**2+1.d0)/(nzero-1.d0)**2.d0-(nzero/(nzero-1.d0)-effect_strain/strainzero)**2)
   Et= Nn5/(strainzero*Nn6)
   Es = effect_stress/ effect_strain

else

effect_stress = stresszero * ((effect_strain/strainzero)**(1.d0/nzero))
Et = (1.d0/nzero)*(effect_strain**(1.d0/nzero-1.d0)*(1.d0/strainzero)**1/nzero)
Es = effect_stress/ effect_strain

end if


stress(1:6) = 0.d0

if (effect_strain == 0.d0) then

stress(1:3) = kzero * volumetric_strain

else

stress(1) = 2.d0/3.d0* effect_stress * deviatoric_strain(1)/effect_strain + kzero * volumetric_strain
stress(2) = 2.d0/3.d0* effect_stress * deviatoric_strain(2)/effect_strain + kzero * volumetric_strain
stress(3) = 2.d0/3.d0* effect_stress * deviatoric_strain(3)/effect_strain + kzero * volumetric_strain
stress(4) = 2.d0/3.d0* effect_stress * deviatoric_strain(4)/effect_strain
stress(5) = 2.d0/3.d0* effect_stress * deviatoric_strain(5)/effect_strain
stress(6) = 2.d0/3.d0* effect_stress * deviatoric_strain(6)/effect_strain

end if


e_dyadic_e = spread (deviatoric_strain,dim=2,ncopies=6) * spread(deviatoric_strain,dim=1,ncopies=6)

Aa = 0.d0
Bb = 0.d0

Aa(1,1) = 2.d0
Aa(2,2) = 2.d0
Aa(3,3) = 2.d0
Aa(4,4) = 1.d0
Aa(5,5) = 1.d0
Aa(6,6) = 1.d0

Bb(1:3,1:3) = 1.d0

if (effect_strain == 0.d0) then

D =  Et/3.d0*Aa + (kzero-2.d0*Et/9.d0)*Bb

else

D = 4.d0/(9.d0*effect_strain**2.d0)*(Et-Es)*e_dyadic_e + Es/3.d0*Aa + (kzero-2.d0*Es/9.d0)*Bb
end if

write(6,*) total_strain
end subroutine hypoelastic_material

