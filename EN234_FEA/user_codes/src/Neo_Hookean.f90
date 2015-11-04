!     Subroutines for Neo_Hookean

!==========================SUBROUTINE Neo_Hookean ==============================
subroutine Neo_Hookean(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties, element_coords, length_coord_array, &                      ! Input variables
    dof_increment, dof_total, length_dof_array, &                                                ! Input variables
    n_state_variables, initial_state_variables, &                                                ! Input variables
    updated_state_variables,element_stiffness,element_residual, fail)      ! Output variables                          ! Output variables
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
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
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
    integer      :: n_points,kint,i
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  B_star(9,length_dof_array)
    real (prec)  ::  G (6,9)
    real (prec)  ::  Sigma(length_dof_array,length_dof_array)
    real (prec)  ::  dNdy(length_coord_array/3,3)
    real (prec)  :: DuDx(3,3)
    real (prec)  :: F(3,3), invert_F(3,3)
    real (prec)  :: J
    real (prec)  :: S(3,length_dof_array/3)
    real (prec)  :: Pvec(length_dof_array), Svec(length_dof_array)
    real (prec)  :: Pmat(length_dof_array, length_dof_array)
    real (prec)  :: Smat(length_dof_array, length_dof_array)
    real (prec)  :: U_new(3,length_dof_array/3)
    real (prec)  :: B_finger(3,3)
    real (prec)  :: tao_stress(6)

    fail = .false.
    
    x = reshape(element_coords,(/3,length_coord_array/3/))
U_new = reshape (dof_total+dof_increment, (/3,length_dof_array/3/))


    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0
    element_stiffness = 0.d0


    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

DuDx(1:3,1:3) = matmul(U_new(1:3,1:n_nodes),dNdx(1:n_nodes,1:3))

F(1:3,1:3) = eye3_D + DuDx(1:3,1:3)

B_finger = matmul (F , transpose(F))

call invert_small(F,invert_F,J)


dNdy(1:n_nodes,1:3) = matmul (dNdx(1:n_nodes,1:3) , invert_F(1:3,1:3))

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


      G(1:6,1:9) = 0.d0
      G(1,1) = 2*B_finger(1,1)
      G(1,4) = 2*B_finger(1,2)
      G(1,6) = 2*B_finger(1,3)
      G(2,2) = 2*B_finger(2,2)
      G(2,5) = 2*B_finger(1,2)
      G(2,8) = 2*B_finger(2,3)
      G(3,3) = 2*B_finger(3,3)
      G(3,7) = 2*B_finger(1,3)
      G(3,9) = 2*B_finger(2,3)
      G(4,1) = 2*B_finger(1,2)
      G(4,2) = 2*B_finger(1,2)
      G(4,4) = 2*B_finger(2,2)
      G(4,5) = 2*B_finger(1,1)
      G(4,6) = 2*B_finger(2,3)
      G(4,8) = 2*B_finger(1,3)
      G(5,1) = 2*B_finger(1,3)
      G(5,3) = 2*B_finger(1,3)
      G(5,4) = 2*B_finger(2,3)
      G(5,6) = 2*B_finger(3,3)
      G(5,7) = 2*B_finger(1,1)
      G(5,9) = 2*B_finger(1,2)
      G(6,2) = 2*B_finger(2,3)
      G(6,3) = 2*B_finger(2,3)
      G(6,5) = 2*B_finger(1,3)
      G(6,7) = 2*B_finger(1,2)
      G(6,8) = 2*B_finger(3,3)
      G(6,9) = 2*B_finger(2,2)



 B_star(1:9,1:3*n_nodes) = 0.d0
 B_star(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
 B_star(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
 B_star(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)
 B_star(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
 B_star(5,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
 B_star(6,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
 B_star(7,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
 B_star(8,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
 B_star(9,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)

  call Neo_Hookean_material(B_finger,J,element_properties,n_properties,stress,D)

 tao_stress(1:6) = J * stress(1:6)

  S = reshape(matmul(transpose(B),tao_stress),(/3,length_dof_array/3/))
  do i= 1, n_nodes
  Pvec = reshape(spread(transpose(dNdx(i:i,1:3)),dim=2,ncopies=n_nodes),(/3*n_nodes/))
  Pmat(3*i-2:3*i,1:3*n_nodes) = spread(Pvec,dim=1,ncopies=3)
  Svec = reshape(spread(S(1:3, i:i),dim=2,ncopies=n_nodes),(/3*n_nodes/))
  Smat(3*i-2:3*i,1:3*n_nodes) = spread(Svec,dim=1,ncopies=3)
  end do

  Sigma = Pmat * transpose (Smat)

  element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),tao_stress)*w(kint)*determinant

  element_stiffness(1:3*n_nodes,1:3*n_nodes) = element_stiffness(1:3*n_nodes,1:3*n_nodes) &
            +  (matmul(transpose(B), matmul (matmul(D,G),B_star))* w(kint)*determinant &
             - Sigma(1:3*n_nodes,1:3*n_nodes)*w(kint)*determinant)

    end do
  
    return

end subroutine Neo_Hookean


!==========================SUBROUTINE el_linelast_3dbasic_dynamic ==============================
!subroutine el_linelast_3dbasic_dynamic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
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
!end subroutine el_linelast_3dbasic_dynamic


!==========================SUBROUTINE fieldvars_Neo_Hookean ==============================
subroutine fieldvars_Neo_Hookean(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords,length_coord_array, &                                ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                      ! Input variables
    n_state_variables, initial_state_variables,updated_state_variables, &                               ! Input variables
    n_field_variables,field_variable_names, &                                                           ! Field variable definition
    nodal_fieldvariables)      ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only : dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only : dNdx => shape_function_spatial_derivatives_3D
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

    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  sdev(6)                           ! Deviatoric stress
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  p, smises                          ! Pressure and Mises stress
    real (prec)  ::  J
    real (prec)  ::  U_new (3,length_coord_array/3)
    real (prec)  ::  DuDx(3,3)
    real (prec)  ::  F(3,3),invert_F(3,3)
    real (prec)  ::  B_finger(3,3)
    real (prec)  ::  dNdy(length_coord_array/3,3)
    real (prec)  ::  tao_stress(6)

    !     Subroutine to compute element contribution to project element integration point data to nodes

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    x = reshape(element_coords,(/3,length_coord_array/3/))
U_new = reshape (dof_total+dof_increment, (/3,length_dof_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    nodal_fieldvariables = 0.d0
	
  
    !     --  Loop over integration points
     do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)


DuDx(1:3,1:3) = matmul(U_new(1:3,1:n_nodes),dNdx(1:n_nodes,1:3))

F(1:3,1:3) = eye3_D + DuDx(1:3,1:3)

B_finger = matmul (F , transpose(F))

call invert_small(F,invert_F,J)


dNdy(1:n_nodes,1:3) = matmul (dNdx(1:n_nodes,1:3) , invert_F(1:3,1:3))



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

call Neo_Hookean_material(B_finger,J,element_properties,n_properties,stress,D)

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
end subroutine fieldvars_Neo_Hookean




subroutine Neo_Hookean_material(B_finger,J,element_properties,n_properties,stress,D)

   use Types
   use ParamIO

   use Element_Utilities, only : invert_small
   implicit none

   integer, intent( in )      :: n_properties
   real (prec), intent( in )  :: element_properties(n_properties)
   real (prec), intent( in )  :: J
   real (prec), intent( in )  :: B_finger(3,3)
   real (prec), intent( out ) :: stress(6)
   real (prec), intent( out ) :: D(6,6)


real (prec) :: invert_B_finger(3,3)
real (prec) :: I_vector(6),B_finger_vector(6),invert_B_finger_vector(6)
real (prec) :: D1(6,6)
real (prec) :: determinant
real (prec) :: I_dyadic_I(6,6)
real (prec) :: B_dyadic_B_inverse(6,6)
real (prec) :: I_dyadic_B_inverse(6,6)
real (prec) :: Bkk
real (prec) :: nu, Kzero


nu = element_properties(1)
Kzero  = element_properties(2)



I_vector(1:6) = 0.d0
I_vector(1:3) = 1.d0

B_finger_vector(1:6) = 0.d0
B_finger_vector(1) = B_finger(1,1)
B_finger_vector(2) = B_finger(2,2)
B_finger_vector(3) = B_finger(3,3)
B_finger_vector(4) = B_finger(1,2)
B_finger_vector(5) = B_finger(1,3)
B_finger_vector(6) = B_finger(2,3)

call invert_small(B_finger,invert_B_finger,determinant)
invert_B_finger_vector(1:6) = 0.d0
invert_B_finger_vector(1) = invert_B_finger(1,1)
invert_B_finger_vector(2) = invert_B_finger(2,2)
invert_B_finger_vector(3) = invert_B_finger(3,3)
invert_B_finger_vector(4) = invert_B_finger(1,2)
invert_B_finger_vector(5) = invert_B_finger(1,3)
invert_B_finger_vector(6) = invert_B_finger(2,3)

I_dyadic_B_inverse = spread (I_vector,dim=2,ncopies=6) * spread(invert_B_finger_vector,dim=1,ncopies=6)
I_dyadic_I = spread (I_vector,dim=2,ncopies=6) * spread(I_vector,dim=1,ncopies=6)
B_dyadic_B_inverse = spread (B_finger_vector,dim=2,ncopies=6) * spread(invert_B_finger_vector,dim=1,ncopies=6)

Bkk = B_finger_vector(1)+B_finger_vector(2)+B_finger_vector(3)


      D1 = 0.d0
      D1(1,1)= 1.d0
      D1(2,2)= 1.d0
      D1(3,3)= 1.d0
      D1(4,4)= 0.5d0
      D1(5,5)= 0.5d0
      D1(6,6)= 0.5d0

D = nu/(J**(2.d0/3.d0))*D1 + nu/(3.d0*J**(2.d0/3.d0))*((1.d0/3.d0)*Bkk*I_dyadic_B_inverse-I_dyadic_I-B_dyadic_B_inverse) &
+Kzero*J*(J-0.5d0)*I_dyadic_B_inverse

  stress(1:6) = 0.d0
  stress(1) = nu/(J**(5.d0/3.d0))*(B_finger_vector(1)-1.d0/3.d0 *Bkk) + kzero*(J-1)
  stress(2) = nu/(J**(5.d0/3.d0))*(B_finger_vector(2)-1.d0/3.d0 *Bkk) + kzero*(J-1)
  stress(3) = nu/(J**(5.d0/3.d0))*(B_finger_vector(3)-1.d0/3.d0 *Bkk) + kzero*(J-1)
  stress(4) = nu/(J**(5.d0/3.d0))* B_finger_vector(4)
  stress(5) = nu/(J**(5.d0/3.d0))* B_finger_vector(5)
  stress(6) = nu/(J**(5.d0/3.d0))* B_finger_vector(6)

  end subroutine Neo_Hookean_material

