
!==========================SUBROUTINE Cahn_Hilliard_2dbasic ==============================
subroutine Cahn_Hilliard_2dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties, element_coords, length_coord_array, &                      ! Input variables
    dof_increment, dof_total, length_dof_array, &                                                ! Input variables
    n_state_variables, initial_state_variables, &                                                ! Input variables
    updated_state_variables,element_stiffness,element_residual, fail)                            ! Output variables
    use Types
    use ParamIO
    use Globals, only: TIME,DTIME  !For a time dependent problem uncomment this line to access the time increment and total time
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_2D
    use Element_Utilities, only : dNdxi => shape_function_derivatives_2D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_2D
    use Element_Utilities, only : xi => integrationpoints_2D, w => integrationweights_2D
    use Element_Utilities, only : dxdxi => jacobian_2D
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
  
    logical, intent( out )        :: fail                                                   ! Set to .true. to force a timestep cutback
    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_stiffness(length_dof_array,length_dof_array)   ! Element stiffness (ROW,COLUMN)
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)
          

    ! Local Variables
    integer      :: n_points,kint

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, 2e12]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s12]
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(2,2), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(2,length_coord_array/2)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  D_coefficient,K_coefficient,Theta_coefficient   ! Material properties
    real (prec)  ::  nu,dnu,C,dc

    fail = .false.
    
    x = reshape(element_coords,(/2,length_coord_array/2/))

    if (n_nodes == 4) n_points = 4
    if (n_nodes == 3) n_points = 1
    if (n_nodes == 6) n_points = 4
    if (n_nodes == 8) n_points = 9

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0
    element_stiffness = 0.d0


    D_coefficient  = element_properties(1)
    K_coefficient  = element_properties(2)
    Theta_coefficient = element_properties(3)


    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:2,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:2,1:n_nodes),dNdxi(1:n_nodes,1:2))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:2) = matmul(dNdxi(1:n_nodes,1:2),dxidx)
        B = 0.d0
        B(1,1:2*n_nodes-1:2) = N(1:n_nodes)
        B(2,2:2*n_nodes:2)   = N(1:n_nodes)
        B(3,1:2*n_nodes-1:2) = dNdx(1:n_nodes,1)
        B(4,1:2*n_nodes-1:2) = dNdx(1:n_nodes,2)
        B(5,2:2*n_nodes:2)   = dNdx(1:n_nodes,1)
        B(6,2:2*n_nodes:2)   = dNdx(1:n_nodes,2)


        strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)

        C  = strain(2)
        dc = dstrain(2)
        nu = strain(1)
        Dnu = dstrain(1)
      
        D = 0.d0
        D(1,1) = 1.d0
        D(1,2) = 1.d0 - 3.d0*(C + dc)**2.d0
        D(2,2) = 1.d0/DTIME
        D(3,5) = -K_coefficient
        D(4,6) = -K_coefficient
        D(5,3) = Theta_coefficient*D_coefficient
        D(6,4) = Theta_coefficient*D_coefficient


        stress = 0.d0
        stress(1) = nu +Dnu- (C+dc)*((C+dc)**2-1)
        stress(2) = dc/ DTIME
        stress(3) = -K_coefficient*(strain(5)+dstrain(5))
        stress(4) = -K_coefficient*(strain(6)+dstrain(6))
        stress(5) =  D_coefficient*(strain(3)+Theta_coefficient*dstrain(3))
        stress(6) =  D_coefficient*(strain(4)+Theta_coefficient*dstrain(4))

       ! stress = matmul (D, strain+dstrain)

        element_residual(1:2*n_nodes) = element_residual(1:2*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant

                element_stiffness(1:2*n_nodes,1:2*n_nodes) = element_stiffness(1:2*n_nodes,1:2*n_nodes) &
                    + matmul(transpose(B(1:6,1:2*n_nodes)),matmul(D,B(1:6,1:2*n_nodes)))*w(kint)*determinant


    end do
  
    return
end subroutine Cahn_Hilliard_2dbasic

