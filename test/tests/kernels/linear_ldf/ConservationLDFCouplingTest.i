[GlobalParams]

    vy = 2.0

    Dxx = 0.01
    Dyy = 0.01

    u_input = 1.0

[] #END GlobalParams

[Problem]

    coord_type = RZ

[] #END Problem

[Mesh]

    type = GeneratedMesh
    dim = 2
    nx = 3
    ny = 10
    xmin = 0.0
    xmax = 0.5
    ymin = 0.0
    ymax = 1.0

[] # END Mesh

[Variables]

    [./u]
        order = SECOND
        family = MONOMIAL
        initial_condition = 0
    [../]
 
	[./v]
		order = SECOND
		family = MONOMIAL
		initial_condition = 0
	[../]


[] #END Variables

[AuxVariables]


[] #END AuxVariables

[ICs]


[] #END ICs

[Kernels]

    [./u_dot]
        type = CoefTimeDerivative
        variable = u
        Coefficient = 0.5
    [../]

    [./u_gadv]
        type = GAdvection
        variable = u

    [../]

    [./u_gdiff]
        type = GAnisotropicDiffusion
        variable = u
    [../]
 
	[./coupled_time]
		type = CoupledCoeffTimeDerivative
		variable = u
		coupled = v
		time_coeff = 0.5
	[../]
 
	[./v_dot]
		type = CoefTimeDerivative
		variable = v
		Coefficient = 1.0
	[../]
 
	[./v_ldf]
		type = CoupledLinearLDF
		variable = v
		coupled = u
		ldf_coef = 1.0
		linear_coef = 1.0
	[../]

[] #END Kernels

[DGKernels]

    [./u_dgadv]
        type = DGAdvection
        variable = u
    [../]

    [./u_dgdiff]
        type = DGAnisotropicDiffusion
        variable = u
    [../]

[] #END DGKernels

[AuxKernels]


[] #END AuxKernels

[BCs]

    [./u_Flux]
        type = DGFluxBC
        variable = u
        boundary = 'top bottom left right'

    [../]


[] #END BCs

[Materials]


[] #END Materials

[Postprocessors]

    [./u_exit]
        type = SideAverageValue
        boundary = 'top'
        variable = u
        execute_on = 'initial timestep_end'
    [../]

    [./u_enter]
        type = SideAverageValue
        boundary = 'bottom'
        variable = u
        execute_on = 'initial timestep_end'
    [../]

    [./u_avg]
        type = ElementAverageValue
        variable = u
        execute_on = 'initial timestep_end'
    [../]
 
	[./v_avg]
		type = ElementAverageValue
		variable = v
		execute_on = 'initial timestep_end'
	[../]

[] #END Postprocessors

[Executioner]

    type = Transient
    scheme = implicit-euler

    # NOTE: The default tolerances are far to strict and cause the program to crawl
    nl_rel_tol = 1e-6
    nl_abs_tol = 1e-6
    nl_rel_step_tol = 1e-10
    nl_abs_step_tol = 1e-10
    l_tol = 1e-6
    l_max_its = 100
    nl_max_its = 100

    solve_type = newton
    line_search = bt    # Options: default shell none basic l2 bt cp
    start_time = 0.0
    end_time = 0.5
    dtmax = 0.1
    petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
    petsc_options_value = 'hypre boomeramg 100'

    [./TimeStepper]
        #Need to write a custom TimeStepper to enforce a maximum allowable dt
        type = ConstantDT
        dt = 0.05
    [../]

[] #END Executioner

[Preconditioning]

[] #END Preconditioning

[Outputs]

    exodus = true
    csv = true
    print_linear_residuals = false

[] #END Outputs
