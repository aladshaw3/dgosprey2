[GlobalParams]

    vy = 2.0
	vx = 0.0

    Dxx = 2.0
    Dyy = 2.0

    u_input = 1.0
 
	dt = 0.01
	sigma = 1   # Penalty value:  NIPG = 0   otherwise, > 0  (between 0.1 and 10)
	epsilon = 1  #  -1 = SIPG   0 = IIPG   1 = NIPG

[] #END GlobalParams

[Problem]

    coord_type = RZ

[] #END Problem

[Mesh]

    type = GeneratedMesh
    dim = 2
    nx = 5
    ny = 40
    xmin = 0.0
    xmax = 2.5
    ymin = 0.0
    ymax = 10.0

[] # END Mesh

[Variables]

    [./u]
        order = FIRST
        family = MONOMIAL
        initial_condition = 0
    [../]
 
	[./v]
		order = FIRST
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
		ldf_coef = 100.0
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
        boundary = 'top bottom'

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
    scheme = bdf2

    # NOTE: The default tolerances are far to strict and cause the program to crawl
    nl_rel_tol = 1e-10
    nl_abs_tol = 1e-4
    l_tol = 1e-8
    l_max_its = 1000
    nl_max_its = 50

    solve_type = pjfnk
    line_search = basic    # Options: default shell none basic l2 bt cp
    start_time = 0.0
    end_time = 10.0
    dtmax = 0.1

    [./TimeStepper]
        #Need to write a custom TimeStepper to enforce a maximum allowable dt
        type = SolutionTimeAdaptiveDT
    [../]

[] #END Executioner

[Preconditioning]
 
	[./smp]
		type = SMP
		full = true
		petsc_options = '-snes_converged_reason'
		petsc_options_iname = '-pc_type -ksp_gmres_restart  -snes_max_funcs'
		petsc_options_value = 'lu 2000 20000'
	[../]

[] #END Preconditioning

[Outputs]

    exodus = true
    csv = true
    print_linear_residuals = false

[] #END Outputs
