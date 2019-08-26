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
    ny = 41
    xmin = 0.0
    xmax = 2.5
    ymin = 0.0
    ymax = 10.1

[] # END Mesh
 
 [MeshModifiers]
 
	[./block0]
		type = SubdomainBoundingBox
		block_id = 0
		top_right = '2.5 5.0 0'
		bottom_left = '0 0 0'
	[../]
 
	[./block1]
		type = SubdomainBoundingBox
		block_id = 1
		top_right = '2.5 5.1 0'
		bottom_left = '0 5.0 0'
	[../]
 
	[./block2]
		type = SubdomainBoundingBox
		block_id = 2
		top_right = '2.5 10.1 0'
		bottom_left = '0 5.1 0'
	[../]
 
	[./interface01]
		type = SideSetsBetweenSubdomains
		depends_on = 'block0 block1'
		master_block = 0
		paired_block = 1
		new_boundary = 'interface01'
	[../]
 
	[./interface12]
		type = SideSetsBetweenSubdomains
		depends_on = 'block1 block2'
		master_block = 1
		paired_block = 2
		new_boundary = 'interface12'
	[../]

[Variables]

    [./u]
		block = '0 1 2'
        order = FIRST
        family = MONOMIAL
        initial_condition = 0
    [../]
 
	[./v]
 		block = '0 2'
		order = CONSTANT
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
 		block = '0 2'
        type = CoefTimeDerivative
        variable = u
        Coefficient = 0.5
    [../]

    [./u_gadv]
 		block = '0 2'
        type = GAdvection
        variable = u

    [../]
 
	[./u_dot_cond]
		block = '1'
		type = CoefTimeDerivative
		variable = u
		Coefficient = 0.5
	[../]
 
	[./u_gadv_cond]
		block = '1'
		type = GAdvection
		variable = u
 
	[../]

    [./u_gdiff]
 		block = '0 2'
        type = GAnisotropicDiffusion
        variable = u
    [../]
 
	[./coupled_time]
 		block = '0 2'
		type = CoupledCoeffTimeDerivative
		variable = u
		coupled = v
		time_coeff = 0.5
	[../]
 
	[./v_dot]
 		block = '0 2'
		type = CoefTimeDerivative
		variable = v
		Coefficient = 1.0
	[../]
 
	[./v_ldf]
 		block = '0 2'
		type = CoupledLinearLDF
		variable = v
		coupled = u
		ldf_coef = 100.0
		linear_coef = 1.0
	[../]

[] #END Kernels

[DGKernels]

    [./u_dgadv]
 		block = '0 1 2'
        type = DGAdvection
        variable = u
    [../]

    [./u_dgdiff]
 		block = '0 2'
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

    [./u_exit_2]
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
 
	[./u_exit_0]
		type = SideAverageValue
		boundary = 'interface01'
		variable = u
		execute_on = 'initial timestep_end'
	[../]
 
	[./u_exit_1]
		type = SideAverageValue
		boundary = 'interface12'
		variable = u
		execute_on = 'initial timestep_end'
	[../]

    [./u_avg_0]
 		block = '0'
        type = ElementAverageValue
        variable = u
        execute_on = 'initial timestep_end'
    [../]
 
	[./v_avg_0]
 		block = '0'
		type = ElementAverageValue
		variable = v
		execute_on = 'initial timestep_end'
	[../]
 
	[./u_avg_1]
		block = '1'
		type = ElementAverageValue
		variable = u
		execute_on = 'initial timestep_end'
	[../]
 
	[./u_avg_2]
		block = '2'
		type = ElementAverageValue
		variable = u
		execute_on = 'initial timestep_end'
	[../]
 
	[./v_avg_2]
		block = '2'
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
