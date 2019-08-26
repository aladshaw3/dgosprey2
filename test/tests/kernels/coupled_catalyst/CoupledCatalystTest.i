[GlobalParams]

    vy = 2.0

    Dxx = 0.01
    Dyy = 0.01

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

    [./w]
        order = SECOND
        family = MONOMIAL
        initial_condition = 0
    [../]

    [./x]
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
        time_coeff = 1.25
    [../]

    [./v_rate]
        type = CoupledConstChemisorption
        variable = v
        main_variable = v
        coupled_gases = 'u'
        gases_stoichiometry = '-1'
        coupled_adsorption = 'v x'
        adsorbed_sites = '2 4'
        adsorbed_stoichiometry = '2 4'
        max_capacity = 5.0
        forward_rate = 0.5
        reverse_rate = 0.0
    [../]

    [./x_rate]
        type = CoupledCatalyst
        variable = x
        main_variable = x
        coupled_gases = ''
        gases_stoichiometry = ''
        coupled_catalysts = 'w'
        catalyst_stoichiometry = '1'
        coupled_adsorption = 'v x'
        adsorbed_sites = '2 4'
        adsorbed_stoichiometry = '2 4'
        max_capacity = 5.0
        forward_rate = 0.5
        reverse_rate = 0.005
    [../]

    [./w_dot]
        type = CoefTimeDerivative
        variable = w
        Coefficient = 0.5
    [../]

    [./w_gadv]
        type = GAdvection
        variable = w
    [../]

    [./w_gdiff]
        type = GAnisotropicDiffusion
        variable = w
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

    [./w_dgadv]
        type = DGAdvection
        variable = w
    [../]

    [./w_dgdiff]
        type = DGAnisotropicDiffusion
        variable = w
    [../]

[] #END DGKernels

[AuxKernels]


[] #END AuxKernels

[BCs]

    [./u_Flux]
        type = DGFluxBC
        variable = u
        u_input = 1.0
        boundary = 'top bottom left right'
    [../]

    [./w_Flux]
        type = DGFluxBC
        variable = w
        u_input = 1.0
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

    [./w_exit]
        type = SideAverageValue
        boundary = 'top'
        variable = w
        execute_on = 'initial timestep_end'
    [../]

    [./w_enter]
        type = SideAverageValue
        boundary = 'bottom'
        variable = w
        execute_on = 'initial timestep_end'
    [../]

    [./w_avg]
        type = ElementAverageValue
        variable = w
        execute_on = 'initial timestep_end'
    [../]

    [./v_avg]
        type = ElementAverageValue
        variable = v
        execute_on = 'initial timestep_end'
    [../]

    [./x_avg]
        type = ElementAverageValue
        variable = x
        execute_on = 'initial timestep_end'
    [../]

[] #END Postprocessors

[Executioner]

    type = Transient
    scheme = implicit-euler

    # NOTE: The default tolerances are far to strict and cause the program to crawl
    nl_rel_tol = 1e-6
    nl_abs_tol = 1e-5
    l_tol = 1e-6
    l_max_its = 100
    nl_max_its = 20

    solve_type = pjfnk
    line_search = bt    # Options: default shell none basic l2 bt cp
    start_time = 0.0
    end_time = 0.25
    dtmax = 0.1

    [./TimeStepper]
        #Need to write a custom TimeStepper to enforce a maximum allowable dt
        type = ConstantDT
        dt = 0.05
    [../]

[] #END Executioner

[Preconditioning]

    [./smp]
        type = SMP
        full = true
        petsc_options = '-snes_converged_reason'
        petsc_options_iname = '-pc_type -sub_pc_type -pc_hypre_type -ksp_gmres_restart  -snes_max_funcs'
        petsc_options_value = 'lu ilu boomeramg 2000 20000'
    [../]

[] #END Preconditioning

[Outputs]

    exodus = true
    csv = true
    print_linear_residuals = false

[] #END Outputs
