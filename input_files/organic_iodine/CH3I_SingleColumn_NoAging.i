[GlobalParams]

    length = 15.0
    pellet_diameter = 0.18
    inner_diameter = 2.0
    flow_rate = 30000
    dt = 0.1   #NOTE: sometimes you need to increase dt for convergence
    sigma = 1   # Penalty value:  NIPG = 0   otherwise, > 0  (between 0.1 and 10)
    epsilon = 1  #  -1 = SIPG   0 = IIPG   1 = NIPG

 [] #END GlobalParams

[Problem]

    coord_type = RZ

 [] #END Problem

[Mesh]

    type = GeneratedMesh
    dim = 2
    nx = 4
    ny = 60
    xmin = 0.0
    xmax = 1.0 #cm
    ymin = 0.0
    ymax = 15.0 #cm

 [] # END Mesh

[Variables]

    [./N2]
        order = FIRST
        family = MONOMIAL
    [../]

    [./CH3I]
        order = FIRST
        family = MONOMIAL
    [../]

    [./C2H6]
        order = FIRST
        family = MONOMIAL
    [../]

    [./column_temp]
        order = FIRST
        family = MONOMIAL
        initial_condition = 423.15
    [../]

    [./AgI]
        order = FIRST
        family = MONOMIAL
        initial_condition = 0.0
    [../]

 [] #END Variables

[AuxVariables]

    [./total_pressure]
        order = CONSTANT
        family = MONOMIAL
        initial_condition = 101.35
    [../]

    [./ambient_temp]
         order = CONSTANT
        family = MONOMIAL
        initial_condition = 423.15
    [../]

    [./wall_temp]
        order = FIRST
        family = MONOMIAL
        initial_condition = 423.15
    [../]

 [] #END AuxVariables

[ICs]

    [./N2_IC]
        type = ConcentrationIC
        variable = N2
        initial_mole_frac = 1.0
        initial_press = 101.35
        initial_temp = 423.15
    [../]

    [./CH3I_IC]
        type = ConcentrationIC
        variable = CH3I
        initial_mole_frac = 0.0
        initial_press = 101.35
        initial_temp = 423.15
    [../]

    [./C2H6_IC]
        type = ConcentrationIC
        variable = C2H6
        initial_mole_frac = 0.0
        initial_press = 101.35
        initial_temp = 423.15
    [../]

 [] #END ICs

[Kernels]

    [./accumN2]
        type = BedMassAccumulation
        variable = N2
    [../]

    [./diffN2]
        type = GColumnMassDispersion
        variable = N2
        index = 0
    [../]

    [./advN2]
        type = GColumnMassAdvection
        variable = N2
    [../]

    [./accumCH3I]
        type = BedMassAccumulation
        variable = CH3I
    [../]

    [./diffCH3I]
        type = GColumnMassDispersion
        variable = CH3I
        index = 1
    [../]

    [./advCH3I]
        type = GColumnMassAdvection
        variable = CH3I
    [../]

    [./accumC2H6]
        type = BedMassAccumulation
        variable = C2H6
    [../]

    [./diffC2H6]
        type = GColumnMassDispersion
        variable = C2H6
        index = 2
    [../]

    [./advC2H6]
        type = GColumnMassAdvection
        variable = C2H6
    [../]

    [./columnAccum]
        type = BedHeatAccumulation
        variable = column_temp
    [../]

    [./columnConduction]
        type = GColumnHeatDispersion
        variable =column_temp
    [../]

    [./columnAdvection]
        type = GColumnHeatAdvection
        variable =column_temp
    [../]
 
#    [./Aerogel_Aging]
#        type = VariableOrderCoupledCatalyst
#        block = '0 2 4'
#        variable = rSAg
#        main_variable = rSAg
#        coupled_gases = 'rSH H2'
#        gases_stoichiometry = '-2 1'
#        gas_order = '2 1'
#        coupled_catalysts = ''
#        catalyst_stoichiometry = ''
#        catalyst_order = ''
#        coupled_adsorption = 'AgI rSAg'
#        adsorbed_sites = '2 2'
#        site_order = '2 2'
#        adsorbed_stoichiometry = '2 2'
#        ads_order = '2 2'
#        max_capacity = 1.1102
#        forward_rate = 0.02
#        reverse_rate = 0.0
#    [../]
#
#    [./rSH_Adsorption]
#        block = '0 2 4'
#        type = CoupledCoeffTimeDerivative
#        variable = rSH
#        coupled = rSAg
#    [../]
#
#    [./H2_Aging_1]
#        block = '0 4'
#        type = CoupledCoeffTimeDerivative
#        variable = H2
#        coupled = rSH
#        time_coeff = 0.15875
#    [../]
#
#    [./H2_Aging_2]
#        block = '2'
#        type = CoupledCoeffTimeDerivative
#        variable = H2
#        coupled = rSH
#        time_coeff = 0.17875
#    [../]

    [./Aerogel_Adsorption]
        type = VariableOrderCoupledCatalyst
        variable = AgI
        main_variable = AgI
        coupled_gases = 'CH3I C2H6'
        gases_stoichiometry = '-2 1'
        gas_order = '2 1'
        coupled_catalysts = ''
        catalyst_stoichiometry = ''
        catalyst_order = ''
        coupled_adsorption = 'AgI'
        adsorbed_sites = '2'
        site_order = '3'
        adsorbed_stoichiometry = '2'
        ads_order = '2'
        max_capacity = 3.28125
        forward_rate = 6.5e7
        reverse_rate = 0.0
    [../]

    [./CH3I_Adsorption]
        type = SolidMassTransfer
        variable = CH3I
        coupled = AgI
    [../]

    [./C2H6_Adsorption]
        type = CoupledCoeffTimeDerivative
        variable = C2H6
        coupled = AgI
        time_coeff = -0.15875
    [../]
 
 
 [] #END Kernels

[DGKernels]

    [./dg_disp_N2]
        type = DGColumnMassDispersion
        variable = N2
        index = 0
    [../]

    [./dg_adv_N2]
        type = DGColumnMassAdvection
        variable = N2
    [../]

    [./dg_disp_CH3I]
        type = DGColumnMassDispersion
        variable = CH3I
        index = 1
    [../]

    [./dg_adv_CH3I]
        type = DGColumnMassAdvection
        variable = CH3I
    [../]

    [./dg_disp_C2H6]
        type = DGColumnMassDispersion
        variable = C2H6
        index = 2
    [../]

    [./dg_adv_C2H6]
        type = DGColumnMassAdvection
        variable = C2H6
    [../]

    [./dg_disp_heat]
        type = DGColumnHeatDispersion
        variable = column_temp
    [../]

    [./dg_adv_heat]
        type = DGColumnHeatAdvection
        variable = column_temp
    [../]

 [] #END DGKernels

[AuxKernels]

    [./column_pressure]
        type = TotalColumnPressure
        variable = total_pressure
        temperature = column_temp
        coupled_gases = 'N2 CH3I C2H6'
        execute_on = 'initial timestep_end'
    [../]

    [./wall_temp_calc]
        type = WallTemperature
        variable = wall_temp
        column_temp = column_temp
        ambient_temp = ambient_temp
        execute_on = 'initial timestep_end'
    [../]

 [] #END AuxKernels

[BCs]

    [./N2_Flux]
        type = DGMassFluxBC
        variable = N2
        boundary = 'top bottom'
        input_temperature = 423.15
        input_pressure = 101.35
        input_molefraction = 0.99999
        index = 0
    [../]

    [./CH3I_Flux]
        type = DGMassFluxBC
        variable = CH3I
        boundary = 'top bottom'
        input_temperature = 423.15
        input_pressure = 101.35
        input_molefraction = 0.00001
        index = 1
    [../]

    [./C2H6_Flux]
        type = DGMassFluxBC
        variable = C2H6
        boundary = 'top bottom'
        input_temperature = 423.15
        input_pressure = 101.35
        input_molefraction = 0.0
        index = 2
    [../]

    [./Heat_Gas_Flux]
        type = DGHeatFluxBC
        variable = column_temp
        boundary = 'top bottom'
        input_temperature = 423.15
    [../]

    [./Heat_Wall_Flux]
        type = DGColumnWallHeatFluxLimitedBC
        variable = column_temp
        boundary = 'right left'
        wall_temp = wall_temp
    [../]
 
 [] #END BCs

[Materials]

    [./BedMaterials]
        type = BedProperties
        outer_diameter = 2.3
        bulk_porosity = 0.365
        wall_density = 8.0
        wall_heat_capacity = 0.5
        wall_heat_trans_coef = 6.12
        extern_heat_trans_coef = 6.12
    [../]

    [./FlowMaterials]
        type = GasFlowProperties
        molecular_weight = '28.016 141.939 30.07'
        comp_heat_capacity = '1.04 2.232 1.28'
        comp_ref_viscosity = '0.0001781 0.00011 0.00008817'
        comp_ref_temp = '300.55 293.15 280'
        comp_Sutherland_const = '111 93 213'
        temperature = column_temp
        total_pressure = total_pressure
        coupled_gases = 'N2 CH3I C2H6'
    [../]

    [./AdsorbentMaterials]
        type = AdsorbentProperties
        binder_fraction = 0.0
        binder_porosity = 0.384
        crystal_radius = 0.0
        macropore_radius = 2.65e-6
        pellet_density = 0.5
        pellet_heat_capacity = 1.2
        ref_diffusion = '0 0 0'
        activation_energy = '0 0 0'
        ref_temperature = '0 0 0'
        affinity = '0 0 0'
        temperature = column_temp
        coupled_gases = 'N2 CH3I C2H6'
    [../]

    [./AdsorbateMaterials]
        type = ThermodynamicProperties
        temperature = column_temp
        total_pressure = total_pressure
        coupled_gases = 'N2 CH3I C2H6'
        number_sites = '0 0 0'
        maximum_capacity = '0 0 0' #mol/kg
        molar_volume = '0 0 0' #cm^3/mol
        enthalpy_site_1 = '0 0 0'
        enthalpy_site_2 = '0 0 0'
        enthalpy_site_3 = '0 0 0'
        enthalpy_site_4 = '0 0 0'
        enthalpy_site_5 = '0 0 0'
        enthalpy_site_6 = '0 0 0'

        entropy_site_1 = '0 0 0'
        entropy_site_2 = '0 0 0'
        entropy_site_3 = '0 0 0'
        entropy_site_4 = '0 0 0'
        entropy_site_5 = '0 0 0'
        entropy_site_6 = '0 0 0'
    [../]

 [] #END Materials

[Postprocessors]

    [./press_exit]
        type = SideAverageValue
        boundary = 'top'
        variable = total_pressure
        execute_on = 'initial timestep_end'
    [../]

    [./AgI_solid]
        type = ElementAverageValue
        variable = AgI
        execute_on = 'initial timestep_end'
    [../]

    [./CH3I_exit]
        type = SideAverageValue
        boundary = 'top'
        variable = CH3I
        execute_on = 'initial timestep_end'
    [../]

    [./C2H6_exit]
        type = SideAverageValue
        boundary = 'top'
        variable = C2H6
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
    l_max_its = 70
    nl_max_its = 20

    solve_type = pjfnk
    line_search = bt    # Options: default none l2 bt basic
    start_time = 0.0
    end_time = 5000.0
    dtmax = 100.0

    [./TimeStepper]
        type = SolutionTimeAdaptiveDT
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

