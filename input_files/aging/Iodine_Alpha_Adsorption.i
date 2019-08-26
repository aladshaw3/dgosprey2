[GlobalParams]

    length = 1.0
    pellet_diameter = 0.16
    inner_diameter = 3.0
    flow_rate = 424100.0
    dt = 0.25    #NOTE: sometimes you need to increase dt for convergence
    sigma = 1   # Penalty value:  NIPG = 0   otherwise, > 0  (between 0.1 and 10)
    epsilon = 1  #  -1 = SIPG   0 = IIPG   1 = NIPG

[] #END GlobalParams

[Problem]

    coord_type = RZ

[] #END Problem

[Mesh]

    type = GeneratedMesh
    dim = 2
    nx = 1
    ny = 1
    xmin = 0.0
    xmax = 1.5 #cm
    ymin = 0.0
    ymax = 1.0 #cm

[] # END Mesh

[Variables]

    [./N2]
        order = FIRST
        family = MONOMIAL
    [../]

    [./O2]
        order = FIRST
        family = MONOMIAL
    [../]

    [./I2]
        order = FIRST
        family = MONOMIAL
    [../]

    [./column_temp]
        order = FIRST
        family = MONOMIAL
        initial_condition = 423.15
    [../]

    [./I2_Alpha_Adsorbed]
        order = FIRST
        family = MONOMIAL
        initial_condition = 0.0
    [../]

    [./I2_Alpha_AdsorbedHeat]
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
        initial_mole_frac = 0.79
        initial_press = 101.35
        initial_temp = 423.15
    [../]

    [./O2_IC]
        type = ConcentrationIC
        variable = O2
        initial_mole_frac = 0.21
        initial_press = 101.35
        initial_temp = 423.15
    [../]

    [./I2_IC]
        type = ConcentrationIC
        variable = I2
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

    [./accumO2]
        type = BedMassAccumulation
        variable = O2
    [../]

    [./diffO2]
        type = GColumnMassDispersion
        variable = O2
        index = 1
    [../]

    [./advO2]
        type = GColumnMassAdvection
        variable = O2
    [../]

    [./accumI2]
        type = BedMassAccumulation
        variable = I2
    [../]

    [./I2_MT]
        type = SolidMassTransfer
        variable = I2
        coupled = I2_Alpha_Adsorbed
    [../]

    [./diffI2]
        type = GColumnMassDispersion
        variable = I2
        index = 2
    [../]

    [./advI2]
        type = GColumnMassAdvection
        variable = I2
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

    [./I2_columnAdsHeat]
        type = SolidHeatTransfer
        variable = column_temp
        coupled = I2_Alpha_AdsorbedHeat
    [../]

    [./I2_adsheat]
        type = HeatofAdsorption
        variable = I2_Alpha_AdsorbedHeat
        coupled = I2_Alpha_Adsorbed
        index = 2
    [../]

    [./I2_alpha_adsorption]
        type = CoupledExtendedLangmuirFunction
        variable = I2_Alpha_Adsorbed
        main_coupled = I2_Alpha_Adsorption
        coupled_list = 'I2'
        langmuir_coeff = '12300.62'
        max_capacity = 0.74803
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

    [./dg_disp_O2]
        type = DGColumnMassDispersion
        variable = O2
        index = 1
    [../]

    [./dg_adv_O2]
        type = DGColumnMassAdvection
        variable = O2
    [../]

    [./dg_disp_I2]
        type = DGColumnMassDispersion
        variable = I2
        index = 2
    [../]

    [./dg_adv_I2]
        type = DGColumnMassAdvection
        variable = I2
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
        coupled_gases = 'N2 O2 I2'
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
        input_molefraction = 0.78985
        index = 0
    [../]

    [./O2_Flux]
        type = DGMassFluxBC
        variable = O2
        boundary = 'top bottom'
        input_temperature = 423.15
        input_pressure = 101.35
        input_molefraction = 0.2101
        index = 1
    [../]

    [./I2_Flux]
        type = DGMassFluxBC
        variable = I2
        boundary = 'top bottom'
        input_temperature = 423.15
        input_pressure = 101.35
        input_molefraction = 0.00005
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
        block = 0
        outer_diameter = 3.30
        bulk_porosity = 0.7255
        wall_density = 8.0
        wall_heat_capacity = 0.5
        wall_heat_trans_coef = 6.12
        extern_heat_trans_coef = 6.12
    [../]

    [./FlowMaterials]
        type = GasFlowProperties
        block = 0
        molecular_weight = '28.016 32 253.8'
        comp_heat_capacity = '1.04 0.919 0.214'
        comp_ref_viscosity = '0.0001781 0.0002018 0.00013283'
        comp_ref_temp = '300.55 292.25 295.496'
        comp_Sutherland_const = '111 127 573.474'
        temperature = column_temp
        total_pressure = total_pressure
        coupled_gases = 'N2 O2 I2'
    [../]

    [./AdsorbentMaterials]
        type = AdsorbentProperties
        block = 0
        binder_fraction = 0.0
        binder_porosity = 0.384
        crystal_radius = 0.0
        macropore_radius = 2.65e-6
        pellet_density = 3.06
        pellet_heat_capacity = 1.2
        ref_diffusion = '0 0 0'
        activation_energy = '0 0 0'
        ref_temperature = '0 0 0'
        affinity = '0 0 0'
        temperature = column_temp
        coupled_gases = 'N2 O2 I2'
    [../]

    [./AdsorbateMaterials]
        type = ThermodynamicProperties
        block = 0
        temperature = column_temp
        total_pressure = total_pressure
        coupled_gases = 'N2 O2 I2'
        number_sites = '0 0 1'
        maximum_capacity = '0 0 0.55' #mol/kg
        molar_volume = '0 0 41.8' #cm^3/mol
        enthalpy_site_1 = '0 0 -8308'
        enthalpy_site_2 = '0 0 0'
        enthalpy_site_3 = '0 0 0'
        enthalpy_site_4 = '0 0 0'
        enthalpy_site_5 = '0 0 0'
        enthalpy_site_6 = '0 0 0'

        entropy_site_1 = '0 0 92.35'
        entropy_site_2 = '0 0 0'
        entropy_site_3 = '0 0 0'
        entropy_site_4 = '0 0 0'
        entropy_site_5 = '0 0 0'
        entropy_site_6 = '0 0 0'
    [../]

[] #END Materials

[Postprocessors]

    [./O2_enter]
        type = SideAverageValue
        boundary = 'bottom'
        variable = O2
        execute_on = 'initial timestep_end'
    [../]

    [./O2_avg_gas]
        type = ElementAverageValue
        variable = O2
        execute_on = 'initial timestep_end'
    [../]

    [./I2_enter]
        type = SideAverageValue
        boundary = 'bottom'
        variable = I2
        execute_on = 'initial timestep_end'
    [../]

    [./I2_avg_gas]
        type = ElementAverageValue
        variable = I2
        execute_on = 'initial timestep_end'
    [../]

    [./I2_exit]
        type = SideAverageValue
        boundary = 'top'
        variable = I2
        execute_on = 'initial timestep_end'
    [../]

    [./press_exit]
        type = SideAverageValue
        boundary = 'top'
        variable = total_pressure
        execute_on = 'initial timestep_end'
    [../]

    [./Alpha_AgI_solid]
        type = ElementAverageValue
        variable = I2_Alpha_Adsorbed
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
    l_max_its = 100
    nl_max_its = 50

    solve_type = pjfnk
    line_search = basic    # Options: default none l2 bt basic
    start_time = 0.0
    end_time = 400.0
    dtmax = 1.0

    [./TimeStepper]
        type = ConstantDT
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

