 [GlobalParams]

    length = 0.1
    pellet_diameter = 0.18
    inner_diameter = 2.5
    flow_rate = 30000
    dt = 0.01   #NOTE: sometimes you need to increase dt for convergence
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
    xmax = 1.25 #cm
    ymin = 0.0
    ymax = 0.1 #cm

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

    [./NO2]
        order = FIRST
        family = MONOMIAL
    [../]

#    [./H2O]
#        order = FIRST
#        family = MONOMIAL
#        initial_condition = 0.0
#    [../]

    [./column_temp]
        order = FIRST
        family = MONOMIAL
        initial_condition = 423.15
    [../]

#    [./Ag_MOR]
#        order = FIRST
#        family = MONOMIAL
#        initial_condition = 0.0
#    [../]
#
#    [./H_MOR]
#        order = FIRST
#        family = MONOMIAL
#        initial_condition = 0.799344
#    [../]

    [./Ag_NO3]
        order = FIRST
        family = MONOMIAL
        initial_condition = 0.0
    [../]

#    [./Ag0]
#        order = FIRST
#        family = MONOMIAL
#        initial_condition = 1.1102
#    [../]


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

    [./NO2_IC]
        type = ConcentrationIC
        variable = NO2
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

    [./accumNO2]
        type = BedMassAccumulation
        variable = NO2
    [../]

    [./diffNO2]
        type = GColumnMassDispersion
        variable = NO2
        index = 2
    [../]

    [./advNO2]
        type = GColumnMassAdvection
        variable = NO2
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

#    [./AgNO3_MT]
#        type = CoefTimeDerivative
#        variable = Ag_NO3
#    [../]
#
#    [./Ag0_MT]
#        type = CoefTimeDerivative
#        variable = Ag0
#    [../]

    [./NO2_Ag_aging]
        type = VariableOrderCoupledCatalyst
        variable = Ag_NO3
        main_variable = Ag_NO3
        coupled_gases = 'NO2 O2'
        gases_stoichiometry = '-2 -1'
        gas_order = '2 1'
        coupled_catalysts = ''
        catalyst_stoichiometry = ''
        catalyst_order = ''
        coupled_adsorption = 'Ag_NO3'
        adsorbed_sites = '2'
        site_order = '2'
        adsorbed_stoichiometry = '2'
        ads_order = '2'
        max_capacity = 1.1102
        forward_rate = 2.3e10
        reverse_rate = 0.0
    [../]

#    [./Ag_Aging_1]
#        type = VariableOrderReac
#        variable = Ag_NO3
#        main_variable = Ag_NO3
#        coupled_species = 'Ag0 NO2 O2'
#        stoichiometry = '-2 -2 -1'
#        order = '2 2 1'
#        main_stoichiometry = 2
#        main_order = 2
#        forward_rate = 1.0e9
#        reverse_rate = 0.0
#    [../]
#
#    [./Ag_Aging_2]
#        type = VariableOrderReac
#        variable = Ag0
#        main_variable = Ag0
#        coupled_species = 'Ag_NO3 NO2 O2'
#        stoichiometry = '-2 2 1'
#        order = '2 2 1'
#        main_stoichiometry = 2
#        main_order = 2
#        forward_rate = 0.0
#        reverse_rate = 1.0e9
#    [../]
#
#    [./Ag_Aging_to_MOR_1]
#        type = VariableOrderReac
#        variable = Ag_MOR
#        main_variable = Ag_MOR
#        coupled_species = 'Ag0 H_MOR O2 H2O'
#        stoichiometry = '-4 -4 -1 2'
#        order = '4 4 1 2'
#        main_stoichiometry = 4
#        main_order = 4
#        forward_rate = 0.0
#        reverse_rate = 0.0
#    [../]
#
#    [./Ag_Aging_to_MOR_2]
#        type = VariableOrderReac
#        variable = Ag0
#        main_variable = Ag0
#        coupled_species = 'H_MOR O2 Ag_MOR H2O'
#        stoichiometry = '4 1 -4 -2'
#        order = '4 4 1 2'
#        main_stoichiometry = 4
#        main_order = 4
#        forward_rate = 0.0
#        reverse_rate = 0.0
#    [../]
#
#    [./Ag_Aging_to_MOR_3]
#        type = VariableOrderReac
#        variable = H_MOR
#        main_variable = H_MOR
#        coupled_species = 'Ag0 O2 Ag_MOR H2O'
#        stoichiometry = '4 1 -4 -2'
#        order = '4 4 1 2'
#        main_stoichiometry = 4
#        main_order = 4
#        forward_rate = 0.0
#        reverse_rate = 0.0
#    [../]
#
#    [./Ag_Aging_to_MOR_4]
#        type = VariableOrderReac
#        variable = H2O
#        main_variable = H2O
#        coupled_species = 'Ag0 H_MOR O2 Ag_MOR'
#        stoichiometry = '4 4 1 -4'
#        order = '4 4 1 4'
#        main_stoichiometry = 2
#        main_order = 2
#        forward_rate = 0.0
#        reverse_rate = 0.0
#    [../]

    [./O2_MT_1]
        type = SolidMassTransfer
        variable = O2
        coupled = Ag_NO3
    [../]

#    [./O2_MT_2]
#        type = SolidMassTransfer
#        variable = O2
#        coupled = Ag_MOR
#    [../]

    [./NO2_MT]
        type = SolidMassTransfer
        variable = NO2
        coupled = Ag_NO3
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

    [./dg_disp_NO2]
        type = DGColumnMassDispersion
        variable = NO2
        index = 2
    [../]

    [./dg_adv_NO2]
        type = DGColumnMassAdvection
        variable = NO2
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
        coupled_gases = 'N2 O2 NO2'
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
        input_molefraction = 0.77
        index = 0
    [../]

    [./O2_Flux]
        type = DGMassFluxBC
        variable = O2
        boundary = 'top bottom'
        input_temperature = 423.15
        input_pressure = 101.35
        input_molefraction = 0.21
        index = 1
    [../]

    [./NO2_Flux]
        type = DGMassFluxBC
        variable = NO2
        boundary = 'top bottom'
        input_temperature = 423.15
        input_pressure = 101.35
        input_molefraction = 0.02
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
        outer_diameter = 2.8
        bulk_porosity = 0.3336
        wall_density = 8.0
        wall_heat_capacity = 0.5
        wall_heat_trans_coef = 6.12
        extern_heat_trans_coef = 6.12
    [../]

    [./FlowMaterials]
        type = GasFlowProperties
        block = 0
        molecular_weight = '28.016 32.0 46.0'
        comp_heat_capacity = '1.04 0.919 0.8066'
        comp_ref_viscosity = '0.0001781 0.0002018 0.00001963'
        comp_ref_temp = '300.55 292.25 400.0'
        comp_Sutherland_const = '111 127 701.226'
        temperature = column_temp
        total_pressure = total_pressure
        coupled_gases = 'N2 O2 NO2'
    [../]

    [./AdsorbentMaterials]
        type = AdsorbentProperties
        block = 0
        binder_fraction = 0.0
        binder_porosity = 0.384
        crystal_radius = 0.0
        macropore_radius = 2.65e-6
        pellet_density = 3.057
        pellet_heat_capacity = 1.2
        ref_diffusion = '0 0 0'
        activation_energy = '0 0 0'
        ref_temperature = '0 0 0'
        affinity = '0 0 0'
        temperature = column_temp
        coupled_gases = 'N2 O2 NO2'
    [../]

    [./AdsorbateMaterials]
        type = ThermodynamicProperties
        block = 0
        temperature = column_temp
        total_pressure = total_pressure
        coupled_gases = 'N2 O2 NO2'
        number_sites = '0 0 2'
        maximum_capacity = '0 0 4.75' #mol/kg
        molar_volume = '0 0 13500.0' #cm^3/mol
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

    [./Ag_NO3_Sol]
        type = ElementAverageValue
        variable = Ag_NO3
        execute_on = 'initial timestep_end'
    [../]

#    [./Ag0_Sol]
#        type = ElementAverageValue
#        variable = Ag0
#        execute_on = 'initial timestep_end'
#    [../]

    [./NO2_ave]
        type = ElementAverageValue
        variable = NO2
        execute_on = 'initial timestep_end'
    [../]

#    [./H2O_column]
#        type = ElementAverageValue
#        variable = H2O
#        execute_on = 'initial timestep_end'
#    [../]

[] #END Postprocessors

[Executioner]

    type = Transient
    scheme = bdf2

# NOTE: The default tolerances are far to strict and cause the program to crawl
    nl_rel_tol = 1e-10
    nl_abs_tol = 1e-4
    l_tol = 1e-8
    l_max_its = 200
    nl_max_its = 80

    solve_type = pjfnk
    line_search = bt    # Options: default none l2 bt basic
    start_time = 0.0
    end_time = 1.0
    dtmax = 0.1

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


