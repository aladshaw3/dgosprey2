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
 
    [./NO]
        order = FIRST
        family = MONOMIAL
    [../]
 
    [./H2O]
        order = FIRST
        family = MONOMIAL
    [../]
 
    [./NO2]
        order = FIRST
        family = MONOMIAL
    [../]
 
    [./column_temp]
        order = FIRST
        family = MONOMIAL
        initial_condition = 423.15
    [../]
 
    [./Ag0]
        order = FIRST
        family = MONOMIAL
        initial_condition = 1.1102
    [../]
 
    [./AgMOR]
        order = FIRST
        family = MONOMIAL
        initial_condition = 0.0
    [../]

    [./HMOR]
        order = FIRST
        family = MONOMIAL
        initial_condition = 1.1102
    [../]
 
#    [./HMOR]
#        order = FIRST
#        family = MONOMIAL
#        initial_condition = 0.88816
#    [../]
 
    [./Ag2O]
        order = FIRST
        family = MONOMIAL
        initial_condition = 0.0
    [../]
 
    [./AgNO3]
        order = FIRST
        family = MONOMIAL
        initial_condition = 0.0
    [../]
 
    [./AgNO3_Res]
        order = FIRST
        family = MONOMIAL
        initial_condition = 0.0
    [../]
 
    [./X_Sites]
        order = FIRST
        family = MONOMIAL
        initial_condition = 0.24561947
    [../]
 
#    [./X_Sites]
#        order = FIRST
#        family = MONOMIAL
#        initial_condition = 0.22204
#    [../]
 
    [./FeMOR]
        order = FIRST
        family = MONOMIAL
        initial_condition = 0.1764377
    [../]
 
    [./FeNO]
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
        initial_mole_frac = 0.99
        initial_press = 101.35
        initial_temp = 423.15
    [../]
 
    [./NO_IC]
        type = ConcentrationIC
        variable = NO
        initial_mole_frac = 0.01
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

    [./H2O_IC]
       type = ConcentrationIC
       variable = H2O
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
 
    [./accumNO]
        type = BedMassAccumulation
        variable = NO
    [../]
 
    [./diffNO]
        type = GColumnMassDispersion
        variable = NO
        index = 1
    [../]
 
    [./advNO]
        type = GColumnMassAdvection
        variable = NO
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
 
      [./accumH2O]
        type = BedMassAccumulation
        variable = H2O
    [../]

    [./diffH2O]
        type = GColumnMassDispersion
        variable = H2O
        index = 3
    [../]

    [./advH2O]
        type = GColumnMassAdvection
        variable = H2O
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
 
#    [./NO_Ag_aging]
#        type = VariableOrderCoupledCatalyst
#        variable = AgNO3
#        main_variable = AgNO3
#        coupled_gases = 'NO N2'
#        gases_stoichiometry = '-3 1'
#        gas_order = '3 1'
#        coupled_catalysts = ''
#        catalyst_stoichiometry = ''
#        catalyst_order = ''
#        coupled_adsorption = 'AgNO3'
#        adsorbed_sites = '1'
#        site_order = '3'
#        adsorbed_stoichiometry = '1'
#        ads_order = '1'
#        max_capacity = 1.1102
#        forward_rate = 3.0e9
#        reverse_rate = 0.0
#    [../]
#
#    [./NO_MT]
#        type = CoupledCoeffTimeDerivative
#        variable = NO
#        coupled = AgNO3
#        time_coeff = 6.1115544
#    [../]
#
#    [./N2_MT]
#        type = CoupledCoeffTimeDerivative
#        variable = N2
#        coupled = Ag2O
#        time_coeff = 1.0185924
#    [../]
 
    [./AgNO3_MT]
        type = CoefTimeDerivative
        variable = AgNO3
    [../]

    [./Ag0_MT]
        type = CoefTimeDerivative
        variable = Ag0
    [../]
 
    [./Ag2O_MT]
        type = CoefTimeDerivative
        variable = Ag2O
    [../]
 
    [./AgMOR_MT]
        type = CoefTimeDerivative
        variable = AgMOR
    [../]
 
    [./AgNO3_Res_MT]
        type = CoefTimeDerivative
        variable = AgNO3_Res
    [../]
 
    [./XMOR_MT]
        type = CoefTimeDerivative
        variable = X_Sites
    [../]
 
    [./HMOR_MT]
        type = CoefTimeDerivative
        variable = HMOR
    [../]
 
    [./H2O_MT]
        type = CoefTimeDerivative
        variable = H2O
    [../]
 
    [./FeMOR_MT]
        type = CoefTimeDerivative
        variable = FeMOR
    [../]
 
    [./FeNO_MT]
        type = CoefTimeDerivative
        variable = FeNO
    [../]

    [./Ag2O_Adsortpion_1]
        type = VariableOrderReac
        variable = Ag2O
        main_variable = Ag2O
        coupled_species = 'Ag0 N2'
        stoichiometry = '-4 1'
        order = '2 1'
        main_stoichiometry = 2
        main_order = 2
        forward_rate = 0.02
        reverse_rate = 0.0
    [../]
 
    [./Ag2O_Adsortpion_2]
        type = VariableOrderReac
        variable = N2
        main_variable = N2
        coupled_species = 'Ag2O Ag0'
        stoichiometry = '2 -4'
        order = '2 2'
        main_stoichiometry = 1
        main_order = 1
        forward_rate = 0.040744
        reverse_rate = 0.0
    [../]
 
    [./Ag2O_Adsortpion_3]
        type = VariableOrderReac
        variable = Ag0
        main_variable = Ag0
        coupled_species = 'Ag2O N2'
        stoichiometry = '-2 -1'
        order = '2 1'
        main_stoichiometry = 4
        main_order = 2
        forward_rate = 0.0
        reverse_rate = 0.02
    [../]
 
    [./AgNO3_Adsortpion_1]
        type = VariableOrderReac
        variable = AgNO3
        main_variable = AgNO3
        coupled_species = 'Ag0 N2 Ag2O'
        stoichiometry = '2 1 -2'
        order = '2 1 2'
        main_stoichiometry = 2
        main_order = 2
        forward_rate = 15.0
        reverse_rate = 0.0
    [../]
 
    [./AgNO3_Adsortpion_2]
        type = VariableOrderReac
        variable = Ag0
        main_variable = Ag0
        coupled_species = 'AgNO3 N2 Ag2O'
        stoichiometry = '2 1 -2'
        order = '2 1 2'
        main_stoichiometry = 2
        main_order = 2
        forward_rate = 15.0
        reverse_rate = 0.0
    [../]
 
    [./AgNO3_Adsortpion_3]
        type = VariableOrderReac
        variable = N2
        main_variable = N2
        coupled_species = 'AgNO3 Ag0 Ag2O'
        stoichiometry = '2 2 -2'
        order = '2 2 2'
        main_stoichiometry = 1
        main_order = 1
        forward_rate = 30.558
        reverse_rate = 0.0
    [../]
 
    [./AgNO3_Adsortpion_4]
        type = VariableOrderReac
        variable = Ag2O
        main_variable = Ag2O
        coupled_species = 'AgNO3 N2 Ag0'
        stoichiometry = '-2 -1 -2'
        order = '2 1 2'
        main_stoichiometry = 2
        main_order = 2
        forward_rate = 0.0
        reverse_rate = 15.0
    [../]
 
    [./AgNO3_to_AgMOR_1]
        type = VariableOrderReac
        variable = AgMOR
        main_variable = AgMOR
        coupled_species = 'HMOR AgNO3 NO2 H2O'
        stoichiometry = '-2 -2 3 1'
        order = '4 2 3 1'
        main_stoichiometry = 2
        main_order = 2
        forward_rate = 10.0
        reverse_rate = 0.0
    [../]
 
    [./NO2_MT]
        type = CoupledCoeffTimeDerivative
        variable = NO2
        coupled = HMOR
        time_coeff = 3.0547776
    [../]
 
    [./H2O_MT2]
        type = CoupledCoeffTimeDerivative
        variable = H2O
        coupled = HMOR
        time_coeff = 1.0182592
    [../]
 
    [./AgNO3_to_AgMOR_2]
        type = VariableOrderReac
        variable = HMOR
        main_variable = HMOR
        coupled_species = 'AgMOR NO2 H2O AgNO3'
        stoichiometry = '-2 -3 -1 2'
        order = '2 3 1 2'
        main_stoichiometry = 2
        main_order = 4
        forward_rate = 0.0
        reverse_rate = 10.0
    [../]
 
    [./AgNO3_to_AgMOR_3]
        type = VariableOrderReac
        variable = AgNO3
        main_variable = AgNO3
        coupled_species = 'AgMOR NO2 H2O HMOR'
        stoichiometry = '-2 -3 -1 2'
        order = '2 3 1 4'
        main_stoichiometry = 2
        main_order = 2
        forward_rate = 0.0
        reverse_rate = 10.0
    [../]
 
    [./AgNO3_Crystal_1]
        type = VariableOrderReac
        variable = AgNO3_Res
        main_variable = AgNO3_Res
        coupled_species = 'AgNO3 X_Sites'
        stoichiometry = '-1 -1'
        order = '1 1'
        main_stoichiometry = 1
        main_order = 1
        forward_rate = 26.5
        reverse_rate = 0.3
    [../]
 
    [./AgNO3_Crystal_2]
        type = VariableOrderReac
        variable = AgNO3
        main_variable = AgNO3
        coupled_species = 'X_Sites AgNO3_Res'
        stoichiometry = '1 -1'
        order = '1 1'
        main_stoichiometry = 1
        main_order = 1
        forward_rate = 0.3
        reverse_rate = 26.5
    [../]
 
    [./AgNO3_Crystal_3]
        type = VariableOrderReac
        variable = X_Sites
        main_variable = X_Sites
        coupled_species = 'AgNO3 AgNO3_Res'
        stoichiometry = '1 -1'
        order = '1 1'
        main_stoichiometry = 1
        main_order = 1
        forward_rate = 0.3
        reverse_rate = 26.5
    [../]
 
    [./FeNitro_1]
        type = VariableOrderReac
        variable = FeNO
        main_variable = FeNO
        coupled_species = 'FeMOR X_Sites'
        stoichiometry = '-1 -2'
        order = '1 2'
        main_stoichiometry = 1
        main_order = 1
        forward_rate = 12.75
        reverse_rate = 0.000125
    [../]
 
    [./FeNitro_2]
        type = VariableOrderReac
        variable = FeMOR
        main_variable = FeMOR
        coupled_species = 'X_Sites FeNO'
        stoichiometry = '2 -1'
        order = '2 1'
        main_stoichiometry = 1
        main_order = 1
        forward_rate = 0.000125
        reverse_rate = 12.75
    [../]
 
    [./FeNitro_3]
        type = VariableOrderReac
        variable = X_Sites
        main_variable = X_Sites
        coupled_species = 'FeMOR FeNO'
        stoichiometry = '1 -1'
        order = '1 1'
        main_stoichiometry = 2
        main_order = 2
        forward_rate = 0.000125
        reverse_rate = 12.75
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
 
    [./dg_disp_NO]
        type = DGColumnMassDispersion
        variable = NO
        index = 1
    [../]
 
    [./dg_adv_NO]
        type = DGColumnMassAdvection
        variable = NO
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
 
    [./dg_disp_H2O]
        type = DGColumnMassDispersion
        variable = H2O
        index = 3
    [../]

    [./dg_adv_H2O]
        type = DGColumnMassAdvection
        variable = H2O
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
        coupled_gases = 'N2 NO NO2 H2O'
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
        input_molefraction = 0.99
        index = 0
    [../]
 
    [./NO_Flux]
        type = DGMassFluxBC
        variable = NO
        boundary = 'top bottom'
        input_temperature = 423.15
        input_pressure = 101.35
        input_molefraction = 0.01
        index = 1
    [../]
 
    [./NO2_Flux]
        type = DGMassFluxBC
        variable = NO2
        boundary = 'top bottom'
        input_temperature = 423.15
        input_pressure = 101.35
        input_molefraction = 0.0
        index = 2
    [../]

    [./H2O_Flux]
        type = DGMassFluxBC
        variable = H2O
        boundary = 'top bottom'
        input_temperature = 423.15
        input_pressure = 101.35
        input_molefraction = 0.0
        index = 3
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
        molecular_weight = '28.016 30.0 46.0 18.0'
        comp_heat_capacity = '1.04 0.996 0.8066 1.97'
        comp_ref_viscosity = '0.0001781 0.000294 0.00001963 0.0001043'
        comp_ref_temp = '300.55 288.0 400.0 298.16'
        comp_Sutherland_const = '111 128 701.226 784.72'
        temperature = column_temp
        total_pressure = total_pressure
        coupled_gases = 'N2 NO NO2 H2O'
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
        ref_diffusion = '0 0 0 0'
        activation_energy = '0 0 0 0'
        ref_temperature = '0 0 0 0'
        affinity = '0 0 0 0'
        temperature = column_temp
        coupled_gases = 'N2 NO NO2 H2O'
    [../]
 
    [./AdsorbateMaterials]
        type = ThermodynamicProperties
        block = 0
        temperature = column_temp
        total_pressure = total_pressure
        coupled_gases = 'N2 NO NO2 H2O'
        number_sites = '0 0 0 0'
        maximum_capacity = '0 0 0 0' #mol/kg
        molar_volume = '0 0 13500.0 13.91' #cm^3/mol
        enthalpy_site_1 = '0 0 0 -47390'
        enthalpy_site_2 = '0 0 0 -106390'
        enthalpy_site_3 = '0 0 0 -135640'
        enthalpy_site_4 = '0 0 0 0'
        enthalpy_site_5 = '0 0 0 0'
        enthalpy_site_6 = '0 0 0 0'
 
        entropy_site_1 = '0 0 0 -50.44'
        entropy_site_2 = '0 0 0 -165.59'
        entropy_site_3 = '0 0 0 -212.03'
        entropy_site_4 = '0 0 0 0'
        entropy_site_5 = '0 0 0 0'
        entropy_site_6 = '0 0 0 0'
    [../]
 
 [] #END Materials
 
[Postprocessors]
 
    [./press_exit]
        type = SideAverageValue
        boundary = 'top'
        variable = total_pressure
        execute_on = 'initial timestep_end'
    [../]
 
    [./Ag2O_Sol]
        type = ElementAverageValue
        variable = Ag2O
        execute_on = 'initial timestep_end'
    [../]
 
    [./AgNO3_Sol]
        type = ElementAverageValue
        variable = AgNO3
        execute_on = 'initial timestep_end'
    [../]
 
    [./AgNO3_Res_Sol]
        type = ElementAverageValue
        variable = AgNO3_Res
        execute_on = 'initial timestep_end'
    [../]
 
    [./Ag0_Sol]
        type = ElementAverageValue
        variable = Ag0
        execute_on = 'initial timestep_end'
    [../]
 
    [./AgMOR_Sol]
        type = ElementAverageValue
        variable = AgMOR
        execute_on = 'initial timestep_end'
    [../]
 
    [./XMOR_Sol]
        type = ElementAverageValue
        variable = X_Sites
        execute_on = 'initial timestep_end'
    [../]
 
    [./FeMOR_Sol]
        type = ElementAverageValue
        variable = FeMOR
        execute_on = 'initial timestep_end'
    [../]
 
    [./FeNO_Sol]
        type = ElementAverageValue
        variable = FeNO
        execute_on = 'initial timestep_end'
    [../]
 
    [./NO2_column]
        type = ElementAverageValue
        variable = 'NO2'
        execute_on = 'initial timestep_end'
    [../]
 
    [./HMOR_Sol]
        type = ElementAverageValue
        variable = HMOR
        execute_on = 'initial timestep_end'
    [../]
 
    [./H2O_column]
        type = ElementAverageValue
        variable = 'H2O'
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
    l_max_its = 200
    nl_max_its = 40
 
    solve_type = pjfnk
    line_search = bt    # Options: default none l2 bt basic
    start_time = 0.0
    end_time = 1345.0
    dtmax = 0.5
 
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
 
 

