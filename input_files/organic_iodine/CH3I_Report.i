[GlobalParams]

    length = 15.0
    pellet_diameter = 0.18
    inner_diameter = 2.0
    flow_rate = 60000
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

[MeshModifiers]

    [./column_1]
        type = SubdomainBoundingBox
        block_id = 0
        top_right = '1.0 2.0 0'
        bottom_left = '0 0 0'
    [../]

    [./conduit_12]
        type = SubdomainBoundingBox
        block_id = 1
        top_right = '1.0 2.5 0'
        bottom_left = '0 2.0 0'
    [../]

    [./column_2]
        type = SubdomainBoundingBox
        block_id = 2
        top_right = '1.0 6.5 0'
        bottom_left = '0 2.5 0'
    [../]

    [./conduit_23]
        type = SubdomainBoundingBox
        block_id = 3
        top_right = '1.0 7.0 0'
        bottom_left = '0 6.5 0'
    [../]

    [./column_3]
        type = SubdomainBoundingBox
        block_id = 4
        top_right = '1.0 15.0 0'
        bottom_left = '0 7.0 0'
    [../]

    [./interface_112]
        type = SideSetsBetweenSubdomains
        depends_on = 'column_1 conduit_12'
        master_block = 0
        paired_block = 1
        new_boundary = 'interface_112'
    [../]

    [./interface_122]
        type = SideSetsBetweenSubdomains
        depends_on = 'conduit_12 column_2'
        master_block = 1
        paired_block = 2
        new_boundary = 'interface_122'
    [../]

    [./interface_223]
        type = SideSetsBetweenSubdomains
        depends_on = 'column_2 conduit_23'
        master_block = 2
        paired_block = 3
        new_boundary = 'interface_223'
    [../]

    [./interface_233]
        type = SideSetsBetweenSubdomains
        depends_on = 'conduit_23 column_3'
        master_block = 3
        paired_block = 4
        new_boundary = 'interface_233'
    [../]

    [./wall_1]
        type = SideSetsAroundSubdomain
        block = '0'
        new_boundary = 'wall_1'
        normal = '1 0 0'
        force_prepare = true
        depends_on = column_1
    [../]

    [./wall_12]
        type = SideSetsAroundSubdomain
        block = '1'
        new_boundary = 'wall_12'
        normal = '1 0 0'
        force_prepare = true
        depends_on = conduit_12
    [../]

    [./wall_2]
        type = SideSetsAroundSubdomain
        block = '2'
        new_boundary = 'wall_2'
        normal = '1 0 0'
        force_prepare = true
        depends_on = column_2
    [../]

    [./wall_23]
        type = SideSetsAroundSubdomain
        block = '3'
        new_boundary = 'wall_23'
        normal = '1 0 0'
        force_prepare = true
        depends_on = conduit_23
    [../]

    [./wall_3]
        type = SideSetsAroundSubdomain
        block = '4'
        new_boundary = 'wall_3'
        normal = '1 0 0'
        force_prepare = true
        depends_on = column_3
    [../]

[] #END MeshModifiers


[Variables]

    [./N2]
        block = '0 1 2 3 4'
        order = FIRST
        family = MONOMIAL
    [../]

    [./H2]
        block = '0 1 2 3 4'
        order = FIRST
        family = MONOMIAL
    [../]

    [./CH3I]
        block = '0 1 2 3 4'
        order = FIRST
        family = MONOMIAL
    [../]

    [./C2H6]
        block = '0 1 2 3 4'
        order = FIRST
        family = MONOMIAL
    [../]

    [./column_temp]
        block = '0 1 2 3 4'
        order = FIRST
        family = MONOMIAL
        initial_condition = 423.15
    [../]

    [./Ag0]
        block = '0 2 4'
        order = FIRST
        family = MONOMIAL
        initial_condition = 3.28125
    [../]

    [./AgI]
        block = '0 2 4'
        order = FIRST
        family = MONOMIAL
        initial_condition = 0.0
    [../]

    [./rSH]
        block = '0 2 4'
        order = FIRST
        family = MONOMIAL
        initial_condition = 3.28125
    [../]

    [./rSAg]
        block = '0 2 4'
        order = FIRST
        family = MONOMIAL
        initial_condition = 0.0
    [../]

 [] #END Variables

[AuxVariables]

    [./total_pressure]
        block = '0 1 2 3 4'
        order = CONSTANT
        family = MONOMIAL
        initial_condition = 101.35
    [../]

    [./ambient_temp]
        block = '0 1 2 3 4'
         order = CONSTANT
        family = MONOMIAL
        initial_condition = 423.15
    [../]

    [./wall_temp]
        block = '0 1 2 3 4'
        order = FIRST
        family = MONOMIAL
        initial_condition = 423.15
    [../]

 [] #END AuxVariables

[ICs]

    [./N2_IC]
        block = '0 1 2 3 4'
        type = ConcentrationIC
        variable = N2
        initial_mole_frac = 1.0
        initial_press = 101.35
        initial_temp = 423.15
    [../]

    [./H2_IC]
        block = '0 1 2 3 4'
        type = ConcentrationIC
        variable = H2
        initial_mole_frac = 0.0
        initial_press = 101.35
        initial_temp = 423.15
    [../]

    [./CH3I_IC]
        block = '0 1 2 3 4'
        type = ConcentrationIC
        variable = CH3I
        initial_mole_frac = 0.0
        initial_press = 101.35
        initial_temp = 423.15
    [../]

    [./C2H6_IC]
        block = '0 1 2 3 4'
        type = ConcentrationIC
        variable = C2H6
        initial_mole_frac = 0.0
        initial_press = 101.35
        initial_temp = 423.15
    [../]

 [] #END ICs

[Kernels]

    [./accumN2]
        block = '0 1 2 3 4'
        type = BedMassAccumulation
        variable = N2
    [../]

    [./diffN2]
        block = '0 1 2 3 4'
        type = GColumnMassDispersion
        variable = N2
        index = 0
    [../]

    [./advN2]
        block = '0 1 2 3 4'
        type = GColumnMassAdvection
        variable = N2
    [../]

    [./accumH2]
        block = '0 1 2 3 4'
        type = BedMassAccumulation
        variable = H2
    [../]

    [./diffH2]
        block = '0 1 2 3 4'
        type = GColumnMassDispersion
        variable = H2
        index = 1
    [../]

    [./advH2]
        block = '0 1 2 3 4'
        type = GColumnMassAdvection
        variable = H2
    [../]

    [./accumCH3I]
        block = '0 1 2 3 4'
        type = BedMassAccumulation
        variable = CH3I
    [../]

    [./diffCH3I]
        block = '0 1 2 3 4'
        type = GColumnMassDispersion
        variable = CH3I
        index = 2
    [../]

    [./advCH3I]
        block = '0 1 2 3 4'
        type = GColumnMassAdvection
        variable = CH3I
    [../]

    [./accumC2H6]
        block = '0 1 2 3 4'
        type = BedMassAccumulation
        variable = C2H6
    [../]

    [./diffC2H6]
        block = '0 1 2 3 4'
        type = GColumnMassDispersion
        variable = C2H6
        index = 3
    [../]

    [./advC2H6]
        block = '0 1 2 3 4'
        type = GColumnMassAdvection
        variable = C2H6
    [../]

    [./columnAccum]
        block = '0 1 2 3 4'
        type = BedHeatAccumulation
        variable = column_temp
    [../]

    [./columnConduction]
        block = '0 2 4'
        type = GColumnHeatDispersion
        variable =column_temp
    [../]

    [./columnAdvection]
        block = '0 1 2 3 4'
        type = GColumnHeatAdvection
        variable =column_temp
    [../]

    [./Ag0_MT]
        block = '0 2 4'
        type = CoefTimeDerivative
        variable = Ag0
    [../]

    [./AgI_MT]
        block = '0 2 4'
        type = CoefTimeDerivative
        variable = AgI
    [../]

    [./rSH_MT]
        block = '0 2 4'
        type = CoefTimeDerivative
        variable = rSH
    [../]

    [./rSAg_MT]
        block = '0 2 4'
        type = CoefTimeDerivative
        variable = rSAg
    [../]

    [./C2H6_MT]
        block = '0 2 4'
        type = CoefTimeDerivative
        variable = C2H6
    [../]

    [./Aging_1]
        block = '0 2 4'
        type = VariableOrderReac
        variable = rSAg
        main_variable = rSAg
        coupled_species = 'rSH Ag0 H2'
        stoichiometry = '-2 -2 1'
        order = '2 5 1'
        main_stoichiometry = 2
        main_order = 2
        forward_rate = 3.05e-7
        reverse_rate = 0.0
    [../]

    [./Aging_2]
        block = '0 2 4'
        type = VariableOrderReac
        variable = Ag0
        main_variable = Ag0
        coupled_species = 'rSH rSAg H2'
        stoichiometry = '2 -2 -1'
        order = '2 2 1'
        main_stoichiometry = 2
        main_order = 5
        forward_rate = 0.0
        reverse_rate = 3.05e-7
    [../]

    [./Aging_3]
        block = '0 2 4'
        type = VariableOrderReac
        variable = rSH
        main_variable = rSH
        coupled_species = 'rSAg Ag0 H2'
        stoichiometry = '-2 2 -1'
        order = '2 5 1'
        main_stoichiometry = 2
        main_order = 2
        forward_rate = 0.0
        reverse_rate = 3.05e-7
    [../]

    [./H2_Aging_1]
        block = '0 4'
        type = CoupledCoeffTimeDerivative
        variable = H2
        coupled = rSH
        time_coeff = 0.15875
    [../]

    [./H2_Aging_2]
        block = '2'
        type = CoupledCoeffTimeDerivative
        variable = H2
        coupled = rSH
        time_coeff = 0.17875
    [../]

    [./Adsorption_1]
        block = '0 2 4'
        type = VariableOrderReac
        variable = AgI
        main_variable = AgI
        coupled_species = 'CH3I Ag0 C2H6'
        stoichiometry = '-2 -2 1'
        order = '2 3 1'
        main_stoichiometry = 2
        main_order = 2
        forward_rate = 6.5e7
        reverse_rate = 0.0
    [../]

    [./Adsorption_2]
        block = '0 2 4'
        type = VariableOrderReac
        variable = Ag0
        main_variable = Ag0
        coupled_species = 'CH3I AgI C2H6'
        stoichiometry = '2 -2 -1'
        order = '2 2 1'
        main_stoichiometry = 2
        main_order = 3
        forward_rate = 0.0
        reverse_rate = 6.5e7
    [../]

    [./Adsorption_3]
        block = '0 4'
        type = VariableOrderReac
        variable = C2H6
        main_variable = C2H6
        coupled_species = 'CH3I Ag0 AgI'
        stoichiometry = '-2 -2 2'
        order = '2 3 2'
        main_stoichiometry = 1
        main_order = 1
        forward_rate = 10318750.0
        reverse_rate = 0.0
    [../]

    [./Adsorption_4]
        block = '2'
        type = VariableOrderReac
        variable = C2H6
        main_variable = C2H6
        coupled_species = 'CH3I Ag0 AgI'
        stoichiometry = '-2 -2 2'
        order = '2 3 2'
        main_stoichiometry = 1
        main_order = 1
        forward_rate = 11618750.0
        reverse_rate = 0.0
    [../]

    [./CH3I_Adsorption]
        block = '0 2 4'
        type = SolidMassTransfer
        variable = CH3I
        coupled = AgI
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
#
#    [./Aerogel_Adsorption]
#        block = '0 2 4'
#        type = VariableOrderCoupledCatalyst
#        variable = AgI
#        main_variable = AgI
#        coupled_gases = 'CH3I C2H6'
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
#    [./CH3I_Adsorption]
#        block = '0 2 4'
#        type = SolidMassTransfer
#        variable = CH3I
#        coupled = AgI
#    [../]
#
#    [./C2H6_Adsorption_1]
#        block = '0 4'
#        type = CoupledCoeffTimeDerivative
#        variable = C2H6
#        coupled = AgI
#        time_coeff = 0.15875
#    [../]
#
#    [./C2H6_Adsorption_2]
#        block = '2'
#        type = CoupledCoeffTimeDerivative
#        variable = C2H6
#        coupled = AgI
#        time_coeff = 0.17875
#    [../]
 
 
 [] #END Kernels

[DGKernels]

    [./dg_disp_N2]
        block = '0 2 4'
        type = DGColumnMassDispersion
        variable = N2
        index = 0
    [../]

    [./dg_adv_N2]
        block = '0 1 2 3 4'
        type = DGColumnMassAdvection
        variable = N2
    [../]

    [./dg_disp_H2]
        block = '0 2 4'
        type = DGColumnMassDispersion
        variable = H2
        index = 1
    [../]

    [./dg_adv_H2]
        block = '0 1 2 3 4'
        type = DGColumnMassAdvection
        variable = H2
    [../]

    [./dg_disp_CH3I]
        block = '0 2 4'
        type = DGColumnMassDispersion
        variable = CH3I
        index = 2
    [../]

    [./dg_adv_CH3I]
        block = '0 1 2 3 4'
        type = DGColumnMassAdvection
        variable = CH3I
    [../]

    [./dg_disp_C2H6]
        block = '0 2 4'
        type = DGColumnMassDispersion
        variable = C2H6
        index = 3
    [../]

    [./dg_adv_C2H6]
        block = '0 1 2 3 4'
        type = DGColumnMassAdvection
        variable = C2H6
    [../]

    [./dg_disp_heat]
        block = '0 2 4'
        type = DGColumnHeatDispersion
        variable = column_temp
    [../]

    [./dg_adv_heat]
        block = '0 1 2 3 4'
        type = DGColumnHeatAdvection
        variable = column_temp
    [../]

 [] #END DGKernels

[AuxKernels]

    [./column_pressure]
        block = '0 1 2 3 4'
        type = TotalColumnPressure
        variable = total_pressure
        temperature = column_temp
        coupled_gases = 'N2 H2 CH3I C2H6'
        execute_on = 'initial timestep_end'
    [../]

    [./wall_temp_calc]
        block = '0 1 2 3 4'
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
        input_molefraction = 0.999995
        index = 0
    [../]

    [./H2_Flux]
        type = DGMassFluxBC
        variable = H2
        boundary = 'top bottom'
        input_temperature = 423.15
        input_pressure = 101.35
        input_molefraction = 0.0
        index = 1
    [../]

    [./CH3I_Flux]
        type = DGMassFluxBC
        variable = CH3I
        boundary = 'top bottom'
        input_temperature = 423.15
        input_pressure = 101.35
        input_molefraction = 0.000005
        index = 2
    [../]

    [./C2H6_Flux]
        type = DGMassFluxBC
        variable = C2H6
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

    [./Heat_Wall_Flux_1]
        type = DGColumnWallHeatFluxLimitedBC
        variable = column_temp
        boundary = 'wall_1 left'
        wall_temp = wall_temp
    [../]

    [./Heat_Wall_Flux_12]
        type = DGColumnWallHeatFluxLimitedBC
        variable = column_temp
        boundary = 'wall_12 left'
        wall_temp = wall_temp
    [../]

    [./Heat_Wall_Flux_2]
        type = DGColumnWallHeatFluxLimitedBC
        variable = column_temp
        boundary = 'wall_2 left'
        wall_temp = wall_temp
    [../]

    [./Heat_Wall_Flux_23]
        type = DGColumnWallHeatFluxLimitedBC
        variable = column_temp
        boundary = 'wall_23 left'
        wall_temp = wall_temp
    [../]

    [./Heat_Wall_Flux_3]
        type = DGColumnWallHeatFluxLimitedBC
        variable = column_temp
        boundary = 'wall_3 left'
        wall_temp = wall_temp
    [../]


 [] #END BCs

[Materials]

    [./BedMaterials_1]
        type = BedProperties
        block = '0 4'
        outer_diameter = 2.3
        bulk_porosity = 0.365
        wall_density = 8.0
        wall_heat_capacity = 0.5
        wall_heat_trans_coef = 6.12
        extern_heat_trans_coef = 6.12
    [../]

    [./BedMaterials_2]
        type = BedProperties
        block = '2'
        outer_diameter = 2.3
        bulk_porosity = 0.285
        wall_density = 8.0
        wall_heat_capacity = 0.5
        wall_heat_trans_coef = 6.12
        extern_heat_trans_coef = 6.12
    [../]

    [./BedMaterials_3]
        type = BedProperties
        block = '1 3'
        outer_diameter = 2.3
        bulk_porosity = 1.0
        wall_density = 8.0
        wall_heat_capacity = 0.5
        wall_heat_trans_coef = 6.12
        extern_heat_trans_coef = 6.12
    [../]

    [./FlowMaterials]
        type = GasFlowProperties
        block = '0 1 2 3 4'
        molecular_weight = '28.016 2.0 141.939 30.07'
        comp_heat_capacity = '1.04 14.307 2.232 1.28'
        comp_ref_viscosity = '0.0001781 0.0000876 0.00011 0.00008817'
        comp_ref_temp = '300.55 293.85 293.15 280'
        comp_Sutherland_const = '111 72 93 213'
        temperature = column_temp
        total_pressure = total_pressure
        coupled_gases = 'N2 H2 CH3I C2H6'
    [../]

    [./AdsorbentMaterials]
        block = '0 1 2 3 4'
        type = AdsorbentProperties
        binder_fraction = 0.0
        binder_porosity = 0.384
        crystal_radius = 0.0
        macropore_radius = 2.65e-6
        pellet_density = 0.5
        pellet_heat_capacity = 1.2
        ref_diffusion = '0 0 0 0'
        activation_energy = '0 0 0 0'
        ref_temperature = '0 0 0 0'
        affinity = '0 0 0 0'
        temperature = column_temp
        coupled_gases = 'N2 H2 CH3I C2H6'
    [../]

    [./AdsorbateMaterials]
        block = '0 1 2 3 4'
        type = ThermodynamicProperties
        temperature = column_temp
        total_pressure = total_pressure
        coupled_gases = 'N2 H2 CH3I C2H6'
        number_sites = '0 0 0 0'
        maximum_capacity = '0 0 0 0' #mol/kg
        molar_volume = '0 0 0 0' #cm^3/mol
        enthalpy_site_1 = '0 0 0 0'
        enthalpy_site_2 = '0 0 0 0'
        enthalpy_site_3 = '0 0 0 0'
        enthalpy_site_4 = '0 0 0 0'
        enthalpy_site_5 = '0 0 0 0'
        enthalpy_site_6 = '0 0 0 0'

        entropy_site_1 = '0 0 0 0'
        entropy_site_2 = '0 0 0 0'
        entropy_site_3 = '0 0 0 0'
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

    [./AgI_solid_1]
        block = 0
        type = ElementAverageValue
        variable = AgI
        execute_on = 'initial timestep_end'
    [../]

    [./AgI_solid_2]
        block = 2
        type = ElementAverageValue
        variable = AgI
        execute_on = 'initial timestep_end'
    [../]

    [./AgI_solid_3]
        block = 4
        type = ElementAverageValue
        variable = AgI
        execute_on = 'initial timestep_end'
    [../]

    [./rSH_solid_1]
        block = 0
        type = ElementAverageValue
        variable = rSH
        execute_on = 'initial timestep_end'
    [../]

    [./rSH_solid_2]
        block = 2
        type = ElementAverageValue
        variable = rSH
        execute_on = 'initial timestep_end'
    [../]

    [./rSH_solid_3]
        block = 4
        type = ElementAverageValue
        variable = rSH
        execute_on = 'initial timestep_end'
    [../]

    [./rSAg_solid_1]
        block = 0
        type = ElementAverageValue
        variable = rSAg
        execute_on = 'initial timestep_end'
    [../]

    [./rSAg_solid_2]
        block = 2
        type = ElementAverageValue
        variable = rSAg
        execute_on = 'initial timestep_end'
    [../]

    [./rSAg_solid_3]
        block = 4
        type = ElementAverageValue
        variable = rSAg
        execute_on = 'initial timestep_end'
    [../]

    [./CH3I_exit_1]
        type = SideAverageValue
        boundary = 'interface_112'
        variable = CH3I
        execute_on = 'initial timestep_end'
    [../]

    [./CH3I_exit_2]
        type = SideAverageValue
        boundary = 'interface_223'
        variable = CH3I
        execute_on = 'initial timestep_end'
    [../]

    [./CH3I_exit_3]
        type = SideAverageValue
        boundary = 'top'
        variable = CH3I
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
    dtmax = 10.0

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
