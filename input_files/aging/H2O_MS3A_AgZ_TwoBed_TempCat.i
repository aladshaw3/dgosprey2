[GlobalParams]

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
	nx = 10
	ny = 80
	xmin = 0.0
	xmax = 3.0 #cm
	ymin = 0.0
	ymax = 40.4 #cm

 [] # END Mesh

[MeshModifiers]

	[./column_1]
		type = SubdomainBoundingBox
		block_id = 0
		top_right = '3.0 12.7 0'
		bottom_left = '0 0 0'
	[../]

	[./conduit_12]
		type = SubdomainBoundingBox
		block_id = 1
		top_right = '3.0 27.7 0'
		bottom_left = '0 12.7 0'
	[../]

	[./column_2]
		type = SubdomainBoundingBox
		block_id = 2
		top_right = '3.0 40.4 0'
		bottom_left = '0 27.7 0'
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

 [] #END MeshModifiers

[Variables]

	[./N2]
		block = '0 1 2'
		order = FIRST
		family = MONOMIAL
	[../]

	[./O2]
		block = '0 1 2'
		order = FIRST
		family = MONOMIAL
	[../]

	[./H2O]
		block = '0 1 2'
		order = FIRST
		family = MONOMIAL
	[../]

	[./I2]
		block = '0 1 2'
		order = FIRST
		family = MONOMIAL
	[../]

	[./column_temp]
		block = '0 1 2'
		order = FIRST
		family = MONOMIAL
	[../]

	[./H2O_Adsorbed]
		block = '0 2'
		order = FIRST
		family = MONOMIAL
		initial_condition = 0.0
	[../]

	[./H2O_AdsorbedHeat]
		block = '0 2'
		order = FIRST
		family = MONOMIAL
		initial_condition = 0.0
	[../]

	[./I2_Adsorbed]
		block = '0 2'
		order = FIRST
		family = MONOMIAL
		initial_condition = 0.0
	[../]

	[./I2_AdsorbedHeat]
		block = '0 2'
		order = FIRST
		family = MONOMIAL
		initial_condition = 0.0
	[../]

	[./AgOH]
		block = '2'
		order = FIRST
		family = MONOMIAL
		initial_condition = 0.0
	[../]

	[./Ag2O]
		block = '2'
		order = FIRST
		family = MONOMIAL
		initial_condition = 0.0
	[../]


 [] #END Variables

[AuxVariables]

	[./total_pressure]
		block = '0 1 2'
		order = FIRST
		family = MONOMIAL
		initial_condition = 101.35
	[../]

	[./ambient_temp]
		block = '0 1 2'
		order = FIRST
		family = MONOMIAL
	[../]

	[./wall_temp_1]
		block = '0'
		order = FIRST
		family = MONOMIAL
		initial_condition = 423.15
	[../]

	[./wall_temp_12]
		block = '1'
		order = FIRST
		family = MONOMIAL
		initial_condition = 423.15
	[../]

	[./wall_temp_2]
		block = '2'
		order = FIRST
		family = MONOMIAL
		initial_condition = 423.15
	[../]

 [] #END AuxVariables

#Consider custom ICs for all variables to set values based on blocks
[ICs]

	[./N2_IC]
		block = '0 1 2'
		type = ConcentrationIC
		variable = N2
		initial_mole_frac = 0.79
		initial_press = 101.35
		initial_temp = 423.15
	[../]

	[./O2_IC]
		block = '0 1 2'
		type = ConcentrationIC
		variable = O2
		initial_mole_frac = 0.21
		initial_press = 101.35
		initial_temp = 423.15
	[../]

	[./H2O_IC]
		block = '0 1 2'
		type = ConcentrationIC
		variable = H2O
		initial_mole_frac = 0.0
		initial_press = 101.35
		initial_temp = 423.15
	[../]

	[./I2_IC]
		block = '0 1 2'
		type = ConcentrationIC
		variable = I2
		initial_mole_frac = 0.0
		initial_press = 101.35
		initial_temp = 423.15
	[../]

	[./amb_ic_1]
		type = ConstantIC
		block = '0'
		variable = ambient_temp
		value = 423.15
	[../]

	[./amb_ic_12]
		type = ConstantIC
		block = '1'
		variable = ambient_temp
		value = 423.15
	[../]

	[./amb_ic_2]
		type = ConstantIC
		block = '2'
		variable = ambient_temp
		value = 423.15
	[../]

	[./col_ic_1]
		type = ConstantIC
		block = '0'
		variable = column_temp
		value = 423.15
	[../]

	[./col_ic_12]
		type = ConstantIC
		block = '1'
		variable = column_temp
		value = 423.15
	[../]

	[./col_ic_2]
		type = ConstantIC
		block = '2'
		variable = column_temp
		value = 423.15
	[../]

 [] #END ICs

[Kernels]

	[./accumN2]
		block = '0 1 2'
		type = BedMassAccumulation
		variable = N2
	[../]

	[./diffN2]
		block = '0 2'
		type = GColumnMassDispersion
		variable = N2
		index = 0
	[../]

	[./advN2]
		block = '0 1 2'
		type = GColumnMassAdvection
		variable = N2
	[../]

	[./accumO2]
		block = '0 1 2'
		type = BedMassAccumulation
		variable = O2
	[../]

	[./diffO2]
		block = '0 2'
		type = GColumnMassDispersion
		variable = O2
		index = 1
	[../]

	[./advO2]
		block = '0 1 2'
		type = GColumnMassAdvection
		variable = O2
	[../]

	[./accumH2O]
		block = '0 1 2'
		type = BedMassAccumulation
		variable = H2O
	[../]

	[./H2O_MT]
		block = '0 2'
		type = SolidMassTransfer
		variable = H2O
		coupled = H2O_Adsorbed
	[../]

	[./diffH2O]
		block = '0 2'
		type = GColumnMassDispersion
		variable = H2O
		index = 2
	[../]

	[./advH2O]
		block = '0 1 2'
		type = GColumnMassAdvection
		variable = H2O
	[../]

	[./accumI2]
		block = '0 1 2'
		type = BedMassAccumulation
		variable = I2
	[../]

#	[./I2_MT]
#		block = '0 2'
#		type = SolidMassTransfer
#		variable = I2
#		coupled = I2_Adsorbed
#	[../]

	[./I2_MT_1]
		block = '0'
		type = SolidMassTransfer
		variable = I2
		coupled = I2_Adsorbed
	[../]

	[./I2_MT_2]
		block = '2'
		type = CoupledCoeffTimeDerivative
		variable = I2
		coupled = I2_Adsorbed
		time_coeff = 0.42    #rho_b * 0.5
	[../]

	[./diffI2]
		block = '0 2'
		type = GColumnMassDispersion
		variable = I2
		index = 3
	[../]

	[./advI2]
		block = '0 1 2'
		type = GColumnMassAdvection
		variable = I2
	[../]


	[./columnAccum]
		block = '0 1 2'
		type = BedHeatAccumulation
		variable = column_temp
	[../]

	[./columnConduction]
		block = '0 2'
		type = GColumnHeatDispersion
		variable =column_temp
	[../]

	[./columnAdvection]
		block = '0 1 2'
		type = GColumnHeatAdvection
		variable =column_temp
	[../]

	[./H2O_columnAdsHeat]
		block = '0 2'
		type = SolidHeatTransfer
		variable = column_temp
		coupled = H2O_AdsorbedHeat
	[../]

	[./I2_columnAdsHeat]
		block = '0 2'
		type = SolidHeatTransfer
		variable = column_temp
		coupled = I2_AdsorbedHeat
	[../]

	[./H2O_adsheat]
		block = '0 2'
		type = HeatofAdsorption
		variable = H2O_AdsorbedHeat
		coupled = H2O_Adsorbed
		index = 2
	[../]

	[./I2_adsheat]
		block = '0 2'
		type = HeatofAdsorption
		variable = I2_AdsorbedHeat
		coupled = I2_Adsorbed
		index = 3
	[../]

    [./H2O_adsorption]
		block = '0 2'
        type = CoupledGSTALDFmodel
        variable = H2O_Adsorbed
        coupled_gas = H2O
        coupled_temp = column_temp
        index = 2
	[../]

	[./I2_adsorption_1]
		block = '0'
		type = CoupledGSTALDFmodel
		variable = I2_Adsorbed
		coupled_gas = I2
		coupled_temp = column_temp
		index = 3
	[../]

	[./I2_adsorption_2]
		block = '2'
		type = CoupledConstChemisorption
		variable = I2_Adsorbed
		main_variable = I2_Adsorbed
		coupled_gases = 'I2'
		gases_stoichiometry = '-1'
		coupled_adsorption = 'I2_Adsorbed AgOH Ag2O'
		adsorbed_sites = '2 4 4'
		adsorbed_stoichiometry = '2 4 2'
		max_capacity = 0.9764
		forward_rate = 75000.0
		reverse_rate = 1.0E-10
	[../]

	[./AgOH_ads]
		block = '2'
		type = CoupledConstChemisorption
		variable = AgOH
		main_variable = AgOH
		coupled_gases = 'O2 H2O'
		gases_stoichiometry = '-1 -2'
		coupled_adsorption = 'I2_Adsorbed AgOH Ag2O'
		adsorbed_sites = '2 4 4'
		adsorbed_stoichiometry = '2 4 2'
		max_capacity = 0.9764
		forward_rate = 1.5E+9
		reverse_rate = 1.48E-1
	[../]

	[./Ag2O_ads]
		block = '2'
		type = CoupledConstChemisorption
		variable = Ag2O
		main_variable = Ag2O
		coupled_gases = 'O2'
		gases_stoichiometry = '-1'
		coupled_adsorption = 'I2_Adsorbed AgOH Ag2O'
		adsorbed_sites = '2 4 4'
		adsorbed_stoichiometry = '2 4 2'
		max_capacity = 0.9764
		forward_rate = 6.18e-2
		reverse_rate = 1.57e-3
	[../]

 [] #END Kernels

[DGKernels]

	[./dg_disp_N2]
		block = '0 2'
		type = DGColumnMassDispersion
		variable = N2
		index = 0
	[../]

	[./dg_adv_N2]
		block = '0 1 2'
		type = DGColumnMassAdvection
		variable = N2
	[../]

	[./dg_disp_O2]
		block = '0 2'
		type = DGColumnMassDispersion
		variable = O2
		index = 1
	[../]

	[./dg_adv_O2]
		block = '0 1 2'
		type = DGColumnMassAdvection
		variable = O2
	[../]

	[./dg_disp_H2O]
		block = '0 2'
		type = DGColumnMassDispersion
		variable = H2O
		index = 2
	[../]

	[./dg_adv_H2O]
		block = '0 1 2'
		type = DGColumnMassAdvection
		variable = H2O
	[../]

	[./dg_disp_I2]
		block = '0 2'
		type = DGColumnMassDispersion
		variable = I2
		index = 3
	[../]

	[./dg_adv_I2]
		block = '0 1 2'
		type = DGColumnMassAdvection
		variable = I2
	[../]

	[./dg_disp_heat]
		block = '0 2'
		type = DGColumnHeatDispersion
		variable = column_temp
	[../]

	[./dg_adv_heat]
		block = '0 1 2'
		type = DGColumnHeatAdvection
		variable = column_temp
	[../]

 [] #END DGKernels

[AuxKernels]

	[./column_pressure]
		block = '0 1 2'
		type = TotalColumnPressure
		variable = total_pressure
		temperature = column_temp
		coupled_gases = 'N2 O2 H2O'
		execute_on = 'initial timestep_end'
	[../]

	[./wall_temp_calc_1]
		block = '0'
		type = WallTemperature
		variable = wall_temp_1
		column_temp = column_temp
		ambient_temp = ambient_temp
		execute_on = 'initial timestep_end'
	[../]

	[./wall_temp_calc_12]
		block = '1'
		type = WallTemperature
		variable = wall_temp_12
		column_temp = column_temp
		ambient_temp = ambient_temp
		execute_on = 'initial timestep_end'
	[../]

	[./wall_temp_calc_2]
		block = '2'
		type = WallTemperature
		variable = wall_temp_2
		column_temp = column_temp
		ambient_temp = ambient_temp
		execute_on = 'initial timestep_end'
	[../]

 [] #END AuxKernels

[BCs]
#NOTE: Temperature at input does not match because we wanted to represent H2O and I2 at -15 oC (at 100 oC) and 10 ppmv, respectively
	[./N2_Flux]
		type = DGMassFluxBC
		variable = N2
		boundary = 'top bottom'
		input_temperature = 423.15
		input_pressure = 101.35
		input_molefraction = 0.788620587
		index = 0
	[../]

	[./O2_Flux]
		type = DGMassFluxBC
		variable = O2
		boundary = 'top bottom'
		input_temperature = 423.15
		input_pressure = 101.35
		input_molefraction = 0.20974
		index = 1
	[../]

	[./H2O_Flux]
		type = DGMassFluxBC
		variable = H2O
		boundary = 'top bottom'
		input_temperature = 423.15
		input_pressure = 101.35
		input_molefraction = 0.001629413
		index = 2
	[../]

	[./I2_Flux]
		type = DGMassFluxBC
		variable = I2
		boundary = 'top bottom'
		input_temperature = 423.15
		input_pressure = 101.35
		input_molefraction = 0.00005
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
		wall_temp = wall_temp_1
	[../]

	[./Heat_Wall_Flux_12]
		type = DGColumnWallHeatFluxLimitedBC
		variable = column_temp
		boundary = 'wall_12 left'
		wall_temp = wall_temp_12
	[../]

	[./Heat_Wall_Flux_2]
		type = DGColumnWallHeatFluxLimitedBC
		variable = column_temp
		boundary = 'wall_2 left'
		wall_temp = wall_temp_2
	[../]

 [] #END BCs

[Materials]

 # Approximate Properties for MS3A bed
	[./BedMaterials_1]
        type = BedProperties
        block = 0
		length = 12.7
		inner_diameter = 3.0
        outer_diameter = 3.30
        bulk_porosity = 0.541
        wall_density = 8.0
        wall_heat_capacity = 0.5
        wall_heat_trans_coef = 6.12
        extern_heat_trans_coef = 6.12
	[../]

	[./FlowMaterials_1]
        type = GasFlowProperties
        block = 0
		flow_rate = 66000.0
		molecular_weight = '28.016 32 18 253.8'
		comp_heat_capacity = '1.04 0.919 1.97 0.214'
		comp_ref_viscosity = '0.0001781 0.0002018 0.0001043 0.00013283'
		comp_ref_temp = '300.55 292.25 298.16 295.496'
		comp_Sutherland_const = '111 127 784.72 573.474'
		temperature = column_temp
		total_pressure = total_pressure
		coupled_gases = 'N2 O2 H2O I2'
	[../]

	[./AdsorbentMaterials_1]
        type = AdsorbentProperties
        block = 0
        binder_fraction = 0.175
        binder_porosity = 0.27
        crystal_radius = 1.5
		pellet_diameter = 0.236
        macropore_radius = 3.5e-6
        pellet_density = 1.69
        pellet_heat_capacity = 1.045
        ref_diffusion = '0 0 0.8814 0'
        activation_energy = '0 0 0 0'
        ref_temperature = '0 0 267.999 0'
        affinity = '0 0 0 0'
        temperature = column_temp
        coupled_gases = 'N2 O2 H2O I2'
	[../]

	[./AdsorbateMaterials_1]
        type = ThermodynamicProperties
        block = 0
        temperature = column_temp
        total_pressure = total_pressure
		coupled_gases = 'N2 O2 H2O I2'
		number_sites = '0 0 4 1'
		maximum_capacity = '0 0 11.67 0.0236' #mol/kg
		molar_volume = '0 0 13.91 41.8' #cm^3/mol
		enthalpy_site_1 = '0 0 -46597.5 -28663'
		enthalpy_site_2 = '0 0 -125024 0'
		enthalpy_site_3 = '0 0 -193619 0'
		enthalpy_site_4 = '0 0 -272228 0'
		enthalpy_site_5 = '0 0 0 0'
		enthalpy_site_6 = '0 0 0 0'

		entropy_site_1 = '0 0 -53.6994 5.15'
		entropy_site_2 = '0 0 -221.073 0'
		entropy_site_3 = '0 0 -356.728 0'
		entropy_site_4 = '0 0 -567.459 0'
		entropy_site_5 = '0 0 0 0'
		entropy_site_6 = '0 0 0 0'
	[../]

# Properties for the conduit between beds
	[./BedMaterials_12]
		type = BedProperties
		block = 1
		length = 15.0
		inner_diameter = 3.0
		outer_diameter = 3.3
		bulk_porosity = 1.0
		wall_density = 12.0
		wall_heat_capacity = 1.5
		wall_heat_trans_coef = 9.12
		extern_heat_trans_coef = 9.12
	[../]

	[./FlowMaterials_12]
		type = GasFlowProperties
		block = 1
		flow_rate = 66000.0
		molecular_weight = '28.016 32 18 253.8'
		comp_heat_capacity = '1.04 0.919 1.97 0.214'
		comp_ref_viscosity = '0.0001781 0.0002018 0.0001043 0.00013283'
		comp_ref_temp = '300.55 292.25 298.16 295.496'
		comp_Sutherland_const = '111 127 784.72 573.474'
		temperature = column_temp
		total_pressure = total_pressure
		coupled_gases = 'N2 O2 H2O I2'
	[../]

	[./AdsorbentMaterials_12]
		type = AdsorbentProperties
		block = 1
		binder_fraction = 0.0
		binder_porosity = 0.0
		crystal_radius = 0.0
		pellet_diameter = 0.0
		macropore_radius = 0.0
		pellet_density = 0.0
		pellet_heat_capacity = 0.0
		ref_diffusion = '0 0 0 0'
		activation_energy = '0 0 0 0'
		ref_temperature = '0 0 0 0'
		affinity = '0 0 0 0'
		temperature = column_temp
		coupled_gases = 'N2 O2 H2O I2'
	[../]

	[./AdsorbateMaterials_12]
		type = ThermodynamicProperties
		block = 1
		temperature = column_temp
		total_pressure = total_pressure
		coupled_gases = 'N2 O2 H2O I2'
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


# Approximate Properties for Ag0Z bed
	[./BedMaterials_2]
        type = BedProperties
        block = 2
		length = 12.7
		inner_diameter = 3.0
        outer_diameter = 3.3
		bulk_porosity = 0.7255
        wall_density = 8.0
        wall_heat_capacity = 0.5
        wall_heat_trans_coef = 6.12
        extern_heat_trans_coef = 6.12
	[../]

	[./FlowMaterials_2]
        type = GasFlowProperties
        block = 2
		flow_rate = 66000.0
		molecular_weight = '28.016 32 18 253.8'
		comp_heat_capacity = '1.04 0.919 1.97 0.214'
		comp_ref_viscosity = '0.0001781 0.0002018 0.0001043 0.00013283'
		comp_ref_temp = '300.55 292.25 298.16 295.496'
		comp_Sutherland_const = '111 127 784.72 573.474'
		temperature = column_temp
		total_pressure = total_pressure
		coupled_gases = 'N2 O2 H2O I2'
	[../]

	[./AdsorbentMaterials_2]
        type = AdsorbentProperties
        block = 2
        binder_fraction = 0.0
        binder_porosity = 0.384
        crystal_radius = 0.0
		pellet_diameter = 0.16
        macropore_radius = 2.6e-6
        pellet_density = 3.06
		pellet_heat_capacity = 1.2
        ref_diffusion = '0 0 0 0'
        activation_energy = '0 0 0 0'
        ref_temperature = '0 0 0 0'
        affinity = '0 0 0 0'
        temperature = column_temp
        coupled_gases = 'N2 O2 H2O I2'
	[../]

	[./AdsorbateMaterials_2]
        type = ThermodynamicProperties
        block = 2
        temperature = column_temp
        total_pressure = total_pressure
		coupled_gases = 'N2 O2 H2O I2'
		number_sites = '0 0 3 1'
		maximum_capacity = '0 0 4.75 0.55' #mol/kg
		molar_volume = '0 0 13.91 41.8' #cm^3/mol
		enthalpy_site_1 = '0 0 -47390 -8308'
		enthalpy_site_2 = '0 0 -106390 0'
		enthalpy_site_3 = '0 0 -135640 0'
		enthalpy_site_4 = '0 0 0 0'
		enthalpy_site_5 = '0 0 0 0'
		enthalpy_site_6 = '0 0 0 0'

		entropy_site_1 = '0 0 -50.44 92.35'
		entropy_site_2 = '0 0 -165.59 0'
		entropy_site_3 = '0 0 -212.03 0'
		entropy_site_4 = '0 0 0 0'
		entropy_site_5 = '0 0 0 0'
		entropy_site_6 = '0 0 0 0'
	[../]


 [] #END Materials

[Postprocessors]

	[./H2O_exit_1]
		type = SideAverageValue
		boundary = 'interface_112'
		variable = H2O
		execute_on = 'initial timestep_end'
	[../]

	[./H2O_exit_2]
		type = SideAverageValue
		boundary = 'top'
		variable = H2O
		execute_on = 'initial timestep_end'
	[../]

	[./I2_exit_1]
		type = SideAverageValue
		boundary = 'interface_112'
		variable = I2
		execute_on = 'initial timestep_end'
	[../]

	[./I2_exit_2]
		type = SideAverageValue
		boundary = 'top'
		variable = I2
		execute_on = 'initial timestep_end'
	[../]

   [./temp_exit_2]
       type = SideAverageValue
       boundary = 'top'
       variable = column_temp
       execute_on = 'initial timestep_end'
   [../]

	[./temp_exit_1]
		type = SideAverageValue
		boundary = 'interface_112'
		variable = column_temp
		execute_on = 'initial timestep_end'
	[../]

	[./wall_temp_1]
		type = SideAverageValue
		boundary = 'wall_1'
		variable = wall_temp_1
		execute_on = 'initial timestep_end'
	[../]

	[./wall_temp_12]
		type = SideAverageValue
		boundary = 'wall_12'
		variable = wall_temp_12
		execute_on = 'initial timestep_end'
	[../]

	[./wall_temp_2]
		type = SideAverageValue
		boundary = 'wall_2'
		variable = wall_temp_2
		execute_on = 'initial timestep_end'
	[../]

	[./H2O_MS3A]
		block = 0
		type = ElementAverageValue
		variable = H2O_Adsorbed
		execute_on = 'initial timestep_end'
	[../]

	[./H2O_Ag0Z]
		block = 2
		type = ElementAverageValue
		variable = H2O_Adsorbed
		execute_on = 'initial timestep_end'
	[../]

	[./I2_MS3A]
		block = 0
		type = ElementAverageValue
		variable = I2_Adsorbed
		execute_on = 'initial timestep_end'
	[../]

	[./2x_AgI]
		block = 2
		type = ElementAverageValue
		variable = I2_Adsorbed
		execute_on = 'initial timestep_end'
	[../]

	[./AgOH]
		block = 2
		type = ElementAverageValue
		variable = AgOH
		execute_on = 'initial timestep_end'
	[../]

	[./Ag2O]
		block = 2
		type = ElementAverageValue
		variable = Ag2O
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
	nl_max_its = 30

	solve_type = pjfnk
	line_search = basic    # Options: default none l2 bt
	start_time = 0.0
	end_time = 500.0
	dtmax = 1.0

	[./TimeStepper]
		type = SolutionTimeAdaptiveDT
	[../]

 [] #END Executioner

# [Adaptivity]
#NOTE: Cannot use adaptivity with StatefulMaterials + DG Methods
# [] #END Adaptivity

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

