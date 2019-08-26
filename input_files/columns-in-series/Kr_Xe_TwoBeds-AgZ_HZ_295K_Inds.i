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
	nx = 5
	ny = 503
	xmin = 0.0
	xmax = 37.25 #cm
	ymin = 0.0
	ymax = 550.0 #cm
 
 [] # END Mesh
 
[MeshModifiers]
 
	[./column_1]
		type = SubdomainBoundingBox
		block_id = 0
		top_right = '37.25 250.0 0'
		bottom_left = '0 0 0'
	[../]
 
	[./conduit_12]
		type = SubdomainBoundingBox
		block_id = 1
		top_right = '37.25 300.0 0'
		bottom_left = '0 250.0 0'
	[../]
 
	[./column_2]
		type = SubdomainBoundingBox
		block_id = 2
		top_right = '37.25 550.0 0'
		bottom_left = '0 300.0 0'
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
 
	[./Kr]
		block = '0 1 2'
		order = FIRST
		family = MONOMIAL
	[../]
 
	[./Xe]
		block = '0 1 2'
		order = FIRST
		family = MONOMIAL
	[../]
 
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
 
	[./column_temp]
		block = '0 1 2'
		order = FIRST
		family = MONOMIAL
	[../]
 
	[./Kr_Adsorbed]
		block = '0 2'
		order = FIRST
		family = MONOMIAL
		initial_condition = 0.0
	[../]
 
	[./Xe_Adsorbed]
		block = '0 2'
		order = FIRST
		family = MONOMIAL
		initial_condition = 0.0
	[../]
 
	[./N2_Adsorbed]
		block = '0 2'
		order = FIRST
		family = MONOMIAL
		initial_condition = 0.0
	[../]
 
	[./Kr_AdsorbedHeat]
		block = '0 2'
		order = FIRST
		family = MONOMIAL
		initial_condition = 0.0
	[../]
 
	[./Xe_AdsorbedHeat]
		block = '0 2'
		order = FIRST
		family = MONOMIAL
		initial_condition = 0.0
	[../]
 
	[./N2_AdsorbedHeat]
		block = '0 2'
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
		initial_condition = 295.15
	[../]
 
	[./wall_temp_12]
		block = '1'
		order = FIRST
		family = MONOMIAL
		initial_condition = 222.15
	[../]
 
	[./wall_temp_2]
		block = '2'
		order = FIRST
		family = MONOMIAL
		initial_condition = 191.15
	[../]
 
 [] #END AuxVariables
 
#Consider custom ICs for all variables to set values based on blocks
[ICs]
 
	[./Kr_IC]
		type = ConcentrationIC
		block = '0 1 2'
		variable = Kr
		initial_mole_frac = 0.0
		initial_press = 101.35
		initial_temp = 295.15
	[../]
 
	[./Xe_IC]
		type = ConcentrationIC
		block = '0 1 2'
		variable = Xe
		initial_mole_frac = 0.0
		initial_press = 101.35
		initial_temp = 295.15
	[../]
 
	[./N2_IC]
		type = ConcentrationIC
		block = '0 1 2'
		variable = N2
		initial_mole_frac = 0.79
		initial_press = 101.35
		initial_temp = 295.15
	[../]
 
	[./O2_IC]
		type = ConcentrationIC
		block = '0 1 2'
		variable = O2
		initial_mole_frac = 0.21
		initial_press = 101.35
		initial_temp = 295.15
	[../]
 
 [./amb_ic_1]
 type = ConstantIC
 block = '0'
 variable = ambient_temp
 value = 295.15
	[../]
 
	[./amb_ic_12]
 type = ConstantIC
 block = '1'
 variable = ambient_temp
 value = 222.15
	[../]
 
	[./amb_ic_2]
 type = ConstantIC
 block = '2'
 variable = ambient_temp
 value = 191.15
	[../]
 
	[./col_ic_1]
 type = ConstantIC
 block = '0'
 variable = column_temp
 value = 295.15
	[../]
 
	[./col_ic_12]
 type = ConstantIC
 block = '1'
 variable = column_temp
 value = 222.15
	[../]
 
	[./col_ic_2]
 type = ConstantIC
 block = '2'
 variable = column_temp
 value = 191.15
	[../]
 
 [] #END ICs
 
[Kernels]
 
 [./accumKr]
	type = BedMassAccumulation
	block = '0 1 2'
	variable = Kr
 [../]
 
[./Kr_MT]
	type = SolidMassTransfer
	block = '0 2'
	variable = Kr
	coupled = Kr_Adsorbed
 [../]
 
[./diffKr]
	type = GColumnMassDispersion
	block = '0 2'
	variable = Kr
	index = 0
 [../]
 
[./advKr]
	type = GColumnMassAdvection
	block = '0 1 2'
	variable = Kr
 [../]
 
[./accumXe]
	type = BedMassAccumulation
	block = '0 1 2'
	variable = Xe
 [../]
 
[./Xe_MT]
	type = SolidMassTransfer
	block = '0 2'
	variable = Xe
	coupled = Xe_Adsorbed
 [../]
 
[./diffXe]
	type = GColumnMassDispersion
	block = '0 2'
	variable = Xe
	index = 1
 [../]
 
[./advXe]
	type = GColumnMassAdvection
	block = '0 1 2'
	variable = Xe
 [../]
 
[./accumN2]
	type = BedMassAccumulation
	block = '0 1 2'
	variable = N2
 [../]
 
 [./N2_MT]
	type = SolidMassTransfer
	block = '0 2'
	variable = N2
	coupled = N2_Adsorbed
 [../]
 
[./diffN2]
	type = GColumnMassDispersion
	block = '0 2'
	variable = N2
	index = 2
 [../]
 
[./advN2]
	type = GColumnMassAdvection
	block = '0 1 2'
	variable = N2
 [../]
 
 [./accumO2]
	type = BedMassAccumulation
	block = '0 1 2'
	variable = O2
 [../]
 
[./diffO2]
	type = GColumnMassDispersion
	block = '0 2'
	variable = O2
	index = 3
 [../]
 
[./advO2]
	type = GColumnMassAdvection
	block = '0 1 2'
	variable = O2
 [../]
 
[./columnAccum]
	type = BedHeatAccumulation
	block = '0 1 2'
	variable = column_temp
 [../]
 
[./columnConduction]
	type = GColumnHeatDispersion
	block = '0 2'
	variable =column_temp
 [../]
 
[./columnAdvection]
	type = GColumnHeatAdvection
	block = '0 1 2'
	variable =column_temp
 [../]
 
[./columnAdsHeat_Kr]
	type = SolidHeatTransfer
	block = '0 2'
	variable = column_temp
	coupled = Kr_AdsorbedHeat
 [../]
 
[./columnAdsHeat_Xe]
	type = SolidHeatTransfer
	block = '0 2'
	variable = column_temp
	coupled = Xe_AdsorbedHeat
 [../]
 
 [./columnAdsHeat_N2]
	type = SolidHeatTransfer
	block = '0 2'
	variable = column_temp
	coupled = N2_AdsorbedHeat
 [../]
 
[./Kr_adsheat]
	type = HeatofAdsorption
	block = '0 2'
	variable = Kr_AdsorbedHeat
	coupled = Kr_Adsorbed
	index = 0
 [../]
 
[./Xe_adsheat]
	type = HeatofAdsorption
	block = '0 2'
	variable = Xe_AdsorbedHeat
	coupled = Xe_Adsorbed
	index = 1
 [../]
 
 [./N2_adsheat]
	type = HeatofAdsorption
	block = '0 2'
	variable = N2_AdsorbedHeat
	coupled = N2_Adsorbed
	index = 2
 [../]
 
[./Kr_adsorption_1]
	type = CoupledExtendedLangmuirLDFModel
 	block = '0'
	variable = Kr_Adsorbed
	main_coupled = Kr
	coupled_temp = column_temp
	coupled_list = 'Kr Xe N2'
	enthalpies = '-31773 -23713 -9495'
	entropies = '-91.8 -22.59 -25.36'
	max_capacity = 1.35
	index = 0
	alpha = 15.0
	beta = 15.0
 [../]
 
[./Xe_adsorption_1]
	type = CoupledExtendedLangmuirLDFModel
 	block = '0'
	variable = Xe_Adsorbed
	main_coupled = Xe
	coupled_temp = column_temp
	coupled_list = 'Kr Xe N2'
	enthalpies = '-31773 -23713 -9495'
	entropies = '-91.8 -22.59 -25.36'
	max_capacity = 1.07
	index = 1
	alpha = 15.0
	beta = 15.0
 [../]
 
[./N2_adsorption_1]
	type = CoupledExtendedLangmuirLDFModel
 	block = '0'
	variable = N2_Adsorbed
	main_coupled = N2
	coupled_temp = column_temp
	coupled_list = 'Kr Xe N2'
	enthalpies = '-31773 -23713 -9495'
	entropies = '-91.8 -22.59 -25.36'
	max_capacity = 0.096
	index = 2
	alpha = 15.0
	beta = 15.0
 [../]
 
	[./Kr_adsorption_2]
		type = CoupledExtendedLangmuirModel
		block = '2'
		variable = Kr_Adsorbed
		main_coupled = Kr
		coupled_temp = column_temp
		coupled_list = 'Kr Xe N2'
 enthalpies = '-13212.6 -22875.8 -10816.1'
 entropies = '2.094 -29.866 -21.4657'
 max_capacity = 1.5965
	[../]
 
	[./Xe_adsorption_2]
		type = CoupledExtendedLangmuirModel
		block = '2'
		variable = Xe_Adsorbed
		main_coupled = Xe
		coupled_temp = column_temp
		coupled_list = 'Kr Xe N2'
 enthalpies = '-13212.6 -22875.8 -10816.1'
 entropies = '2.094 -29.866 -21.4657'
 max_capacity = 1.95
	[../]
 
	[./N2_adsorption_2]
		type = CoupledExtendedLangmuirModel
		block = '2'
		variable = N2_Adsorbed
		main_coupled = N2
		coupled_temp = column_temp
		coupled_list = 'Kr Xe N2'
 enthalpies = '-13212.6 -22875.8 -10816.1'
 entropies = '2.094 -29.866 -21.4657'
 max_capacity = 0.03
	[../]

 
 [] #END Kernels
 
[DGKernels]
 
 [./dg_disp_Kr]
	type = DGColumnMassDispersion
	block = '0 2'
	variable = Kr
	index = 0
 [../]
 
[./dg_adv_Kr]
	type = DGColumnMassAdvection
	block = '0 1 2'
	variable = Kr
 [../]
 
[./dg_disp_Xe]
	type = DGColumnMassDispersion
	block = '0 2'
	variable = Xe
	index = 1
 [../]
 
[./dg_adv_Xe]
	type = DGColumnMassAdvection
	block = '0 1 2'
	variable = Xe
 [../]
 
[./dg_disp_N2]
	type = DGColumnMassDispersion
	block = '0 2'
	variable = N2
	index = 2
 [../]
 
[./dg_adv_N2]
	type = DGColumnMassAdvection
	block = '0 1 2'
	variable = N2
 [../]
 
 [./dg_disp_O2]
	type = DGColumnMassDispersion
	block = '0 2'
	variable = O2
	index = 3
 [../]
 
[./dg_adv_O2]
	type = DGColumnMassAdvection
	block = '0 1 2'
	variable = O2
 [../]
 
[./dg_disp_heat]
	type = DGColumnHeatDispersion
	block = '0 2'
	variable = column_temp
 [../]
 
[./dg_adv_heat]
	type = DGColumnHeatAdvection
	block = '0 1 2'
	variable = column_temp
 [../]

 
 [] #END DGKernels
 
[AuxKernels]
 
	[./column_pressure]
		block = '0 1 2'
		type = TotalColumnPressure
		variable = total_pressure
		temperature = column_temp
		coupled_gases = 'Kr Xe N2 O2'
		execute_on = 'initial timestep_end'
	[../]
 
	[./wall_temp_calc_1]
		block = '0'
		type = WallTemperature
		variable = wall_temp_1
		column_temp = column_temp
		ambient_temp = 295.15
		execute_on = 'initial timestep_end'
	[../]
 
	[./wall_temp_calc_12]
		block = '1'
		type = WallTemperature
		variable = wall_temp_12
		column_temp = column_temp
		ambient_temp = 222.15
		execute_on = 'initial timestep_end'
	[../]
 
	[./wall_temp_calc_2]
		block = '2'
		type = WallTemperature
		variable = wall_temp_2
		column_temp = column_temp
		ambient_temp = 191.15
		execute_on = 'initial timestep_end'
	[../]
 
 [] #END AuxKernels
 
[BCs]
 
 [./Kr_Flux]
	type = DGMassFluxBC
	variable = Kr
	boundary = 'top bottom'
	input_temperature = 295.15
	input_pressure = 101.35
	input_molefraction = 0.000150122
	index = 0
 [../]
 
[./Xe_Flux]
	type = DGMassFluxBC
	variable = Xe
	boundary = 'top bottom'
	input_temperature = 295.15
	input_pressure = 101.35
	input_molefraction = 0.00099807
	index = 1
 [../]
 
[./N2_Flux]
	type = DGMassFluxBC
	variable = N2
	boundary = 'top bottom'
	input_temperature = 295.15
	input_pressure = 101.35
	input_molefraction = 0.7885
	index = 2
 [../]

 [./O2_Flux]
	type = DGMassFluxBC
	variable = O2
	boundary = 'top bottom'
	input_temperature = 295.15
	input_pressure = 101.35
	input_molefraction = 0.21035
	index = 3
 [../]

	[./Heat_Gas_Flux]
		type = DGHeatFluxBC
		variable = column_temp
		boundary = 'top bottom'
		input_temperature = 295.15
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
 
 # Approximate Properties for AgZ bed
	[./BedMaterials_1]
        type = BedProperties
        block = 0
		length = 250.0
		inner_diameter = 74.5
		outer_diameter = 76.5
		bulk_porosity = 0.8772				#not known
		wall_density = 7.7
		wall_heat_capacity = 0.5
		wall_heat_trans_coef = 9.0
		extern_heat_trans_coef = 9.0       #not known
	[../]
 
	[./FlowMaterials_1]
        type = GasFlowProperties
        block = 0
		flow_rate = 2.62e5
		molecular_weight = '83.8 131.29 28.016 32'
		comp_heat_capacity = '0.25 0.16 1.04 0.919'
		comp_ref_viscosity = '0.00023219 0.00021216 0.0001781 0.0002018'
		comp_ref_temp = '273.15 273.15 300.55 292.25'
		comp_Sutherland_const = '266.505 232.746 111 127'
		temperature = column_temp
		total_pressure = total_pressure
		coupled_gases = 'Kr Xe N2 O2'
	[../]
 
	[./AdsorbentMaterials_1]
        type = AdsorbentProperties
        block = 0
		pellet_diameter = 0.045
		binder_fraction = 0.0				#not known
		binder_porosity = 0.134				#not known
		crystal_radius = 0.0				#not known
		macropore_radius = 1.335e-7			#not Known
		pellet_density = 2.519				#not Known
		pellet_heat_capacity = 1.2  		#not known
		ref_diffusion = '0 0 0 0'
		activation_energy = '0 0 0 0'
		ref_temperature = '0 0 0 0'
		affinity = '0 0 0 0'
		temperature = column_temp
		coupled_gases = 'Kr Xe N2 O2'
	[../]
 
	[./AdsorbateMaterials_1]
        type = ThermodynamicProperties
        block = 0
		temperature = column_temp
		total_pressure = total_pressure
		coupled_gases = 'Kr Xe N2 O2'
		number_sites = '1 1 1 0'
		maximum_capacity = '1.35 1.07 0.096 0' #mol/kg
		molar_volume = '20.785 25.412 15.8 0' #cm^3/mol
 
		enthalpy_site_1 = '-31773 -23713 -9495 0'
		enthalpy_site_2 = '0 0 0 0'
		enthalpy_site_3 = '0 0 0 0'
		enthalpy_site_4 = '0 0 0 0'
		enthalpy_site_5 = '0 0 0 0'
		enthalpy_site_6 = '0 0 0 0'
 
		entropy_site_1 = '-91.8 -22.59 -25.36 0'
		entropy_site_2 = '0 0 0 0'
		entropy_site_3 = '0 0 0 0'
		entropy_site_4 = '0 0 0 0'
		entropy_site_5 = '0 0 0 0'
		entropy_site_6 = '0 0 0 0'
	[../]
 
	# Properties for the conduit between beds
	[./BedMaterials_12]
		type = BedProperties
		block = 1
		length = 50.0
		inner_diameter = 74.5
		outer_diameter = 76.5
		bulk_porosity = 1.0				#not known
		wall_density = 7.7
		wall_heat_capacity = 0.5
		wall_heat_trans_coef = 9.0
		extern_heat_trans_coef = 9.0       #not known
	[../]
 
	[./FlowMaterials_12]
		type = GasFlowProperties
		block = 1
		flow_rate = 2.62e5
		molecular_weight = '83.8 131.29 28.016 32'
		comp_heat_capacity = '0.25 0.16 1.04 0.919'
		comp_ref_viscosity = '0.00023219 0.00021216 0.0001781 0.0002018'
		comp_ref_temp = '273.15 273.15 300.55 292.25'
		comp_Sutherland_const = '266.505 232.746 111 127'
		temperature = column_temp
		total_pressure = total_pressure
		coupled_gases = 'Kr Xe N2 O2'
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
		coupled_gases = 'Kr Xe N2 O2'
	[../]
 
	[./AdsorbateMaterials_12]
		type = ThermodynamicProperties
		block = 1
		temperature = column_temp
		total_pressure = total_pressure
		coupled_gases = 'Kr Xe N2 O2'
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
 
 
	# Approximate Properties for HZ bed
	[./BedMaterials_2]
        type = BedProperties
        block = 2
		length = 250.0
		inner_diameter = 74.5
		outer_diameter = 76.5
		bulk_porosity = 0.95				#not known
		wall_density = 7.7
		wall_heat_capacity = 0.5
		wall_heat_trans_coef = 9.0
		extern_heat_trans_coef = 9.0       #not known
	[../]
 
	[./FlowMaterials_2]
        type = GasFlowProperties
        block = 2
		flow_rate = 2.62e5
		molecular_weight = '83.8 131.29 28.016 32'
		comp_heat_capacity = '0.25 0.16 1.04 0.919'
		comp_ref_viscosity = '0.00023219 0.00021216 0.0001781 0.0002018'
		comp_ref_temp = '273.15 273.15 300.55 292.25'
		comp_Sutherland_const = '266.505 232.746 111 127'
		temperature = column_temp
		total_pressure = total_pressure
		coupled_gases = 'Kr Xe N2 O2'
	[../]
 
	[./AdsorbentMaterials_2]
        type = AdsorbentProperties
        block = 2
		pellet_diameter = 0.056
		binder_fraction = 0.0				#not known
		binder_porosity = 0.25				#not known
		crystal_radius = 0.0				#not known
		macropore_radius = 1.945e-7			#not Known
		pellet_density = 2.5				#not Known
		pellet_heat_capacity = 1.2  		#not known
		ref_diffusion = '0 0 0 0'
		activation_energy = '0 0 0 0'
		ref_temperature = '0 0 0 0'
		affinity = '0 0 0 0'
		temperature = column_temp
		coupled_gases = 'Kr Xe N2 O2'
	[../]
 
	[./AdsorbateMaterials_2]
        type = ThermodynamicProperties
        block = 2
		temperature = column_temp
		total_pressure = total_pressure
		coupled_gases = 'Kr Xe N2 O2'
		number_sites = '1 1 1 0'
		maximum_capacity = '1.2 1.94 0.03 0' #mol/kg
		molar_volume = '20.785 25.412 15.8 0' #cm^3/mol
 
		enthalpy_site_1 = '-15758 -37630 -20959 0'
		enthalpy_site_2 = '0 0 0 0'
		enthalpy_site_3 = '0 0 0 0'
		enthalpy_site_4 = '0 0 0 0'
		enthalpy_site_5 = '0 0 0 0'
		enthalpy_site_6 = '0 0 0 0'
 
		entropy_site_1 = '-5.28 -96.9 -62.9 0'
		entropy_site_2 = '0 0 0 0'
		entropy_site_3 = '0 0 0 0'
		entropy_site_4 = '0 0 0 0'
		entropy_site_5 = '0 0 0 0'
		entropy_site_6 = '0 0 0 0'
	[../]
 

 [] #END Materials
 
[Postprocessors]
 
	[./Kr_enter_1]
		type = SideAverageValue
		boundary = 'bottom'
		variable = Kr
		execute_on = 'initial timestep_end'
	[../]
 
	[./Kr_enter_2]
		type = SideAverageValue
		boundary = 'interface_122'
		variable = Kr
		execute_on = 'initial timestep_end'
	[../]
 
	[./Kr_exit_1]
		type = SideAverageValue
		boundary = 'interface_112'
		variable = Kr
		execute_on = 'initial timestep_end'
	[../]
 
	[./Kr_avg_gas_1]
		block = 0
		type = ElementAverageValue
		variable = Kr
		execute_on = 'initial timestep_end'
	[../]
 
	[./Kr_avg_gas_2]
		block = 2
		type = ElementAverageValue
		variable = Kr
		execute_on = 'initial timestep_end'
	[../]
 
	[./Kr_exit_2]
		type = SideAverageValue
		boundary = 'top'
		variable = Kr
		execute_on = 'initial timestep_end'
	[../]
 
 [./N2_enter_1]
 type = SideAverageValue
 boundary = 'bottom'
 variable = N2
 execute_on = 'initial timestep_end'
	[../]
 
	[./N2_enter_2]
 type = SideAverageValue
 boundary = 'interface_122'
 variable = N2
 execute_on = 'initial timestep_end'
	[../]
 
	[./N2_exit_1]
 type = SideAverageValue
 boundary = 'interface_112'
 variable = N2
 execute_on = 'initial timestep_end'
	[../]
 
	[./N2_avg_gas_1]
 block = 0
 type = ElementAverageValue
 variable = N2
 execute_on = 'initial timestep_end'
	[../]
 
	[./N2_avg_gas_2]
 block = 2
 type = ElementAverageValue
 variable = N2
 execute_on = 'initial timestep_end'
	[../]
 
	[./N2_exit_2]
 type = SideAverageValue
 boundary = 'top'
 variable = N2
 execute_on = 'initial timestep_end'
	[../]
 
 [./Xe_enter_1]
	type = SideAverageValue
	boundary = 'bottom'
	variable = Xe
	execute_on = 'initial timestep_end'
	[../]
 
	[./Xe_enter_2]
	type = SideAverageValue
	boundary = 'interface_122'
	variable = Xe
	execute_on = 'initial timestep_end'
	[../]
 
	[./Xe_exit_1]
	type = SideAverageValue
	boundary = 'interface_112'
	variable = Xe
	execute_on = 'initial timestep_end'
	[../]
 
	[./Xe_avg_gas_1]
	block = 0
	type = ElementAverageValue
	variable = Xe
	execute_on = 'initial timestep_end'
	[../]
 
	[./Xe_avg_gas_2]
	block = 2
	type = ElementAverageValue
	variable = Xe
	execute_on = 'initial timestep_end'
	[../]
 
	[./Xe_exit_2]
	type = SideAverageValue
	boundary = 'top'
	variable = Xe
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
 
	[./Kr_solid_1]
		block = 0
		type = ElementAverageValue
		variable = Kr_Adsorbed
		execute_on = 'initial timestep_end'
	[../]
 
	[./Kr_solid_2]
		block = 2
		type = ElementAverageValue
		variable = Kr_Adsorbed
		execute_on = 'initial timestep_end'
	[../]
 
 [./Xe_solid_1]
		block = 0
		type = ElementAverageValue
		variable = Xe_Adsorbed
		execute_on = 'initial timestep_end'
	[../]
 
	[./Xe_solid_2]
		block = 2
		type = ElementAverageValue
		variable = Xe_Adsorbed
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
	end_time = 1000.0
	dtmax = 1.0
	
	[./TimeStepper]
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
