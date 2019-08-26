[GlobalParams]
 
	sigma = 1   # Penalty value:  NIPG = 0   otherwise, > 0
	epsilon = 1  #  -1 = SIPG   0 = IIPG   1 = NIPG
 
	flow_rate = 1.2e5
	length = 50.8
	inner_diameter = 1.7399
	pellet_diameter = 0.06305
	dt = 0.01
 
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
	xmax = 0.86995 #cm
	ymin = 0.0
	ymax = 50.8 #cm
 
 [] # END Mesh
 
[Variables]
 
[./Kr]
	order = FIRST
	family = MONOMIAL
 [../]
 
[./Xe]
	order = FIRST
	family = MONOMIAL
 [../]
 
[./N2]
	order = FIRST
	family = MONOMIAL
 [../]
 
 [./O2]
	order = FIRST
	family = MONOMIAL
 [../]
 
[./column_temp]
	order = FIRST
	family = MONOMIAL
	initial_condition = 253.15
 [../]
 
[./Kr_Adsorbed]
	order = FIRST
	family = MONOMIAL
	initial_condition = 0.0
 [../]
 
[./Xe_Adsorbed]
	order = FIRST
	family = MONOMIAL
	initial_condition = 0.0
 [../]
 
 [./N2_Adsorbed]
	order = FIRST
	family = MONOMIAL
	initial_condition = 0.0
 [../]
 
[./Kr_AdsorbedHeat]
	order = FIRST
	family = MONOMIAL
	initial_condition = 0.0
 [../]
 
[./Xe_AdsorbedHeat]
	order = FIRST
	family = MONOMIAL
	initial_condition = 0.0
 [../]
 
 [./N2_AdsorbedHeat]
	order = FIRST
	family = MONOMIAL
	initial_condition = 0.0
 [../]
 
 [] #END Variables
 
[AuxVariables]
 
[./total_pressure]
	order = FIRST
	family = MONOMIAL
	initial_condition = 101.35
 [../]
 
[./ambient_temp]
	order = FIRST
	family = MONOMIAL
	initial_condition = 253.15
 [../]
 
[./wall_temp]
	order = FIRST
	family = MONOMIAL
	initial_condition = 253.15
 [../]
 
 
 [] #END AuxVariables
 
[ICs]
 
[./Kr_IC]
	type = ConcentrationIC
	variable = Kr
	initial_mole_frac = 0.0
	initial_press = 101.35
	initial_temp = 253.15
 [../]
 
[./Xe_IC]
	type = ConcentrationIC
	variable = Xe
	initial_mole_frac = 0.0
	initial_press = 101.35
	initial_temp = 253.15
 [../]
 
[./N2_IC]
	type = ConcentrationIC
	variable = N2
	initial_mole_frac = 0.8
	initial_press = 101.35
	initial_temp = 253.15
 [../]
 
 [./O2_IC]
	type = ConcentrationIC
	variable = O2
	initial_mole_frac = 0.2
	initial_press = 101.35
	initial_temp = 253.15
 [../]
 
 [] #END ICs
 
[Kernels]
 
[./accumKr]
	type = BedMassAccumulation
	variable = Kr
 [../]
 
[./Kr_MT]
	type = SolidMassTransfer
	variable = Kr
	coupled = Kr_Adsorbed
 [../]
 
[./diffKr]
	type = GColumnMassDispersion
	variable = Kr
	index = 0
 [../]
 
[./advKr]
	type = GColumnMassAdvection
	variable = Kr
 [../]
 
[./accumXe]
	type = BedMassAccumulation
	variable = Xe
 [../]
 
[./Xe_MT]
	type = SolidMassTransfer
	variable = Xe
	coupled = Xe_Adsorbed
 [../]
 
[./diffXe]
	type = GColumnMassDispersion
	variable = Xe
	index = 1
 [../]
 
[./advXe]
	type = GColumnMassAdvection
	variable = Xe
 [../]
 
[./accumN2]
	type = BedMassAccumulation
	variable = N2
 [../]
 
 [./N2_MT]
	type = SolidMassTransfer
	variable = N2
	coupled = N2_Adsorbed
 [../]
 
[./diffN2]
	type = GColumnMassDispersion
	variable = N2
	index = 2
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
	index = 3
 [../]
 
[./advO2]
	type = GColumnMassAdvection
	variable = O2
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
 
[./columnAdsHeat_Kr]
	type = SolidHeatTransfer
	variable = column_temp
	coupled = Kr_AdsorbedHeat
 [../]
 
[./columnAdsHeat_Xe]
	type = SolidHeatTransfer
	variable = column_temp
	coupled = Xe_AdsorbedHeat
 [../]
 
 [./columnAdsHeat_N2]
	type = SolidHeatTransfer
	variable = column_temp
	coupled = N2_AdsorbedHeat
 [../]
 
[./Kr_adsheat]
	type = HeatofAdsorption
	variable = Kr_AdsorbedHeat
	coupled = Kr_Adsorbed
	index = 0
 [../]
 
[./Xe_adsheat]
	type = HeatofAdsorption
	variable = Xe_AdsorbedHeat
	coupled = Xe_Adsorbed
	index = 1
 [../]
 
 [./N2_adsheat]
	type = HeatofAdsorption
	variable = N2_AdsorbedHeat
	coupled = N2_Adsorbed
	index = 2
 [../]
 
[./Kr_adsorption]
	type = CoupledExtendedLangmuirLDFModel
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
 
[./Xe_adsorption]
	type = CoupledExtendedLangmuirLDFModel
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
 
[./N2_adsorption]
	type = CoupledExtendedLangmuirLDFModel
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
 
[] #END Kernels
 
[DGKernels]
 
[./dg_disp_Kr]
	type = DGColumnMassDispersion
	variable = Kr
	index = 0
 [../]
 
[./dg_adv_Kr]
	type = DGColumnMassAdvection
	variable = Kr
 [../]
 
[./dg_disp_Xe]
	type = DGColumnMassDispersion
	variable = Xe
	index = 1
 [../]
 
[./dg_adv_Xe]
	type = DGColumnMassAdvection
	variable = Xe
 [../]
 
[./dg_disp_N2]
	type = DGColumnMassDispersion
	variable = N2
	index = 2
 [../]
 
[./dg_adv_N2]
	type = DGColumnMassAdvection
	variable = N2
 [../]
 
 [./dg_disp_O2]
	type = DGColumnMassDispersion
	variable = O2
	index = 3
 [../]
 
[./dg_adv_O2]
	type = DGColumnMassAdvection
	variable = O2
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
	coupled_gases = 'Kr Xe N2 O2'
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
 
[./Kr_Flux]
	type = DGMassFluxBC
	variable = Kr
	boundary = 'top bottom'
	input_temperature = 253.15
	input_pressure = 101.35
	input_molefraction = 0.000128689
	index = 0
 [../]
 
[./Xe_Flux]
	type = DGMassFluxBC
	variable = Xe
	boundary = 'top bottom'
	input_temperature = 253.15
	input_pressure = 101.35
	input_molefraction = 0.000857934
	index = 1
 [../]
 
[./N2_Flux]
	type = DGMassFluxBC
	variable = N2
	boundary = 'top bottom'
	input_temperature = 253.15
	input_pressure = 101.35
	input_molefraction = 0.67775
	index = 2
 [../]
 
 [./O2_Flux]
	type = DGMassFluxBC
	variable = O2
	boundary = 'top bottom'
	input_temperature = 253.15
	input_pressure = 101.35
	input_molefraction = 0.32126
	index = 3
 [../]
 
[./Heat_Gas_Flux]
	type = DGHeatFluxBC
	variable = column_temp
	boundary = 'top bottom'
	input_temperature = 253.15
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
	outer_diameter = 1.905
	bulk_porosity = 0.691				#not known
	wall_density = 7.7
	wall_heat_capacity = 0.5
	wall_heat_trans_coef = 9.0
	extern_heat_trans_coef = 9.0       #not known
 [../]
 
[./FlowMaterials]
	type = GasFlowProperties
	block = 0
	molecular_weight = '83.8 131.29 28.016 32'
	comp_heat_capacity = '0.25 0.16 1.04 0.919'
	comp_ref_viscosity = '0.00023219 0.00021216 0.0001781 0.0002018'
	comp_ref_temp = '273.15 273.15 300.55 292.25'
	comp_Sutherland_const = '266.505 232.746 111 127'
	temperature = column_temp
	total_pressure = total_pressure
	coupled_gases = 'Kr Xe N2 O2'
 [../]
 
[./AdsorbentMaterials]
	type = AdsorbentProperties
	block = 0
	binder_fraction = 0.0			#not known
	binder_porosity = 0.134			#not known
	crystal_radius = 0.0				#not known
	macropore_radius = 1.335e-7		#not Known
	pellet_density = 1.2				#not Known
	pellet_heat_capacity = 1.2  		#not known
	ref_diffusion = '0 0 0 0'
	activation_energy = '0 0 0 0'
	ref_temperature = '0 0 0 0'
	affinity = '0 0 0 0'
	temperature = column_temp
	coupled_gases = 'Kr Xe N2 O2'
 [../]
 
[./AdsorbateMaterials]
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
 
 [] #END Materials
 
[Postprocessors]
 
[./Kr_exit]
	type = SideAverageValue
	boundary = 'top'
	variable = Kr
	execute_on = 'initial timestep_end'
 [../]
 
[./Xe_exit]
	type = SideAverageValue
	boundary = 'top'
	variable = Xe
	execute_on = 'initial timestep_end'
 [../]
 
 [./N2_exit]
	type = SideAverageValue
	boundary = 'top'
	variable = N2
	execute_on = 'initial timestep_end'
 [../]
 
[./Kr_avg]
	type = ElementAverageValue
	variable = Kr
	execute_on = 'initial timestep_end'
 [../]
 
[./Xe_avg]
	type = ElementAverageValue
	variable = Xe
	execute_on = 'initial timestep_end'
 [../]
 
 [./N2_avg]
	type = ElementAverageValue
	variable = N2
	execute_on = 'initial timestep_end'
 [../]
 
	[./temp_exit]
		type = SideAverageValue
		boundary = 'top'
		variable = column_temp
		execute_on = 'initial timestep_end'
	[../]
 
	[./wall_temp]
		type = SideAverageValue
		boundary = 'right'
		variable = wall_temp
		execute_on = 'initial timestep_end'
	[../]
 
[./Kr_solid]
	type = ElementAverageValue
	variable = Kr_Adsorbed
	execute_on = 'initial timestep_end'
 [../]
 
[./Xe_solid]
	type = ElementAverageValue
	variable = Xe_Adsorbed
	execute_on = 'initial timestep_end'
 [../]
 
 [./N2_solid]
	type = ElementAverageValue
	variable = N2_Adsorbed
	execute_on = 'initial timestep_end'
 [../]
 
 [] #END Postprocessors
 
[Executioner]
 
	type = Transient
	scheme = implicit-euler
 
	# NOTE: The default tolerances are far to strict and cause the program to crawl
	nl_rel_tol = 1e-10
	nl_abs_tol = 1e-4
	l_tol = 1e-6
	l_max_its = 200
	nl_max_its = 30
 
	solve_type = pjfnk
	line_search = basic    # Options: default none basic l2 bt
	start_time = 0.0
	end_time = 5.0
	dtmax = 0.5
 
	[./TimeStepper]
	#Need to write a custom TimeStepper to enforce a maximum allowable dt
#		type = ConstantDT
		type = SolutionTimeAdaptiveDT
	[../]
 
 [] #END Executioner
 
[Preconditioning]
 
[./precond]
	type = SMP
 
	full = true
	petsc_options = '-snes_converged_reason'
	petsc_options_iname = '-pc_type -ksp_gmres_restart -snes_max_funcs'
	petsc_options_value = 'lu 2000 60000'
 [../]
 
 [] #END Preconditioning
 
[Outputs]
 
	exodus = true
	csv = true
	print_linear_residuals = false
 
 [] #END Outputs
