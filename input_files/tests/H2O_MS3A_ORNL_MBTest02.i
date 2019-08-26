[GlobalParams]
	
	length = 250.0
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
 	ny = 20
 	xmin = 0.0
	xmax = 37.25 #cm
 	ymin = 0.0
	ymax = 250.0 #cm

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

	[./H2O]
		order = FIRST
		family = MONOMIAL
	[../]

	[./column_temp]
 		order = FIRST
 		family = MONOMIAL
 		initial_condition = 303.15
	[../]
 
	[./H2O_Adsorbed]
		order = FIRST
		family = MONOMIAL
		initial_condition = 0.0
	[../]
 
	[./H2O_AdsorbedHeat]
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
        initial_condition = 273.15
# 		initial_condition = 303.15
 	[../]
 
	[./wall_temp]
		order = FIRST
		family = MONOMIAL
		initial_condition = 303.15
	[../]

 [] #END AuxVariables

[ICs]

	[./N2_IC]
		type = ConcentrationIC
		variable = N2
		initial_mole_frac = 0.79
		initial_press = 101.35
		initial_temp = 303.15
	[../]

	[./O2_IC]
		type = ConcentrationIC
		variable = O2
		initial_mole_frac = 0.21
 		initial_press = 101.35
 		initial_temp = 303.15
	[../]

	[./H2O_IC]
		type = ConcentrationIC
		variable = H2O
		initial_mole_frac = 0.0
 		initial_press = 101.35
 		initial_temp = 303.15
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

 	[./accumH2O]
 		type = BedMassAccumulation
 		variable = H2O
 	[../]
 
	[./H2O_MT]
		type = SolidMassTransfer
		variable = H2O
		coupled = H2O_Adsorbed
	[../]

	[./diffH2O]
		type = GColumnMassDispersion
		variable = H2O
		index = 2
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
 
	[./H2O_columnAdsHeat]
		type = SolidHeatTransfer
		variable = column_temp
		coupled = H2O_AdsorbedHeat
	[../]

	[./H2O_adsheat]
		type = HeatofAdsorption
		variable = H2O_AdsorbedHeat
		coupled = H2O_Adsorbed
		index = 2
	[../]

	[./H2O_adsorption]
		type = CoupledLinearForcingFunction
		variable = H2O_Adsorbed
		coupled = H2O
		coeff = 16578.78
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

	[./dg_disp_H2O]
		type = DGColumnMassDispersion
		variable = H2O
		index = 2
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
		coupled_gases = 'N2 O2 H2O'
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
 		input_temperature = 303.15
 		input_pressure = 101.35
 		input_molefraction = 0.775
 		index = 0
 	[../]

 	[./O2_Flux]
		type = DGMassFluxBC
 		variable = O2
 		boundary = 'top bottom'
 		input_temperature = 303.15
 		input_pressure = 101.35
 		input_molefraction = 0.21
 		index = 1
 	[../]

 	[./H2O_Flux]
		type = DGMassFluxBC
 		variable = H2O
 		boundary = 'top bottom'
 		input_temperature = 303.15
 		input_pressure = 101.35
 		input_molefraction = 0.015
 		index = 2
 	[../]

	[./Heat_Gas_Flux]
 		type = DGHeatFluxBC
 		variable = column_temp
 		boundary = 'top bottom'
 		input_temperature = 303.15
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
		inner_diameter = 74.5
		outer_diameter = 75.5
		bulk_porosity = 0.585
		wall_density = 8.0
		wall_heat_capacity = 0.5
		wall_heat_trans_coef = 6.12
		extern_heat_trans_coef = 6.12
	[../]

	[./FlowMaterials]
		type = GasFlowProperties
		block = 0
		molecular_weight = '28.016 32 18'
		comp_heat_capacity = '1.04 0.919 1.97'
		comp_ref_viscosity = '0.0001781 0.0002018 0.0001043'
		comp_ref_temp = '300.55 292.25 298.16'
		comp_Sutherland_const = '111 127 784.72'
		flow_rate = 2.62e3
		temperature = column_temp
 		total_pressure = total_pressure
		coupled_gases = 'N2 O2 H2O'
	[../]

	[./AdsorbentMaterials]
		type = AdsorbentProperties
		block = 0
		binder_fraction = 0.175
		binder_porosity = 0.27
		crystal_radius = 1.5
		pellet_diameter = 0.236
		macropore_radius = 3.5e-6
		pellet_density = 1.69
		pellet_heat_capacity = 1.045
		ref_diffusion = '0 0 0.8814'
		activation_energy = '0 0 0'
		ref_temperature = '0 0 267.999'
		affinity = '0 0 0'
		temperature = column_temp
		coupled_gases = 'N2 O2 H2O'
	[../]

	[./AdsorbateMaterials]
		type = ThermodynamicProperties
		block = 0
		temperature = column_temp
		total_pressure = total_pressure
		coupled_gases = 'N2 O2 H2O'
		number_sites = '0 0 4'
		maximum_capacity = '0 0 11.67' #mol/kg 11.67
		molar_volume = '0 0 13.91' #cm^3/mol
		enthalpy_site_1 = '0 0 -46597.5'
		enthalpy_site_2 = '0 0 -125024'
		enthalpy_site_3 = '0 0 -193619'
		enthalpy_site_4 = '0 0 -272228'
		enthalpy_site_5 = '0 0 0'
		enthalpy_site_6 = '0 0 0'

		entropy_site_1 = '0 0 -53.6994'
		entropy_site_2 = '0 0 -221.073'
		entropy_site_3 = '0 0 -356.728'
		entropy_site_4 = '0 0 -567.459'
		entropy_site_5 = '0 0 0'
		entropy_site_6 = '0 0 0'
	[../]
 
#	[./KineticMaterials]
#		type = KineticProperties
#		block = 0
#		dirichlet_bc = false
#		heterogeneous = true
#		surface_diffusion = true
#       coupled_adsorption = H2O_Adsorbed (NOTE: This is causing the error. Need to have all species as adsorbed for kinetics)
#   [../]

 [] #END Materials

[Postprocessors]
 
	[./H2O_enter]
		type = SideAverageValue
		boundary = 'bottom'
		variable = H2O
		execute_on = 'initial timestep_end'
	[../]
 
	[./H2O_avg_gas]
		type = ElementAverageValue
		variable = H2O
		execute_on = 'initial timestep_end'
	[../]

	[./H2O_exit]
		type = SideAverageValue
		boundary = 'top'
		variable = H2O
		execute_on = 'initial timestep_end'
	[../]

	[./temp_exit]
		type = SideAverageValue
		boundary = 'top'
		variable = column_temp
		execute_on = 'initial timestep_end'
	[../]

	[./press_exit]
		type = SideAverageValue
		boundary = 'top'
		variable = total_pressure
		execute_on = 'initial timestep_end'
	[../]

 	[./wall_temp]
 		type = SideAverageValue
 		boundary = 'right'
 		variable = wall_temp
		execute_on = 'initial timestep_end'
 	[../]

	[./H2O_solid]
		type = ElementAverageValue
		variable = H2O_Adsorbed
		execute_on = 'initial timestep_end'
	[../]
 
	[./H2O_heat]
		type = ElementAverageValue
		variable = H2O_AdsorbedHeat
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
    line_search = bt    # Options: default none l2 bt
	start_time = 0.0
	end_time = 24.0

	[./TimeStepper]
		type = SolutionTimeAdaptiveDT
		dt = 0.001
	[../]

 [] #END Executioner

[Preconditioning]
 
	active = 'smp'
 
	[./none]
		type = SMP
		petsc_options = '-snes_converged_reason'
		petsc_options_iname = '-pc_type -ksp_gmres_restart'
		petsc_options_value = 'lu 2000'
	[../]
 
	[./smp]
		type = SMP
		full = true
		petsc_options = '-snes_converged_reason'
		petsc_options_iname = '-pc_type -ksp_gmres_restart  -snes_max_funcs'
		petsc_options_value = 'lu 2000 20000'
	[../]
 
	[./fdp]
		type = FDP
		full = true
		petsc_options = '-snes_converged_reason'
		petsc_options_iname = '-mat_fd_coloring_err -mat_fd_type'
		petsc_options_value = '1e-6 ds'
	[../]
 
	[./pbp]
		type = PBP
		solve_order = 'N2 O2 H2O H2O_Adsorbed H2O_AdsorbedHeat column_temp'
		preconditioner = 'lu lu lu lu lu lu'
		off_diag_row =    'O2 H2O H2O H2O_Adsorbed H2O_Adsorbed H2O_Adsorbed H2O_AdsorbedHeat H2O_AdsorbedHeat H2O_AdsorbedHeat H2O_AdsorbedHeat column_temp column_temp column_temp column_temp column_temp'
		off_diag_column = 'N2 N2  O2  N2           O2           H2O          N2               O2               H2O              H2O_Adsorbed     N2 O2 H2O H2O_Adsorbed H2O_AdsorbedHeat'
		petsc_options = '-snes_converged_reason'
		petsc_options_iname = '-pc_type -ksp_gmres_restart'
		petsc_options_value = 'lu 2000'
	[../]

[] #END Preconditioning

[Outputs]

	exodus = true
	csv = true
	print_linear_residuals = false

 [] #END Outputs
