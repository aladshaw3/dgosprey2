[GlobalParams]

[] #END GlobalParams

[Problem]
 
	coord_type = RZ
 
[] #END Problem

[Mesh]
 
	type = GeneratedMesh
	dim = 2
	nx = 3
	ny = 5
	xmin = 0.0
	xmax = 0.8636 #cm
	ymin = 0.0
	ymax = 22.86 #cm
 
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

	[./He]
		order = FIRST
		family = MONOMIAL
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
		initial_condition = 220.15
	[../]
 
	[./wall_temp]
		order = FIRST
		family = MONOMIAL
		initial_condition = 220.15
	[../]
 
	[./column_temp]
		order = FIRST
		family = MONOMIAL
		initial_condition = 220.15
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

[] #END AuxVariables

[ICs]

	[./Kr_IC]
		type = ConcentrationIC
		variable = Kr
		initial_mole_frac = 0.0
		initial_press = 101.35
		initial_temp = 220.15
	[../]

	[./Xe_IC]
		type = ConcentrationIC
		variable = Xe
		initial_mole_frac = 0.0
		initial_press = 101.35
		initial_temp = 220.15
	[../]

	[./He_IC]
		type = ConcentrationIC
		variable = He
		initial_mole_frac = 1.0
		initial_press = 101.35
		initial_temp = 220.15
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

	[./accumHe]
		type = BedMassAccumulation
		variable = He
	[../]

	[./diffHe]
		type = GColumnMassDispersion
		variable = He
		index = 2
	[../]

	[./advHe]
		type = GColumnMassAdvection
		variable = He
	[../]
 
	[./Kr_adsorption]
		type = CoupledLinearForcingFunction
		variable = Kr_Adsorbed
		coupled = Kr
		coeff = 277.55
	[../]
 
	[./Xe_adsorption]
		type = CoupledLinearForcingFunction
		variable = Xe_Adsorbed
		coupled = Xe
		coeff = 3154.93
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

	[./dg_disp_He]
		type = DGColumnMassDispersion
		variable = He
		index = 2
	[../]

	[./dg_adv_He]
		type = DGColumnMassAdvection
		variable = He
	[../]

[] #END DGKernels

[AuxKernels]

	[./column_pressure]
		type = TotalColumnPressure
		execute_on = 'initial timestep_end'
		variable = total_pressure
		temperature = column_temp
		coupled_gases = 'Kr Xe He'
	[../]

[] #END AuxKernels

[BCs]

	[./Kr_Flux]
		type = DGMassFluxBC
		variable = Kr
		boundary = 'top bottom'
		input_temperature = 220.15
		input_pressure = 101.35
		input_molefraction = 0.000114729
		index = 0
	[../]

	[./Xe_Flux]
		type = DGMassFluxBC
		variable = Xe
		boundary = 'top bottom'
		input_temperature = 220.15
		input_pressure = 101.35
		input_molefraction = 0.000750891
		index = 1
	[../]

	[./He_Flux]
		type = DGMassFluxBC
		variable = He
		boundary = 'top bottom'
		input_temperature = 220.15
		input_pressure = 101.35
		input_molefraction = 0.99913438
		index = 2
	[../]

[] #END BCs

[Materials]

	[./BedMaterials]
		type = BedProperties
		block = 0
		length = 22.86
		inner_diameter = 1.7272
		outer_diameter = 1.905
		bulk_porosity = 0.798				#not known
		wall_density = 7.7
		wall_heat_capacity = 0.5
		wall_heat_trans_coef = 9.0
		extern_heat_trans_coef = 90.0       #not known
	[../]

	[./FlowMaterials]
		type = GasFlowProperties
		block = 0
		molecular_weight = '83.8 131.29 4.0026'
		comp_heat_capacity = '0.25 0.16 5.1916'
		comp_ref_viscosity = '0.00023219 0.00021216 0.0001885'
		comp_ref_temp = '273.15 273.15 273.15'
		comp_Sutherland_const = '266.505 232.746 80.0'
		flow_rate = 2994.06
		temperature = column_temp
		total_pressure = total_pressure
		coupled_gases = 'Kr Xe He'
	[../]

	[./AdsorbentMaterials]
		type = AdsorbentProperties
		block = 0
		binder_fraction = 0.175				#not known
		binder_porosity = 0.27				#not known
		crystal_radius = 1.5				#not known
		pellet_diameter = 0.236				#not known
		macropore_radius = 3.5e-6			#not Known
		pellet_density = 1.69				#not Known
		pellet_heat_capacity = 1.045		#not known
		ref_diffusion = '0 0 0'
		activation_energy = '0 0 0'
		ref_temperature = '0 0 0'
		affinity = '0 0 0'
		temperature = column_temp
		coupled_gases = 'Kr Xe He'
	[../]

	[./AdsorbateMaterials]
		type = ThermodynamicProperties
		block = 0
		temperature = column_temp
		total_pressure = total_pressure
		coupled_gases = 'Kr Xe He'
		number_sites = '2 3 0'
		maximum_capacity = '1.716 1.479 0' #mol/kg
		molar_volume = '20.785 25.412 0' #cm^3/mol
 
		enthalpy_site_1 = '-44696.86 -18455.18 0'
		enthalpy_site_2 = '-65465.52 -35511.74 0'
		enthalpy_site_3 = '0 -53315.13 0'
		enthalpy_site_4 = '0 0 0'
		enthalpy_site_5 = '0 0 0'
		enthalpy_site_6 = '0 0 0'

		entropy_site_1 = '-170.45 -23.25 0'
		entropy_site_2 = '-248.55 -62.45 0'
		entropy_site_3 = '0 -100.10 0'
		entropy_site_4 = '0 0 0'
		entropy_site_5 = '0 0 0'
		entropy_site_6 = '0 0 0'
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

	[./He_exit]
		type = SideAverageValue
		boundary = 'top'
		variable = He
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

	[./Kr_solid]
		type = ElementAverageValue
		variable = Kr_Adsorbed
		execute_on = 'initial timestep_end'
	[../]

	[./Kr_heat]
		type = ElementAverageValue
		variable = Kr_AdsorbedHeat
		execute_on = 'initial timestep_end'
	[../]

	[./Xe_solid]
		type = ElementAverageValue
		variable = Xe_Adsorbed
		execute_on = 'initial timestep_end'
	[../]

	[./Xe_heat]
		type = ElementAverageValue
		variable = Xe_AdsorbedHeat
		execute_on = 'initial timestep_end'
	[../]

[] #END Postprocessors

[Executioner]

	type = Transient
	scheme = implicit-euler

	# NOTE: The default tolerances are far to strict and cause the program to crawl
	nl_rel_tol = 1e-6
	nl_abs_tol = 1e-6
	nl_rel_step_tol = 1e-10
	nl_abs_step_tol = 1e-10
	l_tol = 1e-6
	l_max_its = 100
	nl_max_its = 20

	solve_type = jfnk
	line_search = bt    # Options: default shell none basic l2 bt cp
	start_time = 0.0
	end_time = 0.0002
	dtmax = 0.1
	petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
	petsc_options_value = 'hypre boomeramg 100'

	[./TimeStepper]
		#Need to write a custom TimeStepper to enforce a maximum allowable dt
		type = ConstantDT
#		type = SolutionTimeAdaptiveDT
		dt = 0.0001
	[../]

[] #END Executioner

[Preconditioning]
	
	[./precond]
		type = PBP
		solve_order = 'He Kr Kr_Adsorbed Xe Xe_Adsorbed'
		preconditioner = 'AMG AMG AMG AMG AMG'
	[../]

[] #END Preconditioning

[Outputs]

	exodus = true
	csv = true
	print_linear_residuals = false

[] #END Outputs
