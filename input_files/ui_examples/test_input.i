[GlobalParams]
	initial_dt = 0.01
	length = 22.86
[] #END GlobalParams

[Problem]
	coord_type = RZ
[] #END Problem

[Mesh]
	dim = 2
	nx = 10
	ny = 40
	type = GeneratedMesh
	xmax = 0.8636
	xmin = 0.0
	ymax = 22.86
	ymin = 0.0
[] #END Mesh

[Variables]
	[./He]
		family = MONOMIAL
		order = CONSTANT
	[../]
	[./Kr]
		family = MONOMIAL
		order = CONSTANT
	[../]
	[./Xe]
		family = MONOMIAL
		order = CONSTANT
	[../]
	[./column_temp]
		family = MONOMIAL
		initial_condition = 253.15
		order = CONSTANT
	[../]
	[./wall_temp]
		family = MONOMIAL
		initial_condition = 253.15
		order = CONSTANT
	[../]
[] #END Variables

[AuxVariables]
	[./He_Adsorbed]
		family = MONOMIAL
		initial_condition = 0.0
		order = CONSTANT
	[../]
	[./He_AdsorbedHeat]
		family = MONOMIAL
		initial_condition = 0.0
		order = CONSTANT
	[../]
	[./He_Perturb]
		family = MONOMIAL
		initial_condition = 0.0
		order = CONSTANT
	[../]
	[./Kr_Adsorbed]
		family = MONOMIAL
		initial_condition = 0.0
		order = CONSTANT
	[../]
	[./Kr_AdsorbedHeat]
		family = MONOMIAL
		initial_condition = 0.0
		order = CONSTANT
	[../]
	[./Kr_Perturb]
		family = MONOMIAL
		initial_condition = 0.0
		order = CONSTANT
	[../]
	[./Xe_Adsorbed]
		family = MONOMIAL
		initial_condition = 0.0
		order = CONSTANT
	[../]
	[./Xe_AdsorbedHeat]
		family = MONOMIAL
		initial_condition = 0.0
		order = CONSTANT
	[../]
	[./Xe_Perturb]
		family = MONOMIAL
		initial_condition = 0.0
		order = CONSTANT
	[../]
	[./ambient_temp]
		family = MONOMIAL
		initial_condition = 253.15
		order = CONSTANT
	[../]
	[./total_pressure]
		family = MONOMIAL
		initial_condition = 101.35
		order = CONSTANT
	[../]
[] #END AuxVariables

[ICs]
	[./He_IC]
		initial_mole_frac = 1.0
		initial_press = 101.35
		initial_temp = 253.15
		type = ConcentrationIC
		variable = He
	[../]
	[./Kr_IC]
		initial_mole_frac = 0.0
		initial_press = 101.35
		initial_temp = 253.15
		type = ConcentrationIC
		variable = Kr
	[../]
	[./Xe_IC]
		initial_mole_frac = 0.0
		initial_press = 101.35
		initial_temp = 253.15
		type = ConcentrationIC
		variable = Xe
	[../]
[] #END ICs

[Kernels]
	[./He_accum]
		index = 2
		type = BedMassAccumulation
		variable = He
	[../]
	[./He_adv]
		type = GColumnMassAdvection
		variable = He
	[../]
	[./He_diff]
		index = 2
		type = GColumnMassDispersion
		variable = He
	[../]
	[./He_masstrans]
		solid_conc = He_Adsorbed
		type = AdsorptionMassTransfer
		variable = He
	[../]
	[./Kr_accum]
		index = 0
		type = BedMassAccumulation
		variable = Kr
	[../]
	[./Kr_adv]
		type = GColumnMassAdvection
		variable = Kr
	[../]
	[./Kr_diff]
		index = 0
		type = GColumnMassDispersion
		variable = Kr
	[../]
	[./Kr_masstrans]
		solid_conc = Kr_Adsorbed
		type = AdsorptionMassTransfer
		variable = Kr
	[../]
	[./Xe_accum]
		index = 1
		type = BedMassAccumulation
		variable = Xe
	[../]
	[./Xe_adv]
		type = GColumnMassAdvection
		variable = Xe
	[../]
	[./Xe_diff]
		index = 1
		type = GColumnMassDispersion
		variable = Xe
	[../]
	[./Xe_masstrans]
		solid_conc = Xe_Adsorbed
		type = AdsorptionMassTransfer
		variable = Xe
	[../]
	[./columnAccum]
		type = BedHeatAccumulation
		variable = column_temp
	[../]
	[./columnAdvection]
		type = GColumnHeatAdvection
		variable = column_temp
	[../]
	[./columnConduction]
		type = GColumnHeatDispersion
		variable = column_temp
	[../]
	[./column_AdsHeat]
		solid_heats = 'Kr_AdsorbedHeat Xe_AdsorbedHeat He_AdsorbedHeat'
		type = AdsorptionHeatAccumulation
		variable = column_temp
	[../]
	[./wallAccum]
		type = WallHeatAccumulation
		variable = wall_temp
	[../]
	[./wall_amb_trans]
		coupled = ambient_temp
		type = WallAmbientHeatTransfer
		variable = wall_temp
	[../]
	[./wall_bed_trans]
		coupled = column_temp
		type = BedWallHeatTransfer
		variable = wall_temp
	[../]
[] #END Kernels

[DGKernels]
	[./DGcolumnAdvection]
		type = DGColumnHeatAdvection
		variable = column_temp
	[../]
	[./DGcolumnConduction]
		type = DGColumnHeatDispersion
		variable = column_temp
	[../]
	[./He_DGadv]
		type = DGColumnMassAdvection
		variable = He
	[../]
	[./He_DGdiff]
		index = 2
		type = DGColumnMassDispersion
		variable = He
	[../]
	[./Kr_DGadv]
		type = DGColumnMassAdvection
		variable = Kr
	[../]
	[./Kr_DGdiff]
		index = 0
		type = DGColumnMassDispersion
		variable = Kr
	[../]
	[./Xe_DGadv]
		type = DGColumnMassAdvection
		variable = Xe
	[../]
	[./Xe_DGdiff]
		index = 1
		type = DGColumnMassDispersion
		variable = Xe
	[../]
[] #END DGKernels

[AuxKernels]
	[./He_ads]
		index = 2
		type = Scopsowl_Adsorption
		variable = He_Adsorbed
	[../]
	[./He_ads_heat]
		index = 2
		solid_conc = He_Adsorbed
		type = MAGPIE_AdsorptionHeat
		variable = He_AdsorbedHeat
	[../]
	[./He_pert]
		index = 2
		type = Scopsowl_Perturbation
		variable = He_Perturb
	[../]
	[./Kr_ads]
		index = 0
		type = Scopsowl_Adsorption
		variable = Kr_Adsorbed
	[../]
	[./Kr_ads_heat]
		index = 0
		solid_conc = Kr_Adsorbed
		type = MAGPIE_AdsorptionHeat
		variable = Kr_AdsorbedHeat
	[../]
	[./Kr_pert]
		index = 0
		type = Scopsowl_Perturbation
		variable = Kr_Perturb
	[../]
	[./Xe_ads]
		index = 1
		type = Scopsowl_Adsorption
		variable = Xe_Adsorbed
	[../]
	[./Xe_ads_heat]
		index = 1
		solid_conc = Xe_Adsorbed
		type = MAGPIE_AdsorptionHeat
		variable = Xe_AdsorbedHeat
	[../]
	[./Xe_pert]
		index = 1
		type = Scopsowl_Perturbation
		variable = Xe_Perturb
	[../]
	[./column_pressure]
		coupled_gases = 'Kr Xe He'
		temperature = column_temp
		type = TotalColumnPressure
		variable = total_pressure
	[../]
[] #END AuxKernels

[BCs]
	[./He_Flux]
		boundary = 'top bottom'
		index = 2
		input_molefraction = 0.999005101
		input_pressure = 101.35
		input_temperature = 253.15
		type = DGMassFluxLimitedBC
		variable = He
	[../]
	[./Heat_Gas_Flux]
		boundary = 'top bottom'
		input_temperature = 253.15
		type = DGHeatFluxLimitedBC
		variable = column_temp
	[../]
	[./Heat_Wall_Flux]
		boundary = 'right left'
		type = DGColumnWallHeatFluxLimitedBC
		variable = column_temp
		wall_temp = wall_temp
	[../]
	[./Kr_Flux]
		boundary = 'top bottom'
		index = 0
		input_molefraction = 0.000131792
		input_pressure = 101.35
		input_temperature = 253.15
		type = DGMassFluxLimitedBC
		variable = Kr
	[../]
	[./Xe_Flux]
		boundary = 'top bottom'
		index = 1
		input_molefraction = 0.000863107
		input_pressure = 101.35
		input_temperature = 253.15
		type = DGMassFluxLimitedBC
		variable = Xe
	[../]
[] #END BCs

[Materials]
	[./AdsorbateMaterials]
		block = 0
		coupled_gases = 'Kr Xe He'
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
		maximum_capacity = '1.716 1.479 0'
		molar_volume = '20.785 25.412 0'
		number_sites = '2 3 0'
		temperature = column_temp
		total_pressure = total_pressure
		type = MagpieAdsorbateProperties
	[../]
	[./AdsorbentMaterials]
		activation_energy = '0 0 0'
		affinity = '0 0 0'
		binder_porosity = 0.384
		block = 0
		coupled_gases = 'Kr Xe He'
		macropore_radius = 1.5e-4
		pellet_density = 3.06
		pellet_diameter = 0.16
		pellet_heat_capacity = 1.045
		ref_diffusion = '0 0 0'
		ref_temperature = '0 0 0'
		temperature = column_temp
		type = AdsorbentProperties
	[../]
	[./BedMaterials]
		axial_conductivity = 6.292E-05
		block = 0
		bulk_porosity = 0.885
		coupled_gases = 'Kr Xe He'
		extern_heat_trans_coef = 90.0
		inner_diameter = 1.7272
		outer_diameter = 1.905
		temperature = column_temp
		type = BedProperties
		wall_density = 7.7
		wall_heat_capacity = 0.5
		wall_heat_trans_coef = 9.0
	[../]
	[./FlowMaterials]
		block = 0
		comp_Sutherland_const = '266.505 232.746 80.0'
		comp_heat_capacity = '0.25 0.16 5.1916'
		comp_ref_temp = '273.15 273.15 273.15'
		comp_ref_viscosity = '0.00023219 0.00021216 0.0001885'
		coupled_adsorption = 'Kr_Adsorbed Xe_Adsorbed He_Adsorbed'
		coupled_gases = 'Kr Xe He'
		coupled_perturbation = 'Kr_Perturb Xe_Perturb He_Perturb'
		flow_rate = 2994.06
		molecular_weight = '83.8 131.29 4.0026'
		temperature = column_temp
		total_pressure = total_pressure
		type = FlowProperties
	[../]
	[./KineticMaterials]
		block = 0
		coupled_adsorption = 'Kr_Adsorbed Xe_Adsorbed He_Adsorbed'
		coupled_gases = 'Kr Xe He'
		dirichlet_bc = false
		heterogeneous = false
		macro_length = 0.4
		macro_spheres = false
		surface_diffusion = true
		type = ScopsowlProperties
	[../]
[] #END Materials

[Postprocessors]
	[./He_exit]
		boundary = 'top'
		execute_on = 'initial timestep_end'
		type = SideAverageValue
		variable = He
	[../]
	[./Kr_adsorbed]
		execute_on = 'initial timestep_end'
		type = ElementAverageValue
		variable = Kr_Adsorbed
	[../]
	[./Kr_exit]
		boundary = 'top'
		execute_on = 'initial timestep_end'
		type = SideAverageValue
		variable = Kr
	[../]
	[./Xe_adsorbed]
		execute_on = 'initial timestep_end'
		type = ElementAverageValue
		variable = Xe_Adsorbed
	[../]
	[./Xe_exit]
		boundary = 'top'
		execute_on = 'initial timestep_end'
		type = SideAverageValue
		variable = Xe
	[../]
	[./pressure_exit]
		boundary = 'top'
		execute_on = 'initial timestep_end'
		type = SideAverageValue
		variable = total_pressure
	[../]
	[./temp_exit]
		boundary = 'top'
		execute_on = 'initial timestep_end'
		type = SideAverageValue
		variable = column_temp
	[../]
	[./wall_temp]
		boundary = 'right'
		execute_on = 'initial timestep_end'
		type = SideAverageValue
		variable = wall_temp
	[../]
[] #END Postprocessors

[Executioner]
	dtmax = 1.0
	end_time = 50.0
	l_max_its = 100
	l_tol = 1e-6
	line_search = none
	nl_abs_step_tol = 1e-10
	nl_abs_tol = 1e-6
	nl_max_its = 10
	nl_rel_step_tol = 1e-10
	nl_rel_tol = 1e-6
	petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
	petsc_options_value = 'hypre boomeramg 100'
	scheme = implicit-euler
	solve_type = newton
	start_time = 0.0
	type = Transient
	[./TimeStepper]
		dt = 0.01
		type = SolutionTimeAdaptiveDT
	[../]
[] #END Executioner

[Preconditioning]
[] #END Preconditioning

[Outputs]
	csv = true
	exodus = true
	print_linear_residuals = true
[] #END Outputs

