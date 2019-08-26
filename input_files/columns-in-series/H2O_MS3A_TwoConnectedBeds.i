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
	ny = 40
	xmin = 0.0
	xmax = 1.27 #cm
	ymin = 0.0
	ymax = 12.7 #cm
 
 [] # END Mesh
 
[MeshModifiers]
 
	[./column_one]
		type = SubdomainBoundingBox
		block_id = 0
		top_right = '1.27 6.35 0'
		bottom_left = '0 0 0'
	[../]
 
	[./column_two]
		type = SubdomainBoundingBox
		block_id = 1
		top_right = '1.27 12.7 0'
		bottom_left = '0 6.35 0'
	[../]
 
	[./interface]
		type = SideSetsBetweenSubdomains
		depends_on = 'column_one column_two'
		master_block = 0
		paired_block = 1
		new_boundary = 'interface'
 [../]
 
 [] #END MeshModifiers
 
[Variables]
 
	[./N2]
		block = '0 1'
		order = FIRST
		family = MONOMIAL
	[../]
 
	[./O2]
		block = '0 1'
		order = FIRST
		family = MONOMIAL
	[../]
 
	[./H2O]
		block = '0 1'
		order = FIRST
		family = MONOMIAL
	[../]
 
	[./column_temp]
		block = '0 1'
		order = FIRST
		family = MONOMIAL
		initial_condition = 298.15
	[../]
 
	[./H2O_Adsorbed]
		block = '0 1'
		order = FIRST
		family = MONOMIAL
		initial_condition = 0.0
	[../]
 
	[./H2O_AdsorbedHeat]
		block = '0 1'
		order = FIRST
		family = MONOMIAL
		initial_condition = 0.0
	[../]
 
 
 [] #END Variables
 
[AuxVariables]
 
	[./total_pressure]
		block = '0 1'
		order = FIRST
		family = MONOMIAL
		initial_condition = 101.35
	[../]
 
	[./ambient_temp]
		block = '0 1'
		order = FIRST
		family = MONOMIAL
		initial_condition = 298.15
	[../]
 
	[./wall_temp]
		block = '0 1'
		order = FIRST
		family = MONOMIAL
		initial_condition = 298.15
	[../]
 
 [] #END AuxVariables
 
[ICs]
 
	[./N2_IC]
		block = '0 1'
		type = ConcentrationIC
		variable = N2
		initial_mole_frac = 0.79
		initial_press = 101.35
		initial_temp = 298.15
	[../]
 
	[./O2_IC]
		block = '0 1'
		type = ConcentrationIC
		variable = O2
		initial_mole_frac = 0.21
		initial_press = 101.35
		initial_temp = 298.15
	[../]
 
	[./H2O_IC]
		block = '0 1'
		type = ConcentrationIC
		variable = H2O
		initial_mole_frac = 0.0
		initial_press = 101.35
		initial_temp = 298.15
	[../]
 
 [] #END ICs
 
[Kernels]
 
	[./accumN2]
		block = '0 1'
		type = BedMassAccumulation
		variable = N2
	[../]
 
	[./diffN2]
		block = '0 1'
		type = GColumnMassDispersion
		variable = N2
		index = 0
	[../]
 
	[./advN2]
		block = '0 1'
		type = GColumnMassAdvection
		variable = N2
	[../]
 
	[./accumO2]
		block = '0 1'
		type = BedMassAccumulation
		variable = O2
	[../]
 
	[./diffO2]
		block = '0 1'
		type = GColumnMassDispersion
		variable = O2
		index = 1
	[../]
 
	[./advO2]
		block = '0 1'
		type = GColumnMassAdvection
		variable = O2
	[../]
 
	[./accumH2O]
		block = '0 1'
		type = BedMassAccumulation
		variable = H2O
	[../]
 
	[./H2O_MT]
		block = '0 1'
		type = SolidMassTransfer
		variable = H2O
		coupled = H2O_Adsorbed
	[../]
 
	[./diffH2O]
		block = '0 1'
		type = GColumnMassDispersion
		variable = H2O
		index = 2
	[../]
 
	[./advH2O]
		block = '0 1'
		type = GColumnMassAdvection
		variable = H2O
	[../]
 
	[./columnAccum]
		block = '0 1'
		type = BedHeatAccumulation
		variable = column_temp
	[../]
 
	[./columnConduction]
		block = '0 1'
		type = GColumnHeatDispersion
		variable =column_temp
	[../]
 
	[./columnAdvection]
		block = '0 1'
		type = GColumnHeatAdvection
		variable =column_temp
	[../]
 
	[./H2O_columnAdsHeat]
		block = '0 1'
		type = SolidHeatTransfer
		variable = column_temp
		coupled = H2O_AdsorbedHeat
	[../]
 
	[./H2O_adsheat]
		block = '0 1'
		type = HeatofAdsorption
		variable = H2O_AdsorbedHeat
		coupled = H2O_Adsorbed
		index = 2
	[../]
 
#	[./H2O_adsorption]
#		block = '0 1'
#		type = CoupledLangmuirForcingFunction
#		variable = H2O_Adsorbed
#		coupled = H2O
#		langmuir_coeff = 100.0
#		max_capacity = 1.17
#	[../]
 
#	[./H2O_adsorption]
#		block = '0 1'
#		type = CoupledGSTAisotherm
#		variable = H2O_Adsorbed
#		coupled_gas = H2O
#		coupled_temp = column_temp
#		max_capacity = 11.67
#		num_sites = 4
#		gsta_params = '228357.3949 22688965955 1.93815E+15 1.1268E+18'
#	[../]
 
    [./H2O_adsorption]
		block = '0 1'
        type = CoupledGSTALDFmodel
        variable = H2O_Adsorbed
        coupled_gas = H2O
        coupled_temp = column_temp
        index = 2
	[../]
 
 [] #END Kernels
 
[DGKernels]
 
	[./dg_disp_N2]
		block = '0 1'
		type = DGColumnMassDispersion
		variable = N2
		index = 0
	[../]
 
	[./dg_adv_N2]
		block = '0 1'
		type = DGColumnMassAdvection
		variable = N2
	[../]
 
	[./dg_disp_O2]
		block = '0 1'
		type = DGColumnMassDispersion
		variable = O2
		index = 1
	[../]
 
	[./dg_adv_O2]
		block = '0 1'
		type = DGColumnMassAdvection
		variable = O2
	[../]
 
	[./dg_disp_H2O]
		block = '0 1'
		type = DGColumnMassDispersion
		variable = H2O
		index = 2
	[../]
 
	[./dg_adv_H2O]
		block = '0 1'
		type = DGColumnMassAdvection
		variable = H2O
	[../]
 
	[./dg_disp_heat]
		block = '0 1'
		type = DGColumnHeatDispersion
		variable = column_temp
	[../]
 
	[./dg_adv_heat]
		block = '0 1'
		type = DGColumnHeatAdvection
		variable = column_temp
	[../]
 
 [] #END DGKernels
 
[AuxKernels]
 
	[./column_pressure]
		block = '0 1'
		type = TotalColumnPressure
		variable = total_pressure
		temperature = column_temp
		coupled_gases = 'N2 O2 H2O'
		execute_on = 'initial timestep_end'
	[../]
 
	[./wall_temp_calc]
		block = '0 1'
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
		input_temperature = 298.15
		input_pressure = 101.35
		input_molefraction = 0.78863
		index = 0
	[../]

	[./O2_Flux]
		type = DGMassFluxBC
		variable = O2
		boundary = 'top bottom'
		input_temperature = 298.15
		input_pressure = 101.35
		input_molefraction = 0.20974
		index = 1
	[../]

	[./H2O_Flux]
		type = DGMassFluxBC
		variable = H2O
		boundary = 'top bottom'
		input_temperature = 298.15
		input_pressure = 101.35
		input_molefraction = 0.00163
		index = 2
	[../]

	[./Heat_Gas_Flux]
		type = DGHeatFluxBC
		variable = column_temp
		boundary = 'top bottom'
		input_temperature = 298.15
	[../]
 
	[./Heat_Wall_Flux]
		type = DGColumnWallHeatFluxLimitedBC
		variable = column_temp
		boundary = 'right left'
		wall_temp = wall_temp
	[../]
 
 [] #END BCs
 
[Materials]
 
	[./BedMaterials_One]
        type = BedProperties
        block = 0
		length = 6.35
		inner_diameter = 2.54
        outer_diameter = 2.84
        bulk_porosity = 0.541
        wall_density = 8.0
        wall_heat_capacity = 0.5
        wall_heat_trans_coef = 6.12
        extern_heat_trans_coef = 6.12
	[../]
 
	[./FlowMaterials_One]
        type = GasFlowProperties
        block = 0
		flow_rate = 211680.0
        molecular_weight = '28.016 32 18'
        comp_heat_capacity = '1.04 0.919 1.97'
        comp_ref_viscosity = '0.0001781 0.0002018 0.0001043'
        comp_ref_temp = '300.55 292.25 298.16'
        comp_Sutherland_const = '111 127 784.72'
        temperature = column_temp
        total_pressure = total_pressure
        coupled_gases = 'N2 O2 H2O'
	[../]
 
	[./AdsorbentMaterials_One]
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
 
	[./AdsorbateMaterials_One]
        type = ThermodynamicProperties
        block = 0
        temperature = column_temp
        total_pressure = total_pressure
        coupled_gases = 'N2 O2 H2O'
        number_sites = '0 0 4'
        maximum_capacity = '0 0 11.67' #mol/kg
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
 
 
	[./BedMaterials_Two]
        type = BedProperties
        block = 1
		length = 6.35
		inner_diameter = 2.54
        outer_diameter = 2.84
        bulk_porosity = 0.541
        wall_density = 8.0
        wall_heat_capacity = 0.5
        wall_heat_trans_coef = 6.12
        extern_heat_trans_coef = 6.12
	[../]
 
	[./FlowMaterials_Two]
        type = GasFlowProperties
        block = 1
		flow_rate = 211680.0
        molecular_weight = '28.016 32 18'
        comp_heat_capacity = '1.04 0.919 1.97'
        comp_ref_viscosity = '0.0001781 0.0002018 0.0001043'
        comp_ref_temp = '300.55 292.25 298.16'
        comp_Sutherland_const = '111 127 784.72'
        temperature = column_temp
        total_pressure = total_pressure
        coupled_gases = 'N2 O2 H2O'
	[../]
 
	[./AdsorbentMaterials_Two]
        type = AdsorbentProperties
        block = 1
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
 
	[./AdsorbateMaterials_Two]
        type = ThermodynamicProperties
        block = 1
        temperature = column_temp
        total_pressure = total_pressure
        coupled_gases = 'N2 O2 H2O'
        number_sites = '0 0 4'
        maximum_capacity = '0 0 11.67' #mol/kg
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
 

 [] #END Materials
 
[Postprocessors]
 
	[./H2O_enter]
		type = SideAverageValue
		boundary = 'bottom'
		variable = H2O
		execute_on = 'initial timestep_end'
	[../]
 
	[./H2O_mid]
		type = SideAverageValue
		boundary = 'interface'
		variable = H2O
		execute_on = 'initial timestep_end'
	[../]
 
	[./H2O_avg_gas_one]
		block = 0
		type = ElementAverageValue
		variable = H2O
		execute_on = 'initial timestep_end'
	[../]
 
	[./H2O_avg_gas_two]
		block = 1
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
 
	[./temp_mid]
		type = SideAverageValue
		boundary = 'interface'
		variable = column_temp
		execute_on = 'initial timestep_end'
	[../]
 
   [./wall_temp]
       type = SideAverageValue
       boundary = 'right'
       variable = wall_temp
       execute_on = 'initial timestep_end'
   [../]
 
	[./H2O_solid_one]
		block = 0
		type = ElementAverageValue
		variable = H2O_Adsorbed
		execute_on = 'initial timestep_end'
	[../]
 
	[./H2O_solid_two]
		block = 1
		type = ElementAverageValue
		variable = H2O_Adsorbed
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
	nl_max_its = 50
 
	solve_type = pjfnk
	line_search = basic    # Options: default none l2 bt
	start_time = 0.0
	end_time = 90.0
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
