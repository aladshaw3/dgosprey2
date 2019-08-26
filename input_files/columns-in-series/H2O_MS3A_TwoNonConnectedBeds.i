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
	ny = 60
	xmin = 0.0
	xmax = 1.27 #cm
	ymin = 0.0
	ymax = 30.4 #cm
 
 [] # END Mesh
 
[MeshModifiers]
 
	[./column_1]
		type = SubdomainBoundingBox
		block_id = 0
		top_right = '1.27 12.7 0'
		bottom_left = '0 0 0'
	[../]
 
	[./conduit_12]
		type = SubdomainBoundingBox
		block_id = 1
		top_right = '1.27 17.7 0'
		bottom_left = '0 12.7 0'
	[../]
 
	[./column_2]
		type = SubdomainBoundingBox
		block_id = 2
		top_right = '1.27 30.4 0'
		bottom_left = '0 17.7 0'
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
 
	[./column_temp]
		block = '0 1 2'
		order = FIRST
		family = MONOMIAL
		initial_condition = 298.15
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
		initial_condition = 298.15
	[../]
 
	[./wall_temp_1]
		block = '0'
		order = FIRST
		family = MONOMIAL
		initial_condition = 298.15
	[../]
 
	[./wall_temp_12]
		block = '1'
		order = FIRST
		family = MONOMIAL
		initial_condition = 298.15
	[../]
 
	[./wall_temp_2]
		block = '2'
		order = FIRST
		family = MONOMIAL
		initial_condition = 298.15
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
		initial_temp = 298.15
	[../]
 
	[./O2_IC]
		block = '0 1 2'
		type = ConcentrationIC
		variable = O2
		initial_mole_frac = 0.21
		initial_press = 101.35
		initial_temp = 298.15
	[../]
 
	[./H2O_IC]
		block = '0 1 2'
		type = ConcentrationIC
		variable = H2O
		initial_mole_frac = 0.0
		initial_press = 101.35
		initial_temp = 298.15
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
 
	[./H2O_adsheat]
		block = '0 2'
		type = HeatofAdsorption
		variable = H2O_AdsorbedHeat
		coupled = H2O_Adsorbed
		index = 2
	[../]
 
#	[./H2O_adsorption]
#		block = '0 2'
#		type = CoupledLangmuirForcingFunction
#		variable = H2O_Adsorbed
#		coupled = H2O
#		langmuir_coeff = 100.0
#		max_capacity = 1.17
#	[../]
 
#	[./H2O_adsorption]
#		block = '0 2'
#		type = CoupledGSTAisotherm
#		variable = H2O_Adsorbed
#		coupled_gas = H2O
#		coupled_temp = column_temp
#		max_capacity = 11.67
#		num_sites = 4
#		gsta_params = '228357.3949 22688965955 1.93815E+15 1.1268E+18'
#	[../]
 
    [./H2O_adsorption]
		block = '0 2'
        type = CoupledGSTALDFmodel
        variable = H2O_Adsorbed
        coupled_gas = H2O
        coupled_temp = column_temp
        index = 2
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
		length = 6.35
		inner_diameter = 2.54
        outer_diameter = 2.84
        bulk_porosity = 0.541
        wall_density = 8.0
        wall_heat_capacity = 0.5
        wall_heat_trans_coef = 6.12
        extern_heat_trans_coef = 6.12
	[../]
 
	[./FlowMaterials_1]
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
        ref_diffusion = '0 0 0.8814'
        activation_energy = '0 0 0'
        ref_temperature = '0 0 267.999'
        affinity = '0 0 0'
        temperature = column_temp
        coupled_gases = 'N2 O2 H2O'
	[../]
 
	[./AdsorbateMaterials_1]
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
 
# Properties for the conduit between beds
	[./BedMaterials_12]
		type = BedProperties
		block = 1
		length = 5.0
		inner_diameter = 2.54
		outer_diameter = 2.84
		bulk_porosity = 1.0
		wall_density = 12.0
		wall_heat_capacity = 1.5
		wall_heat_trans_coef = 9.12
		extern_heat_trans_coef = 9.12
	[../]
 
	[./FlowMaterials_12]
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
		ref_diffusion = '0 0 0.0'
		activation_energy = '0 0 0'
		ref_temperature = '0 0 0.0'
		affinity = '0 0 0'
		temperature = column_temp
		coupled_gases = 'N2 O2 H2O'
	[../]
 
	[./AdsorbateMaterials_12]
		type = ThermodynamicProperties
		block = 1
		temperature = column_temp
		total_pressure = total_pressure
		coupled_gases = 'N2 O2 H2O'
		number_sites = '0 0 0'
		maximum_capacity = '0 0 0' #mol/kg
		molar_volume = '0 0 0' #cm^3/mol
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
 
 
# Approximate Properties for Ag0Z bed
	[./BedMaterials_2]
        type = BedProperties
        block = 2
		length = 6.35
		inner_diameter = 2.54
        outer_diameter = 2.84
		bulk_porosity = 0.888  #not known
        wall_density = 8.0
        wall_heat_capacity = 0.5
        wall_heat_trans_coef = 6.12
        extern_heat_trans_coef = 6.12
	[../]
 
	[./FlowMaterials_2]
        type = GasFlowProperties
        block = 2
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
 
	[./AdsorbentMaterials_2]
        type = AdsorbentProperties
        block = 2
        binder_fraction = 0.0
        binder_porosity = 0.384
        crystal_radius = 0.0
		pellet_diameter = 0.16
        macropore_radius = 2.6e-6
        pellet_density = 3.06
		pellet_heat_capacity = 1.5 #not known
        ref_diffusion = '0 0 0.0'
        activation_energy = '0 0 0'
        ref_temperature = '0 0 0.0'
        affinity = '0 0 0'
        temperature = column_temp
        coupled_gases = 'N2 O2 H2O'
	[../]
 
	[./AdsorbateMaterials_2]
        type = ThermodynamicProperties
        block = 2
        temperature = column_temp
        total_pressure = total_pressure
        coupled_gases = 'N2 O2 H2O'
        number_sites = '0 0 3'
        maximum_capacity = '0 0 4.589' #mol/kg
        molar_volume = '0 0 13.91' #cm^3/mol
        enthalpy_site_1 = '0 0 -46991.5'
        enthalpy_site_2 = '0 0 -111956.3'
        enthalpy_site_3 = '0 0 -144189.7'
        enthalpy_site_4 = '0 0 0'
        enthalpy_site_5 = '0 0 0'
        enthalpy_site_6 = '0 0 0'
 
        entropy_site_1 = '0 0 -46.88'
        entropy_site_2 = '0 0 -179.86'
        entropy_site_3 = '0 0 -231.92'
        entropy_site_4 = '0 0 0'
        entropy_site_5 = '0 0 0'
        entropy_site_6 = '0 0 0'
	[../]
 

 [] #END Materials
 
[Postprocessors]
 
	[./H2O_enter_1]
		type = SideAverageValue
		boundary = 'bottom'
		variable = H2O
		execute_on = 'initial timestep_end'
	[../]
 
	[./H2O_enter_2]
		type = SideAverageValue
		boundary = 'interface_122'
		variable = H2O
		execute_on = 'initial timestep_end'
	[../]
 
	[./H2O_exit_1]
		type = SideAverageValue
		boundary = 'interface_112'
		variable = H2O
		execute_on = 'initial timestep_end'
	[../]
 
	[./H2O_avg_gas_1]
		block = 0
		type = ElementAverageValue
		variable = H2O
		execute_on = 'initial timestep_end'
	[../]
 
	[./H2O_avg_gas_12]
		block = 1
		type = ElementAverageValue
		variable = H2O
		execute_on = 'initial timestep_end'
	[../]
 
	[./H2O_avg_gas_2]
		block = 2
		type = ElementAverageValue
		variable = H2O
		execute_on = 'initial timestep_end'
	[../]
 
	[./H2O_exit_2]
		type = SideAverageValue
		boundary = 'top'
		variable = H2O
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
 
	[./H2O_solid_1]
		block = 0
		type = ElementAverageValue
		variable = H2O_Adsorbed
		execute_on = 'initial timestep_end'
	[../]
 
	[./H2O_solid_2]
		block = 2
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
	nl_max_its = 30
 
	solve_type = pjfnk
	line_search = basic    # Options: default none l2 bt
	start_time = 0.0
	end_time = 90.0
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
		petsc_options_iname = '-pc_type -ksp_gmres_restart  -snes_max_funcs'
		petsc_options_value = 'lu 2000 20000'
	[../]
 
 [] #END Preconditioning
 
[Outputs]
 
	exodus = true
	csv = true
	print_linear_residuals = false
 
 [] #END Outputs
