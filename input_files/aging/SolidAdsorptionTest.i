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
        initial_condition = 1.0
    [../]
 
    [./column_temp]
        order = FIRST
        family = MONOMIAL
        initial_condition = 423.15
    [../]
 
    [./Ag2O]
        order = FIRST
        family = MONOMIAL
        initial_condition = 0.3
    [../]
 
    [./AgNO3]
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
        initial_mole_frac = 1.0
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
 
    [./Ag2O_MT]
        type = CoefTimeDerivative
        variable = Ag2O
    [../]
 
    [./NO_Ag_aging]
        type = VariableOrderCoupledSolid
        variable = Ag2O
        main_variable = Ag2O
        main_order = 2
        main_stoichiometry = 2
        coupled_species = ''
        species_stoichiometry = ''
        species_order = ''
        coupled_catalysts = ''
        catalyst_order = ''
        coupled_adsorption = 'AgNO3'
        adsorbed_stoichiometry = '-2'
        adsorbed_order = '2'
        coupled_all_adsorption = 'Ag2O AgNO3'
        adsorbed_sites = '2 1'
        site_order = 1
        max_capacity = 1.1102
        forward_rate = 1.0e-2
        reverse_rate = 1.0e-2
    [../]
 
    [./NO_MT]
        type = CoefTimeDerivative
        variable = NO
    [../]
 
    [./AgNO3_MT]
        type = CoupledCoeffTimeDerivative
        variable = AgNO3
        coupled = Ag2O
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
        coupled_gases = 'N2'
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
        input_molefraction = 1.0
        index = 0
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
        molecular_weight = '28.016'
        comp_heat_capacity = '1.04'
        comp_ref_viscosity = '0.0001781'
        comp_ref_temp = '300.55'
        comp_Sutherland_const = '111'
        temperature = column_temp
        total_pressure = total_pressure
        coupled_gases = 'N2'
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
        ref_diffusion = '0'
        activation_energy = '0'
        ref_temperature = '0'
        affinity = '0'
        temperature = column_temp
        coupled_gases = 'N2'
    [../]
 
    [./AdsorbateMaterials]
        type = ThermodynamicProperties
        block = 0
        temperature = column_temp
        total_pressure = total_pressure
        coupled_gases = 'N2'
        number_sites = '0'
        maximum_capacity = '0' #mol/kg
        molar_volume = '0' #cm^3/mol
        enthalpy_site_1 = '0'
        enthalpy_site_2 = '0'
        enthalpy_site_3 = '0'
        enthalpy_site_4 = '0'
        enthalpy_site_5 = '0'
        enthalpy_site_6 = '0'
 
        entropy_site_1 = '0'
        entropy_site_2 = '0'
        entropy_site_3 = '0'
        entropy_site_4 = '0'
        entropy_site_5 = '0'
        entropy_site_6 = '0'
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
 
    [./NO_column]
        type = ElementAverageValue
        variable = 'NO'
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
    nl_max_its = 80
 
    solve_type = pjfnk
    line_search = bt    # Options: default none l2 bt basic
    start_time = 0.0
    end_time = 1344.0
    dtmax = 100.0
 
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
