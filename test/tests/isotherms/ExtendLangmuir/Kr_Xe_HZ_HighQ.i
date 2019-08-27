 [GlobalParams]

sigma = 1   # Penalty value:  NIPG = 0   otherwise, > 0
epsilon = 1  #  -1 = SIPG   0 = IIPG   1 = NIPG

flow_rate = 1.2e5
length = 50.8
inner_diameter = 1.905
pellet_diameter = 0.056
dt = 0.05

[] #END GlobalParams

[Problem]

coord_type = RZ

[] #END Problem

[Mesh]

type = GeneratedMesh
dim = 2
nx = 3
ny = 8
xmin = 0.0
xmax = 0.9525 #cm
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

[./He]
order = FIRST
family = MONOMIAL
[../]

[./column_temp]
order = FIRST
family = MONOMIAL
initial_condition = 191.15
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
initial_condition = 191.15
[../]

[./wall_temp]
order = FIRST
family = MONOMIAL
initial_condition = 191.15
[../]


[] #END AuxVariables

[ICs]

[./Kr_IC]
type = ConcentrationIC
variable = Kr
initial_mole_frac = 0.0
initial_press = 101.35
initial_temp = 191.15
[../]

[./Xe_IC]
type = ConcentrationIC
variable = Xe
initial_mole_frac = 0.0
initial_press = 101.35
initial_temp = 191.15
[../]

[./He_IC]
type = ConcentrationIC
variable = He
initial_mole_frac = 1.0
initial_press = 101.35
initial_temp = 191.15
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

[./Kr_adsorption]
type = CoupledExtendedLangmuirModel
variable = Kr_Adsorbed
main_coupled = Kr
coupled_temp = column_temp
coupled_list = 'Kr Xe'
enthalpies = '-17306 -22684'
entropies = '-23.4 -23.4'
max_capacity = 1.9375
[../]

[./Xe_adsorption]
type = CoupledExtendedLangmuirModel
variable = Xe_Adsorbed
main_coupled = Xe
coupled_temp = column_temp
coupled_list = 'Kr Xe'
enthalpies = '-17306 -22684'
entropies = '-23.4 -23.4'
max_capacity = 1.3666
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
execute_on = 'initial timestep_end'
variable = total_pressure
temperature = column_temp
coupled_gases = 'Kr Xe He'
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
input_temperature = 191.15
input_pressure = 101.35
input_molefraction = 0.00015
index = 0
[../]

[./Xe_Flux]
type = DGMassFluxBC
variable = Xe
boundary = 'top bottom'
input_temperature = 191.15
input_pressure = 101.35
input_molefraction = 0.001
index = 1
[../]

[./He_Flux]
type = DGMassFluxBC
variable = He
boundary = 'top bottom'
input_temperature = 191.15
input_pressure = 101.35
input_molefraction = 0.99885
index = 2
[../]

[./Heat_Gas_Flux]
type = DGHeatFluxBC
variable = column_temp
boundary = 'top bottom'
input_temperature = 191.15
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
outer_diameter = 2.0828
bulk_porosity = 0.90507				#not known
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
temperature = column_temp
total_pressure = total_pressure
coupled_gases = 'Kr Xe He'
[../]

[./AdsorbentMaterials]
type = AdsorbentProperties
block = 0
binder_fraction = 0.175				#not known
binder_porosity = 0.13				#not known
crystal_radius = 1.5				#not known
macropore_radius = 1.945e-7			#not Known
pellet_density = 2.174				#not Known
pellet_heat_capacity = 1.2  		#not known
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

[] #END Postprocessors

[Executioner]

type = Transient
scheme = bdf2

# NOTE: The default tolerances are far to strict and cause the program to crawl
nl_rel_tol = 1e-8
nl_abs_tol = 1e-4
l_tol = 1e-6
l_max_its = 2000
nl_max_its = 30

solve_type = pjfnk
line_search = basic    # Options: default none basic l2 bt
start_time = 0.0
end_time = 0.4
dtmax = 0.1

[./TimeStepper]
#Need to write a custom TimeStepper to enforce a maximum allowable dt
type = ConstantDT
#type = DGOSPREY_TimeStepper
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
