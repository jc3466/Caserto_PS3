# Problem set 3 CHEME 7770
# Flux.jl provided by Prof. Varner
# Urea cycle in human cells
# Flux Balance Analysis (FBA)

include("Flux.jl")  # Julia script that optimizes flux distribution
using LinearAlgebra
tdouble = 20 # human cell doubling time (hrs), give in problem statement

# stoichiometric array for the metablic network: Urea Cycle
stoichiometric_array = [
-1	0.0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0
-1	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
-1	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0
1	0	0	0	0	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0
1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	0
1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	1	0	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0	0
0	1	-1	0	-1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1
0	0	1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	1	0	0	0	0	0	0	-1	0	0	0	0	0	0	0	0	0	0	0
0	0	0	-1	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	-1	0	0	0	0
0	0	0	0	-1	1	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0
0	0	0	0	1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	-1	0
0	0	0	0	-1	1	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0
0	0	0	0	-1	1	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0
0	0	0	0	1	-1	0	0	0	0	0	0	0	-1	0	0	0	0	0	0	0
0	0	-1	0	1	-1	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0
]


# kcat values provided in problem statement
kcat1 = 203.0 # sec^-1
kcat2 = 34.5 # sec^-1
kcat3 = 249.0 # sec^-1
kcat4 = 88.1 # sec^-1
kcat5 = 13.7 # sec^-1
E = 0.01e-3 # mmol/gDW, given in problem statement

water_fraction = 0.70 # Wet % of total weight in HeLa cell, Bionumbers ID: 100387
mass_cell = 2.3E-9 # mass of HeLa cell  in grams, Bionumbers ID: 103720
vol_cell = 3.7E-12 # volume of Hela cell in Liters, Bionumbers ID: 105879
E1 = (E * (mass_cell / (vol_cell*water_fraction)))  # mmol/L

# Km values from Park et al. units are mmol/L
# Supplementary Data Set 2, values for homo sapiens
Km_v1_asp = 10.0   # aspartate, 1.00e-2 M
Km_v1_ATP = 6.14e-1 # 6.14e-4 M
Km_v3_arg = 2.5 # arginine, 2.50e-3 M
Km_v5_arg = 2.5   # arginine, 2.50e-3 M
Km_v5_NADPH = 5.25e-3  # 5.25e-6 M
# Could not find values for homo sapiens for:
# citrulline, argininosuccinate, carbamoyl phosphate, ornithine


# METABOLITE CONCENTRATIONS from Park et al. paper
# Comprehensive absolute cellular metabolite concentrations, Mammalian iBMK
# Supplementary Text and Figures
# Supplementary Data Set 2
# alL concentrations are in mmol/L
# all values are for homo sapiens, because trying to model urea cycle in human cells

# v1 reaction, intracellular
ATP = 4.67  # 4.67e-3 M
aspartate = 14.9  # 1.49e-2 M
# didn't provide value for citrulline

# v2 reaction, intracellular
# value for homo sapiens not given in park et al.
#argin_succin =         # argininosuccinate

# v3 reaction, intracellular
arginine = 2.55e-1  # 2.55e-4 M

# v4 reaction, intracellular
# values for homo sapiens not given in park et al.
# carb_phosph =
# ornithine =

# v5 reaction, intracellular
# arginine concentration already defined in reaction 3
NADPH = 6.54e-2  # 6.54e-5 M


# Saturation, using substrates (reactants) for each reaction in the metabolic network
# general form: x/(x+Km) for each substrate, from lecture on 2-27-20
# for saturations in which metabolite concentrations were not found in Park et al. ---> Assumed saturation = 1
saturation_1 = (aspartate/(aspartate+Km_v1_asp))*(ATP/(ATP+Km_v1_ATP))
# saturation_2 = (argsuccinate/(argsuccinate+Km_v2_argsuccinate))  # could not find values in Park et al.
saturation_3 = (arginine/(arginine+Km_v3_arg))
# saturation_4 = (ornithine/(ornithine+Km_v4_orn)) # could not find values in Park et al.
saturation_5 = (arginine/(arginine+Km_v5_arg))*(NADPH/(NADPH+Km_v5_NADPH))

# upper bounds using Vmax = kcat*E*saturation for intracellular reactions
# upper bound = 10 mmol/gDW-hr for all external fluxes, b
# need to convert units to mmol/L/s for b upper bound in order to have consistent units in all bound constraints
upper_bound_1 = kcat1 * saturation_1 * E1  # v1 reactions
upper_bound_2 = kcat2 * E1     # v2 reaction
upper_bound_3 = kcat3 * saturation_3 * E1   # v3 reaction
upper_bound_4 = kcat4  * E1   # v4 reaction
upper_bound_5 = kcat5 * saturation_5* E1    # v5 reaction
upper_bound_b = (10 * (mass_cell / (vol_cell * water_fraction)))/ 3600 # mmol/L/s   # external fluxes, upper bound given in problem statement 0<b<10 mmol/gDWl/hr


# dimensions of default_bounds_array are (rxnsx2), 21 fluxes in stoich array so (21x2)
Default_bounds_array = [0.0 upper_bound_1;
 0 upper_bound_2;
 0 upper_bound_3;
 0 upper_bound_4;
 0 upper_bound_5;
 0 upper_bound_5;
 0 upper_bound_b;
 0 upper_bound_b;
 0 upper_bound_b;
 0 upper_bound_b;
 0 upper_bound_b;
 0 upper_bound_b;
 0 upper_bound_b;
 0 upper_bound_b;
 0 upper_bound_b;
 0 upper_bound_b;
 0 upper_bound_b;
 0 upper_bound_b;
 0 upper_bound_b;
 0 upper_bound_b;
 0 upper_bound_b
 ]

Species_bounds_array=zeros(18,2);  # array of zeros, 18 metabolites in stoich array

# maximize rate of Urea production, use the objective coefficient array, Urea leaves the cell so -1, Urea is 10th row of Stoich array
# (21x1) objective coefficient array
Objective_coefficient_array = [0.0; 0; 0; 0; 0; 0; 0; 0; 0; -1; 0; 0; 0; 0; 0; 0; 0; 0; 0;0;0]

optimize = calculate_optimal_flux_distribution(stoichiometric_array, Default_bounds_array, Species_bounds_array, Objective_coefficient_array)
r= optimize[2]
urea_flux=r[10] * 1/(mass_cell / (vol_cell * water_fraction)) * 3600  # need to convert back to mmol/gDW/hr
println("max rate of urea production  = ", urea_flux, " ", "mmol/gDW/hr")
