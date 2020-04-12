# Caserto_PS3
CHEME 7770 PS3

Urea Cycle Reconstruction 

This repository can be used to calculate the maximum rate of urea production in growing human cells. Cell mass, volume, and wet % values
were found for HeLa cells using Bionumbers. 

Flux.jl is used to optimize the flux distribution of the metabolic network, given the bound constraints and stoichiometric array. 
Flux.jl was provided by Professor Varner and is included in this repository. 

Metabolite concentrations and Km values were found in Park et al., and KEGG was used to construct the metabolic network
and the stoichiometric matrix S. Elemental balances were used to confirm that the intracellular reactions were balanced. 

To execute the code, enter the following command in Julia: 
include("PS3_partc.jl)

Maximum rate of urea production was calculated to be 0.83 mmol/gDW-hr. 

Further details regarding the construction of the metabolic network and determination of the bound constraints are provided 
in the pdf file named README. 
