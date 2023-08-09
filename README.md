## Temperature-sensitive enzyme network (TSEN) model package
This is a simulation and analysis package for TSEN models described in Knapp et al. (2023) (https://www.biorxiv.org/content/10.1101/2023.07.22.550177v1). 

The TSEN model describes a closed autocatalytic enzyme reaction network composed of chained, fully temperature-sensitive Michaelis-Menten reactions (https://en.wikipedia.org/wiki/Michaelis%E2%80%93Menten_kinetics), which imports metabolites through the system's cellular envelope to produce a final growth reaction. 

The TSEN scheme is composed of reactions encompassing (1) Import, (2) Production, and (3) Growth

### 1. Generalized linear TSEN model (TSEN_Linear_Generalized)
This TSEN model is a linear reaction network (i.e., without branching) generalized to incorporate any number of intermediate production reactions, but with a single import reaction and growth reaction.
