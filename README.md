## Temperature-sensitive enzyme network (TSEN) model package
This is a simulation and analysis package for TSEN models described in Knapp et al. (2023) (https://www.biorxiv.org/content/10.1101/2023.07.22.550177v1). 

The TSEN model describes a closed autocatalytic enzyme reaction network composed of chained, fully temperature-sensitive Michaelis-Menten reactions (https://en.wikipedia.org/wiki/Michaelis%E2%80%93Menten_kinetics), which imports metabolites through the system's cellular envelope to produce a final growth reaction. 

## TSEN model scheme
The TSEN scheme is composed of reactions encompassing (1) Import, (2) Production, and (3) Growth. The production rate, $k_{prod}$, for each reaction is described by $k_{prod} = \dfrac{k_{cat}(T) ec}{K_M(T)+c}$, where the catalytic rate ($k_{cat}$) and Michaelis-Menten constant ($K_M$) are both temperature sensitive according to an Arrhenius relationship ($e^{-E_a/k_BT}$) (each with a separate activation energy $E_a$). The substrate concentration, $c$, is consumed by an enzyme at concentration $e$. 

In chained Michealis-Menten reactions, the total production rate of an intermediate metabolite $c_i$ is the production rate from the previous reaction, which consumes $c_{i-1}$, minus the next reaction's production rate, which consumes $c_i$.

<p align="center">$\dot{c_i} = k_{prod}(c_{i-1}) - k_{prod}(c_i)$</p>

While this rate of production of $c_i$ is true under fixed volume, the TSEN system is contained within a growing envelope of size $V$, which expands at a growth rate $g = \dot{V}/V$. To account for this dilution

<p align="center">$\dot{c_i} = k_{prod}(c_{i-1}) - k_{prod}(c_i) - c_i g$</p>

The final growth rate reaction is defined as a Michaelis-Menten reaction, but with a constant growth efficiency factor ($\gamma_o$) that captures the final reaction's efficiency at converting the final reactant ($c_N$, network length $N$) to mature cell envelope material.

<p align="center">$g = \gamma_o\dfrac{k_N e_Nc_N}{K_N+c_N}$</p>

## File descriptions
### 1. Generalized linear TSEN model (TSEN_Linear_Generalized)
This TSEN model is a linear reaction network (i.e., without branching) generalized to incorporate any number of intermediate production reactions, but with a single import reaction and growth reaction.

### 2. Production-less linear TSEN model (TSEN_ProductionLess)
This model is similar to the generalized TSEN model (#1), but removes the intermediate production reaction (i.e., is composed only of import and growth reactions). This model is useful as an analytically tractable system, thus allowing deeper physical insight into the system's dynamics.
