## README
This repository contains the code necessary to replicate the analyses presented in the manuscript:
"Rethinking tipping points in spatial ecosystems"

Authors: Swarnendu Banerjee (swarnendubanerjee92@gmail.com), Mara Baudena, Paul Carter, Robbin Bastiaansen, Arjen Doelman, Max Rietkerk

## Study summary
Tipping point theory has garnered substantial attention over the last decades. It predicts abrupt and often irreversible transitions from one ecosystem state to an alternative state. However, typically, ecosystem models that predict tipping neglect spatial dynamics. Recent studies reveal that incorporating spatial dynamics may enable ecosystems to evade tipping predicted by non-spatial models. Here, we use a dryland and a savanna-forest model to synthesize mechanisms by which spatial processes can alter the theory of tipping. We further propose that the underlying drivers of positive feedback leading to alternative stable states may provide insight into the tipping evasion mechanisms most relevant to a specific ecosystem. For instance, while positive feedbacks may arise in drylands from direct self-facilitation, such as enhancing the uptake of a limiting resource, at the savanna-forest boundary, it may arise from mutual inhibition between two ecosystem components. In the former case, ecosystems can evade tipping by forming self-organized patterns, whereas in the latter, the presence of environmental heterogeneity may be required. Our study highlights that deepening our understanding of how ecological feedbacks connect to tipping evasion mechanisms is crucial to formulate better strategies to increase ecosystem resilience.


## Directory and Script Descriptions
Below is a description of the files required to produce the tipping evasion mechanisms discussed in the main text.

## dryland_CoexistenceStates/
- dryland_coexistence_beyondtipping.m: The coexistence states in the dryland exist beyond the tipping point
- dryland_multifront.m: Multiple propagating fronts interacts to form stable coexistence states
  
  Data:
- v_init.mat: Initial stable coexistence state used in dryland_coexistence_beyondtipping.m
- w_init.mat: Initial stable coexistence state used in dryland_coexistence_beyondtipping.m

## dryland_FrontInstability/
- solve_grda-pde_semi_rectangular.m: solves 2D dryland PDE model and plot the solution after each time "tend" for "K" iterations.
  (initialized with vegetation covering half of the landscape)
- solve_grda-pde_circular.m: solves 2D dryland PDE model and plot the solution after each time "tend" for "K" iterations.
  (initialized with vegetation covering a circular patch of the landscape)
- grda_pde_rhs.m: the reaction-diffusion advection PDE is described here.
- Dgrda_pde_rhs.m: the Jacobian of the reaction -diffusion advection PDE is described here
  
  To run the code enter the following in terminal solve_grda-pde_semi_rectangular(tend,K) or solve_grda-pde_circular(tend, K)
  
  Videos:
- SemiRectangular_50_100.avi: Fingering instability and pattern formation when initialized with vegetation covering half of the landscape
- circular_50_100.avi: Fingering instability and pattern formation when initialized with vegetation covering a circular patch

## dryland_GradualInvasion/
- dryland_front.m: solves the 1D dryland PDE model and demonstrate invasion of vegetated state into barren state and vice-versa
- dryland_frontBCfun.m: boundary conditions are described here
- dryland_frontICfun.m: initial condition is described here
- dryland_frontPDEfun.m: dryland model is described here

  To run the code, run dryland_front.m

## dryland_TuringBeforeTipping/ 
- turing_patterns_1D.m: solves the 1D dryland PDE model and demonstrate Turing before Tipping

## savannaForest_CoexistenceStates/
- Maxwell_amplitude.m: plots the Maxwell point and threshold amplitude required for coexistence states as a function of \delta
- savannaForest_Coexistence_pattern.m : solves the 1D savanna-forest model and demonstrate coexistence states in the presence of spatial environmental heterogeneity

## savannaForest_GradualInvasion/
- sav_forest_front.m: solves the 1D savanna-forest PDE model and demonstrate invasion of savanna state into forest state and vice-versa
- sav_forest_frontBCfun.m: boundary conditions are described here
- sav_forest_frontICfun.m: initial condition is described here
- sav_forest_frontPDEfun.m: savanna-forest model is described here

  To run the code, run sav_forest_front.m

## Software Requirements
 - All simulations are performed using `MATLAB 2023b` and `MATLAB 2025a` on Microsoft Windows 11 (x64).
 
