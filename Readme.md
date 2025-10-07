## README
This repository contains the code necessary to replicate the analyses presented in the manuscript:
"Rethinking tipping points in spatial ecosystems"

Authors: Swarnendu Banerjee (swarnendubanerjee92@gmail.com), Mara Baudena, Paul Carter, Robbin Bastiaansen, Arjen Doelman, Max Rietkerk

## Study summary
Tipping point theory has garnered substantial attention over the last decades. It predicts abrupt and often irreversible transitions from one ecosystem state to an alternative state. However, typically, ecosystem models that predict tipping neglect spatial dynamics. Recent studies reveal that incorporating spatial dynamics may enable ecosystems to evade tipping predicted by non-spatial models. Here, we use a dryland and a savanna-forest model to synthesize mechanisms by which spatial processes can alter the theory of tipping. We further propose that the underlying drivers of positive feedback leading to alternative stable states may provide insight into the tipping evasion mechanisms most relevant to a specific ecosystem. For instance, while positive feedbacks may arise in drylands from direct self-facilitation, such as enhancing the uptake of a limiting resource, at the savanna-forest boundary, it may arise from mutual inhibition between two ecosystem components. In the former case, ecosystems can evade tipping by forming self-organized patterns, whereas in the latter, the presence of environmental heterogeneity may be required. Our study highlights that deepening our understanding of how ecological feedbacks connect to tipping evasion mechanisms is crucial to formulate better strategies to increase ecosystem resilience.


## Directory and Script Descriptions
Below is a description of the files required to produce the tipping evasion mechanisms discussed in the main text.

## dryland

## Font instability
- solve_grda-pde.m: solves reaction-diffusion-advection PDE and plot the solution after each time "tend" for "K" iterations.
- grda_pde_rhs.m: the reaction-diffusion advection PDE is described here.
- Dgrda_pde_rhs.m: the Jacobian of the reaction -diffusion advection PDE is described here

  To run the code enter the following in terminal solve_grda-pde(tend,K)

## turing_patterns_1D.m: 5.8848 ; now I have 5.8734

## 
