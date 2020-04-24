# reacting_outbreaks_neighboring_localities
Code for the preprint on "Reacting to outbreaks at neighboring localities"

The repository contains 5 files:

SEIR_j.m: function script that simulates one step of the SEIR ODE model at a locality given a step size. This function is called repetitively in all other scripts

run_meta_population.m: Script for running the networked metapopulation model. Generates Figure 1(left) given the set of parameters. 

test_distancing.m: Script that runs the model for varying values of the social distancing responsiveness parameter. Assumes there is no adopted awareness. Generates Figure 1(Right).

test_w_2.m: Script that tests the effects of varying adopted awareness levels with respect to the response at localities and mobility levels. Generatesd Figures 2 and 3 (Left and Right). 

test_mobility.m: Script that tests the effects of varying mobility on outbreak at Locality 2 given a strong response at the origin and weak response at importing locality. 

