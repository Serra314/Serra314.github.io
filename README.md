# Projects repository of Francesco Serafini

This repository contains the code and documents I have produced in various context.

For now, there is 

1) Multinomial_with_INLA - Tutorial on how to implement Multinomial models using INLA

2) Ranking_earthquake_forecast - The code used to produce all the figures of a paper "Ranking earthquake forecasts using proper scoring rules: Binary events in a low probability environment" than we are about to publish. The preprint version of the paper can be found at https://arxiv.org/abs/2105.12065

3) Toy_examples - A folder containing toy examples on various arguments. LGCP_toy : example on how it is implemented an LGCP model with Inlabru. Specifically, we show how to retireve the same results using the Inlabru function lgcp and the general bru function for fitting Poisson models.

3) Hawkes_process_with_Inlabru - A folder containing code and examples to approximate Hawkes process models using Inlabru. For now, spatial_ETAS_utils.R contains the code to be used in all the other documents (for now it contains function to sample from an ETAS model and to calculate the conditional log-intensity). SpatialETAS_sampling : contains an example on how ETAS sampling works using Italy data, it also shows how to run forecasting experiments on the number of earthquakes in multiple forecasting instances.
