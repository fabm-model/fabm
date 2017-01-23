## History

MedERGOM is based on ERGOM (Neumann, 2000) with some specific adaptations in order to better represent the Mediterranean Sea pelagic biogeochemistry. ERGOM is an initially adequate frame to represent the Mediterranean ecosystem as it includes the co-limitation by the two macronutrients (nitrate and phosphate) described as controlling production in this basin. It is also adequate to represent the two main pathways for food and energy transfer in the Mediterranean, the classical herbivores–carnivores food path usually present in eutrophic regions and the small-sized microbial community more usual in the open-sea regions (Siokou-Frangou et al., 2010).


## Process description
However, and in spite of being a potentially suitable candidate to represent the Mediterranean ecosystem, ERGOM was initially created and further developed to simulate the Baltic Sea, which has some obvious differences with the Mediterranean so not only the values of biological rates needs to be changed but also the functional expression of certain equations (see details in Macias et al., 2014a). The main modifications to the original ERGOM formulation are:

1. Zooplankton losses:

   The two original losses from zooplankton (i.e., the closure of the model) in the original ERGOM formulation (to detritus and to ammonium) where dependent on the square of zooplankton abundance (i.e., were density-dependent). In MedERGOM we adopted a different formulation where physiological rates (the two losses described above) are dependent on the zooplankton abundance (not squared) and we added a new loss term, dependent on the squared of zooplankton representing predation by higher trophic levels, which is typically a density-dependent rate (e.g., Ohman and Hshieh, 2008).
2. Light limitation of primary production:

   The light limitation expression for the three phytoplankton groups were the same in the original ERGOM formulation. We needed to create a group-specific formulation of light limitation for the Mediterranean due to the specific composition of the floristic community in this basin.  We also changed the way the optimal light intensity is computed (in the original ERGOM is dependent on the incident light) and also set different optima for each phytoplankton group.


## References

Macías, D., Stips, A., García-Gorriz, E. (2014) The relevance of deep chlorophyll maxima in the open Mediterranean Sea evaluated through 3D hydrodynamic-biogeochemical coupled simulations. Ecological Modelling, 281, 26-37,  [doi](http://dx.doi.org/10.1016/j.ecolmodel.2014.03.002)

Macías, D., García-Gorriz, E., Piroddi, C.,  Stips, A. (2014) Biogeochemical control of marine productivity in the Mediterranean Sea during the last 50 years. Global Biochemical Cycles, 28(8), 897-907, [doi](http://dx.doi.org/10.1002/2014GB004846)

Piroddi, C., Coll, M., Steenbeek, J., Macías, D., Christensen, V. (2015) Modelling the Mediterranean marine ecosystem as a whole: addressing the challenge of complexity. Marine Ecology Progress Series, 533, 47-65. [doi](http://dx.doi.org/10.3354/meps11387).

Macías, D., Garcia-Gorriz, E., Stips, A. (2015) Productivity changes in the Mediterranean Sea for the twenty-first century in response to changes in the regional atmospheric forcing. Frontiers in Marine Science, 2(79). [doi](http://dx.doi.org/10.3389/fmars.2015.00079)

Liquete, C., Piroddi, C., Macias, D., Druon, JN., Zulian, G. Ecosystem services sustainability in the Mediterranean Sea: assessment of status and trends using multiple modelling approaches. Nature Scientific Reports , 6:34162, [doi](http://dx.doi.org/10.1038/srep34162)
