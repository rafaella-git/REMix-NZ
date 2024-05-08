---
title: What is REMix?
lang: en-US
---

(about_introduction)=

# What is REMix?

REMix stands for "Renewable energy mix for a sustainable energy supply". It is a
framework for setting up linear optimization models written in GAMS. By
framework, we mean a collection of mutually compatible source codes required for
one particular model, which can be combined in a modular manner. In this way
you can reutilize the same modeling concepts (and associated source codes) to
address different content focuses based on a common set of available model
{ref}`features <modeling_concept_label>`.

REMix is developed for studies in the field of energy system modeling. This
means you would usually use it to set up energy system optimization models
(however, there are other applications apart from energy research
conceivable). In particular, these energy system optimization models are often
characterized as bottom-up models in terms of explicitly modeling different
technologies. In addition, these models are resolved on a spatial and a
temporal dimension.

The following figure from {cite}`Cao2021` may give an exemplary impression how
models created with REMix can be generally characterized.

```{figure} /img/ESOM_scope.svg
:align: center
Figure: Characterization of REMix models.
```

## What to model with REMix?

In practical terms, you can model:

1. the **competition between technologies** which may fulfill the same purpose, e.g., power generation, whereas you get also an answer on **when and where** do I need a particular technology
2. any kind of **transportation problem** where the optimal exchange of a commodity between at least two distinct regions shall be determined
3. any kind of **storage problem** where the optimal balancing for something that is produced and consumed at different points in time shall be determined

### Examples for energy system optimization models based on REMix

You may consider two use-cases where an optimization model would be helpful. A Model A may be supposed to analyze the required amount of photovoltaics to reach the European climate targets, whereas another Model B is required for better understanding where in the German power grid the electricity from photovoltaics could be best balanced from controlled charging of battery electric vehicles.

```{figure} /img/aboutREMix.svg
:align: center
Figure: Framework overview.
```

Both, Model A and Model B can be set up with REMix. The main difference between them is the input data used.

### Further examples for existing and conceivable models based on REMix

The availability of appropriate data is the central point you need to care about if you want to set up your own energy system optimization model with REMix. If you have it, there is many real-world use cases you can investigate. To get started you may want to check out some of our previous projects implemented in REMix:

-   [START-EU](https://gitlab.dlr.de/remix/projects/start-eu) a European model to investigate coupling of energy demand sectors (power, heat, transport) and energy carriers (electricity, green synthetic fuels)
-   [NTVR-LCA](https://gitlab.dlr.de/remix/projects/ntvr-lca) a model to analyze the impact of indirect carbon emissions of energy technologies for the transformation of the energy system
-   [YSSP](https://gitlab.dlr.de/remix/projects/yssp) a German model of the power sector for analysis on transmission grid level
-   [IEEE-24bus](https://gitlab.dlr.de/remix/projects/ieee-rts-24-bus-system) a generic power system model based on one of the IEEE reliability test cases
-   [REFuels](https://gitlab.dlr.de/remix/projects/refuels) a site-specific optimization model for production of green kerosene for the aviation sector in Brazil

Some other potential use cases which may be interesting include:

-   A global energy system model to investigate the large-scale transformation of the energy sector, e.g., to better understand the role of "bridging technologies" and impacts of importing liquid hydrogen, methanol and ammonia from other world regions
-   Energy-water nexus modeling to explore the interrelations between energy technologies and the water sector for example with respect to water electrolysis and desalination of seawater
-   A prospective model to analyze the uptake of critical materials for the transformation of the energy system and derive strategies to balance out annual production of materials, recycling and substitutions of technologies