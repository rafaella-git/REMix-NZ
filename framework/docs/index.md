# Welcome to REMix

REMix is addressing research questions in the field of energy system analysis.
The main focus is on the broad techno-economical assessment of possible future
energy system designs and analysis of interactions between technologies.

```{figure} /_static/images/DLR_Logo_REMix_full.svg
:class: only-light
:scale: 200 %
:alt: REMix logo
```

```{figure} /_static/images/DLR_Logo_REMix_full_darkmode.svg
:class: only-dark
:scale: 200 %
:alt: REMix logo
```

This will allow system analysts to inform policy makers and technology
researchers to gain a better understanding of both the system and individual
components.

(remix_key_features)=

## Key Features

To know if REMix is apt for your project take into account these key features:

**Large Models**:
REMix is developed with large models in mind.
This means high spatial and technological resolutions.

**Path Optimization**:
Multi-year analyses are built into the framework.

**Custom accounting approaches**:
The indicator module allows for a very flexible definition of what contributes to the objective functions.

**Flexible modeling**:
There is not a single way of modeling technologies in REMix.
With the {ref}`core modules<modeling_concept_label>` you can find the best way of integrating your modeling needs.

**Multi-criteria optimization**:
Apart from running a cost minimization, also other criteria like ecological or resilience indicators can be taken into account in the objective function.

## Navigation

`````{grid} 4
````{grid-item-card}  About REMix
:link: about_label
:link-type: ref

- What is REMix?
- Feature Overview

````
````{grid-item-card}  Getting Started
:link: getting_started_label
:link-type: ref

- Installation
- Learning REMix
- Tutorials
- Example Applications

````
````{grid-item-card}  Documentation
:link: documentation_label
:link-type: ref

- Modeling Concepts
- API Documentation
- Literature References

````
````{grid-item-card}  Contributing
:link: contributing_label
:link-type: ref

- What to contribute
- How to contribute

````
`````

## Installation

Install with:

```
pip install remix.framework
```

To get git versions:

```
git clone [TODO]
pip install -e ./framework
```

## Recent Publications

-   [Wetzel et al. (2023): "Green energy carriers and energy sovereignty in a climate neutral European energy system"](https://doi.org/10.1016/j.renene.2023.04.015)
-   [Gils et al. (2022): "Model-related outcome differences in power system models with sector coupling - quantification and drivers"](https://doi.org/10.1016/j.rser.2022.112177)[^1]
-   [Gils et al. (2021): "Interaction of hydrogen infrastructures with other sector coupling options towards a zero-emission energy system in Germany"](https://doi.org/10.1016/j.renene.2021.08.016)[^1]
-   [Sasanpour et al. (2021): "Strategic policy targets and the contribution of hydrogen in a 100% renewable European power system"](https://doi.org/10.1016/j.egyr.2021.07.005)[^1]

## Contact

Do not hesitate to ask questions about REMix in the [openmod forum](https://forum.openmod.org/tag/remix).

## License

```{eval-rst}
.. include:: /../LICENSE
```

```{toctree}
:maxdepth: 3
:hidden:

about/index
```

```{toctree}
:maxdepth: 3
:hidden:

getting-started/index
```

```{toctree}
:maxdepth: 3
:hidden:

documentation/index
```

```{toctree}
:maxdepth: 3
:hidden:

contributing/index
```

## Acknowledgments

The preparation of the open source version of REMix was financed by the Helmholtz Association's Energy System Design research programme, which also enables the continuous maintenance of the framework. The methodological and content-related development of REMix was made possible by funding from the German Federal Ministries for Economic Affairs and Climate Protection (BMWK) and for Education and Research (BMBF) as part of the projects UNSEEN (BMWK, FKZ 03EI1004A), Sesame Seed (BMWK, FKZ 03EI1021B), Fahrplan Gaswende (BMWK, FKZ 03EI1030B), ReMoDigital (BMWK, FKZ 03EI1020B), and HINT (BMBF, FKZ 03SF0690) as well as the DLR-internal projects NaGsys and CarnotBat, which were also funded by the Helmholtz Association's Energy System Design research programme. The development of earlier REMix versions, which provided the basis for the published version, was made possible by funding from the projects MuSeKo (BMWK, FKZ 03ET4038B), INTEEVER-II (BMWK, FKZ 03ET4069A), START (BMBF, FKZ 03EK3046D), BEAM-ME (BMWK, FKZ 03ET4023A), INTEEVER (BMWK, FKZ 03ET4020A), Plan-DelyKaD (BMWK, FKZ 0325501), “Lastausgleich” (BMWK, FKZ 0328009), and “Elektromobilitaet” (BMWK, FKZ 0328005A) as well as the Helmholtz Association's Energy System Design research programme and its predecessors.

## Footnotes

[^1]: These papers were still using a non-open legacy version of the REMix framework
