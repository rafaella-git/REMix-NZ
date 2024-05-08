---
title: Commodities
lang: en-US
---

(explanations_commodities_label)=

# Commodities

## General concept
REMix builds on modeling regions. Commodities represent energy carriers or physical goods. For `each commodity` in
`each region` of the model there is `one commodity bus`. The commodity buses are the nodes of the model's mathematical
graph representation. The commodities can be, for example,

- electricity,
- methane,
- CO<sub>2</sub>,
- water.

## Commodity bus and commodity balance
In the model for `each commodity bus` conservation of the respective `commodity` is enforced per time step by employing
a `commodity balance equation`. The equation constrains the sum of all incoming (positive) and outgoing (negative) flows
of the commodity to zero. At the example of electricity, incoming flows could be

- production of electricity through conversion from gas,
- import from outside the model scope,
- import from other model regions through the electricity grid or
- discharging of batteries.

Outgoing flows could be

- conversion of electricity to heat,
- consumption of electricity by exogenously defined demand,
- export to other model regions through the electricity grid or
- charging of batteries.

```{figure} /img/REMix_commodities.svg
:align: center
Figure: Illustration of commodities within a region.
```
