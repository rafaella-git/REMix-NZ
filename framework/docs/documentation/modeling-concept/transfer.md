---
title: Energy transfer
lang: en-US
---

(remix_explanations_energy_transfer)=

# Energy transfer

## General concept

Commodity transfer from one region to another is implemented in REMix as `transfer`.
Links are associated with different parameters such as length and capacity.
The transfer can be realized uni- or bidirectionally.
The amount of a commodity that can be transferred at each point in time depends
on flow limits and the size of the transfer, which can be pre-existing or
newly built.
Additionally, losses during the transfer can be considered.
In this section, the general approach is described and its particularities highlighted.

```{tip}
To model electricity networks REMix also offers linear optimal
power-flow (LOPF).
Under what circumstances you should apply linear optimal power flow instead of a
transshipment model (as explained in the general modeling approach), scroll down
to {ref}`this section <remix_explanations_powerflow>`.
```

## General modeling approach

First, consider a two different model nodes (regions R1 and R2) which have the
same commodity and are connected with a link as shown in the figure
below.

- The link has a specific transfer capacity.
- The energy is transferred in either direction with the `transfer_flows`.
- Energy losses (`transfer_losses`) may occur, e.g. due to voltage
  transformation, electrical resistance or leakage of a pipeline.
- Investments in new links or decommissioning can be done.
- Marginal cost can be accounted for.

```{figure} /img/REMix_transfer_general.svg
:align: center
Figure: Abstract modeling concept for transfer.
```

Each link can be built from smaller segments as shown in the next
figure.
Each segment has a specific transfer_type `T` and a respective length `L`.
For example, energy losses can be associated with the `length` in each segment
and the transfer's `flow`.
Or, investment in the capacity of the link can be specific to the length of
the segments.

```{figure} /img/REMix_transfer_segments.svg
:align: center
Figure: Segments of a transfer link.
```

```{note}
For very simple and generic transfer modeling, which does not take lengths
between nodes into consideration, it is possible to define a link without
segments.
```

Apart from the segments of the link, general transfer parameters can be specified
to characterize the model.
These parameters are mostly related to the link's capacity and the actual flow
of energy.
For example, length-independent energy losses per flow or marginal cost per
flow in the transfer.
In the following sections, the modeling approach for the energy losses and the
accounting are described in detail.

### Energy losses

The total energy loss {math}`\dot{E}_\text{loss}` of a transfer is the
product of the flow {math}`\dot{E}_\text{transfer}` with the total loss
coefficient {math}`c_\text{loss}`.
The loss coefficient itself consists of the sum of all coefficients per length
`coefPerLength` {math}`c_\text{pL}` per link segment (ls) {math}`i`
multiplied with its `length` {math}`D_i` plus the general loss coefficient
`coefPerFlow` of the transfer {math}`c_\text{pF}`:

```{math}
c_\text{loss} = \sum_{i \in \text{ls}} \left( c_{\text{pD,}i} \cdot D_i \right) + c_\text{pF}\\
\dot{E}_\text{loss} = |\dot{E}_\text{transfer}| \cdot c_\text{loss}
```

The figure below illustrates this mathematical formulation:

```{figure} /img/REMix_transfer_losses.svg
:align: center
Figure: Illustration of the energy losses over a transfer link.
```

The losses are equally distributed to both regions. That means that the region
with the starting point of the link will provide power
({math}`\dot{E}_\text{link,R1}`) to the link.

```{math}
\dot{E}_\text{link,R1} = \dot{E}_\text{transfer} + \frac{|\dot{E}_\text{loss}|}{2}
```

And the ending point of the link receives power
({math}`\dot{E}_\text{link,R2}`) from the link.

```{math}
\dot{E}_\text{link,R2} = \dot{E}_\text{transfer} - \frac{|\dot{E}_\text{loss}|}{2}
```

In case of reversed flow direction, the same equations apply. The transfer
flow {math}`\dot{E}_\text{transfer}` will be negative in that case.

````{caution}
If one wants to define an efficiency value {math}`\eta_\text{transfer}` of the
transfer (ratio of energy output to energy input), the value is not easily
connected to the provided coefficients, since the coefficients are relative to
the transfer flow, which is not equal to the energy input or output but
the mean value of both.

```{math}
\dot{E}_\text{transfer} = \frac{|\dot{E}_\text{link,R1}| + |\dot{E}_\text{link,R2}|}{2}\\
\eta_\text{transfer} = \frac{|\dot{E}_\text{transfer}| - \frac{|\dot{E}_\text{loss}|}
{2}}{|\dot{E}_\text{transfer}| + \frac{|\dot{E}_\text{loss}|}{2}}
```

For example, consider a total loss coefficient {math}`c_\text{loss}` of 20 %
given a transfer flow of 100:

```{math}
\begin{align}
\eta_\text{transfer} &= \frac{100 - \frac{100 \cdot 0.2}{2}}{100 + \frac{100 \cdot 0.2}{2}}\\
 &= \frac{90}{110}\\
 &= 81.81 \%
\end{align}
```

If you want to specify your total loss coefficient to result in a specific
efficiency, e.g. of 80 % (or a total loss of 20 % with respect to the input
value), you can calculate that value as follows:

```{math}
\begin{align}
c_\text{loss} &= \frac{2 \cdot \left( 1 - \eta\right)}{1 + \eta}\\
 &= \frac{0.4}{1.8}\\
 &= 22.22 \%
\end{align}
```

````

(remix_explanations_powerflow)=

## Modeling electrical power flow in REMix

Actual real-world power transmission systems often consist of both HVAC and HVDC
transmission links.
To model power flows in HVDC links or in spatially strongly aggregated HVAC
systems (e.g., power markets) it is usually sufficient to model electrical power
flow with the generic transfer model.
Modeling power flows in (roughly aggregated) electricity transmission networks
with that approach is only sufficient if the magnitude and distribution of power
flows can be fully controlled.

In case of spatially resolved HVAC grids, the power transmission capacities
should rather represent the grid transfer capabilities (GTCs), a physical
quantity derived from the maximal current that a transmission line is able to
conduct.
Therefore, additional constraints become necessary and REMix offers to perform
LOPF (Linear Optimal Power Flow):

LOPF summarizes modeling approaches, where the distribution of active power
flows is defined by distribution factors.
These factors are not part of the optimization.
They are usually pre-processed either by a linearization of pre-executed AC
power-flow simulations or based on electrical properties of the transmission
links.
Accordingly, you may provide additional technical parameters, such as the
reactance per length for specified transmission technologies (i.e., conductor
types).

The following flow-chart might help you to decide which modeling approach to use
for electrical grid modeling:

```{figure} /img/REMix_power_flow.svg
:align: center
Figure: Modeling concept for power flow.
```
