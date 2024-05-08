## Part d: bonus tasks

The bonus tasks are about analysing the effects of the integration of the lithium-ion batteries on the optimal energy system. We will also attempt a little sensitivity analysis to find out how costs for storage units influence the optimal system setup.

### (1) Analyze the effects of integration of storages

-   Renewable energy share (fraction of renewable generation to overall annual electricity demand)
-   Carbon emission (in kt of CO2)
-   Total system cost (in Mio EUR)

|                    | PV and Wind + no storages | PV only + storages | Wind only + storages | PV and Wind + storages |
| ------------------ | ------------------------- | ------------------ | -------------------- | ---------------------- |
| Renewable share    |                           |                    |                      |                        |
| CO2 emissions      |                           |                    |                      |                        |
| Total system costs |                           |                    |                      |                        |

### (2) Analyze the cost sensitivity of the storage reservoirs

The main cost driver for lithium-ion batteries is the storage capacity itself, since the remaining power electronics is quite cheap. There are other storage types such as vanadium-redox-flow batteries which show a different behavior with expensive membranes for the charging and discharging but very cheap reservoir costs.

Try reducing the specific reservoir costs of the lithium-ion storage (with the current restrictions we allow up to 30\*8 = 240 GWh of reservoir capacity). See, what effect this has on the indicators you assessed in task 1.

### (3) Adding a second storage technology

Instead of modifying the techno-economical parameters of the lithium-ion batteries, we can also see what the optimal system looks like if we let the system choose between two different storage options. For that purpose try adding a second storage technology called VRFB (vanadium redox flow battery) to the model. You can assume the following data:

```
Life time: 25 years
Charging / discharging efficiency: 91 %
Investment charging unit: 200 million EUR/GW
Amortization time: 20 years
Interest rate: 6 %
Fixed cost storage unit: 6 million EUR per year and unit
Investment storage unit: 35 million EUR/GWh
Amortization time: 20 years
Interest rate: 6 %
Fixed cost storage unit: 1.05 million EUR per year and unit
```

After the integration, run the model again and see the effects on the indicators assessed in task 1.
