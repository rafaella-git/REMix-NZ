## Part d: bonus tasks

The bonus tasks are about analyzing the effects of the integration of electric transmission on the optimal energy system. We will also run a sensitivity analysis to find out how costs for storage units influence the optimal system setup.

### (1) Analyze the effects of the integration of a transmission grid on:

-   Renewable energy share (fraction of renewable generation to overall annual electricity demand)
-   Carbon emission (in kt of CO2)
-   Total system cost (in Mio EUR)

|                    | PV and Wind   | PV and Wind | PV and Wind + no storage | PV and Wind + storage |
|                    | + no storage  | + storage   | + transmission grid      | + transmission grid   |
|--------------------|---------------|-------------|--------------------------|-----------------------|
| Renewable share    |               |             |                          |                       |
| CO2 emissions      |               |             |                          |                       |
| Total system costs |               |             |                          |                       |

-   you have to adjust the storage and transmission upper limits for this task

### (2) Evaluate the influence on your system by additional flexibility from hydro power plants from another node

For this task we need a new model node NO. NO should have just one converter technology (hydro) with an upper limit of 10 GW (see accounting indicators below). NO is connected to DE with a transmission line (upper limit 10 GW). Please create a time series for the hydro converter with just random values, so that the production from hydro is not always = 100 % of the installed capacity.

```
Hydro
Life time: 40 years
Amortization time: 40 years
Interest rate: 6 %
Investment cost: 1665 million EUR/GW
Amortization time: 40 years
Interest rate: 6 %
```

Please evaluate the influence of this additional flexibility on your system via the indicators from task 1. Also have a look into the new line flows and line capacities between the model nodes. What has changed there?
