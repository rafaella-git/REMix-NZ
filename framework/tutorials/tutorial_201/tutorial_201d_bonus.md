## Part d: bonus tasks

The bonus tasks are about analyzing the effects of the integration of electric vehicles (EVs) on the optimal energy system. We will also run a sensitivity analysis to find out how costs for stationary and mobile storage units influences the optimal system setup.

### (1) Changed scenario for the stationary storage technology

To see how the contribution of electric vehicles might change, we can also see what the optimal system looks like in a scenario in which the stationary battery storage system becomes way more expensive.
For this purpose try modifying the storage technology in the model. You can assume the following data:

```
Life time: 25 years
Charging / discharging efficiency: 95 %
Investment charging unit: 60 million EUR/GW
Amortization time: 25 years
Interest rate: 6 %
Fixed cost storage unit: 120 million EUR per year and unit
Amortization time: 25 years
Interest rate: 6 %
```

After the modification, run the model again and see the effects on the indicators.
For a sample solution to the bonus tasks, you can refer to the dedicated part in the REMix documentation.

## (2) Analyze the effects of the integration of storages and/or EVs

- Renewable energy share (fraction of renewable generation to overall annual electricity demand)
- Carbon emission (in kt of CO2)
- Total system cost (in Mio EUR)

|                    | PV and Wind   | PV and Wind  | PV and Wind  | PV and Wind     |
|                    | + no storage  | + storage    | + EVs        | + storage + EVs |
|--------------------|---------------|-----------------------------|-----------------|
| Renewable share    |               |              |              |                 |
| CO2 emissions      |               |              |              |                 |
| Total system costs |               |              |              |                 |

### (3) Analyze the effects of controlled charging for electric vehicles

To see how the positive contribution of controlled charging for EVs, try to change the share of the fleet that can be controlled and reasses how the resulting system operation changes.
