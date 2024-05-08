* SPDX-FileCopyrightText: Copyright (c) 2023-2024 German Aerospace Center (DLR)
* SPDX-License-Identifier: BSD-3-Clause

* // # core_storage
* // The equations in this file describe the storages in the model.

* // ## Variables
* // {special_table_storage_variables}
$offListing
$onEmpty

* ==== declaration of variables ====
integer variables
  storage_unitsTotal_MIP(nodesModel,years,storage_techs,vintage
    ) "Total integer number of active storage units in the system"

positive variables
** // OUTPUT: storage_unitsBuild | OEO_00000350:quantity value
* // ### storage_unitsBuild
* // Title: Storage units build
  storage_unitsBuild(nodesModel,years,storage_techs,vintage
    ) "Number of storage units built"
** // OUTPUT: storage_unitsDecom | OEO_00000350:quantity value
* // ### storage_unitsDecom
* // Title: Storage units decomissioned
  storage_unitsDecom(nodesModel,years,storage_techs,vintage
    ) "Number of storage units decommissioned"
  storage_unitsTotal(nodesModel,years,storage_techs,vintage
    ) "Total number of active storage units in the system"

  storage_level(timeModel,nodesModel,years,storage_techs,vintage,commodity
    ) "Storage level per commodity"
  storage_charge(timeModel,nodesModel,years,storage_techs,vintage,commodity
    ) "Positive storage level changes"
  storage_discharge(timeModel,nodesModel,years,storage_techs,vintage,commodity
    ) "Negative storage level changes"
** // OUTPUT: storage_unitsStateTracker | OEO_00000350:quantity value
* // ### storage_unitsStateTracker
* // Title: Storage units state tracker
  storage_unitsStateTracker(nodesModel,years,yearsCom,storage_techs,vintage,degradation_states
    ) "Degradation state tracking based on unit age and usage"
** // OUTPUT: storage_unitsStateTrackerDecom | OEO_00000350:quantity value
* // ### storage_unitsStateTrackerDecom
* // Title: Storage units decomissioning tracker
  storage_unitsStateTrackerDecom(nodesModel,years,yearsCom,storage_techs,vintage,degradation_states
    ) "Number of decommissioned units having been used up to a certain number of cycles."
** // OUTPUT: storage_levelPerAge | OEO_00000350:quantity value
* // ### storage_levelPerAge
* // Title: Storage levels per unit type and age
  storage_levelPerAge(timeModel,nodesModel,years,yearsCom,storage_techs,vintage,commodity
    ) "Storage level per unit of a certain degradation_state, and commodity"
** // OUTPUT: storage_chargePerAge | OEO_00000350:quantity value
* // ### storage_chargePerAge
* // Title: Storage charges per unit type and age
  storage_chargePerAge(timeModel,nodesModel,years,yearsCom,storage_techs,vintage,commodity
    ) "Commodity charged into the units of a certain degradation state,Influx"

  storage_losses(timeModel,nodesModel,years,storage_techs,vintage,commodity
    ) "Amount of stored commodities lost due to self discharge"
** // OUTPUT: storage_unitsSoC | OEO_00000350:quantity value
* // ### storage_unitsSoC
* // Title: Storage units per state of charge
  storage_unitsSoC(timeModel,nodesModel,years,storage_techs,vintage,soc_statesIn
    ) "Units at certain SoC."

SOS1 variables
  storage_unitsSoC_activeRange(timeModel,nodesModel,years,storage_techs,vintage,soc_statesIn
    ) "Allow SoCs in specific range to be used."
  storage_unitsStateTracker_activeRange(nodesModel,years,yearsCom,storage_techs,vintage,degradation_states
    ) "Allow only one range of degradation states per age group and year to enforce sequential degradation state usage."
;

* The table headers in the parameter input sets represent the following:
* "Title | Description| Constraints | Type |  Units | Ontology"

* // ## Input Files
** // INPUT: storage_reservoirParam | IAO:0000100:data set
* // ### storage_reservoirParam
* // Title: Storage Reservoir Parameters
* // Description: Storage reservoir parameters describe the expansion boundary conditions for storage.
* // {table_storage_reservoirParam}
set pc_storage_reservoirParam
    /
    unitsBuild          "Units Build | Exogenously given units to be built in a given year | minimum:0 | number | {units} | OEO_00010257:power capacity",
    unitsLowerLimit     "Units Lower Limit | Lower limit on total number of units for all vintage classes | minimum:0;lessThan:['unitsBuild','unitsUpperLimit'] | number | {units} | OEO_00000104:constraint",
    unitsUpperLimit     "Units Upper Limit | Upper limit on total number of units for all vintage classes | minimum:0;moreThan:['unitsBuild'];default:'inf' | number | {units} | OEO_00000104:constraint",
    noExpansion         "No Expansion | Prevent expansion beyond exogenous capacity expansion | | boolean | {none} | OEO_00000339:program parameter"
    /;
table storage_reservoirParamIn(nodesData,years,storage_techs,pc_storage_reservoirParam)
$onDelim
$if exist "%instancedir%/storage_reservoirparam.csv" $include "%instancedir%/storage_reservoirparam.csv"
$offDelim
$if not exist "%instancedir%/storage_reservoirparam.csv" $log "No bounds for storage units included"
;

parameter storage_reservoirParam(nodesModel,years,storage_techs,pc_storage_reservoirParam);
$batinclude %aggregateNodes% storage_reservoirParam(nodesModel,years,storage_techs,pc_storage_reservoirParam) storage_reservoirParamIn(nodesData,years,storage_techs,pc_storage_reservoirParam) sum
$batinclude %aggregateNodes% storage_reservoirParam(nodesModel,years,storage_techs,"noExpansion")            storage_reservoirParamIn(nodesData,years,storage_techs,"noExpansion")            smin

** // INPUT: storage_techParam | IAO:0000100:data set
* // ### storage_techParam
* // Title: Storage Technology Parameters
* // Description: Storage technology parameters describe the operational properties of the storage.
* // {table_storage_techParam}
set pc_storage_techParam
    /
    lifeTime                     "Technical lifetime | Technical lifetime of the unit for calculating decommissioning | minimum:0 | integer | {span} | OEO_00000339:operational life time"
    freeDecom                    "Free decommissioning | Allow decommissioning of the unit before the end of the technical life time | | boolean | {none} | OEO_00000339:program parameter"
    mipUnits                     "MIP units | Model the units of the technology as integer values | | boolean | {none} | OEO_00000339:program parameter"
    levelLowerLimit              "Storage level lower limit | Factor per unit for the lower limit of the storage level, overwritten if respective pc_storage_levelProfile is set | minimum:0;maximum:1;lessThan:['levelUpperLimit'] | number | {none} | OEO_00000104:constraint"
    levelUpperLimit              "Storage level upper limit | Factor per unit for the upper limit of the storage level, overwritten if respective pc_storage_levelProfile is set | minimum:0;maximum:1 | number | {none} | OEO_00000104:constraint"
    annualDegradation            "Annual Degradation | Yearly reduction of upper limit due to storage unit degradation | minimum:0;maximum:1;lessThan:['levelUpperLimit'] | number | {none} | OEO_00000339:program parameter"
    usageDegradation             "Usage Degradation | Indicate whether the storage unit has charge-based degradation states | | boolean | {none} | OEO_00000339:program parameter"
    sequentialDegradationStates  "Sequential Degradation States | The unit states of degradation have to be used in sequential order and can only involve two neighboring states of degradation at any time | | boolean | {none} | OEO_00000339:program parameter"
    sequentialSoC                "Sequential SoC | The unit states of charge have to be used in sequential order and can only involve two neighboring SoCs at any time | | boolean | {none} | OEO_00000339:program parameter"
    maxCRate                     "Max Charging Rate | Maximum charging rate in unit capacities per hour | minimum:0 | number | {rate} | OEO_00030019:process attribute"
    maxERate                     "Max Discharging Rate | Maximum discharging rate in unit capacities per hour | minimum:0 | number | {rate} | OEO_00030019:process attribute"
    chargingLoss                 "Relative charging loss | Charging loss per unit provided to the storage | minimum:0 | number | {rate} | OEO_00030019:process attribute"
    dischargingLoss              "Relative discharging loss | Discharging loss per unit taken from the storage | minimum:0 | number | {rate} | OEO_00030019:process attribute"

    /;
table storage_techParam(storage_techs,vintage,pc_storage_techParam)
$onDelim
$if exist "%instancedir%/storage_techparam.csv" $include "%instancedir%/storage_techparam.csv"
$offDelim
$if not exist "%instancedir%/storage_techparam.csv" $log "No technology parameters for storage units included"
;
set storage_hasDegradation(storage_techs,vintage);
storage_hasDegradation(storage_techs,vintage)
  $(storage_techParam(storage_techs,vintage,"annualDegradation")
    or storage_techParam(storage_techs,vintage,"usageDegradation"))
  = yes;

** // INPUT: storage_sizeParam | IAO:0000100:data set
* // ### storage_sizeParam
* // Title: Storage Size Parameters
* // Description: Storage size parameters set the size and self discharge rate of storage per unit and time step.
* // {table_storage_sizeParam}
set pc_storage_sizeParam
    /
    size                "Size coefficient | Coefficient for reservoir size per unit | minimum:0 | number | {capacity_per_unit} | OEO_00230000:storage capacity"
    selfdischarge       "Self discharge rate | Rate of self discharge per time step relative to the current state of charge | maximum:0 | number | {rate} | MISSING_TERM:self-discharge"
    selfdischargeAbs    "Absolute self discharge | Absolute rate of self discharge per time step | maximum:0 | number | {rate} | MISSING_TERM:self-discharge"
    /;
table storage_sizeParam(storage_techs,vintage,commodity,pc_storage_sizeParam)
$onDelim
$if exist "%instancedir%/storage_sizeparam.csv" $include "%instancedir%/storage_sizeparam.csv"
$offDelim
$if not exist "%instancedir%/storage_sizeparam.csv" $log "No reservoir capacities and losses for storage units included"
;
$onempty
** // INPUT: storage_SoCParam | IAO:0000100:data set
* // ### storage_SoCParam
* // Title: Storage State of Charge Parameters
* // Description: Storage state-of-charge (SoC) parameters describe SOC-dependent storage characteristics.
* // {table_storage_SoCParam}
set pc_storage_SoCParam
    /
    SoC             "State of charge | Nominal SoC of this state | | number | {none} | MISSING_TERM:State of charge"
    cRate           "Maximum charge rate | Maximum charge rate at current SoC | minimum:0 | number | {rate} | OEO_00030019:process atribute"
    eRate           "Maximum discharge rate | Maximum discharge rate at current SoC | minimum:0 | number | {rate} | OEO_00030019:process atribute"
    selfdischarge   "Self discharge (SoC) | Self-discharge rate at current SoC | maximum:0 | number | {none} | MISSING_TERM:self-discharge"
    /;
table storage_socparam(storage_techs,vintage,soc_statesIn,pc_storage_SoCParam)
$ondelim
$if exist "%instancedir%/storage_socparam.csv" $include "%instancedir%/storage_socparam.csv"
$offdelim
$if not exist "%instancedir%/storage_socparam.csv" $log "No reservoir state of charge dependent information included"
;

set soc_states(soc_statesIn);
soc_states(soc_statesIn)$(smax((storage_techs,vintage), storage_SoCParam(storage_techs,vintage,soc_statesIn,"cRate")) > 0
                          or smax((storage_techs,vintage), storage_SoCParam(storage_techs,vintage,soc_statesIn,"eRate")) > 0) = yes;

** // INPUT: storage_degradationParam | IAO:0000100:data set
* // ### storage_degradationParam
* // Title: Storage Degradation Parameters
* // Description: Storage degradation parameters describe the degradation of the storage.
* // {table_storage_degradationParam}
set pc_storage_degradationParam
    /
    remainingCapacity   "Remaining capacity coefficient | Coefficient of remaining reservoir capacity per unit in given degradation state | | number | {capacity_per_unit} | OEO_00000339:program parameter"
    maxFullCycles       "Maximum full charge cycles | Maximum number of full charge cycles a unit can remain in the given state | | number | {none} | OEO_00000339:program parameter"
    /;
table storage_degradationParam(storage_techs,vintage,degradation_states,pc_storage_degradationParam)
$ondelim
$if exist "%instancedir%/storage_degradationparam.csv" $include "%instancedir%/storage_degradationparam.csv"
$offdelim
$if not exist "%instancedir%/storage_degradationparam.csv" $log "No reservoir degradation states included"
;

* Set default new state parameters if no usable states are present
storage_degradationParam(storage_techs,vintage,"new","remainingCapacity")
  $storage_hasDegradation(storage_techs,vintage) = 1.0;
storage_degradationParam(storage_techs,vintage,"new","maxFullCycles")
  $storage_hasDegradation(storage_techs,vintage) = 0;
storage_degradationParam(storage_techs,vintage,"new","maxFullCycles")
  $(storage_hasDegradation(storage_techs,vintage)
      and smax(degradation_states, storage_degradationParam(storage_techs,vintage,degradation_states,"maxFullCycles")) = 0)
  = inf;

set storage_usedDegradation(storage_techs,vintage,degradation_states);
option storage_usedDegradation < storage_degradationParam;

** // INPUT: storage_levelProfile | IAO:0000100:data set
* // ### storage_levelProfile
* // Title: Storage Level Profiles
* // Description: Storage level profiles set the limits of storage activities at an hourly level.
* // {table_storage_levelProfile}
set pc_storage_levelProfile
    /
    lower               "Lower storage level profile | Lower profile for storage level, overwrites levelLowerLimit set in pc_storage_techParam |"
    upper               "Upper storage level profile | Upper profile for storage level, overwrites levelUpperLimit set in pc_storage_techParam |"
    fixed               "Fixed storage level profile | Fixed profile for storage level, overwrites level limits set in pc_storage_techParam |"
    /;
table storage_levelProfileLoad(nodesData,years,storage_techs,pc_storage_levelProfile,timeData)
$onDelim
$if exist "%instancedir%/storage_levelprofile.csv" $include "%instancedir%/storage_levelprofile.csv"
$if exist "%instancedir%/storage_levelprofile.csv" $log "Storage level profiles given. Warning: will overwrite level limits for respective technologies."
$offDelim
$if not exist "%instancedir%/storage_levelprofile.csv" $log "No profile for storage level limits included"
;
parameter storage_levelProfileIn(timeData,nodesData,years,storage_techs,pc_storage_levelProfile);
option storage_levelProfileIn < storage_levelProfileLoad;
option clear = storage_levelProfileLoad;

$offEmpty
$onListing


* === SoC mappings ===
set storage_usedTechSoCState(storage_techs,vintage,soc_statesIn);
storage_usedTechSoCState(storage_techs,vintage,soc_states(soc_statesIn))
    $(storage_SoCParam(storage_techs,vintage,soc_states,"eRate") > 0
      and storage_SoCParam(storage_techs,vintage,soc_states,"cRate") > 0)
    = yes;

set storage_validSoCRange(storage_techs,vintage);
storage_validSoCRange(storage_techs,vintage)
    = smax(soc_states$storage_usedTechSoCState(storage_techs,vintage,soc_states), storage_SoCParam(storage_techs,vintage,soc_states,"SoC")) = 1
      and smin(soc_states$storage_usedTechSoCState(storage_techs,vintage,soc_states), storage_SoCParam(storage_techs,vintage,soc_states,"SoC")) = 0;

* === Translate SoC parameters to ordered SoC-set ===
* set ordered_socs / soc1*soc20 /;
* alias(soc_states, soc_states_c);
* parameter storage_SoCParam_ordered(storage_techs,vintage,ordered_socs,pc_storage_SoCParam);
* storage_SoCParam_ordered(storage_techs,vintage,ordered_socs,pc_storage_SoCParam)
*     $storage_validSoCRange(storage_techs,vintage)
*     = sum(soc_states$(ord(ordered_socs) = sum(soc_states_c$(storage_SoCParam(storage_techs,vintage,soc_states_c,"SoC") <= storage_SoCParam(storage_techs,vintage,soc_states,"SoC")
*                                                            and storage_usedTechSoCState(storage_techs,vintage,soc_states_c)), 1)
*                       and storage_usedTechSoCState(storage_techs,vintage,soc_states)),
*           storage_SoCParam(storage_techs,vintage,soc_states,pc_storage_SoCParam));
*
* set storage_usedTechSoCStateOrdered(storage_techs,vintage,ordered_socs);
* storage_usedTechSoCStateOrdered(storage_techs,vintage,ordered_socs)
*     $(storage_SoCParam_ordered(storage_techs,vintage,ordered_socs,"eRate") > 0
*       and storage_SoCParam_ordered(storage_techs,vintage,ordered_socs,"cRate") > 0)
*     = yes;

parameter storage_bigM(storage_techs,vintage);
storage_bigM(storage_techs,vintage) = smax((nodesModel,years), storage_reservoirParam(nodesModel,years,storage_techs,"unitsUpperLimit"));
storage_bigM(storage_techs,vintage)$(storage_bigM(storage_techs,vintage) = INF) = 10000;
storage_bigM(storage_techs,vintage)$(storage_bigM(storage_techs,vintage) = 0) = 10000;

* === calculate the number of years represented by single year to calc ===
alias(yearsToCalc, yearsToCalc_a)
parameter representedYears(years);
representedYears(years) = 0.5 *(smin(yearsToCalc$(yearsToCalc.val > years.val or yearsToCalc.val = smax(yearsToCalc_a, yearsToCalc_a.val)), yearsToCalc.val)
                                - smax(yearsToCalc$(yearsToCalc.val < years.val or yearsToCalc.val = smin(yearsToCalc_a, yearsToCalc_a.val)), yearsToCalc.val));

* ==== calculation of mappings ====

* Technologies with a lifeTime > 0 are available
set storage_availTech(nodesModel,years,storage_techs,vintage);
storage_availTech(nodesModel,years,storage_techs,vintage)
    $(vintage.val = smax(vintage_a$(vintage_a.val <= years.val
        and storage_techParam(storage_techs,vintage_a,"lifeTime") > 0), vintage_a.val)) = yes;

* Technologies to optimize become unavailable if they have an unitsUpperLimit of 0
storage_availTech(nodesModel,years,storage_techs,vintage)
    $(yearsToCalc(years) and storage_reservoirParam(nodesModel,years,storage_techs,"unitsUpperLimit") = 0 ) = no;

* Technologies already built become unavailable if they have an unitsBuild of 0
storage_availTech(nodesModel,years,storage_techs,vintage)
    $( ( not yearsToCalc(years)) and storage_reservoirParam(nodesModel,years,storage_techs,"unitsBuild") = 0 ) = no;

* Used technologies are available technologies over their technical lifeTime
set storage_usedTech(nodesModel,years,storage_techs,vintage);
storage_usedTech(nodesModel,years,storage_techs,vintage)
    $(vintage.val <= years.val
        and years.val < smax(years_a$storage_availTech(nodesModel,years_a,storage_techs,vintage),
                                            years_a.val + storage_techParam(storage_techs,vintage,"lifeTime"))
        ) = yes;

* Technologies have to be decomissioned in the interval of first avail + lifetime to last avail + lifetime
set storage_decomTech(nodesModel,years,storage_techs,vintage);
storage_decomTech(nodesModel,years,storage_techs,vintage)
  $(sum(years_a$storage_usedTech(nodesModel,years_a,storage_techs,vintage), 1)
    and sum(yearsToCalc
      $(sameas(years, yearsToCalc)
        and yearsToCalc.val >= smin(years_a$storage_availTech(nodesModel,years_a,storage_techs,vintage), years_a.val) + storage_techParam(storage_techs,vintage,"lifeTime")
        and yearsToCalc.val <= smax(years_a$storage_availTech(nodesModel,years_a,storage_techs,vintage), years_a.val) + storage_techParam(storage_techs,vintage,"lifeTime")), 1))
  = yes;

* Extend the decom frame to the year after the last year of usedTech
storage_decomTech(nodesModel,yearsToCalc,storage_techs,vintage)
  $(storage_usedTech(nodesModel,yearsToCalc-1,storage_techs,vintage)
    and storage_decomTech(nodesModel,yearsToCalc-1,storage_techs,vintage))
  = yes;

* Mapping for used commodities
set storage_usedCom(storage_techs,vintage,commodity);
option storage_usedCom < storage_sizeParam;

set storage_usedTechCom(nodesModel,years,storage_techs,vintage,commodity);
storage_usedTechCom(nodesModel,years,storage_techs,vintage,commodity)
    $(storage_usedTech(nodesModel,years,storage_techs,vintage)
        and storage_usedCom(storage_techs,vintage,commodity))
    = yes;


* // ## Load units from gdx file
$iftheni.cfgdx not %capsfromgdx%==None
storage_reservoirParam(nodesModel,years,storage_techs,"unitsBuild")
  = sum((accNodesModel,accYears,vintage)
          $(sameas(accNodesModel,nodesModel) and sameas(accYears,years)),
        storage_units(accNodesModel,accYears,storage_techs,vintage,"build"));
storage_reservoirParam(nodesModel,years,storage_techs,"unitsUpperLimit")
  = sum((accNodesModel,accYears,vintage)
          $(sameas(accNodesModel,nodesModel) and sameas(accYears,years)),
        storage_units(accNodesModel,accYears,storage_techs,vintage,"total"));
storage_reservoirParam(nodesModel,years,storage_techs,"unitsLowerLimit")
  = storage_reservoirParam(nodesModel,years,storage_techs,"unitsUpperLimit");
storage_reservoirParam(nodesModel,years,storage_techs,"noExpansion")
  = yes;
$endif.cfgdx

* ==== aggregation of profiles ====
* derive upper and lower profiles then aggregate
set storage_level_hasProfileIn(nodesData,years,storage_techs,pc_storage_levelProfile);
option storage_level_hasProfileIn < storage_levelProfileIn;

set storage_level_hasProfile(nodesModel,years,storage_techs,pc_storage_levelProfile);
storage_level_hasProfile(nodesModelToCalc,yearsToCalc,storage_techs,pc_storage_levelProfile)
    = sum(nodesData$aggregateNodesModel(nodesData,nodesModelToCalc),
            storage_level_hasProfileIn(nodesData,yearsToCalc,storage_techs,pc_storage_levelProfile));

* weighted average for a real sum of units upper limits, plain average if one region has an upper limit of inf
parameter storage_levelProfile(timeModel,nodesModel,years,storage_techs,vintage,pc_storage_levelProfile);
storage_levelProfile(timeModelToCalc,storage_usedTech(nodesModelToCalc,yearsToCalc,storage_techs,vintage),"upper")
    = storage_techParam(storage_techs,vintage,"levelUpperLimit");
storage_levelProfile(timeModelToCalc,storage_usedTech(nodesModelToCalc,yearsToCalc,storage_techs,vintage),"lower")
    = storage_techParam(storage_techs,vintage,"levelLowerLimit");

set storage_finiteUnitLimit(nodesModelToCalc,yearsToCalc,storage_techs);
storage_finiteUnitLimit(nodesModelToCalc,yearsToCalc,storage_techs) = sum(nodesData$aggregateNodesModel(nodesData,nodesModelToCalc), storage_reservoirParamIn(nodesData,yearsToCalc,storage_techs,"unitsUpperLimit")) > 0
                                                 and sum(nodesData$aggregateNodesModel(nodesData,nodesModelToCalc), storage_reservoirParamIn(nodesData,yearsToCalc,storage_techs,"unitsUpperLimit")) < inf;

storage_levelProfile(timeModelToCalc,storage_usedTech(nodesModelToCalc,yearsToCalc,storage_techs,vintage),pc_storage_levelProfile)
    $( storage_level_hasProfile(nodesModelToCalc,yearsToCalc,storage_techs,pc_storage_levelProfile)
        and storage_finiteUnitLimit(nodesModelToCalc,yearsToCalc,storage_techs))
    = sum(nodesData$aggregateNodesModel(nodesData,nodesModelToCalc),
            sum(timeData$timeMapping(timeData,timeModelToCalc),
                  storage_levelProfileIn(timeData,nodesData,yearsToCalc,storage_techs,pc_storage_levelProfile))
              / timeLength(timeModelToCalc)
            * storage_reservoirParamIn(nodesData,yearsToCalc,storage_techs,"unitsUpperLimit"))
    / sum(nodesData$aggregateNodesModel(nodesData,nodesModelToCalc),
            storage_reservoirParamIn(nodesData,yearsToCalc,storage_techs,"unitsUpperLimit"));

storage_levelProfile(timeModelToCalc,storage_usedTech(nodesModelToCalc,yearsToCalc,storage_techs,vintage),pc_storage_levelProfile)
    $( storage_level_hasProfile(nodesModelToCalc,yearsToCalc,storage_techs,pc_storage_levelProfile)
        and sum(nodesData$aggregateNodesModel(nodesData,nodesModelToCalc),
                    storage_reservoirParamIn(nodesData,yearsToCalc,storage_techs,"unitsUpperLimit")) = inf )
    = sum(nodesData$(aggregateNodesModel(nodesData,nodesModelToCalc)
                    and storage_reservoirParamIn(nodesData,yearsToCalc,storage_techs,"unitsUpperLimit") = inf ),
            sum(timeData$timeMapping(timeData,timeModelToCalc),
                  storage_levelProfileIn(timeData,nodesData,yearsToCalc,storage_techs,pc_storage_levelProfile))
              / timeLength(timeModelToCalc))
    / sum(nodesData$(aggregateNodesModel(nodesData,nodesModelToCalc)
                    and storage_reservoirParamIn(nodesData,yearsToCalc,storage_techs,"unitsUpperLimit") = inf ), 1);

* for fixed profiles overwrite upper and lower profile
storage_levelProfile(timeModelToCalc,nodesModelToCalc,yearsToCalc,storage_techs,vintage,"lower")
    $storage_level_hasProfile(nodesModelToCalc,yearsToCalc,storage_techs,"fixed")
    = storage_levelProfile(timeModelToCalc,nodesModelToCalc,yearsToCalc,storage_techs,vintage,"fixed");

storage_levelProfile(timeModelToCalc,nodesModelToCalc,yearsToCalc,storage_techs,vintage,"upper")
    $storage_level_hasProfile(nodesModelToCalc,yearsToCalc,storage_techs,"fixed")
    = storage_levelProfile(timeModelToCalc,nodesModelToCalc,yearsToCalc,storage_techs,vintage,"fixed");


* ==== parameter modifications ====
storage_reservoirParam(nodesModel,years,storage_techs,"unitsLowerLimit")
        $sum(vintage, storage_techParam(storage_techs,vintage,"mipUnits"))
    = floor(storage_reservoirParam(nodesModel,years,storage_techs,"unitsLowerLimit"));
storage_reservoirParam(nodesModel,years,storage_techs,"unitsUpperLimit")
        $sum(vintage, storage_techParam(storage_techs,vintage,"mipUnits"))
    = ceil(storage_reservoirParam(nodesModel,years,storage_techs,"unitsUpperLimit"));


* ==== definition of variables ====

* Initialise variables for unitsBuild
storage_unitsBuild.l(nodesModel,years,storage_techs,vintage)
    $storage_availTech(nodesModel,years,storage_techs,vintage)
    = storage_reservoirParam(nodesModel,years,storage_techs,"unitsBuild");
storage_unitsBuild.lo(nodesModel,yearsToCalc,storage_techs,vintage)
    $storage_availTech(nodesModel,yearsToCalc,storage_techs,vintage)
    = storage_unitsBuild.l(nodesModel,yearsToCalc,storage_techs,vintage);
storage_unitsBuild.fx(nodesModel,years,storage_techs,vintage)
    $storage_reservoirParam(nodesModel,years,storage_techs,"noExpansion")
    = storage_unitsBuild.l(nodesModel,years,storage_techs,vintage);

* Initialise variables for unitsDecom
storage_unitsDecom.l(nodesModel,years,storage_techs,vintage)
    $(storage_decomTech(nodesModel,years,storage_techs,vintage)
      and years.val < sum(yearsToCalc$(ord(yearsToCalc) = 1), yearsToCalc.val))
    = sum((years_a,years_aa)$(sameas(years-1, years_aa)
                      and years_a.val > years_aa.val - storage_techParam(storage_techs,vintage,'lifeTime')
                      and years_a.val <= years.val - storage_techParam(storage_techs,vintage,'lifeTime')
                      and storage_availTech(nodesModel,years_a,storage_techs,vintage)),
        storage_unitsBuild.l(nodesModel,years_a,storage_techs,vintage));

storage_unitsDecom.l(nodesModel,yearsToCalc,storage_techs,vintage)
  $storage_decomTech(nodesModel,yearsToCalc,storage_techs,vintage)
  = sum(years$
        (years.val < sum(yearsToCalc_a$(ord(yearsToCalc_a) = 1), yearsToCalc_a.val)
          and storage_availTech(nodesModel,years,storage_techs,vintage)
          and years.val > sum(years_a$sameas(years_a, yearsToCalc-1), years_a.val) - storage_techParam(storage_techs,vintage,'lifeTime')
          and years.val <= yearsToCalc.val - storage_techParam(storage_techs,vintage,'lifeTime')),
      storage_unitsBuild.l(nodesModel,years,storage_techs,vintage))
    + sum(yearsToCalc_a$
        (yearsToCalc_a.val < sum(yearsToCalc_aa$(ord(yearsToCalc_aa) > 1), yearsToCalc_a.val)
          and storage_availTech(nodesModel,yearsToCalc_a,storage_techs,vintage)
          and yearsToCalc_a.val > sum(years_a$sameas(years_a, yearsToCalc-1), years_a.val) - storage_techParam(storage_techs,vintage,'lifeTime')
          and yearsToCalc_a.val <= yearsToCalc.val - storage_techParam(storage_techs,vintage,'lifeTime')),
      storage_unitsBuild.l(nodesModel,yearsToCalc_a,storage_techs,vintage));
      ;

storage_unitsDecom.lo(nodesModel,yearsToCalc,storage_techs,vintage)
    $(storage_usedTech(nodesModel,yearsToCalc,storage_techs,vintage)
      and not storage_techParam(storage_techs,vintage,"freeDecom"))
    = storage_unitsDecom.l(nodesModel,yearsToCalc,storage_techs,vintage)

* Calculate planned unit expansion
parameter storage_unitsPlanned(nodesModel,years,storage_techs,vintage);
storage_unitsPlanned(nodesModel,years,storage_techs,vintage) = 0;
loop(years,
  storage_unitsPlanned(nodesModel,years,storage_techs,vintage)
    =
    storage_unitsPlanned(nodesModel,years-1,storage_techs,vintage)
        $storage_usedTech(nodesModel,years-1,storage_techs,vintage)
    + storage_unitsBuild.l(nodesModel,years,storage_techs,vintage)
        $storage_availTech(nodesModel,years,storage_techs,vintage)
    - storage_unitsDecom.l(nodesModel,years,storage_techs,vintage)
        $storage_usedTech(nodesModel,years,storage_techs,vintage);
);

* Set initial state for planned units
storage_unitsTotal.l(nodesModel,years,storage_techs,vintage)
  = storage_unitsPlanned(nodesModel,years,storage_techs,vintage)

* Calculate if planned unit expansion is bounded by upper and lower limits
set storage_unitBoundsFixed(nodesModel,years,storage_techs);
storage_unitBoundsFixed(nodesModel,years,storage_techs)
  $(sum(vintage$storage_usedTech(nodesModel,years,storage_techs,vintage),
        storage_unitsPlanned(nodesModel,years,storage_techs,vintage))
    = storage_reservoirParam(nodesModel,years,storage_techs,"unitsUpperLimit")
  and sum(vintage$storage_usedTech(nodesModel,years,storage_techs,vintage),
        storage_unitsPlanned(nodesModel,years,storage_techs,vintage))
    = storage_reservoirParam(nodesModel,years,storage_techs,"unitsLowerLimit"))
  = yes;

* Fix unitsBuild, unitsDecom, unitsTotal if levels are predetermined by upper and lower limits
storage_unitsBuild.fx(nodesModel,years,storage_techs,vintage)
  $(storage_availTech(nodesModel,years,storage_techs,vintage)
    and storage_unitBoundsFixed(nodesModel,years,storage_techs))
  = storage_unitsBuild.l(nodesModel,years,storage_techs,vintage);
storage_unitsDecom.fx(nodesModel,years,storage_techs,vintage)
  $(storage_usedTech(nodesModel,years,storage_techs,vintage)
    and storage_unitBoundsFixed(nodesModel,years,storage_techs))
  = storage_unitsDecom.l(nodesModel,years,storage_techs,vintage);
storage_unitsTotal.fx(nodesModel,years,storage_techs,vintage)
  $(storage_usedTech(nodesModel,years,storage_techs,vintage)
    and storage_unitBoundsFixed(nodesModel,years,storage_techs))
  = storage_unitsTotal.l(nodesModel,years,storage_techs,vintage);

storage_unitsTotal_MIP.up(nodesModel,years,storage_techs,vintage)
    $(storage_usedTech(nodesModel,years,storage_techs,vintage)
      and storage_techParam(storage_techs,vintage,"mipUnits") = 1)
    = storage_reservoirParam(nodesModel,years,storage_techs,"unitsUpperLimit");

* Add parameter for fixing capacities during myopic runs
parameter storage_unitsDelta_upper(nodesModel,years,storage_techs);
parameter storage_unitsDelta_lower(nodesModel,years,storage_techs);


* ==== declaration of equations ====

equations
  Eq_storage_unitsBalance(nodesModel,years,storage_techs,vintage
    ) "Ensures the units balance over the planning period."
  Eq_storage_unitsFixedDecom(nodesModel,years,storage_techs,vintage
    ) "Restricts the fixed decommissioning of units over the planning period."
  Eq_storage_unitsFreeDecom(nodesModel,years,storage_techs,vintage
    ) "Restricts the free decommissioning of units over the planning period."
  Eq_storage_unitsUpperLimit(nodesModel,years,storage_techs
    ) "Upper bound for the total number of units."
  Eq_storage_unitsLowerLimit(nodesModel,years,storage_techs
    ) "Lower bound for the total number of units."
  Eq_storage_unitsTotalMIP(nodesModel,years,storage_techs,vintage
    ) "Fixes the total number of units to the corresponding integer variable."

  Eq_storage_levelUpperLimit(timeModel,nodesModel,years,storage_techs,vintage,commodity
    ) "Upper bound for the total number of units."
  Eq_storage_levelUpperLimit_degradation(timeModel,nodesModel,years,storage_techs,vintage,commodity
    ) "Upper bound for the total number of units if the storage technology accounts for degradation."
  Eq_storage_levelLowerLimit(timeModel,nodesModel,years,storage_techs,vintage,commodity
    ) "Lower bound for the total number of units."
  Eq_storage_losses(timeModel,nodesModel,years,storage_techs,vintage,commodity
    ) "Stored commodities lost due to self discharge."
  Eq_storage_unitsBalanceStates(nodesModel,years,storage_techs,vintage
    ) "Currently available units must have been built at some point in time and have some state of degradation."
  Eq_storage_unitsUpperLimitPerState(nodesModel,years,yearsCom,storage_techs,vintage
    ) "There cannot be more units from a certain year than the amount built in that year."
  Eq_storage_unitsStatesNoRecovery(nodesModel,years,yearsCom,storage_techs,vintage,degradation_states
    ) "There is no recovery."
  Eq_storage_cRateLimit(timeModel,nodesModel,years,storage_techs,vintage,commodity
    ) "Limit charging rate based on unit capacity."
  Eq_storage_eRateLimit(timeModel,nodesModel,years,storage_techs,vintage,commodity
    ) "Limit discharging rate based on unit capacity."
  Eq_storage_charge(timeModel,nodesModel,years,storage_techs,vintage,commodity
    ) "Positive change in state of charge is influx."
  Eq_storage_discharge(timeModel,nodesModel,years,storage_techs,vintage,commodity
    ) "Negative change in state of charge is outflux."

  Eq_storage_levelStateSum(timeModel,nodesModel,years,storage_techs,vintage,commodity
    ) "Sum of degradation class storage levels is the global level."
  Eq_storage_levelUpperLimitPerAge(timeModel,nodesModel,years,yearsCom,storage_techs,vintage,commodity
    ) "Upper bound for units of a certain state."
  Eq_storage_chargingPerAge(timeModel,nodesModel,years,yearsCom,storage_techs,vintage,commodity
    ) "Positive change in state of charge is influx for each degradation class."
  Eq_storage_chargeBasedDegradationDistribution(nodesModel,years,yearsCom,storage_techs,vintage,commodity
    ) "Determine degradation states of storage units."
  Eq_storage_unitsDecomStateSum(nodesModel,years,storage_techs,vintage
    ) "All decomminsioned units must have a degradation state."
  Eq_storage_unitsDegradation(nodesModel,years,yearsCom,storage_techs,vintage,degradation_states
    ) "Only allow one range of degradation states per year and age group."
  Eq_storage_unitsDegradation_onlyOneRange(nodesModel,years,yearsCom,storage_techs,vintage
    ) "Limit usable degradation states to one range, i.e., two neighboring states."

  Eq_storage_unitsSoC_sum(timeModel,nodesModel,years,storage_techs,vintage
    ) "Every unit has a state of charge."
  Eq_storage_levelSoC(timeModel,nodesModel,years,storage_techs,vintage,commodity
    ) "Every unit has a state of charge."
  Eq_storage_unitsSoC(timeModel,nodesModel,years,storage_techs,vintage,soc_statesIn
    ) "Number of units in particular SoC."
  Eq_storage_unitsSoC_onlyOneRange(timeModel,nodesModel,years,storage_techs,vintage
    ) "Limit usable SoC states to one range, i.e., two neighboring states."

  Eq_storage_cRateLimit_SoC(timeModel,nodesModel,years,storage_techs,vintage,commodity
    ) "Limit charging rate based on unit capacity."
  Eq_storage_eRateLimit_SoC(timeModel,nodesModel,years,storage_techs,vintage,commodity
    ) "Limit discharging rate based on unit capacity."
  ;

* ==== equation definition ====
* // ## Equations
* // ### Storage Units Balance
* // Ensures that the total units are consistent with the built and decommissioned units.
* // {Eq_storage_unitsBalance}
Eq_storage_unitsBalance(nodesModelSel,yearsSel,storage_techs,vintage)
    $((storage_usedTech(nodesModelSel,yearsSel,storage_techs,vintage)
          or sum(years$sameas(years,yearsSel), storage_usedTech(nodesModelSel,years-1,storage_techs,vintage)))
        and not storage_unitBoundsFixed(nodesModelSel,yearsSel,storage_techs))
    ..
    storage_unitsTotal(nodesModelSel,yearsSel,storage_techs,vintage)
    =e=
        sum(yearsToCalc$(ord(yearsToCalc) = 1 and sameas(yearsToCalc, yearsSel)),
      sum(years$sameas(years, yearsToCalc),
        storage_unitsTotal(nodesModelSel,years-1,storage_techs,vintage)
          $storage_usedTech(nodesModelSel,years-1,storage_techs,vintage)))
    + sum((yearsToCalc)$(ord(yearsToCalc) > 1 and sameas(yearsToCalc, yearsSel)),
      storage_unitsTotal(nodesModelSel,yearsToCalc-1,storage_techs,vintage)
        $storage_usedTech(nodesModelSel,yearsToCalc-1,storage_techs,vintage))
    + storage_unitsBuild(nodesModelSel,yearsSel,storage_techs,vintage)
        $storage_availTech(nodesModelSel,yearsSel,storage_techs,vintage)
    - storage_unitsDecom(nodesModelSel,yearsSel,storage_techs,vintage)
        $storage_usedTech(nodesModelSel,yearsSel,storage_techs,vintage);

* // ### Storage Units Fixed Decommission
* // Restricts the fixed decommissioning of storage units over the planning period.
* // {Eq_storage_unitsFixedDecom}
Eq_storage_unitsFixedDecom(nodesModelSel,yearsSel,storage_techs,vintage)
    $(storage_decomTech(nodesModelSel,yearsSel,storage_techs,vintage)
        and not storage_techParam(storage_techs,vintage,"freeDecom") = 1
        and not storage_unitBoundsFixed(nodesModelSel,yearsSel,storage_techs))
    ..
    storage_unitsDecom(nodesModelSel,yearsSel,storage_techs,vintage)
    =e=
    sum(years$
        (storage_availTech(nodesModelSel,years,storage_techs,vintage)
          and years.val > sum(yearsToCalc$sameas(yearsToCalc+1, yearsSel), yearsToCalc.val) - storage_techParam(storage_techs,vintage,'lifeTime')
          and years.val <= yearsSel.val - storage_techParam(storage_techs,vintage,'lifeTime')),
      storage_unitsBuild(nodesModelSel,years,storage_techs,vintage));

* // ### Storage Units Free Decomission
* // Restricts the free decommissioning of storage units over the planning period.
* // {Eq_storage_unitsFreeDecom}
Eq_storage_unitsFreeDecom(nodesModelSel,yearsSel,storage_techs,vintage)
    $((storage_decomTech(nodesModelSel,yearsSel,storage_techs,vintage)
        or storage_usedTech(nodesModelSel,yearsSel,storage_techs,vintage))
        and storage_techParam(storage_techs,vintage,"freeDecom") = 1)
    ..
    sum(years$
          ((storage_decomTech(nodesModelSel,years,storage_techs,vintage)
            or storage_usedTech(nodesModelSel,years,storage_techs,vintage))
            and years.val < sum(yearsToCalc_a$(ord(yearsToCalc_a)=1), yearsToCalc_a.val)),
        storage_unitsDecom(nodesModelSel,years,storage_techs,vintage))
    + sum(yearsToCalc$
          ((storage_decomTech(nodesModelSel,yearsToCalc,storage_techs,vintage)
            or storage_usedTech(nodesModelSel,yearsToCalc,storage_techs,vintage))
            and yearsToCalc.val >= sum(yearsToCalc_a$(ord(yearsToCalc_a)=1), yearsToCalc_a.val)
            and yearsToCalc.val <= yearsSel.val),
        storage_unitsDecom(nodesModelSel,yearsToCalc,storage_techs,vintage))
    =g=
    sum(years$
          (storage_availTech(nodesModelSel,years,storage_techs,vintage)
            and years.val < sum(yearsToCalc$(ord(yearsToCalc)=1), yearsToCalc.val) - storage_techParam(storage_techs,vintage,'lifeTime')),
        storage_unitsBuild(nodesModelSel,years,storage_techs,vintage))
    + sum(yearsToCalc$
          (storage_availTech(nodesModelSel,yearsToCalc,storage_techs,vintage)
            and yearsToCalc.val >= sum(yearsToCalc_a$(ord(yearsToCalc_a)=1), yearsToCalc_a.val)
            and yearsToCalc.val <= yearsSel.val - storage_techParam(storage_techs,vintage,'lifeTime')),
        storage_unitsBuild(nodesModelSel,yearsToCalc,storage_techs,vintage));

* // ### Storage Units Lower Limit
* // Lower bound for the total number of storage units.
* // {Eq_storage_unitsLowerLimit}
Eq_storage_unitsLowerLimit(nodesModelSel,yearsSel,storage_techs)
    $(storage_reservoirParam(nodesModelSel,yearsSel,storage_techs,'unitsLowerLimit') > 0 )
    ..
    sum(vintage$storage_usedTech(nodesModelSel,yearsSel,storage_techs,vintage),
        storage_unitsTotal(nodesModelSel,yearsSel,storage_techs,vintage))
    =g=
    storage_reservoirParam(nodesModelSel,yearsSel,storage_techs,"unitsLowerLimit");

* // ### Storage Units Upper Limit
* // Upper bound for the total number of storage units.
* // {Eq_storage_unitsUpperLimit}
Eq_storage_unitsUpperLimit(nodesModelSel,yearsSel,storage_techs)
    $(storage_reservoirParam(nodesModelSel,yearsSel,storage_techs,'unitsUpperLimit') >= 0
        and storage_reservoirParam(nodesModelSel,yearsSel,storage_techs,'unitsUpperLimit') < +inf )
    ..
    sum(vintage$storage_usedTech(nodesModelSel,yearsSel,storage_techs,vintage),
        storage_unitsTotal(nodesModelSel,yearsSel,storage_techs,vintage))
    =l=
    storage_reservoirParam(nodesModelSel,yearsSel,storage_techs,"unitsUpperLimit");

* // ### Storage Units Total MIP
* // Fixes the total number of storage units to the corresponding integer variable.
* // {Eq_storage_unitsTotalMIP}
Eq_storage_unitsTotalMIP(nodesModelSel,yearsSel,storage_techs,vintage)
    $( storage_usedTech(nodesModelSel,yearsSel,storage_techs,vintage)
        and storage_techParam(storage_techs,vintage,"mipUnits") = 1 )
    ..
    storage_unitsTotal(nodesModelSel,yearsSel,storage_techs,vintage)
    =e=
    storage_unitsTotal_MIP(nodesModelSel,yearsSel,storage_techs,vintage);

* // ### Storage Level Lower Limit
* // Lower limit on the storage level.
* // {Eq_storage_levelLowerLimit}
Eq_storage_levelLowerLimit(timeModelSel,nodesModelSel,yearsSel,storage_techs,vintage,commodity)
    $(storage_usedTechCom(nodesModelSel,yearsSel,storage_techs,vintage,commodity)
        and storage_levelProfile(timeModelSel,nodesModelSel,yearsSel,storage_techs,vintage,"lower") > 0)
    ..
    storage_level(timeModelSel,nodesModelSel,yearsSel,storage_techs,vintage,commodity)
    =g=
    storage_levelProfile(timeModelSel,nodesModelSel,yearsSel,storage_techs,vintage,"lower")
    * storage_sizeParam(storage_techs,vintage,commodity,"size")
    * storage_unitsTotal(nodesModelSel,yearsSel,storage_techs,vintage);

* // ### Storage Level Upper Limit
* // Upper limit on the storage level.
* // {Eq_storage_levelUpperLimit}
Eq_storage_levelUpperLimit(timeModelSel,nodesModelSel,yearsSel,storage_techs,vintage,commodity)
    $(storage_usedTechCom(nodesModelSel,yearsSel,storage_techs,vintage,commodity)
        and not storage_hasDegradation(storage_techs,vintage)
        and storage_levelProfile(timeModelSel,nodesModelSel,yearsSel,storage_techs,vintage,"upper") >= 0)
    ..
    storage_level(timeModelSel,nodesModelSel,yearsSel,storage_techs,vintage,commodity)
    =l=
    storage_levelProfile(timeModelSel,nodesModelSel,yearsSel,storage_techs,vintage,"upper")
    * storage_sizeParam(storage_techs,vintage,commodity,"size")
    * storage_unitsTotal(nodesModelSel,yearsSel,storage_techs,vintage);

* // ### Storage Level Upper Limit (degradation)
* // Upper limit on the storage level if the storage technology accounts for degradation.
* // {Eq_storage_levelUpperLimit_degradation}
Eq_storage_levelUpperLimit_degradation(timeModelSel,nodesModelSel,yearsSel,storage_techs,vintage,commodity)
    $(storage_usedTechCom(nodesModelSel,yearsSel,storage_techs,vintage,commodity)
        and storage_hasDegradation(storage_techs,vintage)
        and storage_levelProfile(timeModelSel,nodesModelSel,yearsSel,storage_techs,vintage,"upper") >= 0)
    ..
    storage_level(timeModelSel,nodesModelSel,yearsSel,storage_techs,vintage,commodity)
    =l=
    storage_levelProfile(timeModelSel,nodesModelSel,yearsSel,storage_techs,vintage,"upper")
    * storage_sizeParam(storage_techs,vintage,commodity,"size")
    * sum((degradation_states,yearsCom)$storage_usedDegradation(storage_techs,vintage,degradation_states),
            (storage_degradationParam(storage_techs,vintage,degradation_states,"remainingCapacity")
              - (yearsSel.val - yearsCom.val) * storage_techParam(storage_techs,vintage,"annualDegradation"))
              * storage_unitsStateTracker(nodesModelSel,yearsSel,yearsCom,storage_techs,vintage,degradation_states));

* // ### Storage Unit States Sum (degradation)
* // The number of units in all degradation state must match the total unit number.
* // {Eq_storage_unitsBalanceStates}
Eq_storage_unitsBalanceStates(nodesModelSel,yearsSel,storage_techs,vintage)
    $(storage_usedTech(nodesModelSel,yearsSel,storage_techs,vintage)
      and storage_hasDegradation(storage_techs,vintage))
    ..
    storage_unitsTotal(nodesModelSel,yearsSel,storage_techs,vintage)
    =e=
    sum((yearsCom,degradation_states)$storage_usedDegradation(storage_techs,vintage,degradation_states),
          storage_unitsStateTracker(nodesModelSel,yearsSel,yearsCom,storage_techs,vintage,degradation_states));

* // ### Storage Unit States Upper Limit(degradation)
* // The states of all commissioned units in one year (unitsBuilt) must either still be tracked or have been decommissioned.
* // {Eq_storage_unitsUpperLimitPerState}
Eq_storage_unitsUpperLimitPerState(nodesModelSel,yearsSel,yearsCom,storage_techs,vintage)
    $(storage_usedTech(nodesModelSel,yearsSel,storage_techs,vintage)
        and storage_hasDegradation(storage_techs,vintage))
    ..
    sum(degradation_states$storage_usedDegradation(storage_techs,vintage,degradation_states),
          storage_unitsStateTracker(nodesModelSel,yearsSel,yearsCom,storage_techs,vintage,degradation_states)
          + sum(years$(years.val <= yearsSel.val),
                storage_unitsStateTrackerDecom(nodesModelSel,years,yearsCom,storage_techs,vintage,degradation_states)))
    =e=
    storage_unitsBuild(nodesModelSel,yearsCom,storage_techs,vintage)$(yearsCom.val <= yearsSel.val);

* // ### Storage Unit States Progression (degradation)
* // Unit recovery by reassigning the storage cycles to other units is disabled. Therefore, the number of units in a particular degradation state can only increase if the number in a less degraded state is decreased by at least the same amount.
* // {Eq_storage_unitsStatesNoRecovery}
alias(degradation_states, dc_states);
Eq_storage_unitsStatesNoRecovery(nodesModelSel,yearsSel(yearsToCalc),yearsCom,storage_techs,vintage,degradation_states)
    $(storage_usedTech(nodesModelSel,yearsSel,storage_techs,vintage)
      and yearsToCalc.val > yearsCom.val
      and storage_usedDegradation(storage_techs,vintage,degradation_states))
    ..
    sum(dc_states$(storage_usedDegradation(storage_techs,vintage,dc_states)
              and storage_degradationParam(storage_techs,vintage,dc_states,"maxFullCycles") <= storage_degradationParam(storage_techs,vintage,degradation_states,"maxFullCycles")),
            storage_unitsStateTracker(nodesModelSel,yearsToCalc,yearsCom,storage_techs,vintage,dc_states)
            + sum(years$(years.val <= yearsToCalc.val), storage_unitsStateTrackerDecom(nodesModelSel,years,yearsCom,storage_techs,vintage,dc_states)))
    =l=
    sum(dc_states$(storage_usedDegradation(storage_techs,vintage,dc_states)
                    and storage_degradationParam(storage_techs,vintage,dc_states,"maxFullCycles") <= storage_degradationParam(storage_techs,vintage,degradation_states,"maxFullCycles")),
            storage_unitsStateTracker(nodesModelSel,yearsToCalc-1,yearsCom,storage_techs,vintage,dc_states)
            + sum(years$(years.val < yearsToCalc.val), storage_unitsStateTrackerDecom(nodesModelSel,years,yearsCom,storage_techs,vintage,dc_states)));

* // ### Storage Losses
* // Accumulation of storage losses.
* // {Eq_storage_losses}
Eq_storage_losses(timeModelSel(timeModelToCalc),nodesModelSel,yearsSel,storage_techs,vintage,commodity)
    $storage_usedTechCom(nodesModelSel,yearsSel,storage_techs,vintage,commodity)
    ..
    storage_losses(timeModelSel,nodesModelSel,yearsSel,storage_techs,vintage,commodity)
    =e=
    - storage_level(timeModelToCalc--1,nodesModelSel,yearsSel,storage_techs,vintage,commodity)
      * storage_sizeParam(storage_techs,vintage,commodity,"selfdischarge")
    + storage_sizeParam(storage_techs,vintage,commodity,"selfdischargeAbs")
    + (storage_techParam(storage_techs,vintage,"chargingLoss")
        /(1 - storage_techParam(storage_techs,vintage,"chargingLoss")))
      * storage_charge(timeModelSel,nodesModelSel,yearsSel,storage_techs,vintage,commodity)
        $(storage_techParam(storage_techs,vintage,"chargingLoss") > 0)
    + storage_techParam(storage_techs,vintage,"dischargingLoss")
      * storage_discharge(timeModelSel,nodesModelSel,yearsSel,storage_techs,vintage,commodity)
        $(storage_techParam(storage_techs,vintage,"dischargingLoss") > 0)
    - storage_sizeParam(storage_techs,vintage,commodity,"size")
      * sum(soc_states$(storage_usedTechSoCState(storage_techs,vintage,soc_states)
                        and storage_validSoCRange(storage_techs,vintage)),
            storage_unitsSoC(timeModelToCalc--1,nodesModelSel,yearsSel,storage_techs,vintage,soc_states)
            * storage_SoCParam(storage_techs,vintage,soc_states,"SoC")
            * storage_SoCParam(storage_techs,vintage,soc_states,"selfdischarge"));

* // ### C-Rate Limit
* // The increase in storage level per time step, i.e., the charging rate, is limited relative to the storage capacity.
* // {Eq_storage_cRateLimit}
Eq_storage_cRateLimit(timeModelSel(timeModelToCalc),nodesModelSel,yearsSel,storage_techs,vintage,commodity)
    $(storage_usedTechCom(nodesModelSel,yearsSel,storage_techs,vintage,commodity)
      and storage_techParam(storage_techs,vintage,"maxCRate") > 0)
    ..
    storage_level(timeModelSel,nodesModelSel,yearsSel,storage_techs,vintage,commodity)
    - storage_level(timeModelToCalc--1,nodesModelSel,yearsSel,storage_techs,vintage,commodity)
    =l=
    storage_techParam(storage_techs,vintage,"maxCRate")
    * storage_sizeParam(storage_techs,vintage,commodity,"size")
    * storage_unitsTotal(nodesModelSel,yearsSel,storage_techs,vintage);

* // ### E-Rate Limit
* // The decrease in storage level per time step, i.e., the discharging rate, is limited relative to the storage capacity.
* // {Eq_storage_eRateLimit}
Eq_storage_eRateLimit(timeModelSel(timeModelToCalc),nodesModelSel,yearsSel,storage_techs,vintage,commodity)
    $(storage_usedTechCom(nodesModelSel,yearsSel,storage_techs,vintage,commodity)
      and storage_techParam(storage_techs,vintage,"maxERate") > 0)
    ..
    storage_level(timeModelSel,nodesModelSel,yearsSel,storage_techs,vintage,commodity)
    - storage_level(timeModelToCalc--1,nodesModelSel,yearsSel,storage_techs,vintage,commodity)
    =g=
    - storage_techParam(storage_techs,vintage,"maxERate")
    * storage_sizeParam(storage_techs,vintage,commodity,"size")
    * storage_unitsTotal(nodesModelSel,yearsSel,storage_techs,vintage);

* // ### Storage Charging
* // Increases in storage levels are accounted as charging amounts.
* // {Eq_storage_charge}
Eq_storage_charge(timeModelSel(timeModelToCalc),nodesModelSel,yearsSel,storage_techs,vintage,commodity)
    $(storage_usedTechCom(nodesModelSel,yearsSel,storage_techs,vintage,commodity)
      and storage_techParam(storage_techs,vintage,"chargingLoss") > 0)
    ..
    storage_charge(timeModelSel,nodesModelSel,yearsSel,storage_techs,vintage,commodity)
    =g=
    storage_level(timeModelToCalc,nodesModelSel,yearsSel,storage_techs,vintage,commodity)
    - storage_level(timeModelToCalc--1,nodesModelSel,yearsSel,storage_techs,vintage,commodity);

* // ### Storage Discharging
* // Decreases in storage levels are accounted as discharging amounts.
* // {Eq_storage_discharge}
Eq_storage_discharge(timeModelSel(timeModelToCalc),nodesModelSel,yearsSel,storage_techs,vintage,commodity)
    $(storage_usedTechCom(nodesModelSel,yearsSel,storage_techs,vintage,commodity)
      and storage_techParam(storage_techs,vintage,"dischargingLoss") > 0)
    ..
    storage_discharge(timeModelSel,nodesModelSel,yearsSel,storage_techs,vintage,commodity)
    =g=
    storage_level(timeModelToCalc--1,nodesModelSel,yearsSel,storage_techs,vintage,commodity)
    - storage_level(timeModelToCalc,nodesModelSel,yearsSel,storage_techs,vintage,commodity);

* // ### Storage Level Sum (degradation)
* // The storage level is accounted individually by commissioning year in the case of degradation in order to prohibit the model to assign pre-existing storage cycles to newly built storage reservoirs and thereby avoiding degradation. All those storage levels represent the total storage level.
* // {Eq_storage_levelStateSum}
Eq_storage_levelStateSum(timeModelSel,nodesModelSel,yearsSel,storage_techs,vintage,commodity)
    $(storage_usedTechCom(nodesModelSel,yearsSel,storage_techs,vintage,commodity)
      and storage_techParam(storage_techs,vintage,"usageDegradation"))
    ..
    storage_level(timeModelSel,nodesModelSel,yearsSel,storage_techs,vintage,commodity)
    =e=
    sum(yearsCom, storage_levelPerAge(timeModelSel,nodesModelSel,yearsSel,yearsCom,storage_techs,vintage,commodity));

* // ### Storage Level Upper Limit per Age Group (degradation)
* // Upper limit on the storage level per commissioning year if the storage technology accounts for degradation.
* // {Eq_storage_levelUpperLimitPerAge}
Eq_storage_levelUpperLimitPerAge(timeModelSel,nodesModelSel,yearsSel,yearsCom,storage_techs,vintage,commodity)
    $(storage_usedTechCom(nodesModelSel,yearsSel,storage_techs,vintage,commodity)
      and storage_techParam(storage_techs,vintage,"usageDegradation"))
    ..
    storage_levelPerAge(timeModelSel,nodesModelSel,yearsSel,yearsCom,storage_techs,vintage,commodity)
    =l=
    storage_techParam(storage_techs,vintage,"levelUpperLimit")
    * storage_sizeParam(storage_techs,vintage,commodity,"size")
    * sum(degradation_states,
            (storage_degradationParam(storage_techs,vintage,degradation_states,"remainingCapacity")
              - (yearsSel.val - yearsCom.val) * storage_techParam(storage_techs,vintage,"annualDegradation"))
              * storage_unitsStateTracker(nodesModelSel,yearsSel,yearsCom,storage_techs,vintage,degradation_states));

* // ### Storage Cycle Distribution to Degradation States per Commissioning Year (degradation)
* // The charging amounts are converted to equivalent full cycles which then must be represented by an adequate distribution of degradation states.
* // {Eq_storage_chargeBasedDegradationDistribution}
alias(years, prev_years);
Eq_storage_chargeBasedDegradationDistribution(nodesModelSel,yearsSel,yearsCom,storage_techs,vintage,commodity)
    $(storage_usedTechCom(nodesModelSel,yearsSel,storage_techs,vintage,commodity)
      and storage_techParam(storage_techs,vintage,"usageDegradation"))
    ..
    sum(prev_years$(prev_years.val <= yearsSel.val), representedYears(prev_years)
        * sum(timeModel, storage_chargePerAge(timeModel,nodesModelSel,prev_years,yearsCom,storage_techs,vintage,commodity)))
    =l=
    storage_sizeParam(storage_techs,vintage,commodity,"size")
    * sum(dc_states$storage_usedDegradation(storage_techs,vintage,dc_states),
        storage_degradationParam(storage_techs,vintage,dc_states,"maxFullCycles")
        * (storage_unitsStateTracker(nodesModelSel,yearsSel,yearsCom,storage_techs,vintage,dc_states)
           + sum(prev_years$(prev_years.val <= yearsSel.val), storage_unitsStateTrackerDecom(nodesModelSel,prev_years,yearsCom,storage_techs,vintage,dc_states))));

* // ### Storage Charging per Commissioning Year (degradation)
* // Increases in storage levels are accounted as charging amounts.
* // {Eq_storage_chargingPerAge}
Eq_storage_chargingPerAge(timeModelSel(timeModelToCalc),nodesModelSel,yearsSel,yearsCom,storage_techs,vintage,commodity)
    $(storage_usedTechCom(nodesModelSel,yearsSel,storage_techs,vintage,commodity)
      and storage_techParam(storage_techs,vintage,"usageDegradation"))
    ..
    storage_chargePerAge(timeModelSel,nodesModelSel,yearsSel,yearsCom,storage_techs,vintage,commodity)
    =g=
    storage_levelPerAge(timeModelToCalc,nodesModelSel,yearsSel,yearsCom,storage_techs,vintage,commodity)
    - storage_levelPerAge(timeModelToCalc--1,nodesModelSel,yearsSel,yearsCom,storage_techs,vintage,commodity);

* // ### Storage Unit Decommissioning States Sum (degradation)
* // The number of decommissioned units in all degradation state must match the total decommissioning unit number.
* // {Eq_storage_unitsDecomStateSum}
Eq_storage_unitsDecomStateSum(nodesModelSel,yearsSel,storage_techs,vintage)
    $(storage_techParam(storage_techs,vintage,"usageDegradation"))
    ..
    storage_unitsDecom(nodesModelSel,yearsSel,storage_techs,vintage)
    =e=
    sum((yearsCom,degradation_states)$storage_usedDegradation(storage_techs,vintage,degradation_states),
        storage_unitsStateTrackerDecom(nodesModelSel,yearsSel,yearsCom,storage_techs,vintage,degradation_states));

* // ### Storage Unit Sequential Degradation (degradation)
* // Only active degradation ranges, i.e., two neighboring degradation states, can be used.
* // {Eq_storage_unitsDegradation}
alias(degradation_states, degradation_states_a, degradation_states_b);
Eq_storage_unitsDegradation(nodesModelSel,yearsSel,yearsCom,storage_techs,vintage,degradation_states)
    $(storage_techParam(storage_techs,vintage,"usageDegradation")
      and storage_techParam(storage_techs,vintage,"sequentialDegradationStates"))
    ..
    sum(degradation_states_a$((storage_degradationParam(storage_techs,vintage,degradation_states_a,"maxFullCycles") = storage_degradationParam(storage_techs,vintage,degradation_states,"maxFullCycles")
                      or storage_degradationParam(storage_techs,vintage,degradation_states_a,"maxFullCycles")
                                       = smax(degradation_states_b$(storage_degradationParam(storage_techs,vintage,degradation_states_b,"maxFullCycles") < storage_degradationParam(storage_techs,vintage,degradation_states,"maxFullCycles")
                                                            and storage_usedDegradation(storage_techs,vintage,degradation_states_b)),
                                                            storage_degradationParam(storage_techs,vintage,degradation_states_b,"maxFullCycles")))
                      and storage_usedDegradation(storage_techs,vintage,degradation_states_a)),
        storage_unitsStateTracker_activeRange(nodesModelSel,yearsSel,yearsCom,storage_techs,vintage,degradation_states_a))
    * storage_bigM(storage_techs,vintage)
    =g=
    storage_unitsStateTracker(nodesModelSel,yearsSel,yearsCom,storage_techs,vintage,degradation_states);

* // ### Storage Unit Sequential Degradation Range (degradation)
* // Only one degradation range can be active.
* // {Eq_storage_unitsDegradation_onlyOneRange}
Eq_storage_unitsDegradation_onlyOneRange(nodesModelSel,yearsSel,yearsCom,storage_techs,vintage)
    $(storage_techParam(storage_techs,vintage,"usageDegradation") and storage_techParam(storage_techs,vintage,"sequentialDegradationStates"))
    ..
    sum(degradation_states$storage_usedDegradation(storage_techs,vintage,degradation_states),
          storage_unitsStateTracker_activeRange(nodesModelSel,yearsSel,yearsCom,storage_techs,vintage,degradation_states))
    =e= 1;

* // ### Storage Unit Sequential State of Charge (SoC)
* // Only active state of charge ranges, i.e., two neighboring states of charge, can be used.
* // {Eq_storage_unitsSoC}
alias(soc_states, soc_states_a, soc_states_b);
Eq_storage_unitsSoC(timeModelSel(timeModelToCalc),nodesModelSel,yearsSel,storage_techs,vintage,soc_states)
    $(storage_validSoCRange(storage_techs,vintage) and storage_techParam(storage_techs,vintage,"sequentialSoC"))
    ..
    sum(soc_states_a$((storage_SoCParam(storage_techs,vintage,soc_states_a,"SoC") = storage_SoCParam(storage_techs,vintage,soc_states,"SoC")
                      or storage_SoCParam(storage_techs,vintage,soc_states_a,"SoC")
                                       = smax(soc_states_b$(storage_SoCParam(storage_techs,vintage,soc_states_b,"SoC") < storage_SoCParam(storage_techs,vintage,soc_states,"SoC")
                                                            and storage_usedTechSoCState(storage_techs,vintage,soc_states_b)),
                                                            storage_SoCParam(storage_techs,vintage,soc_states_b,"SoC")))
                      and storage_usedTechSoCState(storage_techs,vintage,soc_states_a)),
        storage_unitsSoC_activeRange(timeModelSel,nodesModelSel,yearsSel,storage_techs,vintage,soc_states_a))
    * storage_bigM(storage_techs,vintage)
    =g=
    storage_unitsSoC(timeModelSel,nodesModelSel,yearsSel,storage_techs,vintage,soc_states);

* // ### Storage Unit Sequential State of Charge Range (SoC)
* // Only one state of charge range can be active.
* // {Eq_storage_unitsSoC_onlyOneRange}
Eq_storage_unitsSoC_onlyOneRange(timeModelSel(timeModelToCalc),nodesModelSel,yearsSel,storage_techs,vintage)
    $(storage_validSoCRange(storage_techs,vintage) and storage_techParam(storage_techs,vintage,"sequentialSoC"))
    ..
    sum(soc_states$storage_usedTechSoCState(storage_techs,vintage,soc_states),
          storage_unitsSoC_activeRange(timeModelSel,nodesModelSel,yearsSel,storage_techs,vintage,soc_states)) =e= 1;

* // ### Storage Unit State of Charge Sum (SoC)
* // Each storage unit must have one state of charge.
* // {Eq_storage_unitsSoC_sum}
Eq_storage_unitsSoC_sum(timeModelSel(timeModelToCalc),nodesModelSel,yearsSel,storage_techs,vintage)
    $storage_validSoCRange(storage_techs,vintage)
    ..
    storage_unitsTotal(nodesModelSel,yearsSel,storage_techs,vintage)
    =e=
    sum(soc_states$storage_usedTechSoCState(storage_techs,vintage,soc_states),
          storage_unitsSoC(timeModelSel,nodesModelSel,yearsSel,storage_techs,vintage,soc_states));

* // ### Storage Level State of Charge Sum (SoC)
* // The total storage level must be represented by units in their specific states of charge.
* // {Eq_storage_levelSoC}
Eq_storage_levelSoC(timeModelSel(timeModelToCalc),nodesModelSel,yearsSel,storage_techs,vintage,commodity)
    $(storage_validSoCRange(storage_techs,vintage) and storage_usedTechCom(nodesModelSel,yearsSel,storage_techs,vintage,commodity))
    ..
    storage_level(timeModelSel,nodesModelSel,yearsSel,storage_techs,vintage,commodity)
    =e=
    storage_sizeParam(storage_techs,vintage,commodity,"size")
    * sum(soc_states$storage_usedTechSoCState(storage_techs,vintage,soc_states),
            storage_unitsSoC(timeModelSel,nodesModelSel,yearsSel,storage_techs,vintage,soc_states)
            * storage_SoCParam(storage_techs,vintage,soc_states,"SoC"));

* // ### C-Rate Limit (SoC)
* // The increase in storage level per time step, i.e., the charging rate, is limited relative to the storage capacity. The coefficients can vary between states of charge.
* // {Eq_storage_cRateLimit_SoC}
Eq_storage_cRateLimit_SoC(timeModelSel(timeModelToCalc),nodesModelSel,yearsSel,storage_techs,vintage,commodity)
    $(storage_usedTechCom(nodesModelSel,yearsSel,storage_techs,vintage,commodity)
      and storage_validSoCRange(storage_techs,vintage)
      and sum(soc_states, storage_SoCParam(storage_techs,vintage,soc_states,"cRate")) < inf)
    ..
    storage_level(timeModelSel,nodesModelSel,yearsSel,storage_techs,vintage,commodity)
    - storage_level(timeModelToCalc--1,nodesModelSel,yearsSel,storage_techs,vintage,commodity)
    =l=
    storage_sizeParam(storage_techs,vintage,commodity,"size")
    * sum(soc_states$storage_usedTechSoCState(storage_techs,vintage,soc_states),
                    storage_SoCParam(storage_techs,vintage,soc_states,"cRate")
                    * storage_unitsSoC(timeModelSel,nodesModelSel,yearsSel,storage_techs,vintage,soc_states));

* // ### E-Rate Limit (SoC)
* // The decrease in storage level per time step, i.e., the discharging rate, is limited relative to the storage capacity. The coefficients can vary between states of charge.
* // {Eq_storage_cRateLimit_SoC}
Eq_storage_eRateLimit_SoC(timeModelSel(timeModelToCalc),nodesModelSel,yearsSel,storage_techs,vintage,commodity)
    $(storage_usedTechCom(nodesModelSel,yearsSel,storage_techs,vintage,commodity)
      and storage_validSoCRange(storage_techs,vintage)
      and sum(soc_states, storage_SoCParam(storage_techs,vintage,soc_states,"eRate")) < inf)
    ..
    storage_level(timeModelSel,nodesModelSel,yearsSel,storage_techs,vintage,commodity)
    - storage_level(timeModelToCalc--1,nodesModelSel,yearsSel,storage_techs,vintage,commodity)
    =g=
    - storage_sizeParam(storage_techs,vintage,commodity,"size")
    * sum(soc_states$storage_usedTechSoCState(storage_techs,vintage,soc_states),
                    storage_SoCParam(storage_techs,vintage,soc_states,"eRate")
                    * storage_unitsSoC(timeModelSel,nodesModelSel,yearsSel,storage_techs,vintage,soc_states));


* ==== model definition ====

Model M_storage
/
  Eq_storage_unitsBalance
  Eq_storage_unitsFixedDecom
  Eq_storage_unitsFreeDecom
  Eq_storage_unitsLowerLimit
  Eq_storage_unitsUpperLimit
  Eq_storage_levelUpperLimit_degradation
  Eq_storage_unitsTotalMIP
  Eq_storage_levelLowerLimit
  Eq_storage_levelUpperLimit
$iftheni.pips %method%==pips
$else.pips
  Eq_storage_losses
$endif.pips
  Eq_storage_unitsBalanceStates
  Eq_storage_unitsUpperLimitPerState
  Eq_storage_unitsStatesNoRecovery
  Eq_storage_cRateLimit
  Eq_storage_eRateLimit
  Eq_storage_charge
  Eq_storage_discharge
  Eq_storage_levelStateSum
  Eq_storage_levelUpperLimitPerAge
  Eq_storage_chargeBasedDegradationDistribution
  Eq_storage_chargingPerAge
  Eq_storage_unitsDecomStateSum
  Eq_storage_unitsDegradation
  Eq_storage_unitsDegradation_onlyOneRange

  Eq_storage_unitsSoC
  Eq_storage_unitsSoC_sum
  Eq_storage_unitsSoC_onlyOneRange
  Eq_storage_levelSoC
  Eq_storage_cRateLimit_SoC
  Eq_storage_eRateLimit_SoC
/;
