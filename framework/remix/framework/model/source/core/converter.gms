* SPDX-FileCopyrightText: Copyright (c) 2023-2024 German Aerospace Center (DLR)
* SPDX-License-Identifier: BSD-3-Clause

* // # core_converter
* // The equations in this file describe the converters in the model.

* // ## Advanced features
* //
* // You can find a more detailed explanation of the %curly_open%ref%curly_close%`MIP features <explanations_mip_label>`.
** // There is also a more detailed explanation on the modeling of outages at the %curly_open%ref%curly_close%`resilience section <explanations_resilience_label>`.
* //
* // In the following table you can see which modifications you have to make to your data to use the different features.
* //
* // | Feature | Modification |
* // | ------ | ------ |
* // | MIP expansion | The converter_tech_Parameter "mipUnits" has to be set to 1 in converter_techParam.dat for the technologies and years where descrete expansion is wanted |
* // | Minimum up/downtime | The converter_tech_Parameters "minUptime" and/or "minDowntime" need to be set to the according amount of minimum time steps.Furthermore, the converter_tech_Parameter "mipDispatch" has to be set to 1 in converter_techParam.dat
* // for the technologies and years where a minimum up- or downtime should be applied so that the units can be shut down. "mipUnits" will automatically be set to 1 |
* // | Partial load | The converter_tech_Parameter "mipDispatch" has to be set to 1 in converter_techParam.dat for the technologies and years where you want to make use of the partial load feature. You also have to set a value to
* // the converter_coefficient_parameters "minLoad" and/or "maxLoad" for the according technologies, vintages, activities and commodities. You can then for example set different "coefficient" values for each activity in converter_coefficient.dat |
** // | Outages | You have to include the converter_outageProfile.dat-file and while executing the model you have to set --method=resilience |

* // ## %curly_open%ref%curly_close%`sets <remix_model_sets_label>`
* //
* // ### set_converter_techs.dat
* // In this file all converter technologies are listed that can be used within your model.
* //
* // ### set_activities.dat
* // All modes with which the converter technologies can be operated need to be added here.
* // Examples: Charge (for the charging mode of a storage), Powergen (for the mode in which a converter produces electricity), Curtailment (for the mode in which a renewable power plant is shut down), …
* // If you would like to use the partial load feature, all activity modes that should be available for this feature need to be defined here as well.
* //
* // ### set_commodities.dat
* // In this file all commodities are listed that can be used and converted by the converter technologies.
* // Examples: Electricity, Coal, Biomass, Hydrogen, …

* // ## Variables
* // {special_table_converter_variables}
$offListing
$onEmpty

* ==== declaration of variables ====
integer variables
  converter_unitsTotal_MIP(nodesModel,years,converter_techs,vintage
    ) "Total integer number of converter units in the system"
  converter_unitsOnline_MIP(timeModel,nodesModel,years,converter_techs,vintage
    ) "Number of active converter units in the system"
** // OUTPUT: converter_unitsUsingActivity_MIP | OEO_00000350:quantity value
* // ### converter_unitsUsingActivity_MIP
* // Title: MIP Converter units using activity
  converter_unitsUsingActivity_MIP(timeModel,nodesModel,years,converter_techs,vintage,activity
    ) "Number of units actively using a specific activity"
  ;

positive variables
  converter_unitsBuild(nodesModel,years,converter_techs,vintage
    ) "Number of converter units built"
  converter_unitsDecom(nodesModel,years,converter_techs,vintage
    ) "Number of converter units decommissioned"
  converter_unitsTotal(nodesModel,years,converter_techs,vintage
    ) "Total number of active converter units in the system"
** // OUTPUT: converter_unitsOnline | OEO_00000350:quantity value
* // ### converter_unitsOnline
* // Title: Converter units online
  converter_unitsOnline(timeModel,nodesModel,years,converter_techs,vintage
    ) "Units that are available for activity in each time step"
** // OUTPUT: converter_activity | OEO_00000350:quantity value
* // ### converter_activity
* // Title: Converter activity
  converter_activity(timeModel,nodesModel,years,converter_techs,vintage,activity
    ) "Activation level per unit and activity as a factor, 1 corresponds to 1 unit under full load"
  converter_rampPos(timeModel,nodesModel,years,converter_techs,vintage
    ) "Positive ramping of the activity level in units"
  converter_rampNeg(timeModel,nodesModel,years,converter_techs,vintage
    ) "Negative ramping of the activity level in units"
** // OUTPUT: converter_unitStartups | OEO_00000350:quantity value
* // ### converter_unitStartups
* // Title: Converter unit startups
  converter_unitStartups(timeModel,nodesModel,years,converter_techs,vintage
    ) "Number of units coming online"
  converter_unitShutdowns(timeModel,nodesModel,years,converter_techs,vintage
    ) "Number of units going offline"
  ;

* The table headers in the parameter input sets represent the following:
* "Title | Description| Constraints | Type |  Units | Ontology"

* // ## Input Files
** // INPUT: converter_capacityParam | IAO:0000100:data set
* // ### converter_capacityParam
* // Title: Converter Capacity Parameters
* // Description: Capacity parameters describe the expansion boundary conditions.
* // {table_converter_capacityParam}
set pc_converter_capacityParam
    /
    unitsBuild          "Units Build | Exogenously given units to be built in a given year. | minimum:0;lessThan:['unitsUpperLimit'] | number | {capacity} | OEO_00010257:power capacity"
    unitsLowerLimit     "Units Lower Limit | Lower limit on total units for all vintage classes. | minimum:0;lessThan:['unitsBuild','unitsUpperLimit'];required:True | number | {capacity} | OEO_00000104:constraint"
    unitsUpperLimit     "Units Upper Limit | Upper limit on total number of units of all vintage classes. | minimum:0;required:True;default:'inf' | number | {capacity} | OEO_00000104:constraint"
    noExpansion         "No Expansion | Prevent expansion beyond exogenous capacity expansion. | | boolean | {none} | OEO_00000339:program parameter"
    /;
table converter_capacityParamIn(nodesData,years,converter_techs,pc_converter_capacityParam)
$onDelim
$if exist "%instancedir%/converter_capacityparam.csv" $include "%instancedir%/converter_capacityparam.csv"
$offDelim
$if not exist "%instancedir%/converter_capacityparam.csv" $log "No bounds for converter units included"
;
* // ```%curly_open%note%curly_close%
* // unitsLowerLimit are the minimum number of units that need to be expanded by the optimizer.
* // The difference between the unitsUpperLimit and the unitsLowerLimit indicates the number of units that can be optimized.
* // unitsBuild can be used for path optimization.
* // ```
parameter converter_capacityParam(nodesModel,years,converter_techs,pc_converter_capacityParam);
$batinclude %aggregateNodes% converter_capacityParam(nodesModel,years,converter_techs,pc_converter_capacityParam) converter_capacityParamIn(nodesData,years,converter_techs,pc_converter_capacityParam) sum
$batinclude %aggregateNodes% converter_capacityParam(nodesModel,years,converter_techs,"noExpansion")              converter_capacityParamIn(nodesData,years,converter_techs,"noExpansion")              smin

** // INPUT: converter_techParam | IAO:0000100:data set
* // ### converter_techParam
* // Title: Converter Technology Parameters
* // Description: Technology parameters describe the operational properties of the converter.
* // {table_converter_techParam}
set pc_converter_techParam
    /
    lifeTime                "Technology Lifetime | Technical life time of the unit for calculating decommissioning. | minimum:0;required:True | integer | {span} | OEO_00000339:operational life time"
    freeDecom               "Free Decommissioning | Allow decommissioning of the unit before the end of the technical life time. | | boolean | {none} | OEO_00000339:program parameter"
    mipUnits                "MIP Units | Model the units of the technology as integer values. | | boolean | {none} | OEO_00000339:program parameter"
    mipDispatch             "MIP Dispatch | Track online state of units as integer values. | | boolean | {none} | OEO_00000339:program parameter"
    activityLowerLimit      "Activity Lower Limit | Lower limit for the sum of all activity factors, overwritten if respective pc_converter_activityProfile is set. | minimum:0;lessThan:['activityUpperLimit'] | number | {none} | OEO_00000104:constraint"
    activityUpperLimit      "Activity Upper Limit | Upper limit for the sum of all activity factors, overwritten if respective pc_converter_activityProfile is set. | minimum:0;default:1 | number | {none} | OEO_00000104:constraint"
    activityRampLimit       "Activity Ramping Limit | Maximum change of total activity per unit and time step. | minimum:0;lessThan:['activityUpperLimit'] | number | {none} | OEO_00000104:constraint"
    startupLimit            "Startup Limit | Number of startups allowed per unit. | minimum:0 | integer | {startups} | OEO_00000104:constraint"
    minUptime               "Minimum Uptime | Minimum number of simulation intervals a unit is required to remain online after startup. | minimum:0 | integer | {interval} | OEO_00000339:program parameter"
    minDowntime             "Minimum Downtime | Minimum number of simulation intervals a unit is required to remain offline after shutdown. | minimum:0 | integer | {interval} | OEO_00000339:program parameter"
    mipDetailedPartialLoad  "MIP Detailed Partial Load | Units in a certain mode can access all activities which do not require a stricter mode.| | boolean | {none} | OEO_00000339:program parameter"
    /;
table converter_techParam(converter_techs,vintage,pc_converter_techParam)
$onDelim
$if exist "%instancedir%/converter_techparam.csv" $include "%instancedir%/converter_techparam.csv"
$offDelim
$if not exist "%instancedir%/converter_techparam.csv" $log "No technology parameters for converter units included"
;
* // ```%curly_open%note%curly_close%
* //  A life time of zero means that the technology is not available.
* // ```
* // ```%curly_open%warning%curly_close%
* //  mipDetailedPartialLoad: One unit can use several activities. This increases computational complexity.
* // ```
* //
** // INPUT: converter_activityProfile | IAO:0000100:data set
* // ### converter_activityProfile
* // Title: Converter Activity Profiles
* // Description: Activity profiles set the limits of converter activities per time step.
* // {table_converter_activityProfile}
set pc_converter_activityProfile
    /
    lower               "Lower profile | Lower profile for converter activity, overwrites activityLowerLimit set in pc_converter_techParam. |"
    upper               "Upper profile | Upper profile for converter activity, overwrites activityUpperLimit set in pc_converter_techParam. |"
    fixed               "Fixed profile | Fixed profile for converter activity, overwrites activity limits set in pc_converter_techParam. |"
    /;
table converter_activityProfileLoad(nodesData,years,converter_techs,pc_converter_activityProfile,timeData)
$onDelim
$if exist "%instancedir%/converter_activityprofile.csv" $include "%instancedir%/converter_activityprofile.csv"
$if exist "%instancedir%/converter_activityprofile.csv" $log "Converter activity profiles given. Warning: will overwrite activity limits for respective technologies"
$offDelim
$if not exist "%instancedir%/converter_activityprofile.csv" $log "No profile for converter activity limits included."
;
parameter converter_activityProfileIn(timeData,nodesData,years,converter_techs,pc_converter_activityProfile);
option converter_activityProfileIn < converter_activityProfileLoad;
option clear = converter_activityProfileLoad;

** // INPUT: converter_coefficient | IAO:0000100:data set
* // ### converter_coefficient
* // Title: Converter Coefficients
* // Description: Converter coefficients are used to indicate the rates at which the commodities are transformed by the activities.
* // {table_converter_coefficient}
set pc_converter_coefficient
    /
    coefficient         "Converter Coefficients | Coefficients for used or provided commodities per activity and unit in full load, overwritten by converter_coefficientProfile if given. | required:True | number | {rate} | OEO_00030019:process attribute"
    constant            "Constant Commodity Production | Constant commodity production or usage which is independent of the unit load. | | number | {none} | OEO_00030019:process attribute"
    minLoad             "Minimum Load | Minimum required total unit load to gain access to this activity. | | number | {rate} | OEO_00000339:program parameter"
    maxLoad             "Maximum Load | Maximum allowed total unit load to gain access to this activity. | | number | {rate} | OEO_00000339:program parameter"
    /;
table converter_coefficient(converter_techs,vintage,activity,commodity,pc_converter_coefficient)
$ondelim
$if exist "%instancedir%/converter_coefficient.csv" $include "%instancedir%/converter_coefficient.csv"
$if not exist "%instancedir%/converter_coefficient.csv" $if exist "%instancedir%/converter_activityparam.csv" $include "%instancedir%/converter_activityparam.csv"
$if not exist "%instancedir%/converter_coefficient.csv" $if exist "%instancedir%/converter_activityparam.csv" $log "Project is still using the old file name converter_activityParam.dat, please change it to converter_coefficient.dat in the future."
$offdelim
$if not exist "%instancedir%/converter_coefficient.csv" $if not exist "%instancedir%/converter_activityparam.csv" $log "No coefficients for converter activities included"
;
* // ```%curly_open%note%curly_close%
* //  In converter coeffiecients, negative values indicate inputs, positive values outputs.
* // ```
* //
** // INPUT: converter_coefficientprofile | IAO:0000100:data set
* // ### converter_coefficientprofile
* // Title: Converter Coefficients Profiles
* // Description: Converter coefficient profiles are used to assign coefficient values to individual time steps.
table converter_coefficientProfileLoad(nodesData,years,converter_techs,vintage,activity,commodity,timeData)
$ondelim
$if exist "%instancedir%/converter_coefficientprofile.csv" $include "%instancedir%/converter_coefficientprofile.csv"
$if exist "%instancedir%/converter_coefficientprofile.csv" $log "Converter coefficient profiles given. Warning: will overwrite coefficients for respective technologies"
$offdelim
$if not exist "%instancedir%/converter_coefficientprofile.csv" $log "No time-dependent coefficients for converter activities included"
;
parameter converter_coefficientProfileIn(timeData,nodesData,years,converter_techs,vintage,activity,commodity);
option converter_coefficientProfileIn < converter_coefficientProfileLoad;

$offEmpty
$onListing

* ==== loading units from gdx file ====

* Load units from gdx file
$iftheni.cfgdx not %capsfromgdx%==None
converter_capacityParam(nodesModel,years,converter_techs,"unitsBuild")
  = sum((accNodesModel,accYears,vintage)
          $(sameas(accNodesModel,nodesModel) and sameas(accYears,years)),
        converter_units(accNodesModel,accYears,converter_techs,vintage,"build"));
converter_capacityParam(nodesModel,years,converter_techs,"unitsUpperLimit")
  = sum((accNodesModel,accYears,vintage)
          $(sameas(accNodesModel,nodesModel) and sameas(accYears,years)),
        converter_units(accNodesModel,accYears,converter_techs,vintage,"total"));
converter_capacityParam(nodesModel,years,converter_techs,"unitsLowerLimit")
  = converter_capacityParam(nodesModel,years,converter_techs,"unitsUpperLimit");
converter_capacityParam(nodesModel,years,converter_techs,"noExpansion")
  = yes;
$endif.cfgdx

* ==== calculation of mappings ====

* Technologies with a lifeTime > 0 are available
set converter_availTech(nodesModel,years,converter_techs,vintage);
converter_availTech(nodesModel,years,converter_techs,vintage)
    $(vintage.val = smax(vintage_a$(vintage_a.val <= years.val
        and converter_techParam(converter_techs,vintage_a,"lifeTime") > 0), vintage_a.val)) = yes;

* Technologies to optimize become unavailable if they have an unitsUpperLimit of 0
converter_availTech(nodesModel,years,converter_techs,vintage)
    $(yearsToCalc(years) and converter_capacityParam(nodesModel,years,converter_techs,"unitsUpperLimit") = 0 ) = no;

* Technologies already built become unavailable if they have an unitsBuild of 0
converter_availTech(nodesModel,years,converter_techs,vintage)
    $( ( not yearsToCalc(years)) and converter_capacityParam(nodesModel,years,converter_techs,"unitsBuild") = 0 ) = no;

* Used technologies are available technologies over their technical lifeTime
set converter_usedTech(nodesModel,years,converter_techs,vintage);
converter_usedTech(nodesModel,years,converter_techs,vintage)
    $(vintage.val <= years.val
        and years.val < smax(years_a$converter_availTech(nodesModel,years_a,converter_techs,vintage),
                                years_a.val + converter_techParam(converter_techs,vintage,"lifeTime"))
        ) = yes;

* Technologies have to be decomissioned in the interval of first avail + lifetime to last avail + lifetime
set converter_decomTech(nodesModel,years,converter_techs,vintage);
converter_decomTech(nodesModel,years,converter_techs,vintage)
  $(sum(years_a$converter_usedTech(nodesModel,years_a,converter_techs,vintage), 1)
    and sum(yearsToCalc
      $(sameas(years, yearsToCalc)
        and yearsToCalc.val >= smin(years_a$converter_availTech(nodesModel,years_a,converter_techs,vintage), years_a.val) + converter_techParam(converter_techs,vintage,"lifeTime")
        and yearsToCalc.val <= smax(years_a$converter_availTech(nodesModel,years_a,converter_techs,vintage), years_a.val) + converter_techParam(converter_techs,vintage,"lifeTime")), 1))
  = yes;

* Extend the decom frame to the year after the last year of usedTech
converter_decomTech(nodesModel,yearsToCalc,converter_techs,vintage)
  $(converter_usedTech(nodesModel,yearsToCalc-1,converter_techs,vintage)
    and converter_decomTech(nodesModel,yearsToCalc-1,converter_techs,vintage))
  = yes;

* Mapping for used activities and commodities
set converter_usedActCom(converter_techs,vintage,activity,commodity);
option converter_usedActCom < converter_coefficient;

set converter_usedAct(converter_techs,vintage,activity);
option converter_usedAct < converter_usedActCom;

set converter_usedCom(converter_techs,vintage,commodity);
option converter_usedCom < converter_usedActCom;

set converter_usedTechAct(nodesModel,years,converter_techs,vintage,activity);
converter_usedTechAct(nodesModel,years,converter_techs,vintage,activity)
    $(converter_usedTech(nodesModel,years,converter_techs,vintage)
        and converter_usedAct(converter_techs,vintage,activity))
    = yes;

set converter_useRampPos(nodesModel,years,converter_techs,vintage);
converter_useRampPos(nodesModelToCalc,yearsToCalc,converter_techs,vintage)
  $(sum(indicator$accounting_converterStartup(indicator,nodesModelToCalc,converter_techs,vintage,"perRamp"), 1)
    or sum(indicator$accounting_converterStartup(indicator,nodesModelToCalc,converter_techs,vintage,"perRampPos"), 1))
  = 1;

set converter_useRampNeg(nodesModel,years,converter_techs,vintage);
converter_useRampNeg(nodesModelToCalc,yearsToCalc,converter_techs,vintage)
  $(sum(indicator$accounting_converterStartup(indicator,nodesModelToCalc,converter_techs,vintage,"perRamp"), 1)
    or sum(indicator$accounting_converterStartup(indicator,nodesModelToCalc,converter_techs,vintage,"perRampNeg"), 1))
  = 1;


* ==== aggregation of profiles ====

* derive upper and lower profiles then aggregate
set converter_activity_hasProfileIn(nodesData,years,converter_techs,pc_converter_activityProfile);
option converter_activity_hasProfileIn < converter_activityProfileIn;

set converter_activity_hasProfile(nodesModel,years,converter_techs,pc_converter_activityProfile);
converter_activity_hasProfile(nodesModelToCalc,yearsToCalc,converter_techs,pc_converter_activityProfile)
    = sum(nodesData$aggregateNodesModel(nodesData,nodesModelToCalc),
            converter_activity_hasProfileIn(nodesData,yearsToCalc,converter_techs,pc_converter_activityProfile));

set converter_coefficient_hasProfileIn(nodesData,years,converter_techs,vintage,activity,commodity);
option converter_coefficient_hasProfileIn < converter_coefficientProfileIn;

set converter_coefficient_hasProfile(nodesModel,years,converter_techs,vintage,activity,commodity);
converter_coefficient_hasProfile(nodesModelToCalc,yearsToCalc,converter_techs,vintage,activity,commodity)
    = sum(nodesData$aggregateNodesModel(nodesData,nodesModelToCalc),
            converter_coefficient_hasProfileIn(nodesData,yearsToCalc,converter_techs,vintage,activity,commodity));

* aagregate time dimension
parameter converter_activityProfileIn_aggTime(timeModel,nodesData,yearsToCalc,converter_techs,pc_converter_activityProfile);
converter_activityProfileIn_aggTime(timeModelToCalc,nodesData,yearsToCalc,converter_techs,pc_converter_activityProfile)
  $converter_activity_hasProfileIn(nodesData,yearsToCalc,converter_techs,pc_converter_activityProfile)
  = sum(timeData$timeMapping(timeData,timeModelToCalc),
          converter_activityProfileIn(timeData,nodesData,yearsToCalc,converter_techs,pc_converter_activityProfile)
          / timeLength(timeModelToCalc));
option clear = converter_activityProfileIn;

* weighted average for a real sum of units upper limits, plain average if one region has an upper limit of inf
parameter converter_activityProfile(timeModel,nodesModel,years,converter_techs,vintage,pc_converter_activityProfile);
converter_activityProfile(timeModelToCalc,converter_usedTech(nodesModelToCalc,yearsToCalc,converter_techs,vintage),"upper")
    = converter_techParam(converter_techs,vintage,"activityUpperLimit");
converter_activityProfile(timeModelToCalc,converter_usedTech(nodesModelToCalc,yearsToCalc,converter_techs,vintage),"lower")
    = converter_techParam(converter_techs,vintage,"activityLowerLimit");


set converter_finiteUnitLimit(nodesModelToCalc,yearsToCalc,converter_techs);
converter_finiteUnitLimit(nodesModelToCalc,yearsToCalc,converter_techs)
    = sum(nodesData$aggregateNodesModel(nodesData,nodesModelToCalc), converter_capacityParamIn(nodesData,yearsToCalc,converter_techs,"unitsUpperLimit")) > 0
        and sum(nodesData$aggregateNodesModel(nodesData,nodesModelToCalc), converter_capacityParamIn(nodesData,yearsToCalc,converter_techs,"unitsUpperLimit")) < inf;

set converter_infiniteUnitLimit(nodesModelToCalc,yearsToCalc,converter_techs);
converter_infiniteUnitLimit(nodesModelToCalc,yearsToCalc,converter_techs)
    = sum(nodesData$aggregateNodesModel(nodesData,nodesModelToCalc), converter_capacityParamIn(nodesData,yearsToCalc,converter_techs,"unitsUpperLimit")) = inf;

converter_activityProfile(timeModelToCalc,converter_usedTech(nodesModelToCalc,yearsToCalc,converter_techs,vintage),pc_converter_activityProfile)
    $(converter_activity_hasProfile(nodesModelToCalc,yearsToCalc,converter_techs,pc_converter_activityProfile)
        and converter_finiteUnitLimit(nodesModelToCalc,yearsToCalc,converter_techs))
    = sum(nodesData$(aggregateNodesModel(nodesData,nodesModelToCalc)
                and converter_capacityParamIn(nodesData,yearsToCalc,converter_techs,"unitsUpperLimit") < inf ),
              converter_activityProfileIn_aggTime(timeModelToCalc,nodesData,yearsToCalc,converter_techs,pc_converter_activityProfile)
              * converter_capacityParamIn(nodesData,yearsToCalc,converter_techs,"unitsUpperLimit"))
    / sum(nodesData$aggregateNodesModel(nodesData,nodesModelToCalc),
            converter_capacityParamIn(nodesData,yearsToCalc,converter_techs,"unitsUpperLimit"));

converter_activityProfile(timeModelToCalc,converter_usedTech(nodesModelToCalc,yearsToCalc,converter_techs,vintage),pc_converter_activityProfile)
    $(converter_activity_hasProfile(nodesModelToCalc,yearsToCalc,converter_techs,pc_converter_activityProfile)
        and converter_infiniteUnitLimit(nodesModelToCalc,yearsToCalc,converter_techs))
    = sum(nodesData$(aggregateNodesModel(nodesData,nodesModelToCalc)
                and converter_capacityParamIn(nodesData,yearsToCalc,converter_techs,"unitsUpperLimit") = inf),
              converter_activityProfileIn_aggTime(timeModelToCalc,nodesData,yearsToCalc,converter_techs,pc_converter_activityProfile))
    / sum(nodesData$(aggregateNodesModel(nodesData,nodesModelToCalc)
                    and converter_capacityParamIn(nodesData,yearsToCalc,converter_techs,"unitsUpperLimit") = inf ),
            1);
option clear = converter_activityProfileIn_aggTime;

* for fixed profiles overwrite upper and lower profile
converter_activityProfile(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage,"lower")
    $converter_activity_hasProfile(nodesModelToCalc,yearsToCalc,converter_techs,"fixed")
    = converter_activityProfile(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage,"fixed");

converter_activityProfile(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage,"upper")
    $converter_activity_hasProfile(nodesModelToCalc,yearsToCalc,converter_techs,"fixed")
    = converter_activityProfile(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage,"fixed");

$ifthene.roundcoefs %roundcoefs%==1
converter_activityProfile(timeModelToCalc,converter_usedTech(nodesModelToCalc,yearsToCalc,converter_techs,vintage),pc_converter_activityProfile)
  = round(converter_activityProfile(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage,pc_converter_activityProfile), 3);
converter_activityProfile(timeModelToCalc,converter_usedTech(nodesModelToCalc,yearsToCalc,converter_techs,vintage),pc_converter_activityProfile)
  $(converter_activityProfile(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage,pc_converter_activityProfile) < 1e-3)
  = 0;
$endif.roundcoefs

* weighted average for a real sum of units upper limits, plain average if one region has an upper limit of inf
parameter converter_coefficientProfile(timeModel,nodesModel,years,converter_techs,vintage,activity,commodity);
converter_coefficientProfile(timeModelToCalc,converter_usedTechAct(nodesModelToCalc,yearsToCalc,converter_techs,vintage,activity),commodity)
    = converter_coefficient(converter_techs,vintage,activity,commodity,"coefficient");

converter_coefficientProfile(timeModelToCalc,converter_usedTechAct(nodesModelToCalc,yearsToCalc,converter_techs,vintage,activity),commodity)
    $( converter_coefficient_hasProfile(nodesModelToCalc,yearsToCalc,converter_techs,vintage,activity,commodity)
        and converter_finiteUnitLimit(nodesModelToCalc,yearsToCalc,converter_techs) )
    = sum(nodesData$aggregateNodesModel(nodesData,nodesModelToCalc),
          sum(timeData$timeMapping(timeData,timeModelToCalc),
                  converter_coefficientProfileIn(timeData,nodesData,yearsToCalc,converter_techs,vintage,activity,commodity))
              / timeLength(timeModelToCalc)
            * converter_capacityParamIn(nodesData,yearsToCalc,converter_techs,"unitsUpperLimit") )
    / sum(nodesData$aggregateNodesModel(nodesData,nodesModelToCalc),
            converter_capacityParamIn(nodesData,yearsToCalc,converter_techs,"unitsUpperLimit"));

converter_coefficientProfile(timeModelToCalc,converter_usedTechAct(nodesModelToCalc,yearsToCalc,converter_techs,vintage,activity),commodity)
    $( converter_coefficient_hasProfile(nodesModelToCalc,yearsToCalc,converter_techs,vintage,activity,commodity)
        and converter_infiniteUnitLimit(nodesModelToCalc,yearsToCalc,converter_techs))
    = sum(nodesData$(aggregateNodesModel(nodesData,nodesModelToCalc)
                    and converter_capacityParamIn(nodesData,yearsToCalc,converter_techs,"unitsUpperLimit") = inf ),
          sum(timeData$timeMapping(timeData,timeModelToCalc),
                  converter_coefficientProfileIn(timeData,nodesData,yearsToCalc,converter_techs,vintage,activity,commodity))
              / timeLength(timeModelToCalc))
    / sum(nodesData$(aggregateNodesModel(nodesData,nodesModelToCalc)
                    and converter_capacityParamIn(nodesData,yearsToCalc,converter_techs,"unitsUpperLimit") = inf ), 1);

$ifthene.roundcoefs %roundcoefs%==1
converter_coefficientProfile(timeModelToCalc,converter_usedTechAct(nodesModelToCalc,yearsToCalc,converter_techs,vintage,activity),commodity)
  = round(converter_coefficientProfile(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage,activity,commodity), 3);
converter_coefficientProfile(timeModelToCalc,converter_usedTechAct(nodesModelToCalc,yearsToCalc,converter_techs,vintage,activity),commodity)
  $(converter_coefficientProfile(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage,activity,commodity) < 1e-3
        and converter_coefficientProfile(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage,activity,commodity) > -1e-3)
  = 0;
$endif.roundcoefs

* ==== prepare partial load behavior parameters ====

* set disabled max load requirements to 1 to avoid excessive logical expressions later on
alias(commodity, com);
converter_coefficient(converter_techs,vintage,activity,commodity,"maxLoad")
    $((converter_coefficient(converter_techs,vintage,activity,commodity,"maxLoad") = 0)
      and converter_usedAct(converter_techs,vintage,activity))
    = 1;

* indicate that at least one activity of a technology makes use of partial load behavior
set converter_hasMaxLoad(converter_techs,vintage);
converter_hasMaxLoad(converter_techs,vintage)
  = smin((activity,commodity)$converter_usedAct(converter_techs,vintage,activity),
            converter_coefficient(converter_techs,vintage,activity,commodity,"maxLoad")) < 1;

set converter_hasMinLoad(converter_techs,vintage);
converter_hasMinLoad(converter_techs,vintage)
  = smax((activity,commodity)$converter_usedAct(converter_techs,vintage,activity),
            converter_coefficient(converter_techs,vintage,activity,commodity,"minLoad")) > 0;

set converter_hasConstantFluxInActivity(converter_techs,vintage);
converter_hasConstantFluxInActivity(converter_techs,vintage)
  = sum((activity,commodity)$converter_usedAct(converter_techs,vintage,activity),
                        abs(converter_coefficient(converter_techs,vintage,activity,commodity,"constant"))) > 0;

* setting up a requirements parameter to avoid having to cycle through commodities in equations
set pc_converter_activityRequirements
    /
    minLoad             "Lower profile for converter activity"
    maxLoad             "Upper profile for converter activity"
    /;
parameter converter_activityRequirements(converter_techs,vintage,activity,pc_converter_activityRequirements);
converter_activityRequirements(converter_techs,vintage,activity,"minLoad")
  = smax(commodity, converter_coefficient(converter_techs,vintage,activity,commodity,"minLoad"));
converter_activityRequirements(converter_techs,vintage,activity,"maxLoad")
  = smin(commodity, converter_coefficient(converter_techs,vintage,activity,commodity,"maxLoad"));

* ==== activate MIP units for MIP dispatch or partial load technologies ====

set converter_hasMinUptime(converter_techs,vintage);
converter_hasMinUptime(converter_techs,vintage)
  $(converter_techParam(converter_techs,vintage,"minUptime")
    and converter_techParam(converter_techs,vintage,"mipDispatch"))
  = yes;

set converter_hasMinDowntime(converter_techs,vintage);
converter_hasMinDowntime(converter_techs,vintage)
  $(converter_techParam(converter_techs,vintage,"minDowntime")
    and converter_techParam(converter_techs,vintage,"mipDispatch"))
  = yes;

* require integer unit counts if online state is to be tracked
converter_techParam(converter_techs,vintage,"mipUnits")
    $(converter_techParam(converter_techs,vintage,"mipDispatch")
      or converter_hasMinLoad(converter_techs, vintage)
      or converter_hasMaxLoad(converter_techs, vintage)
      or converter_hasConstantFluxInActivity(converter_techs,vintage))
    = 1;

* ==== floor mip converter units to integer values ====
converter_capacityParam(nodesModelToCalc,yearsToCalc,converter_techs,"unitsLowerLimit")
    $sum(vintage, converter_techParam(converter_techs,vintage,"mipUnits"))
    = floor(converter_capacityParam(nodesModelToCalc,yearsToCalc,converter_techs,"unitsLowerLimit"));
converter_capacityParam(nodesModelToCalc,yearsToCalc,converter_techs,"unitsUpperLimit")
    $sum(vintage, converter_techParam(converter_techs,vintage,"mipUnits"))
    = ceil(converter_capacityParam(nodesModelToCalc,yearsToCalc,converter_techs,"unitsUpperLimit"));

* ==== scalars and sets for minUptime and minDowntime ====

scalar converter_maxUptimeReq;
converter_maxUptimeReq = smax((converter_techs,vintage), converter_techParam(converter_techs,vintage,"minUptime"));
set uptimeSearchRange(timeModelToCalc);
uptimeSearchRange(timeModelToCalc) = ord(timeModelToCalc) <= converter_maxUptimeReq;

scalar converter_maxDowntimeReq;
converter_maxDowntimeReq = smax((converter_techs,vintage), converter_techParam(converter_techs,vintage,"minDowntime"));
set downtimeSearchRange(timeModelToCalc);
downtimeSearchRange(timeModelToCalc) = ord(timeModelToCalc) <= converter_maxDowntimeReq;

* ==== definition of variables ====

* Initialise variables for unitsBuild
converter_unitsBuild.l(nodesModel,years,converter_techs,vintage)
    $converter_availTech(nodesModel,years,converter_techs,vintage)
    = converter_capacityParam(nodesModel,years,converter_techs,"unitsBuild");
converter_unitsBuild.lo(nodesModel,yearsToCalc,converter_techs,vintage)
    $converter_availTech(nodesModel,yearsToCalc,converter_techs,vintage)
    = converter_unitsBuild.l(nodesModel,yearsToCalc,converter_techs,vintage);
converter_unitsBuild.fx(nodesModel,years,converter_techs,vintage)
    $converter_capacityParam(nodesModel,years,converter_techs,"noExpansion")
    = converter_unitsBuild.l(nodesModel,years,converter_techs,vintage);

* Initialise variables for unitsDecom
converter_unitsDecom.l(nodesModel,years,converter_techs,vintage)
  $(converter_decomTech(nodesModel,years,converter_techs,vintage)
    and years.val < sum(yearsToCalc$(ord(yearsToCalc) = 1), yearsToCalc.val))
  = sum((years_a,years_aa)$(sameas(years-1, years_aa)
                      and years_a.val > years_aa.val - converter_techParam(converter_techs,vintage,'lifeTime')
                      and years_a.val <= years.val - converter_techParam(converter_techs,vintage,'lifeTime')
                      and converter_availTech(nodesModel,years_a,converter_techs,vintage)),
        converter_unitsBuild.l(nodesModel,years_a,converter_techs,vintage));

converter_unitsDecom.l(nodesModel,yearsToCalc,converter_techs,vintage)
  $converter_decomTech(nodesModel,yearsToCalc,converter_techs,vintage)
  = sum(years$
        (years.val < sum(yearsToCalc_a$(ord(yearsToCalc_a) = 1), yearsToCalc_a.val)
          and converter_availTech(nodesModel,years,converter_techs,vintage)
          and years.val > sum(years_a$sameas(years_a, yearsToCalc-1), years_a.val) - converter_techParam(converter_techs,vintage,'lifeTime')
          and years.val <= yearsToCalc.val - converter_techParam(converter_techs,vintage,'lifeTime')),
      converter_unitsBuild.l(nodesModel,years,converter_techs,vintage))
    + sum(yearsToCalc_a$
        (yearsToCalc_a.val < sum(yearsToCalc_aa$(ord(yearsToCalc_aa) > 1), yearsToCalc_a.val)
          and converter_availTech(nodesModel,yearsToCalc_a,converter_techs,vintage)
          and yearsToCalc_a.val > sum(years_a$sameas(years_a, yearsToCalc-1), years_a.val) - converter_techParam(converter_techs,vintage,'lifeTime')
          and yearsToCalc_a.val <= yearsToCalc.val - converter_techParam(converter_techs,vintage,'lifeTime')),
      converter_unitsBuild.l(nodesModel,yearsToCalc_a,converter_techs,vintage));
      ;

converter_unitsDecom.lo(nodesModel,yearsToCalc,converter_techs,vintage)
    $(converter_usedTech(nodesModel,yearsToCalc,converter_techs,vintage)
        and not converter_techParam(converter_techs,vintage,"freeDecom"))
    = converter_unitsDecom.l(nodesModel,yearsToCalc,converter_techs,vintage)

* Calculate planned unit expansion
parameter converter_unitsPlanned(nodesModel,years,converter_techs,vintage);
converter_unitsPlanned(nodesModel,years,converter_techs,vintage) = 0;
loop(years,
  converter_unitsPlanned(nodesModel,years,converter_techs,vintage)
    =
    converter_unitsPlanned(nodesModel,years-1,converter_techs,vintage)
        $converter_usedTech(nodesModel,years-1,converter_techs,vintage)
    + converter_unitsBuild.l(nodesModel,years,converter_techs,vintage)
        $converter_availTech(nodesModel,years,converter_techs,vintage)
    - converter_unitsDecom.l(nodesModel,years,converter_techs,vintage)
        $converter_usedTech(nodesModel,years,converter_techs,vintage);
);

* Set initial state for planned units
converter_unitsTotal.l(nodesModel,years,converter_techs,vintage)
  = converter_unitsPlanned(nodesModel,years,converter_techs,vintage)

* Calculate if planned unit expansion is bounded by upper and lower limits
set converter_unitBoundsFixed(nodesModel,years,converter_techs);
converter_unitBoundsFixed(nodesModel,years,converter_techs)
  $(sum(vintage$converter_usedTech(nodesModel,years,converter_techs,vintage),
        converter_unitsPlanned(nodesModel,years,converter_techs,vintage))
    = converter_capacityParam(nodesModel,years,converter_techs,"unitsUpperLimit")
  and sum(vintage$converter_usedTech(nodesModel,years,converter_techs,vintage),
        converter_unitsPlanned(nodesModel,years,converter_techs,vintage))
    = converter_capacityParam(nodesModel,years,converter_techs,"unitsLowerLimit"))
  = yes;

* Fix unitsBuild, unitsDecom, unitsTotal if levels are predetermined by upper and lower limits
converter_unitsBuild.fx(nodesModel,years,converter_techs,vintage)
  $(converter_availTech(nodesModel,years,converter_techs,vintage)
    and converter_unitBoundsFixed(nodesModel,years,converter_techs))
  = converter_unitsBuild.l(nodesModel,years,converter_techs,vintage);
converter_unitsDecom.fx(nodesModel,years,converter_techs,vintage)
  $(converter_usedTech(nodesModel,years,converter_techs,vintage)
    and converter_unitBoundsFixed(nodesModel,years,converter_techs))
  = converter_unitsDecom.l(nodesModel,years,converter_techs,vintage);
converter_unitsTotal.fx(nodesModel,years,converter_techs,vintage)
  $(converter_usedTech(nodesModel,years,converter_techs,vintage)
    and converter_unitBoundsFixed(nodesModel,years,converter_techs))
  = converter_unitsTotal.l(nodesModel,years,converter_techs,vintage);

converter_unitsOnline_MIP.up(timeModelToCalc,nodesModel,years,converter_techs,vintage)
    $(converter_usedTech(nodesModel,years,converter_techs,vintage)
      and converter_techParam(converter_techs,vintage,"mipDispatch") = 1)
    = converter_capacityParam(nodesModel,years,converter_techs,"unitsUpperLimit");

converter_unitsTotal_MIP.up(nodesModel,years,converter_techs,vintage)
    $(converter_usedTech(nodesModel,years,converter_techs,vintage)
      and converter_techParam(converter_techs,vintage,"mipUnits") = 1)
    = converter_capacityParam(nodesModel,years,converter_techs,"unitsUpperLimit");

converter_unitsUsingActivity_MIP.up(timeModelToCalc,nodesModel,years,converter_techs,vintage,activity)
    $(converter_usedTech(nodesModel,years,converter_techs,vintage)
      and (converter_hasMinLoad(converter_techs, vintage)
            or converter_hasMaxLoad(converter_techs, vintage)))
    = converter_capacityParam(nodesModel,years,converter_techs,"unitsUpperLimit");

* Add parameter for fixing capacities during myopic runs
parameter converter_unitsDelta_upper(nodesModel,years,converter_techs);
parameter converter_unitsDelta_lower(nodesModel,years,converter_techs,vintage);

* ==== declaration of equations ====

equations
  Eq_converter_unitsBalance(nodesModel,years,converter_techs,vintage
    ) "Ensures the units balance over the planning period."
  Eq_converter_unitsFixedDecom(nodesModel,years,converter_techs,vintage
    ) "Restricts the fixed decommissioning of units over the planning period."
  Eq_converter_unitsFreeDecom(nodesModel,years,converter_techs,vintage
    ) "Restricts the free decommissioning of units over the planning period."
  Eq_converter_unitsUpperLimit(nodesModel,years,converter_techs
    ) "Upper bound for the total number of units."
  Eq_converter_unitsLowerLimit(nodesModel,years,converter_techs
    ) "Lower bound for the total number of units."
  Eq_converter_unitsFixedLimit(nodesModel,years,converter_techs
    ) "Fixed bound for the total number of units."
  Eq_converter_unitsTotalMIP(nodesModel,years,converter_techs,vintage
    ) "Fixes the total number of units to the corresponding integer variable."
  Eq_converter_unitsOnlineMIP(timeModel,nodesModel,years,converter_techs,vintage
    ) "Fixes the number of online units to the corresponding integer variable."

  Eq_converter_activityLowerLimit(timeModel,nodesModel,years,converter_techs,vintage
    ) "Lower limit on the activity."
  Eq_converter_activityUpperLimit(timeModel,nodesModel,years,converter_techs,vintage
    ) "Upper limit on the activity."

  Eq_converter_activityFixedLimit(timeModel,nodesModel,years,converter_techs,vintage
    ) "Fixed limit on the activity."
  Eq_converter_rampPos(timeModel,nodesModel,years,converter_techs,vintage
    ) "Positive ramping of unit activity."
  Eq_converter_rampNeg(timeModel,nodesModel,years,converter_techs,vintage
    ) "Negative ramping of unit activity."
  Eq_converter_rampLimit(timeModel,nodesModel,years,converter_techs,vintage
    ) "Restrict ramping up of unit activity."

  Eq_converter_unitsOnline(timeModel,nodesModel,years,converter_techs,vintage
    ) "Set online units to total number of operational units."
  Eq_converter_unitsOnlineUC(timeModel,nodesModel,years,converter_techs,vintage
    ) "Allow shutting down units."
  Eq_converter_activityStartups(timeModel,nodesModel,years,converter_techs,vintage
    ) "Variable counting the number of unit startups."
  Eq_converter_activityShutdowns(timeModel,nodesModel,years,converter_techs,vintage
    ) "Variable tracking the number of unit shutdowns."
  Eq_converter_limitStartups(nodesModel,years,converter_techs,vintage
    ) "Limit the number of startup cycles a unit can perform."
  Eq_converter_minUptime(timeModel,nodesModel,years,converter_techs,vintage
    ) "Require recently started units to remain online for their respective minimum uptime."
  Eq_converter_minDowntime(timeModel,nodesModel,years,converter_techs,vintage
    ) "Require recently shut down units to remain offline for their respective minimum downtime."
  Eq_converter_activityUpperLimitDetailedPartLoadMinReq(timeModel,nodesModel,years,converter_techs,vintage,activity
    ) "Limit activity coefficients of activities with a given load requirement or stricter to the number of units in such modes."
  Eq_converter_activityUpperLimitDetailedPartLoadMaxReq(timeModel,nodesModel,years,converter_techs,vintage,activity
    ) "Limit activity coefficients of activities with a given load requirement or stricter to the number of units in such modes."
  Eq_converter_activityLowerLimitDetailedPartLoadMinReq(timeModel,nodesModel,years,converter_techs,vintage,activity
    ) "Enforce sufficient activity coefficients of activities to justify all active modes."
  Eq_converter_noOnlineIdle(timeModel,nodesModel,years,converter_techs,vintage
    ) "Prevent the circumvention of requirements by keeping units online without using any modes."
  Eq_converter_noOnlineIdleDetailedPartLoad(timeModel,nodesModel,years,converter_techs,vintage
    ) "Prevent the circumvention of requirements by keeping units online without using any modes. Allow more than one mode activation per unit."
  Eq_converter_activityUpperLimitPartLoad(timeModel,nodesModel,years,converter_techs,vintage,activity
    ) "Limit usage of an activity to corresponding units."
  Eq_converter_activityLowerLimitPartLoad(timeModel,nodesModel,years,converter_techs,vintage,activity
    ) "Limit usage of an activity to corresponding units."
  Eq_converter_activityModeLimit(timeModel,nodesModel,years,converter_techs,vintage,activity
    ) "Limit the number of units in one particular mode to the number of operational units."
  ;

* ==== equation definition ====
* // ## Equations
* // ### Converter Units Balance
* // Ensures that the total units are consistent with the built and decommissioned units.
* // {Eq_converter_unitsBalance}
Eq_converter_unitsBalance(nodesModelSel,yearsSel,converter_techs,vintage)
    $((converter_usedTech(nodesModelSel,yearsSel,converter_techs,vintage)
          or sum(years$sameas(years,yearsSel), converter_usedTech(nodesModelSel,years-1,converter_techs,vintage)))
        and not converter_unitBoundsFixed(nodesModelSel,yearsSel,converter_techs))
    ..
    converter_unitsTotal(nodesModelSel,yearsSel,converter_techs,vintage)
    =e=
    sum(yearsToCalc$(ord(yearsToCalc) = 1 and sameas(yearsToCalc, yearsSel)),
      sum(years$sameas(years, yearsToCalc),
        converter_unitsTotal(nodesModelSel,years-1,converter_techs,vintage)
          $converter_usedTech(nodesModelSel,years-1,converter_techs,vintage)))
    + sum((yearsToCalc)$(ord(yearsToCalc) > 1 and sameas(yearsToCalc, yearsSel)),
      converter_unitsTotal(nodesModelSel,yearsToCalc-1,converter_techs,vintage)
        $converter_usedTech(nodesModelSel,yearsToCalc-1,converter_techs,vintage))
    + converter_unitsBuild(nodesModelSel,yearsSel,converter_techs,vintage)
        $converter_availTech(nodesModelSel,yearsSel,converter_techs,vintage)
    - converter_unitsDecom(nodesModelSel,yearsSel,converter_techs,vintage)
        $converter_usedTech(nodesModelSel,yearsSel,converter_techs,vintage);

* // ### Converter Units Fixed Decommission
* // Restricts the fixed decommissioning of units over the planning period.
* // {Eq_converter_unitsFixedDecom}
Eq_converter_unitsFixedDecom(nodesModelSel,yearsSel,converter_techs,vintage)
    $(converter_decomTech(nodesModelSel,yearsSel,converter_techs,vintage)
        and not converter_techParam(converter_techs,vintage,"freeDecom") = 1
        and not converter_unitBoundsFixed(nodesModelSel,yearsSel,converter_techs))
    ..
    converter_unitsDecom(nodesModelSel,yearsSel,converter_techs,vintage)
    + converter_unitsDelta_lower(nodesModelSel,yearsSel,converter_techs,vintage)
    =e=
    sum(years$
        (converter_availTech(nodesModelSel,years,converter_techs,vintage)
          and years.val > sum(yearsToCalc$sameas(yearsToCalc+1, yearsSel), yearsToCalc.val) - converter_techParam(converter_techs,vintage,'lifeTime')
          and years.val <= yearsSel.val - converter_techParam(converter_techs,vintage,'lifeTime')),
      converter_unitsBuild(nodesModelSel,years,converter_techs,vintage));

* // ### Converter Units Free Decommission
* // Restricts the free decommissioning of units over the planning period.
* // {Eq_converter_unitsFreeDecom}
Eq_converter_unitsFreeDecom(nodesModelSel,yearsSel,converter_techs,vintage)
    $((converter_decomTech(nodesModelSel,yearsSel,converter_techs,vintage)
        or converter_usedTech(nodesModelSel,yearsSel,converter_techs,vintage))
        and converter_techParam(converter_techs,vintage,"freeDecom") = 1)
    ..
    sum(years$
          ((converter_decomTech(nodesModelSel,years,converter_techs,vintage)
            or converter_usedTech(nodesModelSel,years,converter_techs,vintage))
            and years.val < sum(yearsToCalc_a$(ord(yearsToCalc_a)=1), yearsToCalc_a.val)),
        converter_unitsDecom(nodesModelSel,years,converter_techs,vintage))
    + sum(yearsToCalc$
          ((converter_decomTech(nodesModelSel,yearsToCalc,converter_techs,vintage)
            or converter_usedTech(nodesModelSel,yearsToCalc,converter_techs,vintage))
            and yearsToCalc.val >= sum(yearsToCalc_a$(ord(yearsToCalc_a)=1), yearsToCalc_a.val)
            and yearsToCalc.val <= yearsSel.val),
        converter_unitsDecom(nodesModelSel,yearsToCalc,converter_techs,vintage))
    =g=
    sum(years$
          (converter_availTech(nodesModelSel,years,converter_techs,vintage)
            and years.val < sum(yearsToCalc$(ord(yearsToCalc)=1), yearsToCalc.val) - converter_techParam(converter_techs,vintage,'lifeTime')),
        converter_unitsBuild(nodesModelSel,years,converter_techs,vintage))
    + sum(yearsToCalc$
          (converter_availTech(nodesModelSel,yearsToCalc,converter_techs,vintage)
            and yearsToCalc.val >= sum(yearsToCalc_a$(ord(yearsToCalc_a)=1), yearsToCalc_a.val)
            and yearsToCalc.val <= yearsSel.val - converter_techParam(converter_techs,vintage,'lifeTime')),
        converter_unitsBuild(nodesModelSel,yearsToCalc,converter_techs,vintage));

* // ### Converter Units Upper Limit
* // Upper bound for the total number of units.
* // {Eq_converter_unitsUpperLimit}
Eq_converter_unitsUpperLimit(nodesModelSel,yearsSel,converter_techs)
    $(converter_capacityParam(nodesModelSel,yearsSel,converter_techs,'unitsUpperLimit') >= 0
        and converter_capacityParam(nodesModelSel,yearsSel,converter_techs,'unitsUpperLimit') < +inf
        and converter_capacityParam(nodesModelSel,yearsSel,converter_techs,'unitsUpperLimit')
            <> converter_capacityParam(nodesModelSel,yearsSel,converter_techs,'unitsLowerLimit'))
    ..
    sum(vintage$converter_usedTech(nodesModelSel,yearsSel,converter_techs,vintage),
        converter_unitsTotal(nodesModelSel,yearsSel,converter_techs,vintage))
    =l=
    converter_capacityParam(nodesModelSel,yearsSel,converter_techs,"unitsUpperLimit");

* // ### Converter Units Lower Limit
* // Lower bound for the total number of units.
* // {Eq_converter_unitsLowerLimit}
Eq_converter_unitsLowerLimit(nodesModelSel,yearsSel,converter_techs)
    $(converter_capacityParam(nodesModelSel,yearsSel,converter_techs,'unitsLowerLimit') > 0
        and converter_capacityParam(nodesModelSel,yearsSel,converter_techs,'unitsLowerLimit')
            <> converter_capacityParam(nodesModelSel,yearsSel,converter_techs,'unitsUpperLimit'))
    ..
    sum(vintage$converter_usedTech(nodesModelSel,yearsSel,converter_techs,vintage),
        converter_unitsTotal(nodesModelSel,yearsSel,converter_techs,vintage))
    =g=
    converter_capacityParam(nodesModelSel,yearsSel,converter_techs,"unitsLowerLimit");

Eq_converter_unitsFixedLimit(nodesModelSel,yearsSel,converter_techs)
    $(converter_capacityParam(nodesModelSel,yearsSel,converter_techs,'unitsLowerLimit')
        = converter_capacityParam(nodesModelSel,yearsSel,converter_techs,'unitsUpperLimit'))
    ..
    sum(vintage$converter_usedTech(nodesModelSel,yearsSel,converter_techs,vintage),
        converter_unitsTotal(nodesModelSel,yearsSel,converter_techs,vintage))
    =e=
    converter_capacityParam(nodesModelSel,yearsSel,converter_techs,"unitsUpperLimit");

* // ### Converter Units Total MIP
* // Fixes the total number of units to the corresponding integer variable.
* // {Eq_converter_unitsTotalMIP}
Eq_converter_unitsTotalMIP(nodesModelSel,yearsSel,converter_techs,vintage)
    $( converter_usedTech(nodesModelSel,yearsSel,converter_techs,vintage)
        and converter_techParam(converter_techs,vintage,"mipUnits") = 1 )
    ..
    converter_unitsTotal(nodesModelSel,yearsSel,converter_techs,vintage)
    =e=
    converter_unitsTotal_MIP(nodesModelSel,yearsSel,converter_techs,vintage);

* // ### Converter Units Online MIP
* // Fixes the number of online units to the corresponding integer variable.
* // {Eq_converter_unitsOnlineMIP}
Eq_converter_unitsOnlineMIP(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage)
    $( converter_usedTech(nodesModelSel,yearsSel,converter_techs,vintage)
        and converter_techParam(converter_techs,vintage,"mipDispatch") = 1 )
    ..
    converter_unitsOnline(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage)
    =e=
    converter_unitsOnline_MIP(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage);

* // ### Converter Activity Lower Limit
* // Lower limit on the activity.
* // {Eq_converter_activityLowerLimit}
Eq_converter_activityLowerLimit(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage)
    $(converter_usedTech(nodesModelSel,yearsSel,converter_techs,vintage)
        and converter_activityProfile(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,"lower") > 0
        and converter_activityProfile(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,"lower")
             <> converter_activityProfile(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,"upper")
        and not converter_hasMinLoad(converter_techs, vintage))
    ..
    sum(activity$converter_usedAct(converter_techs,vintage,activity),
        converter_activity(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,activity)
    )
    =g=
    converter_activityProfile(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,"lower")
$iftheni.pips %method%==pips
    * converter_unitsTotal(nodesModelSel,yearsSel,converter_techs,vintage);
$else.pips
    * converter_unitsOnline(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage);
$endif.pips

* // ### Converter Activity Upper Limit
* // Upper limit on the activity.
* // {Eq_converter_activityUpperLimit}
Eq_converter_activityUpperLimit(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage)
    $(converter_usedTech(nodesModelSel,yearsSel,converter_techs,vintage)
        and converter_activityProfile(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,"upper") >= 0
        and converter_activityProfile(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,"upper")
             <> converter_activityProfile(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,"lower"))
    ..
    sum(activity$converter_usedAct(converter_techs,vintage,activity),
        converter_activity(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,activity)
    )
    =l=
    converter_activityProfile(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,"upper")
$iftheni.pips %method%==pips
    * converter_unitsTotal(nodesModelSel,yearsSel,converter_techs,vintage);
$else.pips
    * converter_unitsOnline(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage);
$endif.pips

Eq_converter_activityFixedLimit(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage)
    $(converter_usedTech(nodesModelSel,yearsSel,converter_techs,vintage)
        and converter_activityProfile(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,"lower")
             = converter_activityProfile(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,"upper"))
    ..
    sum(activity$converter_usedAct(converter_techs,vintage,activity),
        converter_activity(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,activity)
    )
    =e=
    converter_activityProfile(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,"upper")
$iftheni.pips %method%==pips
    * converter_unitsTotal(nodesModelSel,yearsSel,converter_techs,vintage);
$else.pips
    * converter_unitsOnline(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage);
$endif.pips

* // ### Converter Positive Ramping
* // Positive ramping of unit activity.
* // {Eq_converter_rampPos}
Eq_converter_rampPos(timeModelSel(timeModelToCalc),nodesModelSel,yearsSel,converter_techs,vintage)
  $( converter_usedTech(nodesModelSel,yearsSel,converter_techs,vintage)
      and converter_useRampPos(nodesModelSel,yearsSel,converter_techs,vintage))
    ..
    converter_rampPos(timeModelToCalc,nodesModelSel,yearsSel,converter_techs,vintage)
    =g=
    sum(activity$converter_usedAct(converter_techs,vintage,activity),
            converter_activity(timeModelToCalc,nodesModelSel,yearsSel,converter_techs,vintage,activity)
            - converter_activity(timeModelToCalc--1,nodesModelSel,yearsSel,converter_techs,vintage,activity));

* // ### Converter Negative Ramping
* // Negative ramping of unit activity.
* // {Eq_converter_rampNeg}
Eq_converter_rampNeg(timeModelSel(timeModelToCalc),nodesModelSel,yearsSel,converter_techs,vintage)
  $( converter_usedTech(nodesModelSel,yearsSel,converter_techs,vintage)
      and converter_useRampNeg(nodesModelSel,yearsSel,converter_techs,vintage))
    ..
    converter_rampNeg(timeModelToCalc,nodesModelSel,yearsSel,converter_techs,vintage)
    =g=
    - sum(activity$converter_usedAct(converter_techs,vintage,activity),
            converter_activity(timeModelToCalc,nodesModelSel,yearsSel,converter_techs,vintage,activity)
            - converter_activity(timeModelToCalc--1,nodesModelSel,yearsSel,converter_techs,vintage,activity));

* // ### Converter Ramping Limit
* // Restrict ramping up of unit activity.
* // {Eq_converter_rampLimit}
Eq_converter_rampLimit(timeModelSel(timeModelToCalc),nodesModelSel,yearsSel,converter_techs,vintage)
  $( converter_usedTech(nodesModelSel,yearsSel,converter_techs,vintage)
        and converter_techParam(converter_techs,vintage,"activityRampLimit") > 0)
    ..
    sum(activity$converter_usedAct(converter_techs,vintage,activity),
            converter_activity(timeModelToCalc,nodesModelSel,yearsSel,converter_techs,vintage,activity)
            - converter_activity(timeModelToCalc--1,nodesModelSel,yearsSel,converter_techs,vintage,activity))
  =l=
  converter_techParam(converter_techs,vintage,"activityRampLimit")
$iftheni.pips %method%==pips
    * converter_unitsTotal(nodesModelSel,yearsSel,converter_techs,vintage);
$else.pips
    * converter_unitsOnline(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage);
$endif.pips


$iftheni.pips %method%==pips
$else.pips

* // ### Converter MIP Units Online
* // Restrict ramping up of unit activity.
* // {Eq_converter_unitsOnline}
Eq_converter_unitsOnline(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage)
    $( converter_usedTech(nodesModelSel,yearsSel,converter_techs,vintage)
        and converter_techParam(converter_techs,vintage,"mipDispatch") = 0 )
  ..
  converter_unitsOnline(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage)
  =e=
  converter_unitsTotal(nodesModelSel,yearsSel,converter_techs,vintage)
    ;

* // ### Converter MIP Units Shutting Down
* // Allow shutting down units.
* // {Eq_converter_unitsOnlineUC}
Eq_converter_unitsOnlineUC(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage)
    $( converter_usedTech(nodesModelSel,yearsSel,converter_techs,vintage)
        and converter_techParam(converter_techs,vintage,"mipDispatch") = 1 )
  ..
  converter_unitsOnline(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage)
  =l=
  converter_unitsTotal(nodesModelSel,yearsSel,converter_techs,vintage)
    ;
$endif.pips

* // ### Converter Activity Startups
* // Variable counting the number of unit startups.
* // {Eq_converter_activityStartups}
Eq_converter_activityStartups(timeModelSel(timeModelToCalc),nodesModelSel,yearsSel,converter_techs,vintage)
  $( converter_usedTech(nodesModelSel,yearsSel,converter_techs,vintage)
        and converter_techParam(converter_techs,vintage,"mipDispatch") = 1 )
    ..
    converter_unitStartups(timeModelToCalc,nodesModelSel,yearsSel,converter_techs,vintage)
    =g=
	converter_unitsOnline(timeModelToCalc,nodesModelSel,yearsSel,converter_techs,vintage)
    - converter_unitsOnline(timeModelToCalc--1,nodesModelSel,yearsSel,converter_techs,vintage);

* // ### Converter Activity Shutdowns
* // Variable tracking the number of unit shutdowns.
* // {Eq_converter_activityShutdowns}
Eq_converter_activityShutdowns(timeModelSel(timeModelToCalc),nodesModelSel,yearsSel,converter_techs,vintage)
  $( converter_usedTech(nodesModelSel,yearsSel,converter_techs,vintage)
        and converter_techParam(converter_techs,vintage,"mipDispatch") = 1
        and converter_techParam(converter_techs,vintage,"minDowntime") > 0)
    ..
    converter_unitShutdowns(timeModelToCalc,nodesModelSel,yearsSel,converter_techs,vintage)
    =g=
	converter_unitsOnline(timeModelToCalc--1,nodesModelSel,yearsSel,converter_techs,vintage)
    - converter_unitsOnline(timeModelToCalc,nodesModelSel,yearsSel,converter_techs,vintage);

* // ### Converter Activity Startup Limit
* // Limit the number of startup cycles a unit can perform.
* // {Eq_converter_limitStartups}
Eq_converter_limitStartups(nodesModelSel,yearsSel,converter_techs,vintage)
  $( converter_usedTech(nodesModelSel,yearsSel,converter_techs,vintage)
        and converter_techParam(converter_techs,vintage,"mipDispatch") = 1
    and converter_techParam(converter_techs,vintage,"startupLimit") > 0)
    ..
    sum(timeModelToCalc, converter_unitStartups(timeModelToCalc,nodesModelSel,yearsSel,converter_techs,vintage))
  =l=
  converter_techParam(converter_techs,vintage,"startupLimit")
  * converter_unitsTotal(nodesModelSel,yearsSel,converter_techs,vintage);

alias(timeModelToCalc,ttc);

* // ### Converter Units Minimum Uptime
* // Require recently started units to remain online for their respective minimum uptime.
* // {Eq_converter_minUptime}
Eq_converter_minUptime(timeModelSel(timeModelToCalc),nodesModelSel,yearsSel,converter_techs,vintage)
  $( converter_usedTech(nodesModelSel,yearsSel,converter_techs,vintage)
        and converter_techParam(converter_techs,vintage,"mipDispatch") = 1
    and converter_techParam(converter_techs,vintage,"minUptime") > 0)
    ..
    sum(ttc$[uptimeSearchRange(ttc) and ord(ttc)<=converter_techParam(converter_techs,vintage,"minUptime")],
        converter_unitStartups(ttc+[ord(timeModelToCalc)-converter_techParam(converter_techs,vintage,"minUptime")],
            nodesModelSel,yearsSel,converter_techs,vintage))
	=l=
	converter_unitsOnline(timeModelToCalc,nodesModelSel,yearsSel,converter_techs,vintage);

* // ### Converter Units Minimum Downtime
* // Require recently shut down units to remain offline for their respective minimum downtime.
* // {Eq_converter_minDowntime}
Eq_converter_minDowntime(timeModelSel(timeModelToCalc),nodesModelSel,yearsSel,converter_techs,vintage)
  $( converter_usedTech(nodesModelSel,yearsSel,converter_techs,vintage)
        and converter_techParam(converter_techs,vintage,"mipDispatch") = 1
    and converter_techParam(converter_techs,vintage,"minDowntime") > 0)
    ..
    sum(ttc$[downtimeSearchRange(ttc) and ord(ttc)<=converter_techParam(converter_techs,vintage,"minDowntime")],
        converter_unitShutdowns(ttc+[ord(timeModelToCalc)-converter_techParam(converter_techs,vintage,"minDowntime")],
            nodesModelSel,yearsSel,converter_techs,vintage))
	=l=
	converter_unitsTotal(nodesModelSel,yearsSel,converter_techs,vintage)
    - converter_unitsOnline(timeModelToCalc,nodesModelSel,yearsSel,converter_techs,vintage);

* This equation is meant to cause the units to activate particular mode counters to gain access to the corresponding activities.
* Activities are allowed to be used on units with stricter activity requirements but not the other way around.
alias(activity, act);

* // ### Converter MIP Units Activity Upper Limit Minimum Required Load
* // Limit activity coefficients of activities with a given load requirement or stricter to the number of units in such modes.
* // {Eq_converter_activityUpperLimitDetailedPartLoadMinReq}
Eq_converter_activityUpperLimitDetailedPartLoadMinReq(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,activity)
  $( converter_usedTechAct(nodesModelSel,yearsSel,converter_techs,vintage,activity)
        and converter_hasMinLoad(converter_techs, vintage)
        and converter_techParam(converter_techs,vintage,"mipDetailedPartialLoad"))
    ..
    sum(act$(converter_activityRequirements(converter_techs,vintage,act,"minLoad")
                >= converter_activityRequirements(converter_techs,vintage,activity,"minLoad")
             and converter_usedAct(converter_techs,vintage,act)),
    converter_activity(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,act))
  =l=
  converter_activityProfile(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,"upper")
    * sum(act$(converter_activityRequirements(converter_techs,vintage,act,"minLoad")
                >= converter_activityRequirements(converter_techs,vintage,activity,"minLoad")
               and converter_usedAct(converter_techs,vintage,act)),
    converter_unitsUsingActivity_MIP(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,act));

* // ### Converter MIP Units Activity Upper Limit Maximum Required Load
* // Limit activity coefficients of activities with a given load requirement or stricter to the number of units in such modes.
* // {Eq_converter_activityUpperLimitDetailedPartLoadMaxReq}
Eq_converter_activityUpperLimitDetailedPartLoadMaxReq(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,activity)
  $( converter_usedTechAct(nodesModelSel,yearsSel,converter_techs,vintage,activity)
        and converter_hasMaxLoad(converter_techs,vintage)
        and converter_techParam(converter_techs,vintage,"mipDetailedPartialLoad"))
    ..
    sum(act$(converter_activityRequirements(converter_techs,vintage,act,"maxLoad")
                <= converter_activityRequirements(converter_techs,vintage,activity,"maxLoad")
             and converter_usedAct(converter_techs,vintage,act)),
    converter_activity(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,act))
  =l=
    sum(act$(converter_activityRequirements(converter_techs,vintage,act,"maxLoad")
                <= converter_activityRequirements(converter_techs,vintage,activity,"maxLoad")
             and converter_usedAct(converter_techs,vintage,act)),
    min(converter_activityRequirements(converter_techs,vintage,act,"maxLoad"),
          converter_activityProfile(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,"upper"))
          * converter_unitsUsingActivity_MIP(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,act));

* This equation is meant to enforce the lower limit requirements of activities in use.
* Loads produced by activities with less strict requirements can contribute to the minimum load requirement of strict activities but not the other way around,
* because these activities can only run on units in the respectively strict activation state.

* // ### Converter MIP Units Activity Lower Limit Minimum Required Load
* // Enforce sufficient activity coefficients of activities to justify all active modes.
* // {Eq_converter_activityLowerLimitDetailedPartLoadMinReq}
Eq_converter_activityLowerLimitDetailedPartLoadMinReq(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,activity)
  $( converter_usedTechAct(nodesModelSel,yearsSel,converter_techs,vintage,activity)
        and converter_hasMinLoad(converter_techs, vintage)
        and converter_techParam(converter_techs,vintage,"mipDetailedPartialLoad"))
    ..
    sum(act$(converter_activityRequirements(converter_techs,vintage,act,"minLoad")
                <= converter_activityRequirements(converter_techs,vintage,activity,"minLoad")
             and converter_usedAct(converter_techs,vintage,act)),
    converter_activity(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,act))
  =g=
    sum(act$(converter_activityRequirements(converter_techs,vintage,act,"minLoad")
                <= converter_activityRequirements(converter_techs,vintage,activity,"minLoad")
             and converter_usedAct(converter_techs,vintage,act)),
    max(converter_activityRequirements(converter_techs,vintage,act,"minLoad"),
          converter_activityProfile(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,"lower"))
          * converter_unitsUsingActivity_MIP(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,act));


* Simple one-activity-per-unit partial load equations
* // ### Converter MIP Units Activity Upper Limit Partial Load Balance
* // Enforce MIP units partial load upper limit.
* // {Eq_converter_activityUpperLimitPartLoad}
Eq_converter_activityUpperLimitPartLoad(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,activity)
  $( converter_usedTechAct(nodesModelSel,yearsSel,converter_techs,vintage,activity)
        and (converter_hasMinLoad(converter_techs, vintage)
              or converter_hasMaxLoad(converter_techs, vintage)
              or converter_hasConstantFluxInActivity(converter_techs, vintage))
        and not converter_techParam(converter_techs,vintage,"mipDetailedPartialLoad"))
    ..
    converter_activity(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,activity)
  =l=
  min(converter_activityRequirements(converter_techs,vintage,activity,"maxLoad"),
        converter_activityProfile(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,"upper"))
        * converter_unitsUsingActivity_MIP(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,activity);

* // ### Converter MIP Units Activity Lower Limit Partial Load Balance
* // Enforce MIP units partial load lower limit.
* // {Eq_converter_activityLowerLimitPartLoad}
Eq_converter_activityLowerLimitPartLoad(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,activity)
  $( converter_usedTechAct(nodesModelSel,yearsSel,converter_techs,vintage,activity)
        and (converter_hasMinLoad(converter_techs, vintage)
              or converter_hasMaxLoad(converter_techs, vintage))
        and not converter_techParam(converter_techs,vintage,"mipDetailedPartialLoad"))
    ..
    converter_activity(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,activity)
  =g=
  max(converter_activityRequirements(converter_techs,vintage,activity,"minLoad"),
        converter_activityProfile(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,"lower"))
        * converter_unitsUsingActivity_MIP(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,activity);

* // ### Converter MIP Units Idle Online Units
* // Counts idle online units at every time step
* // {Eq_converter_noOnlineIdle}
Eq_converter_noOnlineIdle(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage)
  $( converter_usedTech(nodesModelSel,yearsSel,converter_techs,vintage)
        and (converter_hasMinLoad(converter_techs, vintage)
            or converter_hasMaxLoad(converter_techs, vintage)
            or converter_hasConstantFluxInActivity(converter_techs,vintage))
        and not converter_techParam(converter_techs,vintage,"mipDetailedPartialLoad"))
    ..
    converter_unitsOnline(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage)
    =e=
    sum(activity$converter_usedAct(converter_techs,vintage,activity),
          converter_unitsUsingActivity_MIP(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,activity));

* // ### Converter MIP Units Idle Online Units Partial load
* // Counts idle online units at every time step
* // {Eq_converter_noOnlineIdleDetailedPartLoad}
Eq_converter_noOnlineIdleDetailedPartLoad(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage)
  $( converter_usedTech(nodesModelSel,yearsSel,converter_techs,vintage)
        and (converter_hasMinLoad(converter_techs, vintage)
                or converter_hasMaxLoad(converter_techs, vintage))
        and converter_techParam(converter_techs,vintage,"mipDetailedPartialLoad"))
    ..
    converter_unitsOnline(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage)
    =l=
    sum(activity$converter_usedAct(converter_techs,vintage,activity),
          converter_unitsUsingActivity_MIP(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,activity));

* // ### Converter MIP activity model limit
* // Converter activity model limit
* // {Eq_converter_activityModeLimit}
Eq_converter_activityModeLimit(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,activity)
  $( converter_usedTechAct(nodesModelSel,yearsSel,converter_techs,vintage,activity)
        and (converter_hasMinLoad(converter_techs, vintage)
              or converter_hasMaxLoad(converter_techs, vintage))
        and converter_techParam(converter_techs,vintage,"mipDetailedPartialLoad"))
    ..
    converter_unitsOnline(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage)
    =g=
    converter_unitsUsingActivity_MIP(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,activity);

* ==== model definition ====

Model M_converter
/
  Eq_converter_unitsBalance
  Eq_converter_unitsFixedDecom
  Eq_converter_unitsFreeDecom
  Eq_converter_unitsUpperLimit
  Eq_converter_unitsLowerLimit
  Eq_converter_unitsFixedLimit
  Eq_converter_unitsTotalMIP
  Eq_converter_unitsOnlineMIP
  Eq_converter_activityUpperLimit
  Eq_converter_activityLowerLimit
  Eq_converter_activityFixedLimit
  Eq_converter_rampPos
  Eq_converter_rampNeg
  Eq_converter_rampLimit
$iftheni.pips %method%==pips
$else.pips
  Eq_converter_unitsOnline
  Eq_converter_unitsOnlineUC
$endif.pips
  Eq_converter_activityStartups
  Eq_converter_activityShutdowns
  Eq_converter_limitStartups
  Eq_converter_minUptime
  Eq_converter_minDowntime
  Eq_converter_activityUpperLimitDetailedPartLoadMinReq
  Eq_converter_activityUpperLimitDetailedPartLoadMaxReq
  Eq_converter_activityLowerLimitDetailedPartLoadMinReq
  Eq_converter_activityUpperLimitPartLoad
  Eq_converter_activityLowerLimitPartLoad
  Eq_converter_noOnlineIdle
  Eq_converter_noOnlineIdleDetailedPartLoad
  Eq_converter_activityModeLimit
/;
