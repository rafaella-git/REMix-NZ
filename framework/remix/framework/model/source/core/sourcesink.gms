* SPDX-FileCopyrightText: Copyright (c) 2023-2024 German Aerospace Center (DLR)
* SPDX-License-Identifier: BSD-3-Clause

* // # core_sourcesink
* // The equations in this file describe the sources and sinks in the model.

* // ## Variables
* // {special_table_sourcesink_variables}
$offListing
$onEmpty

* ==== declaration of variables ====
variables
  sourcesink_flow(timeModel,nodesModel,years,sourcesink_techs,commodity
    ) "Flows from sources across the system boundary"
;
* // ## Input Files
** // INPUT: sourcesink_config | IAO:0000100:data set
* // ### sourcesink_config
* // Title: Sources and Sinks Configuration
* // Description: Configuration for sources and sinks with upper, fixed, and lower profiles and with upper, fixed, and lower sums given as bool values.
* // {table_sourcesink_config}
set pc_sourcesink_config
    /
    usesFixedProfile            "Uses fixed profile timeseries | Fixed profile for flows from source / to sink | | boolean | {none}| OEO_00000339:program parameter"
    usesLowerProfile            "Uses lower profile timeseries | Lower bound profile for flows from source  / to sink | | boolean | {none} | OEO_00000339:program parameter"
    usesUpperProfile            "Uses upper profile timeseries | Upper bound profile for flows from source / to sink | | boolean | {none} | OEO_00000339:program parameter"
    usesFixedSum                "Uses fixed sum | Fixed annual sum for flows from source / to sink | | boolean | {none} | OEO_00000339:program parameter"
    usesLowerSum                "Uses lower sum | Lower bound annual sum for flows from source / to sink | | boolean | {none} | OEO_00000339:program parameter"
    usesUpperSum                "Uses upper sum | Upper bound annual sum for flows from source / to sink | | boolean | {none} | OEO_00000339:program parameter"
    scaleFixProfileWithFixSum   "Scale fixed profile with sum | Given fixed profile will be preprocessed to fit with given fixed annual sum | | boolean | {none} | OEO_00000339:program parameter"
    scaleLowProfileWithLowSum   "Scale lower profile with sum | Given lower profile will be preprocessed to fit with given lower annual sum | | boolean | {none} | OEO_00000339:program parameter"
    scaleUpProfileWithUpSum     "Scale upper profile with sum | Given upper profile will be preprocessed to fit with given upper annual sum | | boolean | {none} | OEO_00000339:program parameter"

    /;
table sourcesink_configIn(nodesData,years,sourcesink_techs,commodity,pc_sourcesink_config)
$onDelim
$if exist "%instancedir%/sourcesink_config.csv" $include "%instancedir%/sourcesink_config.csv"
$offDelim
$if not exist "%instancedir%/sourcesink_config.csv" $log "No config for sources and sinks included"
;

parameter sourcesink_config(nodesModel,years,sourcesink_techs,commodity,pc_sourcesink_config);
$batinclude %aggregateNodes% sourcesink_config(nodesModel,years,sourcesink_techs,commodity,pc_sourcesink_config) sourcesink_configIn(nodesData,years,sourcesink_techs,commodity,pc_sourcesink_config) smax

* The table headers in the parameter input sets represent the following:
* "Title | Description| Constraints | Type |  Units | Ontology"

** // INPUT: sourcesink_profile | IAO:0000100:data set
* // ### sourcesink_profile
* // Title: Sources and Sinks Profiles
* // Description: Parameter for the source and sink profiles given as either upper, fixed or lower profile.
* // {table_sourcesink_profile}
set pc_sourcesink_profile
    /
    fixed          "Fixed profile | Fixed value for flows from source / to sink |"
    lower          "Lower profile | Lower bound for flows from source / to sink |"
    upper          "Upper profile | Upper bound for flows from source / to sink |"
    /;
table sourcesink_profileLoad(nodesData,years,sourcesink_techs,commodity,pc_sourcesink_profile,timeData)
$onDelim
$if exist "%instancedir%/sourcesink_profile.csv" $include "%instancedir%/sourcesink_profile.csv"
$if not exist "%instancedir%/sourcesink_profile.csv" $if exist "%instancedir%/sourcesink_timeseries.csv" $include "%instancedir%/sourcesink_timeseries.csv"
$if not exist "%instancedir%/sourcesink_profile.csv" $if exist "%instancedir%/sourcesink_timeseries.csv" $log "Project is still using the old file name sourcesink_timeseries.dat, please change it to sourcesink_profile.dat in the future."
$offdelim
$if not exist "%instancedir%/sourcesink_profile.csv" $if not exist "%instancedir%/sourcesink_timeseries.csv" $log "No profile for sources and sinks included."
;

parameter sourcesink_profileIn(timeData,nodesData,years,sourcesink_techs,commodity,pc_sourcesink_profile);
option sourcesink_profileIn < sourcesink_profileLoad;
option clear = sourcesink_profileLoad;

** // INPUT: sourcesink_annualSum | IAO:0000100:data set
* // ### sourcesink_annualSum
* // Title: Sources and Sinks Annual Sums
* // Description: Value for the upper, fixed or lower annual sums for the sources and sinks.
* // {table_sourcesink_annualSum}
set pc_sourcesink_annualSum
    /
    fixed          "Fixed annual sum | Fixed annual sum for flows from source / to sink | excludes:['upper','lower'] | number | {rate} | OEO_00140056:flow potential"
    lower          "Lower annual sum | Lower bound annual sum for flows from source / to sink | lessThan:['upper'] | number | {rate} | OEO_00140056:flow potential"
    upper          "Upper annual sum | Upper bound annual sum for flows from source / to sink | default:'inf' | number | {rate} | OEO_00140056:flow potential"
    /;
table sourcesink_annualSumIn(nodesData,years,sourcesink_techs,commodity,pc_sourcesink_annualSum)
$onDelim
$if exist "%instancedir%/sourcesink_annualsum.csv" $include "%instancedir%/sourcesink_annualsum.csv"
$offDelim
$if not exist "%instancedir%/sourcesink_annualsum.csv" $log "No annual sums for sources and sinks included"
;
parameter sourcesink_annualSum(nodesModel,years,sourcesink_techs,commodity,pc_sourcesink_annualSum);
$batinclude %aggregateNodes% sourcesink_annualSum(nodesModel,years,sourcesink_techs,commodity,pc_sourcesink_annualSum) sourcesink_annualSumIn(nodesData,years,sourcesink_techs,commodity,pc_sourcesink_annualSum) sum

$offEmpty
$onListing

* // ## Aggregation of profiles
* // Profiles are rescaled based on the annual sum.

set sourcesink_usesProfileScaling(nodesData,years,sourcesink_techs,commodity,pc_sourcesink_profile);
sourcesink_usesProfileScaling(nodesData,years,sourcesink_techs,commodity,"fixed")
    $sourcesink_configIn(nodesData,years,sourcesink_techs,commodity,"scaleFixProfileWithFixSum") = yes;
sourcesink_usesProfileScaling(nodesData,years,sourcesink_techs,commodity,"lower")
    $sourcesink_configIn(nodesData,years,sourcesink_techs,commodity,"scaleLowProfileWithLowSum") = yes;
sourcesink_usesProfileScaling(nodesData,years,sourcesink_techs,commodity,"upper")
    $sourcesink_configIn(nodesData,years,sourcesink_techs,commodity,"scaleUpProfileWithUpSum") = yes;

parameter sourcesink_ProfileSum(nodesData,years,sourcesink_techs,commodity,pc_sourcesink_profile);
sourcesink_ProfileSum(nodesData,years,sourcesink_techs,commodity,pc_sourcesink_profile)
    $sourcesink_usesProfileScaling(nodesData,years,sourcesink_techs,commodity,pc_sourcesink_profile)
    = sum(timeData, sourcesink_profileIn(timeData,nodesData,years,sourcesink_techs,commodity,pc_sourcesink_profile));

parameter sourcesink_ProfileAbsSum(nodesData,years,sourcesink_techs,commodity,pc_sourcesink_profile);
sourcesink_ProfileAbsSum(nodesData,years,sourcesink_techs,commodity,pc_sourcesink_profile)
    $sourcesink_usesProfileScaling(nodesData,years,sourcesink_techs,commodity,pc_sourcesink_profile)
    = sum(timeData, abs(sourcesink_profileIn(timeData,nodesData,years,sourcesink_techs,commodity,pc_sourcesink_profile)));

set sourcesink_ProfileScaleError(nodesData,years,sourcesink_techs,commodity,pc_sourcesink_profile);
sourcesink_ProfileScaleError(nodesData,years,sourcesink_techs,commodity,pc_sourcesink_profile)
    $(sourcesink_usesProfileScaling(nodesData,years,sourcesink_techs,commodity,pc_sourcesink_profile)
        and sourcesink_ProfileSum(nodesData,years,sourcesink_techs,commodity,pc_sourcesink_profile) = 0
        and sourcesink_ProfileAbsSum(nodesData,years,sourcesink_techs,commodity,pc_sourcesink_profile) > 0)
    = yes;

abort$sum((nodesData,years,sourcesink_techs,commodity,pc_sourcesink_profile),
            sourcesink_ProfileScaleError(nodesData,years,sourcesink_techs,commodity,pc_sourcesink_profile)) "One or more profiles cannot be rescaled as their annual sum equals zero!"

* calculate scaling factor for each type
parameter sourcesink_scalingFactor(nodesData,years,sourcesink_techs,commodity,pc_sourcesink_profile);
sourcesink_scalingFactor(nodesData,years,sourcesink_techs,commodity,pc_sourcesink_profile)
    $(sourcesink_ProfileSum(nodesData,years,sourcesink_techs,commodity,pc_sourcesink_profile)
        and sourcesink_usesProfileScaling(nodesData,years,sourcesink_techs,commodity,pc_sourcesink_profile))
    = sum(pc_sourcesink_annualSum$sameas(pc_sourcesink_profile,pc_sourcesink_annualSum),
        sourcesink_annualSumIn(nodesData,years,sourcesink_techs,commodity,pc_sourcesink_annualSum)
            / sourcesink_ProfileSum(nodesData,years,sourcesink_techs,commodity,pc_sourcesink_profile));

* rescale profiles
sourcesink_profileIn(timeData,nodesData,years,sourcesink_techs,commodity,pc_sourcesink_profile)
    $sourcesink_usesProfileScaling(nodesData,years,sourcesink_techs,commodity,pc_sourcesink_profile)
    = sourcesink_profileIn(timeData,nodesData,years,sourcesink_techs,commodity,pc_sourcesink_profile)
        * sourcesink_scalingFactor(nodesData,years,sourcesink_techs,commodity,pc_sourcesink_profile)

* // ## Aggregation of profiles
* // Profiles are aggregated based on the mapping from data nodes to model nodes.

* ==== aggregation of profiles ====
set sourcesink_usesProfileIn(nodesData,years,sourcesink_techs,commodity,pc_sourcesink_profile);
option sourcesink_usesProfileIn < sourcesink_profileIn

set sourcesink_usesProfile(nodesModel,years,sourcesink_techs,commodity,pc_sourcesink_profile);
sourcesink_usesProfile(nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity,pc_sourcesink_profile)
    $sum(nodesData$(aggregateNodesModel(nodesData,nodesModelToCalc)
        and sourcesink_usesProfileIn(nodesData,yearsToCalc,sourcesink_techs,commodity,pc_sourcesink_profile)), 1)
    = yes;

* aagregate time dimension
parameter sourcesink_profileIn_aggTime(timeModel,nodesData,yearsToCalc,sourcesink_techs,commodity,pc_sourcesink_profile);
sourcesink_profileIn_aggTime(timeModelToCalc,nodesData,yearsToCalc,sourcesink_techs,commodity,pc_sourcesink_profile)
  $sourcesink_usesProfileIn(nodesData,yearsToCalc,sourcesink_techs,commodity,pc_sourcesink_profile)
  = sum(timeData$timeMapping(timeData,timeModelToCalc),
          sourcesink_profileIn(timeData,nodesData,yearsToCalc,sourcesink_techs,commodity,pc_sourcesink_profile)
          / timeLength(timeModelToCalc));
option clear = sourcesink_profileIn;

* sum up absolute profiles
parameter sourcesink_profile(timeModel,nodesModel,years,sourcesink_techs,commodity,pc_sourcesink_profile);
sourcesink_profile(timeModelToCalc,nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity,pc_sourcesink_profile)
    $sourcesink_usesProfile(nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity,pc_sourcesink_profile)
    = sum(nodesData$aggregateNodesModel(nodesData,nodesModelToCalc),
            sourcesink_profileIn_aggTime(timeModelToCalc,nodesData,yearsToCalc,sourcesink_techs,commodity,pc_sourcesink_profile));
option clear = sourcesink_profileIn_aggTime;

set sourcesink_enabled(nodesModel,years,sourcesink_techs,commodity);
option sourcesink_enabled < sourcesink_config;

$ifthene.roundcoefs %roundcoefs%==1
sourcesink_profile(timeModelToCalc,nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity,pc_sourcesink_profile)
  = round(sourcesink_profile(timeModelToCalc,nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity,pc_sourcesink_profile), 3);

sourcesink_profile(timeModelToCalc,nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity,pc_sourcesink_profile)
  $(sourcesink_profile(timeModelToCalc,nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity,pc_sourcesink_profile) < 1e-3
        and sourcesink_profile(timeModelToCalc,nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity,pc_sourcesink_profile) > -1e-3)
  = 0;
$endif.roundcoefs

* ==== declaration of variables ====
* // ## Bounding of variables
* // Source-sink variables with either a lower, fixed, or upper profile are bounded to their respective profiles given by the input parameters.

sourcesink_flow.lo(timeModelToCalc,nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity)
    $( sourcesink_enabled(nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity)
        and sourcesink_config(nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity,"usesLowerProfile") = 1
        and sourcesink_profile(timeModelToCalc,nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity,"lower") > -inf )
    = sourcesink_profile(timeModelToCalc,nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity,"lower");

sourcesink_flow.up(timeModelToCalc,nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity)
    $( sourcesink_enabled(nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity)
        and sourcesink_config(nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity,"usesUpperProfile") = 1
        and sourcesink_profile(timeModelToCalc,nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity,"upper") < inf )
    = sourcesink_profile(timeModelToCalc,nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity,"upper");

sourcesink_flow.fx(timeModelToCalc,nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity)
    $( sourcesink_enabled(nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity)
        and sourcesink_config(nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity,"usesFixedProfile") = 1 )
    = sourcesink_profile(timeModelToCalc,nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity,"fixed");

option clear = sourcesink_profile;


* ==== equation definition ====
* // ## Equations
equations
Eq_sourcesink_useFixedSum(nodesModel,years,sourcesink_techs,commodity
    ) "Fixes the annual sum from sources / to sinks"
Eq_sourcesink_useLowerSum(nodesModel,years,sourcesink_techs,commodity
    ) "Limits the lower sum from sources / to sinks"
Eq_sourcesink_useUpperSum(nodesModel,years,sourcesink_techs,commodity
    ) "Limits the upper sum from sources / to sinks"
  ;

* // ### Fixed annual sums for sources and sinks
* // Ensures that the sources and sinks with annual sums given as fixed limits are balanced.
* // {Eq_sourcesink_useFixedSum}
Eq_sourcesink_useFixedSum(nodesModelSel,yearsSel,sourcesink_techs,commodity)
    $( sourcesink_enabled(nodesModelSel,yearsSel,sourcesink_techs,commodity)
        and sourcesink_config(nodesModelSel,yearsSel,sourcesink_techs,commodity,"usesFixedSum") = 1 )
    ..
    sum(timeModelSel,
        sourcesink_flow(timeModelSel,nodesModelSel,yearsSel,sourcesink_techs,commodity)
        * timeLength(timeModelSel))
    =e=
    sourcesink_annualSum(nodesModelSel,yearsSel,sourcesink_techs,commodity,"fixed")
    * timefrac
    ;

* // ### Lower annual sums for sources and sinks
* // Ensures that the sources and sinks with annual sums given as lower limits are balanced.
* // {Eq_sourcesink_useLowerSum}
Eq_sourcesink_useLowerSum(nodesModelSel,yearsSel,sourcesink_techs,commodity)
    $( sourcesink_enabled(nodesModelSel,yearsSel,sourcesink_techs,commodity)
        and sourcesink_config(nodesModelSel,yearsSel,sourcesink_techs,commodity,"usesLowerSum") = 1
        and sourcesink_annualSum(nodesModelSel,yearsSel,sourcesink_techs,commodity,"lower") > -inf )
    ..
    sum(timeModelSel,
        sourcesink_flow(timeModelSel,nodesModelSel,yearsSel,sourcesink_techs,commodity)
        * timeLength(timeModelSel))
    =g=
    sourcesink_annualSum(nodesModelSel,yearsSel,sourcesink_techs,commodity,"lower")
    * timefrac
    ;

* // ### Upper annual sums for sources and sinks
* // Ensures that the sources and sinks with annual sums given as upper limits are balanced.
* // {Eq_sourcesink_useUpperSum}
Eq_sourcesink_useUpperSum(nodesModelSel,yearsSel,sourcesink_techs,commodity)
    $( sourcesink_enabled(nodesModelSel,yearsSel,sourcesink_techs,commodity)
        and sourcesink_config(nodesModelSel,yearsSel,sourcesink_techs,commodity,"usesUpperSum") = 1
        and sourcesink_annualSum(nodesModelSel,yearsSel,sourcesink_techs,commodity,"upper") < inf )
    ..
    sum(timeModelSel,
        sourcesink_flow(timeModelSel,nodesModelSel,yearsSel,sourcesink_techs,commodity)
        * timeLength(timeModelSel))
    =l=
    sourcesink_annualSum(nodesModelSel,yearsSel,sourcesink_techs,commodity,"upper")
    * timefrac
    ;


* ==== model definition ====

Model M_sourcesink
/
  Eq_sourcesink_useFixedSum
  Eq_sourcesink_useLowerSum
  Eq_sourcesink_useUpperSum
/;
