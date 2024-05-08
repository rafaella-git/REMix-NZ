* SPDX-FileCopyrightText: Copyright (c) 2023-2024 German Aerospace Center (DLR)
* SPDX-License-Identifier: BSD-3-Clause

* // # sets
* // Sets are the indices of variables and parameters.

* // ## Reference
* // {special_table_sets}

* // ## Input Files
* // {special_table_set_input_files}
$onEmpty
$offListing

* ==== model time steps ====

** // SET: timeData | Time data | OEO_00140043:time stamp | set_timedata.csv
set timeData   "Time data describes the time steps of the input data." / t0001*t8760 /;

$ifthen.timemap not exist "%instancedir%/map_aggregatetimemodel.csv"

$if not set timeres        $setglobal timeres              1
$eval timeModelSteps ceil(card(timeData)/%timeres%)

** // SET: timeModel | Time model | OEO_00140043:time stamp | set_timemodel.csv
set timeModel  "Time model describes the time steps of the model run." / tm1*tm%timeModelSteps% /;

set timeHelper "Aggregation helper of time mapping operation"          / th1*th%timeres% /;

$ife %DEBUG%=1             $setglobal timestart            8
$ife %DEBUG%=1             $setglobal timeend             16
$ife %DEBUG%=2             $setglobal timestart            1
$ife %DEBUG%=2             $setglobal timeend             48

$if not set timestart      $setglobal timestart            1
$if not set timeend        $setglobal timeend              %timeModelSteps%

abort$(%timeend% > card(timeModel) or %timeend% < 1)
      "Specified end time must be within the time range 1 - %timeModelSteps%."
abort$(%timestart% > card(timeModel) or %timestart% < 1)
      "Specified start time must be within the time range 1 - %timeModelSteps%."
abort$(%timestart% > %timeend%)
      "Specified start time must be smaller or equal to end time."

set timeModelToCalc(timeModel) "Time steps to be calculated"  / tm%timestart%*tm%timeend% /;
$onVerbatim

set timeMappingHelper(timeData, timeModel, timeHelper) / #timeData : (#timeModel . #timeHelper) /;

set timeMapping(timeData, timeModel);
option timeMapping < timeMappingHelper;
$offVerbatim

$else.timemap
set timeModel(*) /
$onDelim
$if exist "%instancedir%/set_timemodel.csv" $include "%instancedir%/set_timemodel.csv"
$offDelim
$if not exist "%instancedir%/set_timemodel.csv" $log "No set elements for timeModel included."
/;

set timeModelToCalc(timeModel) /
$onDelim
$if exist "%instancedir%/set_timemodel.csv" $include "%instancedir%/set_timemodel.csv"
$offDelim
/;

** // MAP: timeData | timeModel | map_aggregatetimemodel.csv
* // ### map_aggregatetimemodel
* // Title: Map Aggregate Time Model
* // Description: Aggregate time model describes the mapping between model time and data time.
set timeMapping(timeData, timeModel) /
$onDelim
$if exist "%instancedir%/map_aggregatetimemodel.csv" $include "%instancedir%/map_aggregatetimemodel.csv"
$offDelim
$if not exist "%instancedir%/map_aggregatetimemodel.csv" $log "No set elements for timeData to timeModel mapping included, automatic mapping."
/;
$offDelim

$endif.timemap
alias(timeModel, timeModel_a);

parameter timeLength(timeModel);
timeLength(timeModel) = sum(timeData$timeMapping(timeData,timeModel), 1);

set timeModelSel(timeModel) "Selected timesteps"


$if not set aggregatelinks   $setglobal aggregatelinks    1


* ==== model regions ====
* Developer note: Nodes model comprises of OEO_00020035 (considered region i.e. Europe) OEO_00020032 (study region i.e. Germany)  OEO_00020036 (interacting region i.e. France) and OEO_00020034:study subregion (study subregion i.e. Baden-Württemberg)
** // SET: nodesModel | Nodes model | BFO_0000006:spatial region | set_nodesmodel.csv
set nodesModel(*) "Nodes model describes the model regions." /
$onDelim
$if exist "%instancedir%/set_nodesmodel.csv" $include "%instancedir%/set_nodesmodel.csv"
$offDelim
$if not exist "%instancedir%/set_nodesmodel.csv" $log "No set elements for nodesModel included - you need to specify regions in order to run the model!"
/;
alias (nodesModel,nodesModel_a,nodesModel_start,nodesModel_end)
* Developer note: Nodes data can consist of OEO_00020036 and OEO_00020034:study subregion  almost never from OEO_00020032 and OEO_00020035 unless they are the same, so one model region models.
** // SET: nodesData | Nodes data | BFO_0000006:spatial region | set_nodesdata.csv
set nodesData(*) "Nodes data describes the model data nodes." /
$onDelim
$if exist "%instancedir%/set_nodesdata.csv" $include "%instancedir%/set_nodesdata.csv"
$offDelim
$if not exist "%instancedir%/set_nodesdata.csv" $log "No set elements for nodesData included - you need to specify regions in order to run the model!"
/;

** // SET: nodesModelToCalc | Nodes Model Selected | OEO_00020034:study subregion | set_nodesmodelsel.csv
set nodesModelToCalc(nodesModel) "Nodes model to calculate describes the model data nodes to be accounted for." /
$onDelim
$if exist "%instancedir%/set_nodesmodelsel.csv" $include "%instancedir%/set_nodesmodelsel.csv"
$offDelim
$if not exist "%instancedir%/set_nodesmodelsel.csv" $log "No set elements for nodesModelSel included, using all elements in nodesModel."
$offDelim
$if not exist "%instancedir%/set_nodesmodelsel.csv" #nodesModel
/;
alias (nodesModelToCalc,nodesModelToCalc_a,nodesModelToCalc_start,nodesModelToCalc_end)

set nodesModelSel(nodesModel) "Selected model nodes."
alias (nodesModelSel,nodesModelSel_a)

** // MAP: nodesData | nodesModel | map_aggregatenodesmodel.csv
* // ### map_aggregatenodesmodel
* // Title: Map Aggregate Nodes Model
* // Description: Aggregate nodes model describes the mapping between model nodes and data nodes.
set aggregateNodesModel(nodesData,nodesModel)  /
$onDelim
$if exist "%instancedir%/map_aggregatenodesmodel.csv" $include "%instancedir%/map_aggregatenodesmodel.csv"
$offDelim
$if not exist "%instancedir%/map_aggregatenodesmodel.csv" $log "No set elements for nodesData to nodesModel mapping included, automatic mapping based on same names."
/;
$offDelim
$if not exist "%instancedir%/map_aggregatenodesmodel.csv" aggregateNodesModel(nodesData,nodesModel)$(sameas(nodesData,nodesModel)) = yes;
set map_nodesModel(nodesModel,nodesData) "Map between link and data nodes";
option map_nodesModel < aggregateNodesModel;

* validate all nodesModelSel are a subset of aggregateNodesModel
set validate_nodesModelToCalc(nodesModel);
validate_nodesModelToCalc(nodesModel)
    $nodesModelToCalc(nodesModel) = yes;
validate_nodesModelToCalc(nodesModel)
    $sum(nodesData$aggregateNodesModel(nodesData,nodesModel), 1) = no;

display$sum(nodesModel$validate_nodesModelToCalc(nodesModel), 1) "Not all nodesModel to be calculated have been mapped from nodesData. Check the map_aggregatenodesmodel for the following entries:"
display$sum(nodesModel$validate_nodesModelToCalc(nodesModel), 1) validate_nodesModelToCalc;
abort$sum(nodesModel$validate_nodesModelToCalc(nodesModel), 1) "Error encountered during mapping of modelNodes. See the lst file for a detailed description."

* ==== model years ====
* enforce chronological ordering of years by adding dummy years to UEL
set dummy_year_ordering          / 1900*2200 /;

* Year dimension for technology vintage classes
** // SAME: vintage | years
set vintage(*) "Vintage describes the technologies vintage classes." /
$onDelim
$if exist "%instancedir%/set_years.csv" $include "%instancedir%/set_years.csv"
$offDelim
$if not exist "%instancedir%/set_years.csv" $log "No set elements for vintage included - you need to specify years in order to run the model!"
/;
alias (vintage,vintage_a);
* Developer note: years contains both OEO_00020097:scenario year (Scenario Year) and OEO_00020098 (Scenario Horizon) but the data input comprises only of OEO_00020097:scenario year
* Developer note: There is currently no concept for vintage years in the ontology, how do we represent the year of operation start?
* Year dimension for the optimization
** // SET: years | Years | OEO_00020097:scenario year | set_years.csv
set years(*) "Years describes the modelling years." /
$onDelim
$if exist "%instancedir%/set_years.csv" $include "%instancedir%/set_years.csv"
$offDelim
$if not exist "%instancedir%/set_years.csv" $log "No set elements for years included - you need to specify years in order to run the model!"
/;
* Developer note: Missing concept for study year.
** // SET: yearsToCalc | Years to calculate | OEO_00020097:scenario year | set_yearssel.csv
set yearsToCalc(years) "Years to calculate describes the subset of years the model calculates." /
$onDelim
$if exist "%instancedir%/set_yearssel.csv" $include "%instancedir%/set_yearssel.csv"
$offDelim
$if not exist "%instancedir%/set_yearssel.csv" $log "No set elements for yearsSel included - you need to specify years in order to run the model!"
/;
alias (years,years_a,years_aa);
** // SAME: yearsCom | years
alias(years, yearsCom);
alias (yearsToCalc,yearsToCalc_a,yearsToCalc_aa);
set yearsToFix(years);
set yearsSel(years) "Years sel is used to determine years to calc" ;

parameter yearsLen(yearsToCalc);
yearsLen(yearsToCalc)
    = smin(yearsToCalc_a$(yearsToCalc_a.val > yearsToCalc.val), yearsToCalc_a.val - yearsToCalc.val);

* ==== activities ====
* Activities are processes, both are very wide concepts.
** // SET: activity | Activities | BFO_0000015:process | set_activities.csv
set activity(*) "Activity describes the storage and converter activities." /
$onDelim
$if exist "%instancedir%/set_activities.csv" $include "%instancedir%/set_activities.csv"
$offDelim
$if not exist "%instancedir%/set_activities.csv" $log "No set elements for activities included."
/;

* ==== grid segments ====
* Developer notes: AC network segments which are frequency coupled. These are groups of links not single links. Not final, very tricky!!
** // SET: gridSegments | Grid segments | BFO_0000026:one-dimensional spatial region | set_gridsegments.csv
set gridSegments(*) "Grid segments" /
$onDelim
$if exist "%instancedir%/set_gridsegments.csv" $include "%instancedir%/set_gridsegments.csv"
$offDelim
$if not exist "%instancedir%/set_gridsegments.csv" $log "No set elements for grid segments included."
/;


* ==== transfer links ====
** // SET: linksData | Links data | BFO_0000026:one-dimensional spatial region | set_linksdata.csv
set linksData(*) "Links data describes the transfer links in the data." /
$onDelim
$if exist "%instancedir%/set_linksdata.csv" $include "%instancedir%/set_linksdata.csv"
$offDelim
$if not exist "%instancedir%/set_linksdata.csv" $log "No set elements for linksData included."
/;

alias (linksData, linksData_a)

$eval maxlinks card(linksData)
** // SET: linksModel | Links Model| OEO_00000255:grid component link | set_linksmodel.csv
set linksModel "Links model describes the transfer links of the model. " / link1*link%maxlinks% /;

$eval maxcycle card(linksModel)
set cycles / c1*c%maxcycle% /;


* ==== link types ====
* Developer note: Link types as portions of matter is very rough. There is no concept for "terrain" in the ontology. Proably we will need to add something in that sense.
** // SET: link_types | Link types | OEO_00000331:portion of matter | set_link_types.csv
set link_types(*) "Link types" /
$onDelim
$if exist "%instancedir%/set_link_types.csv" $include "%instancedir%/set_link_types.csv"
$offDelim
$if not exist "%instancedir%/set_link_types.csv" $log "No set elements for link types included."
/;


* ==== commodities ====
* Developer note: Commodities in REMix also include Energy, which is not the case for the OEO. We probably need to assign them independently at set declaration.
** // SET: commodity | Commodities | OEO_00020067:commodity | set_commodities.csv
set commodity(*) "Commodity describes the energy carriers (e.g. electricity, heat, fuel)." /
$onDelim
$if exist "%instancedir%/set_commodities.csv" $include "%instancedir%/set_commodities.csv"
$offDelim
$if not exist "%instancedir%/set_commodities.csv" $log "No set elements for commodities included."
/;
alias (commodity,commodity_a)

* ==== indicators ====
* Developer note: Indicators are rather wildcards, the OEO does not work well with abstract objects. The most general assignment is quantity value.
** // SET: indicator | Indicators | OEO_00000350:quantity value | set_indicators.csv
set indicator(*) "Indicator describes techno-economic parameters for units, links and activities, which are balanced across temporal and spatial scales." /
$onDelim
$if exist "%instancedir%/set_indicators.csv" $include "%instancedir%/set_indicators.csv"
$offDelim
$if not exist "%instancedir%/set_indicators.csv" $log "No set elements for indicators included."
/;
alias (indicator,indicator_a,indicator_aa)

** // SET: techs | All technologies | OEO_00020102:energy transformation unit | set_techs.csv
* "Scenario indexes help to differentiate scenarios." /
* ==== technologies ====
$onMulti
set techs(*) "Techs describes the converter, sources and sinks, storage and transfer technologies." /
$if exist "%instancedir%/set_converter_techs.csv" $include "%instancedir%/set_converter_techs.csv"
/
set techs(*) /
$if exist "%instancedir%/set_storage_techs.csv" $include "%instancedir%/set_storage_techs.csv"
/
set techs(*) /
$if exist "%instancedir%/set_transfer_techs.csv" $include "%instancedir%/set_transfer_techs.csv"
/
set techs(*) /
$if exist "%instancedir%/set_sourcesink_techs.csv" $include "%instancedir%/set_sourcesink_techs.csv"
/;
$offMulti
** // SET: converter_techs | Converter technologies | OEO_00020102:energy transformation unit | set_converter_techs.csv
set converter_techs(techs) "Converter techs describes the converter technologies in the model." /
$if exist "%instancedir%/set_converter_techs.csv" $include "%instancedir%/set_converter_techs.csv"
$if not exist "%instancedir%/set_converter_techs.csv" $log "No set elements for converter technologies included."
/;
** // SET: storage_techs | Storage technologies | OEO_00000159:energy storage object | set_storage_techs.csv
set storage_techs(techs) "Storage techs describes the storage technologies in the model." /
$if exist "%instancedir%/set_storage_techs.csv" $include "%instancedir%/set_storage_techs.csv"
$if not exist "%instancedir%/set_storage_techs.csv" $log "No set elements for storage technologies included."
/;
* Developer Note: Transfer techs are grid components, there is no class that contains both links and gas pipelines.
** // SET: transfer_techs | Transfer technologies | OEO_00020006:grid component | set_transfer_techs.csv
set transfer_techs(techs) "Transfer techs describes the transfer technologies in the model." /
$if exist "%instancedir%/set_transfer_techs.csv" $include "%instancedir%/set_transfer_techs.csv"
$if not exist "%instancedir%/set_transfer_techs.csv" $log "No set elements for transfer technologies included."
/;
** // SET: sourcesink_techs | Source and Sink technologies | OEO_00020102:energy transformation unit | set_sourcesink_techs.csv
set sourcesink_techs(techs) "Sources and sinks describes the source and sink technologies in the model." /
$if exist "%instancedir%/set_sourcesink_techs.csv" $include "%instancedir%/set_sourcesink_techs.csv"
$if not exist "%instancedir%/set_sourcesink_techs.csv" $log "No set elements for source and sink technologies included."
/;

* Developer note: Degradation states and vintaging. We are missing concepts to represent age of technologies in the system.
* ==== degradation states (storage) ====
** // SET: degradation_states | Degradation states | MISSING_TERM:Degradation state | set_degradation_states.csv
set degradation_states "storage degradation unit states" / new "Default state. Nominal operation at specification capacity without any degradation." /
$onmulti
Set degradation_states "user definded unit states" /
$ondelim
$if exist "%instancedir%/set_degradation_states.csv" $include "%instancedir%/set_degradation_states.csv"
$offdelim
$if not exist "%instancedir%/set_degradation_states.csv" $log "No set elements for storage degradation states included."
/;
$offmulti

* ==== states of charge (storage) ====
** // SET: soc_statesIn | States of charge | MISSING_TERM:state of charge | set_soc.csv
set soc_statesIn(*) "States of charge."
/
$onDelim
$if exist "%instancedir%/set_soc.csv" $include "%instancedir%/set_soc.csv"
$offDelim
$if not exist "%instancedir%/set_soc.csv" $log "No set elements for SoC included."
/;

* ==== transfer links aggregation ====
** // INPUT: transfer_linkStartEnd | IAO:0000100:data set
* // ### transfer_linkStartEnd
* // Title: Transfer Link Start to End
* // Description: Definition of the start and end model nodes of each link.
* // {table_transfer_startEnd}
set pc_transfer_linkStartEnd
    /
    start               "Transfer link start | Region where the transfer connection starts | | boolean | {none} | OEO_00000339:program parameter"
    end                 "Transfer link end | Region where the transfer connection ends | | boolean | {none} | OEO_00000339:program parameter"
    /;
table transfer_linkStartEndLoad(linksData,nodesData,pc_transfer_linkStartEnd)
$onDelim
$if exist "%instancedir%/transfer_linkstartend.csv" $include "%instancedir%/transfer_linkstartend.csv"
$offDelim
$if not exist "%instancedir%/transfer_linkstartend.csv" $log "No start-end mapping for transfer links included - this deactivates all links!"
;
parameter transfer_linkStartEndIn(nodesData,linksData,pc_transfer_linkStartEnd);
option transfer_linkStartEndIn < transfer_linkStartEndLoad;

set check_linkStartEnd(linksData);
check_linkStartEnd(linksData)
    $( sum(nodesData$transfer_linkStartEndIn(nodesData,linksData,"start"), 1) <> 1
        or sum(nodesData$transfer_linkStartEndIn(nodesData,linksData,"end"), 1) <> 1 )
    = yes;

abort$(sum(check_linkStartEnd(linksData), 1) > 0) "Transfer: One or more links do not have exactly one starting and one ending node"

parameter transfer_linkStartEnd(nodesModel,linksData,pc_transfer_linkStartEnd);
$batinclude %aggregateNodes% transfer_linkStartEnd(nodesModel,linksData,pc_transfer_linkStartEnd) transfer_linkStartEndIn(nodesData,linksData,pc_transfer_linkStartEnd) sum


parameter transfer_incidenceData(nodesModel,linksData);
transfer_incidenceData(nodesModel,linksData)
    $( nodesModelToCalc(nodesModel)
        and (transfer_linkStartEnd(nodesModel,linksData,"start")
                xor transfer_linkStartEnd(nodesModel,linksData,"end"))
        and transfer_linkStartEnd(nodesModel,linksData,"end"))
    = 1;

transfer_incidenceData(nodesModel,linksData)
    $( nodesModelToCalc(nodesModel)
        and (transfer_linkStartEnd(nodesModel,linksData,"start")
                xor transfer_linkStartEnd(nodesModel,linksData,"end"))
        and transfer_linkStartEnd(nodesModel,linksData,"start"))
    = -1;

* 1: do not flip links, -1: flip links based on nodesmodel ordering
parameter transfer_incidenceData_flip(linksData);
transfer_incidenceData_flip(linksData) = 1;
transfer_incidenceData_flip(linksData)
    $sum((nodesModel,nodesModel_a)
            $(ord(nodesModel) > ord(nodesModel_a)
                and transfer_incidenceData(nodesModel,linksData) < 0
                and transfer_incidenceData(nodesModel_a,linksData) > 0), 1)
    = -1;

$onVerbatim
$ifthene.aggregateLinks %aggregatelinks%=1
$offVerbatim

parameter transfer_incidenceData_t(linksData,nodesModel);
option transfer_incidenceData_t < transfer_incidenceData;

* find all links connected to each model node
set map_linksDataToNodes(linksData,nodesModel,nodesModel);
map_linksDataToNodes(linksData,nodesModel,nodesModel_a)
    $(transfer_incidenceData_t(linksData,nodesModel) > 0
        and transfer_incidenceData_t(linksData,nodesModel_a) < 0 
        and nodesModelToCalc(nodesModel)
        and nodesModelToCalc(nodesModel_a))
    = yes;

* first calculate the min ordinate for each node node combination then assign the corresponding link model
parameter linkOrdNodesNodes(nodesModel,nodesModel);
linkOrdNodesNodes(nodesModel,nodesModel_a)
    $sum(linksData$(map_linksDataToNodes(linksData,nodesModel,nodesModel_a)
                        or map_linksDataToNodes(linksData,nodesModel_a,nodesModel)), 1)
    = smin(linksData$(map_linksDataToNodes(linksData,nodesModel,nodesModel_a)
                        or map_linksDataToNodes(linksData,nodesModel_a,nodesModel)), ord(linksData));

set map_linksModelToNodes(linksModel,nodesModel,nodesModel);
map_linksModelToNodes(linksModel,nodesModel,nodesModel_a)
    $(linkOrdNodesNodes(nodesModel,nodesModel_a)
        and (ord(linksModel) = linkOrdNodesNodes(nodesModel,nodesModel_a)))
    = yes;


set links_aggregateTemp(linksModel,linksData,nodesModel,nodesModel);
links_aggregateTemp(linksModel,linksData,nodesModel,nodesModel_a)
    $(map_linksDataToNodes(linksData,nodesModel,nodesModel_a)
        and map_linksModelToNodes(linksModel,nodesModel,nodesModel_a))
    = yes;

set links_aggregate(linksModel,linksData)  "Associates identifiers with ther respective data links" ;
option links_aggregate < links_aggregateTemp;

$onVerbatim
$else.aggregateLinks
$offVerbatim
set links_aggregate(linksModel,linksData)  "Associates identifiers with ther respective data links" ;
links_aggregate(linksModel,linksData)
    $(ord(linksModel) = ord(linksData))
    = yes;

$onVerbatim
$endif.aggregateLinks
$offVerbatim
$onOrder
alias (map_linksModel, links_aggregate);

set aggregateLinksModel(linksData,linksModel);
option aggregateLinksModel < links_aggregate;

parameter transfer_incidenceModel(nodesModel,linksModel);
transfer_incidenceModel(nodesModelToCalc,linksModel)
    = sum(linksData$links_aggregate(linksModel,linksData),
            transfer_incidenceData_flip(linksData)
            * transfer_incidenceData(nodesModelToCalc,linksData));

transfer_incidenceModel(nodesModelToCalc,linksModel)
    $(transfer_incidenceModel(nodesModelToCalc,linksModel) < 0)
    = -1;

transfer_incidenceModel(nodesModelToCalc,linksModel)
    $(transfer_incidenceModel(nodesModelToCalc,linksModel) > 0)
    = 1;

set linksModelToCalc(linksModel) "Links in which the model will calculate indicators." ;
option linksModelToCalc < transfer_incidenceModel;

$onMulti
* ==== sets and mappings for aggregations ====
** // SET: accNodes | Accounting nodes | OEO_00020034:study subregion | set_accnodes.csv
set accNodesLoad(*) "Accounting nodes are the model accounting regions."
/
  global
$onDelim
$if exist "%instancedir%/set_accnodes.csv" $include "%instancedir%/set_accnodes.csv"
$offDelim
$if not exist "%instancedir%/set_accnodes.csv" $log "No set elements for accounting nodes included."
/;
** // UNION: accNodesData | Accounting and data nodes | BFO_0000006:spatial region | set_accnodesdata.csv
set accNodesData "Accounting nodes data describes the nodes to be accounted for in the data."
/
  #accNodesLoad
  #nodesData
/;
** // SAME: accNodesModel | nodesModel
set accNodesModel "Accounting nodes model describes the nodes to be accounted for in the model."
/
  #accNodesLoad
  #nodesModel
/;
alias(accNodesModel,accNodesModel_a,accNodesModel_aa,accNodesModel_aggregation);

set accNodes(accNodesModel) "Accounting nodes describes the nodes to be accounted for in the model."
/
  #accNodesLoad
/;
** // SAME: accLinksData | linksData
set accLinksData
/
  global
  #linksData
/;

set accLinksModel
/
  global
  #linksModel
/;

set accLinks(accLinksModel)
/
  global
/;
** // SAME: accYears | years
set accYears "Accounting years are used to calculate indicators."
/
  horizon
  #years
/;
alias(accYears,accYears_a,accYears_aa);
set accYearsSel(accYears) "Selected years to calculate in the model run." ;

set accYearsToFix(accYears);
set map_accNodes(accNodesModel,accNodesModel_aggregation);
set map_accNodesToCalc(accNodesModel,nodesModel);
set map_accNodesPostCalc(accNodesModel,nodesModel);
set map_accLinks(accLinksModel,accLinksModel);
set map_accLinksToCalc(accLinksModel,linksModel);
set map_accLinksPostCalc(accLinksModel,linksModel);
set map_accYears(accYears,accYears);
set map_accYearsToCalc(accYears,years);
set map_accYearsPostCalc(accYears,years);

$offMulti
** // MAP: nodesData | accNodes | map_accnodes.csv
* // ### map_accNodes
* // Title: Map Accounting Nodes
* // Description: Map accounting nodes describes the mapping between model accounting nodes and data nodes.
set map_accNodesLoad(nodesData,accNodesModel)
/
$onDelim
$if exist "%instancedir%/map_accnodes.csv" $include "%instancedir%/map_accnodes.csv"
$offDelim
$if not exist "%instancedir%/map_accnodes.csv" $log "No set elements for nodesData to nodesAcc mapping included."
/;
$onListing
$offEmpty

map_accNodes(accNodesModel,accNodesModel_aggregation)
    $sum((nodesData, nodesModelToCalc)
            $( sameas(nodesModelToCalc,accNodesModel)
                and map_accNodesLoad(nodesData,accNodesModel_aggregation)
                and aggregateNodesModel(nodesData,nodesModelToCalc)), 1) = yes;

map_accNodes(accNodesModel,"global")$sum(nodesModelToCalc$sameas(accNodesModel,nodesModelToCalc), 1) = yes;
map_accNodes(accNodesModel,"global")$sum(accNodes$sameas(accNodesModel,accNodes), 1) = yes;
map_accNodes(accNodesModel,accNodesModel) = yes;
map_accLinks(accLinksModel,"global")$sum(linksModelToCalc$sameas(accLinksModel,linksModelToCalc), 1) = yes;
map_accLinks(accLinksModel,accLinksModel) = yes;
map_accYears(accYears,"horizon")$sum(yearsToCalc$sameas(accYears,yearsToCalc), 1)  = yes;
map_accYears(accYears,accYears) = yes;

map_accNodesToCalc(accNodesModel,nodesModelToCalc)
    $sum((nodesData)
            $( map_accNodesLoad(nodesData,accNodesModel)
                and aggregateNodesModel(nodesData,nodesModelToCalc)), 1) = yes;

map_accNodesToCalc("global",nodesModelToCalc) = yes;
map_accNodesToCalc(accNodesModel,nodesModelToCalc)$sameas(accNodesModel,nodesModelToCalc) = yes;
map_accLinksToCalc("global",linksModelToCalc) = yes;
map_accLinksToCalc(accLinksModel,linksModelToCalc)$sameas(accLinksModel,linksModelToCalc) = yes;
map_accYearsToCalc("horizon",yearsToCalc) = yes;
map_accYearsToCalc(accYears,yearsToCalc)$sameas(accYears,yearsToCalc) = yes;

* Duplicate sets for postcalc reporting without horizon and accYears
map_accNodesPostCalc(accNodesModel,nodesModelToCalc) = map_accNodesToCalc(accNodesModel,nodesModelToCalc);
map_accLinksPostCalc(accLinksModel,linksModelToCalc) = map_accLinksToCalc(accLinksModel,linksModelToCalc);
map_accYearsPostCalc(accYears,yearsToCalc)$sameas(accYears,yearsToCalc) = yes;
map_accYearsPostCalc(accYears,years)$(sameas(accYears,years) and years.val < sum(yearsToCalc$(ord(yearsToCalc) = 1), yearsToCalc.val)) = yes;

set map_nodesAccounting(accNodesModel_aggregation,accNodesModel) "Map accounting nodes describes the mapping between model accounting nodes and data nodes." ;
option map_nodesAccounting < map_accNodes;

* ==== Generic sets  ====
* These set names are repeated across different profile parameters.
* they are not strictly related to each other and are not input data but
* are needed for validation.
** // PROFILE: profileTypes | Profile Types | OEO_00140056:flow potential | set_profiletypes.csv

** // SET: scenario | Scenario | OEO_00000364:scenario | set_scenarios.csv
* "Scenario indexes help to differentiate scenarios." /

** // SET: capType | Capacity Types | OEO_00030019:balance process attribute | set_captypes.csv
* "Capacity types differentiate investment decisions." /

** // SET: balanceType | Balance Types | OEO_00030019:balance process attribute | set_balancetypes.csv
* "Balance types differentiate types of balances." /
