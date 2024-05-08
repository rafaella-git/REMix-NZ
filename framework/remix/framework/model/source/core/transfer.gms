* SPDX-FileCopyrightText: Copyright (c) 2023-2024 German Aerospace Center (DLR)
* SPDX-License-Identifier: BSD-3-Clause

* // # core_transfer
* // The equations in this file describe the commodity transfer in the model.

* // ## Variables
* // {special_table_transfer_variables}
$onVerbatim
$if not set opfmethod        $setglobal opfmethod         angle

$iftheni.opfmethod %opfmethod%==kirchhoff
$log Transfer: Using Kirchhoff DC optimal power flow (for AC links)
$elseifi.opfmethod %opfmethod%==angle
$log Transfer: Using voltage angle-based DC optimal power flow (for AC links)
$elseifi.opfmethod %opfmethod%==disable
$log Transfer: Using DC Transfer instead of DC optimal power flow (for AC links)
$else.opfmethod
$abort Transfer: No valid DC optimal power flow method
$endif.opfmethod
$offVerbatim

$offListing
$onEmpty


* ==== declaration of variables ====
integer variables
  transfer_linksTotal_MIP(linksModel,years,transfer_techs,vintage
    ) "Total integer number of active links in the system"

positive variables
  transfer_linksBuild(linksModel,years,transfer_techs,vintage
    ) "Number of transfer links built"
  transfer_linksDecom(linksModel,years,transfer_techs,vintage
    ) "Number of transfer links decommissioned"
  transfer_linksTotal(linksModel,years,transfer_techs,vintage
    ) "Total number of active transfer links in the system"

  transfer_flowAlong(timeModel,linksModel,years,transfer_techs,vintage
    ) "Flow along the transfer link for each time step"
  transfer_flowAgainst(timeModel,linksModel,years,transfer_techs,vintage
    ) "Flow against the transfer link for each time step (for bidirectional links)"

variables
  transfer_dcopf_voltageAngle(timeModel,nodesModel,years,gridSegments
    ) "Net import / export to the grid segment modeled as DC optimal power flow"
;

* The table headers in the parameter input sets represent the following:
* "Title | Description | Constraints | Type |  Units | Ontology"

* // ## Input Files
** // INPUT: transfer_linksParam | IAO:0000100:data set
* // ### transfer_linksParam
* // Title: Transfer Transfer Parameters
* // Description: Transfer parameters describe the building and decommissioning of transfer links.
* // {table_transfer_linksParam}
set pc_transfer_linksParam
    /
    linksBuild          "Links build | Exogenously given links to be built in a given year | minimum:0;lessThan:['linksUpperLimit'] | number | {capacity} | OEO_00010257:power capacity"
    linksLowerLimit     "Links lower limit | Lower limit on total links for all vintage classes | minimum:0;lessThan:['linksBuild','linksUpperLimit'] | number | {capacity} | OEO_00000104:constraint"
    linksUpperLimit     "Links upper limit | Upper limit on total number of links of all vintage classes | minimum:0;default:'inf' | number | {capacity} | OEO_00000104:constraint"
    linksDelta          "Links delta | Maximum allowed delta of links per year | minimum:0 | number | {capacity} | OEO_00000104:constraint"
    limitFlows          "Limit the flows in either direction | If 1, apply flowAlongLimit or flowAgainstLimit or respective profiles |  | boolean | {none} | OEO_00000339:program parameter"
    flowAlongLimit      "Flow along limit | Upper limit for the flow along a link relative to its capacity | minimum:0 | number | {none} | OEO_00000104:constraint"
    flowAgainstLimit    "Flow against limit | Upper limit for the flow against a link relative to its capacity | minimum:0 | number | {none} | OEO_00000104:constraint"
    noExpansion         "Prevent expansion | Prevent expansion beyond exogenous capacity expansion | | boolean | {none} | OEO_00000339:program parameter"
    circuits            "Circuits | Number of electrical circuits used for DC-opf reactance calculation | minimum:0 | integer | {circuits} | MISSING_TERM:circuits"
    /;
table transfer_linksParamIn(linksData,years,transfer_techs,pc_transfer_linksParam)
$onDelim
$if exist "%instancedir%/transfer_linksparam.csv" $include "%instancedir%/transfer_linksparam.csv"
$offDelim
$if not exist "%instancedir%/transfer_linksparam.csv" $log "No bounds for links included"
;

** // INPUT: transfer_flowProfile | IAO:0000100:data set
* // ### transfer_flowProfile
* // Title: Link-Specific Flow Profiles
* // Description: Profile for transfer upper link flows along and against an edge.
* // {table_transfer_flowProfile}
set pc_transfer_flowProfile
    /
    along               "Flow along | Profile for the upper limit of a flow along a link relative to its capacity |"
    against             "Flow against | Profile for the upper limit of a flow against a link relative to its capacity |"
    /;

table transfer_flowProfileLoad(linksData,years,transfer_techs,pc_transfer_flowProfile,timeData)
$onDelim
$if exist "%instancedir%/transfer_flowprofile.csv" $include "%instancedir%/transfer_flowprofile.csv"
$offDelim
$if not exist "%instancedir%/transfer_flowprofile.csv" $log "No profile for transfer flow limits included."
;
parameter transfer_flowProfileIn(timeData,linksData,years,transfer_techs,pc_transfer_flowProfile);
option transfer_flowProfileIn < transfer_flowProfileLoad;
option clear = transfer_flowProfileLoad;

** // INPUT: transfer_techParam | IAO:0000100:data set
* // ### transfer_techParam
* // Title: Transfer Technology Parameters
* // Description: Transfer tech parameters describe the operational properties of the links.
* // {table_transfer_techParam}
set pc_transfer_techParam
    /
    lifeTime            "Technical lifetime | Technical life time of the unit for calculating decommissioning | minimum:0 | integer | {span} | OEO_00000339:operational life time"
    freeDecom           "Free decommissioning | Allow decommissioning of the unit before the end of the technical life time | | boolean | {none} | OEO_00000339:program parameter"
    mipLinks            "MIP links | Model the links of the technology as integer values | | boolean | {none} | OEO_00000339:program parameter"
    flowUpperLimit      "Flow upper limit | Upper limit for the flow factor per  | minimum:0;default:1 | number | {none} | OEO_00000104:constraint"
    /;
table transfer_techParam(transfer_techs,vintage,pc_transfer_techParam)
$onDelim
$if exist "%instancedir%/transfer_techparam.csv" $include "%instancedir%/transfer_techparam.csv"
$offDelim
$if not exist "%instancedir%/transfer_techparam.csv" $log "No technology parameters for links included"
;

** // INPUT: transfer_coefficient | IAO:0000100:data set
* // ### transfer_coefficient
* // Title: Transfer Coefficients
* // Description: Transfer coefficients describe the flow properties of the links.
* // {table_transfer_coefficient}
set pc_transfer_coefficient
    /
    coefficient         "Coefficient | Coefficient for maximum commodity flow | minimum:0 | number | {rate} | MISSING_TERM:transmission coefficient"
    /;
table transfer_coefficient(transfer_techs,vintage,commodity,pc_transfer_coefficient)
$ondelim
$if exist "%instancedir%/transfer_coefficient.csv" $include "%instancedir%/transfer_coefficient.csv"
$offdelim
$if not exist "%instancedir%/transfer_coefficient.csv" $log "No coefficients for maximum commodity flow included"
;

** // INPUT: transfer_coefPerFlow | IAO:0000100:data set
* // ### transfer_coefPerFlow
* // Title: Transfer Coefficients per Flow
* // Description: Transfer coefficients per flow describe the effect of flow events on links.
* // {table_transfer_coefPerFlow}
set pc_transfer_coefPerFlow
    /
    coefPerFlow         "Coefficient per flow | Coefficients for losses and gains per flow | default:0 | number | {rate} | MISSING_TERM:transmission coefficient"
    /;
table transfer_coefPerFlow(transfer_techs,vintage,commodity,pc_transfer_coefPerFlow)
$ondelim
$if exist "%instancedir%/transfer_coefperflow.csv" $include "%instancedir%/transfer_coefperflow.csv"
$offdelim
$if not exist "%instancedir%/transfer_coefperflow.csv" $log "No coefficients for losses and gains per flow included"
;

** // INPUT: transfer_coefPerLength | IAO:0000100:data set
* // ### transfer_coefPerLength
* // Title: Transfer Coefficients per Distance
* // Description: Transfer coefficients per length describe the effect of flow events on links as function of the length.
* // {table_transfer_coefPerLength}
set pc_transfer_coefPerLength
    /
    coefPerLength     "Coefficient per length | Coefficients for losses and gains per flow times length | default:0 | number | {rate} | OEO_00030019:process attribute"
    /;
table transfer_coefPerLength(transfer_techs,vintage,commodity,link_types,pc_transfer_coefPerLength)
$ondelim
$if exist "%instancedir%/transfer_coefperlength.csv" $include "%instancedir%/transfer_coefperlength.csv"
$offdelim
$if not exist "%instancedir%/transfer_coefperlength.csv" $log "No coefficients for losses and gains per flow and length included"
;

** // INPUT: transfer_reactPerLength | IAO:0000100:data set
* // ### transfer_reactPerLength
* // Title: Transfer Reactance per Distance
* // Description: Transfer reactance per length can be given for electrical power links.
* // {table_transfer_reactPerLength}
set pc_transfer_reactPerLength
    /
    reactPerLength    "Reactance per length | Electrical reactance per length | minimum:0 | number | {rate} | OEO_00030019:process attribute"
    /;
table transfer_reactPerLength(transfer_techs,vintage,link_types,pc_transfer_reactPerLength)
$ondelim
$if exist "%instancedir%/transfer_reactperlength.csv" $include "%instancedir%/transfer_reactperlength.csv"
$offdelim
$if not exist "%instancedir%/transfer_reactperlength.csv" $log "No electrical reactance per length included"
;

** // INPUT: transfer_lengthParam | IAO:0000100:data set
* // ### transfer_lengthParam
* // Title: Transfer Length Parameters
* // Description: Transfer length parameters define the length of transmission links.
* // {table_transfer_lengthParam}
set pc_transfer_lengthParam
    /
    length            "Length | Length of a link connecting two nodes | minimum:0 | number | {length} | BFO_0000026:one-dimensional spatial region"
    /;
table transfer_lengthParamIn(linksData,link_types,pc_transfer_lengthParam)
$onDelim
$if exist "%instancedir%/transfer_lengthparam.csv" $include "%instancedir%/transfer_lengthparam.csv"
$offDelim
$if not exist "%instancedir%/transfer_lengthparam.csv" $log "No lengths for links included"
;

** // INPUT: transfer_gridSegments | IAO:0000100:data set
* // ### transfer_gridSegments
* // Title: Transfer Grid Segments
* // Description: Used to specify calculation method for certain grid segments.
* // {table_transfer_gridSegments}
set pc_transfer_gridSegments
    /
    useDCopf            "Use DC optimal power flow | Use DC optimal power flow for the grid segment with all contained link and transfer-technology combinations | | boolean | {none}| OEO_00000339:program parameter"
    /;
table transfer_gridSegmentsLoad(gridSegments,linksData,transfer_techs,pc_transfer_gridSegments)
$onDelim
$if exist "%instancedir%/transfer_gridsegments.csv" $include "%instancedir%/transfer_gridsegments.csv"
$offDelim
$if not exist "%instancedir%/transfer_gridsegments.csv" $log "No grid segments for links included"
;
parameter transfer_gridSegmentsIn(linksData,transfer_techs,gridSegments,pc_transfer_gridSegments);
option transfer_gridSegmentsIn < transfer_gridSegmentsLoad;

$offEmpty
$onListing

parameter transfer_linksParam(linksModel,years,transfer_techs,pc_transfer_linksParam);
transfer_linksParam(linksModelToCalc,years,transfer_techs,pc_transfer_linksParam)
    = sum(linksData$links_aggregate(linksModelToCalc,linksData), transfer_linksParamIn(linksData,years,transfer_techs,pc_transfer_linksParam));

* // ## Load links from gdx file
$iftheni.cfgdx not %capsfromgdx%==None
transfer_linksParam(linksModel,years,transfer_techs,"linksBuild")
  = sum(vintage,
        transfer_links(%selscen%linksModel,years,transfer_techs,vintage,"build"));
transfer_linksParam(linksModel,years,transfer_techs,"linksUpperLimit")
  = sum(vintage,
        transfer_links(%selscen%linksModel,years,transfer_techs,vintage,"total"));
transfer_linksParam(linksModel,years,transfer_techs,"linksLowerLimit")
  = transfer_linksParam(linksModel,years,transfer_techs,"linksUpperLimit");
transfer_linksParam(linksModel,years,transfer_techs,"noExpansion")
  = yes;
$endif.cfgdx


set transfer_hasflowProfileIn(linksData,years,transfer_techs,pc_transfer_flowProfile);
option transfer_hasflowProfileIn < transfer_flowProfileIn;

* === modify transfer_flowProfileIn vector to fill with default values, in case profile is not specified ===
transfer_linksParamIn(linksData,years,transfer_techs,"flowAlongLimit")
    $(not transfer_linksParamIn(linksData,years,transfer_techs,"limitFlows"))
    = 1;

transfer_linksParamIn(linksData,years,transfer_techs,"flowAgainstLimit")
    $(not transfer_linksParamIn(linksData,years,transfer_techs,"limitFlows"))
    = 1;

transfer_flowProfileIn(timeData,linksData,years,transfer_techs,"along")
    $(not transfer_hasflowProfileIn(linksData,years,transfer_techs,"along"))
    = transfer_linksParamIn(linksData,years,transfer_techs,"flowAlongLimit");

transfer_flowProfileIn(timeData,linksData,years,transfer_techs,"against")
    $(not transfer_hasflowProfileIn(linksData,years,transfer_techs,"against"))
    = transfer_linksParamIn(linksData,years,transfer_techs,"flowAgainstLimit");

* aggregate time dimension
parameter transfer_flowProfileIn_aggTime(timeModel,linksData,yearsToCalc,transfer_techs,pc_transfer_flowProfile);
transfer_flowProfileIn_aggTime(timeModelToCalc,linksData,yearsToCalc,transfer_techs,pc_transfer_flowProfile)
    $transfer_hasflowProfileIn(linksData,yearsToCalc,transfer_techs,pc_transfer_flowProfile)
    = sum(timeData$timeMapping(timeData,timeModelToCalc),
          transfer_flowProfileIn(timeData,linksData,yearsToCalc,transfer_techs,pc_transfer_flowProfile)
          / timeLength(timeModelToCalc));
option clear = transfer_flowProfileIn;

transfer_flowProfileIn_aggTime(timeModelToCalc,linksData,yearsToCalc,transfer_techs,pc_transfer_flowProfile)
    $(not transfer_hasflowProfileIn(linksData,yearsToCalc,transfer_techs,pc_transfer_flowProfile) and
      not transfer_linksParamIn(linksData,yearsToCalc,transfer_techs,"limitFlows"))
    = 1;

transfer_flowProfileIn_aggTime(timeModelToCalc,linksData,yearsToCalc,transfer_techs,"along")
    $(not transfer_hasflowProfileIn(linksData,yearsToCalc,transfer_techs,"along") and
      transfer_linksParamIn(linksData,yearsToCalc,transfer_techs,"limitFlows"))
    = transfer_linksParamIn(linksData,yearsToCalc,transfer_techs,"flowAlongLimit");

transfer_flowProfileIn_aggTime(timeModelToCalc,linksData,yearsToCalc,transfer_techs,"against")
    $(not transfer_hasflowProfileIn(linksData,yearsToCalc,transfer_techs,"against") and 
      transfer_linksParamIn(linksData,yearsToCalc,transfer_techs,"limitFlows"))
    = transfer_linksParamIn(linksData,yearsToCalc,transfer_techs,"flowAgainstLimit");

set transfer_finiteLinkLimit(linksModel,years,transfer_techs);
transfer_finiteLinkLimit(linksModelToCalc,yearsToCalc,transfer_techs)
    = sum(linksData$links_aggregate(linksModelToCalc,linksData), transfer_linksParamIn(linksData,yearsToCalc,transfer_techs,"linksUpperLimit")) > 0
        and sum(linksData$links_aggregate(linksModelToCalc,linksData), transfer_linksParamIn(linksData,yearsToCalc,transfer_techs,"linksUpperLimit")) < inf;

set transfer_infiniteLinkLimit(linksModel,years,transfer_techs);
transfer_infiniteLinkLimit(linksModelToCalc,yearsToCalc,transfer_techs)
    = sum(linksData$links_aggregate(linksModelToCalc,linksData), transfer_linksParamIn(linksData,yearsToCalc,transfer_techs,"linksUpperLimit")) = inf;
$offEmpty

parameter transfer_lengthParam(linksModel,link_types,pc_transfer_lengthParam);
transfer_lengthParam(linksModelToCalc,link_types,pc_transfer_lengthParam)
    = sum(linksData$links_aggregate(linksModelToCalc,linksData), transfer_lengthParamIn(linksData,link_types,pc_transfer_lengthParam))
        / sum(linksData$links_aggregate(linksModelToCalc,linksData), 1);

parameter transfer_gridSegments(linksModel,transfer_techs,gridSegments,pc_transfer_gridSegments);
transfer_gridSegments(linksModelToCalc,transfer_techs,gridSegments,pc_transfer_gridSegments)
    = smax(linksData$links_aggregate(linksModelToCalc,linksData), transfer_gridSegmentsIn(linksData,transfer_techs,gridSegments,pc_transfer_gridSegments));

parameter transfer_dcopf_Xtech(linksModel,years,transfer_techs,vintage,gridSegments);
transfer_dcopf_Xtech(linksModelToCalc,yearsToCalc,transfer_techs,vintage,gridSegments)
    $(transfer_gridSegments(linksModelToCalc,transfer_techs,gridSegments,"useDCopf")
        and transfer_linksParam(linksModelToCalc,yearsToCalc,transfer_techs,"circuits") > 0
        and sum(link_types, transfer_reactPerLength(transfer_techs,vintage,link_types,"reactPerLength")) > 0 )
    = (1 / sum(linksData
                $links_aggregate(linksModelToCalc,linksData),
            1 / ( sum(link_types,
                        transfer_lengthParamIn(linksData,link_types,"length")
                        * transfer_reactPerLength(transfer_techs,vintage,link_types,"reactPerLength"))
                    / transfer_linksParam(linksModelToCalc,yearsToCalc,transfer_techs,"circuits"))));

set transfer_hasflowProfile(linksModel,years,transfer_techs,pc_transfer_flowProfile);
transfer_hasflowProfile(linksModelToCalc,yearsToCalc,transfer_techs,pc_transfer_flowProfile)
    = sum(linksData$links_aggregate(linksModelToCalc,linksData),
            transfer_hasflowProfileIn(linksData,yearsToCalc,transfer_techs,pc_transfer_flowProfile));

* ==== parameter modifications ====
transfer_linksParam(linksModelToCalc,yearsToCalc,transfer_techs,"linksLowerLimit")
    $sum(vintage, transfer_techParam(transfer_techs,vintage,"mipLinks"))
    = floor(transfer_linksParam(linksModelToCalc,yearsToCalc,transfer_techs,"linksLowerLimit"));
transfer_linksParam(linksModelToCalc,yearsToCalc,transfer_techs,"linksUpperLimit")
    $sum(vintage, transfer_techParam(transfer_techs,vintage,"mipLinks"))
    = ceil(transfer_linksParam(linksModelToCalc,yearsToCalc,transfer_techs,"linksUpperLimit"));


* ==== calculation of mappings ====

* Technologies with a lifeTime > 0 are available
set transfer_availTech(linksModel,years,transfer_techs,vintage);
transfer_availTech(linksModel,years,transfer_techs,vintage)
    $(vintage.val = smax(vintage_a$(vintage_a.val <= years.val
        and transfer_techParam(transfer_techs,vintage_a,"lifeTime") > 0), vintage_a.val)) = yes;

* Technologies to optimize become unavailable if they have an linksUpperLimit of 0
transfer_availTech(linksModelToCalc,years,transfer_techs,vintage)
    $(yearsToCalc(years) and transfer_linksParam(linksModelToCalc,years,transfer_techs,"linksUpperLimit") = 0 ) = no;

* Technologies already built become unavailable if they have an linksBuild of 0
transfer_availTech(linksModelToCalc,years,transfer_techs,vintage)
    $( ( not yearsToCalc(years)) and transfer_linksParam(linksModelToCalc,years,transfer_techs,"linksBuild") = 0 ) = no;

* Used technologies are available technologies over their technical lifeTime
set transfer_usedTech(linksModel,years,transfer_techs,vintage);
transfer_usedTech(linksModelToCalc,years,transfer_techs,vintage)
    $(vintage.val <= years.val
        and years.val < smax(years_a$transfer_availTech(linksModelToCalc,years_a,transfer_techs,vintage),
                                years_a.val + transfer_techParam(transfer_techs,vintage,"lifeTime"))
        ) = yes;

* Technologies have to be decomissioned in the interval of first avail + lifetime to last avail + lifetime
set transfer_decomTech(linksModel,years,transfer_techs,vintage);
transfer_decomTech(linksModel,years,transfer_techs,vintage)
  $(sum(years_a$transfer_usedTech(linksModel,years_a,transfer_techs,vintage), 1)
    and sum(yearsToCalc
      $(sameas(years, yearsToCalc)
        and yearsToCalc.val >= smin(years_a$transfer_availTech(linksModel,years_a,transfer_techs,vintage), years_a.val) + transfer_techParam(transfer_techs,vintage,"lifeTime")
        and yearsToCalc.val <= smax(years_a$transfer_availTech(linksModel,years_a,transfer_techs,vintage), years_a.val) + transfer_techParam(transfer_techs,vintage,"lifeTime")), 1))
  = yes;

* Extend the decom frame to the year after the last year of usedTech
transfer_decomTech(linksModel,yearsToCalc,transfer_techs,vintage)
  $(transfer_usedTech(linksModel,yearsToCalc-1,transfer_techs,vintage)
    and transfer_decomTech(linksModel,yearsToCalc-1,transfer_techs,vintage))
  = yes;

parameter transfer_flowProfile(timeModel,linksModel,years,transfer_techs,vintage,pc_transfer_flowProfile);

transfer_flowProfile(timeModelToCalc,linksModelToCalc,yearsToCalc,transfer_techs,vintage,"along")
    $transfer_finiteLinkLimit(linksModelToCalc,yearsToCalc,transfer_techs)
    = (sum(linksData$(links_aggregate(linksModelToCalc,linksData)
                     and transfer_incidenceData_flip(linksData) = 1),
            transfer_flowProfileIn_aggTime(timeModelToCalc,linksData,yearsToCalc,transfer_techs,"along")
            * transfer_linksParamIn(linksData,yearsToCalc,transfer_techs,"linksUpperLimit"))
        + sum(linksData$(links_aggregate(linksModelToCalc,linksData)
                        and transfer_incidenceData_flip(linksData) = -1),
            transfer_flowProfileIn_aggTime(timeModelToCalc,linksData,yearsToCalc,transfer_techs,"against")
            * transfer_linksParamIn(linksData,yearsToCalc,transfer_techs,"linksUpperLimit")))
    / sum(linksData$links_aggregate(linksModelToCalc,linksData),
            transfer_linksParamIn(linksData,yearsToCalc,transfer_techs,"linksUpperLimit"));

transfer_flowProfile(timeModelToCalc,linksModelToCalc,yearsToCalc,transfer_techs,vintage,"against")
    $transfer_finiteLinkLimit(linksModelToCalc,yearsToCalc,transfer_techs)
    = (sum(linksData$(links_aggregate(linksModelToCalc,linksData)
                     and transfer_incidenceData_flip(linksData) = 1),
            transfer_flowProfileIn_aggTime(timeModelToCalc,linksData,yearsToCalc,transfer_techs,"against")
            * transfer_linksParamIn(linksData,yearsToCalc,transfer_techs,"linksUpperLimit"))
        + sum(linksData$(links_aggregate(linksModelToCalc,linksData)
                        and transfer_incidenceData_flip(linksData) = -1),
            transfer_flowProfileIn_aggTime(timeModelToCalc,linksData,yearsToCalc,transfer_techs,"along")
            * transfer_linksParamIn(linksData,yearsToCalc,transfer_techs,"linksUpperLimit")))
    / sum(linksData$links_aggregate(linksModelToCalc,linksData),
            transfer_linksParamIn(linksData,yearsToCalc,transfer_techs,"linksUpperLimit"));

transfer_flowProfile(timeModelToCalc,linksModelToCalc,yearsToCalc,transfer_techs,vintage,"along")
    $transfer_infiniteLinkLimit(linksModelToCalc,yearsToCalc,transfer_techs)
    = (sum(linksData$(links_aggregate(linksModelToCalc,linksData)
                     and transfer_incidenceData_flip(linksData) = 1
                     and transfer_linksParamIn(linksData,yearsToCalc,transfer_techs,"linksUpperLimit") = inf ),
            transfer_flowProfileIn_aggTime(timeModelToCalc,linksData,yearsToCalc,transfer_techs,"along"))
        + sum(linksData$(links_aggregate(linksModelToCalc,linksData)
                        and transfer_incidenceData_flip(linksData) = -1
                        and transfer_linksParamIn(linksData,yearsToCalc,transfer_techs,"linksUpperLimit") = inf ),
            transfer_flowProfileIn_aggTime(timeModelToCalc,linksData,yearsToCalc,transfer_techs,"against")))
    / sum(linksData$(links_aggregate(linksModelToCalc,linksData)
                    and transfer_linksParamIn(linksData,yearsToCalc,transfer_techs,"linksUpperLimit") = inf ),
            1);

transfer_flowProfile(timeModelToCalc,linksModelToCalc,yearsToCalc,transfer_techs,vintage,"against")
    $transfer_infiniteLinkLimit(linksModelToCalc,yearsToCalc,transfer_techs)
    = (sum(linksData$(links_aggregate(linksModelToCalc,linksData)
                     and transfer_incidenceData_flip(linksData) = 1
                     and transfer_linksParamIn(linksData,yearsToCalc,transfer_techs,"linksUpperLimit") = inf ),
            transfer_flowProfileIn_aggTime(timeModelToCalc,linksData,yearsToCalc,transfer_techs,"against"))
        + sum(linksData$(links_aggregate(linksModelToCalc,linksData)
                        and transfer_incidenceData_flip(linksData) = -1
                        and transfer_linksParamIn(linksData,yearsToCalc,transfer_techs,"linksUpperLimit") = inf ),
            transfer_flowProfileIn_aggTime(timeModelToCalc,linksData,yearsToCalc,transfer_techs,"along")))
    / sum(linksData$(links_aggregate(linksModelToCalc,linksData)
                    and transfer_linksParamIn(linksData,yearsToCalc,transfer_techs,"linksUpperLimit") = inf ),
            1);
option clear = transfer_flowProfileIn_aggTime;

* Orientation of links affects the flowAlongLimit/flowAgainstLimit

transfer_linksParam(linksModelToCalc,years,transfer_techs,"flowAlongLimit")
    $transfer_finiteLinkLimit(linksModelToCalc,years,transfer_techs)
    = (sum(linksData$(links_aggregate(linksModelToCalc,linksData)
                     and transfer_incidenceData_flip(linksData) = 1),
            transfer_linksParamIn(linksData,years,transfer_techs,"flowAlongLimit")
            * transfer_linksParamIn(linksData,years,transfer_techs,"linksUpperLimit"))
        + sum(linksData$(links_aggregate(linksModelToCalc,linksData)
                        and transfer_incidenceData_flip(linksData) = -1),
            transfer_linksParamIn(linksData,years,transfer_techs,"flowAgainstLimit")
            * transfer_linksParamIn(linksData,years,transfer_techs,"linksUpperLimit")))
    / sum(linksData$links_aggregate(linksModelToCalc,linksData),
            transfer_linksParamIn(linksData,years,transfer_techs,"linksUpperLimit"));

transfer_linksParam(linksModelToCalc,years,transfer_techs,"flowAgainstLimit")
    $transfer_finiteLinkLimit(linksModelToCalc,years,transfer_techs)
    = (sum(linksData$(links_aggregate(linksModelToCalc,linksData)
                     and transfer_incidenceData_flip(linksData) = -1),
            transfer_linksParamIn(linksData,years,transfer_techs,"flowAlongLimit")
            * transfer_linksParamIn(linksData,years,transfer_techs,"linksUpperLimit"))
        + sum(linksData$(links_aggregate(linksModelToCalc,linksData)
                        and transfer_incidenceData_flip(linksData) = 1),
            transfer_linksParamIn(linksData,years,transfer_techs,"flowAgainstLimit")
            * transfer_linksParamIn(linksData,years,transfer_techs,"linksUpperLimit")))
    / sum(linksData$links_aggregate(linksModelToCalc,linksData),
            transfer_linksParamIn(linksData,years,transfer_techs,"linksUpperLimit"));

transfer_linksParam(linksModelToCalc,years,transfer_techs,"flowAlongLimit")
    $transfer_infiniteLinkLimit(linksModelToCalc,years,transfer_techs)
    = (sum(linksData$(links_aggregate(linksModelToCalc,linksData)
                     and transfer_incidenceData_flip(linksData) = 1
                     and transfer_linksParamIn(linksData,years,transfer_techs,"linksUpperLimit") = inf),
            transfer_linksParamIn(linksData,years,transfer_techs,"flowAlongLimit"))
        + sum(linksData$(links_aggregate(linksModelToCalc,linksData)
                        and transfer_incidenceData_flip(linksData) = -1
                        and transfer_linksParamIn(linksData,years,transfer_techs,"linksUpperLimit") = inf),
            transfer_linksParamIn(linksData,years,transfer_techs,"flowAgainstLimit")))
    / sum(linksData$(links_aggregate(linksModelToCalc,linksData)
                    and transfer_linksParamIn(linksData,years,transfer_techs,"linksUpperLimit") = inf ),
            1);

transfer_linksParam(linksModelToCalc,years,transfer_techs,"flowAgainstLimit")
    $transfer_infiniteLinkLimit(linksModelToCalc,years,transfer_techs)
    = (sum(linksData$(links_aggregate(linksModelToCalc,linksData)
                     and transfer_incidenceData_flip(linksData) = -1
                     and transfer_linksParamIn(linksData,years,transfer_techs,"linksUpperLimit") = inf),
            transfer_linksParamIn(linksData,years,transfer_techs,"flowAlongLimit"))
        + sum(linksData$(links_aggregate(linksModelToCalc,linksData)
                        and transfer_incidenceData_flip(linksData) = 1
                        and transfer_linksParamIn(linksData,years,transfer_techs,"linksUpperLimit") = inf),
            transfer_linksParamIn(linksData,years,transfer_techs,"flowAgainstLimit")))
    / sum(linksData$(links_aggregate(linksModelToCalc,linksData)
                    and transfer_linksParamIn(linksData,years,transfer_techs,"linksUpperLimit") = inf ),
            1);

* Mapping for grid segments using DC optimal power flow
set gridSegments_dcopf(linksModel,transfer_techs,gridSegments);
gridSegments_dcopf(linksModelToCalc,transfer_techs,gridSegments)
    $transfer_gridSegments(linksModelToCalc,transfer_techs,gridSegments,"useDCopf")
$iftheni.opfmethod %opfmethod%==disable
    = no;
$else.opfmethod
    = yes;
$endif.opfmethod

* Ensure each grid segment uses exactly one commodity
parameter checkGridSegmentCommodities(gridSegments);
checkGridSegmentCommodities(gridSegments)
    = sum (commodity$(sum((linksModelToCalc,transfer_techs,vintage)
                            $( transfer_coefficient(transfer_techs,vintage,commodity,"coefficient") > 0
                                and gridSegments_dcopf(linksModelToCalc,transfer_techs,gridSegments)) , 1) > 0), 1);

parameter transfer_incidenceSegments(nodesModel,linksModel,years,gridSegments);
transfer_incidenceSegments(nodesModelToCalc,linksModelToCalc,yearsToCalc,gridSegments)
    $(sum((transfer_techs,vintage)$(transfer_usedTech(linksModelToCalc,yearsToCalc,transfer_techs,vintage)
                                        and gridSegments_dcopf(linksModelToCalc,transfer_techs,gridSegments)
                                        and transfer_incidenceModel(nodesModelToCalc,linksModelToCalc) <> 0), 1) > 0)
    = transfer_incidenceModel(nodesModelToCalc,linksModelToCalc);


$onVerbatim
$iftheni.opfmethod %opfmethod%==kirchhoff
$offVerbatim

* Generate cycles based on the network for the Kirchhoff approach
parameter transfer_Cmat(cycles,linksModel,years,gridSegments);
embeddedCode Python:

    try:
        from remix.framework.tools.transfer import build_cmat
    except ImportError:
        from sys import path
        model_dir = Path(r'%sourcedir% '.strip()).parents[3].as_posix()
        if model_dir not in path:
            path.append(model_dir)
        from remix.framework.tools.transfer import build_cmat

    inc = gams.get("transfer_incidenceSegments")
    C_mat = build_cmat(inc)
    gams.set("transfer_Cmat", C_mat)
endEmbeddedCode transfer_Cmat

$onVerbatim
$endif.opfmethod
$offVerbatim

scalar transfer_enableMIP;
transfer_enableMIP = sum(transfer_usedTech(linksModelToCalc,yearsToCalc,transfer_techs,vintage)$transfer_techParam(transfer_techs,vintage,"mipLinks"), 1 );


* ==== definition of variables ====

* Initialise variables for linksBuild
transfer_linksBuild.l(linksModel,years,transfer_techs,vintage)
    $transfer_availTech(linksModel,years,transfer_techs,vintage)
    = transfer_linksParam(linksModel,years,transfer_techs,"linksBuild");
transfer_linksBuild.lo(linksModel,yearsToCalc,transfer_techs,vintage)
    $transfer_availTech(linksModel,yearsToCalc,transfer_techs,vintage)
    = transfer_linksBuild.l(linksModel,yearsToCalc,transfer_techs,vintage);
transfer_linksBuild.fx(linksModel,years,transfer_techs,vintage)
    $transfer_linksParam(linksModel,years,transfer_techs,"noExpansion")
    = transfer_linksBuild.l(linksModel,years,transfer_techs,vintage);

* Initialise variables for linksDecom
transfer_linksDecom.l(linksModel,years,transfer_techs,vintage)
    $(transfer_decomTech(linksModel,years,transfer_techs,vintage)
      and years.val < sum(yearsToCalc$(ord(yearsToCalc) = 1), yearsToCalc.val))
    = sum((years_a,years_aa)$(sameas(years-1, years_aa)
                      and years_a.val > years_aa.val - transfer_techParam(transfer_techs,vintage,'lifeTime')
                      and years_a.val <= years.val - transfer_techParam(transfer_techs,vintage,'lifeTime')
                      and transfer_availTech(linksModel,years_a,transfer_techs,vintage)),
        transfer_linksBuild.l(linksModel,years_a,transfer_techs,vintage));

transfer_linksDecom.l(linksModel,yearsToCalc,transfer_techs,vintage)
  $transfer_decomTech(linksModel,yearsToCalc,transfer_techs,vintage)
  = sum(years$
        (years.val < sum(yearsToCalc_a$(ord(yearsToCalc_a) = 1), yearsToCalc_a.val)
          and transfer_availTech(linksModel,years,transfer_techs,vintage)
          and years.val > sum(years_a$sameas(years_a, yearsToCalc-1), years_a.val) - transfer_techParam(transfer_techs,vintage,'lifeTime')
          and years.val <= yearsToCalc.val - transfer_techParam(transfer_techs,vintage,'lifeTime')),
      transfer_linksBuild.l(linksModel,years,transfer_techs,vintage))
    + sum(yearsToCalc_a$
        (yearsToCalc_a.val < sum(yearsToCalc_aa$(ord(yearsToCalc_aa) > 1), yearsToCalc_a.val)
          and transfer_availTech(linksModel,yearsToCalc_a,transfer_techs,vintage)
          and yearsToCalc_a.val > sum(years_a$sameas(years_a, yearsToCalc-1), years_a.val) - transfer_techParam(transfer_techs,vintage,'lifeTime')
          and yearsToCalc_a.val <= yearsToCalc.val - transfer_techParam(transfer_techs,vintage,'lifeTime')),
      transfer_linksBuild.l(linksModel,yearsToCalc_a,transfer_techs,vintage));
      ;

transfer_linksDecom.lo(linksModel,yearsToCalc,transfer_techs,vintage)
    $(transfer_usedTech(linksModel,yearsToCalc,transfer_techs,vintage)
        and not transfer_techParam(transfer_techs,vintage,"freeDecom"))
    = transfer_linksDecom.l(linksModel,yearsToCalc,transfer_techs,vintage)

* Calculate planned transfer links expansion
parameter transfer_linksPlanned(linksModel,years,transfer_techs,vintage);
transfer_linksPlanned(linksModel,years,transfer_techs,vintage) = 0;
loop(years,
  transfer_linksPlanned(linksModel,years,transfer_techs,vintage)
    =
    transfer_linksPlanned(linksModel,years-1,transfer_techs,vintage)
        $transfer_usedTech(linksModel,years-1,transfer_techs,vintage)
    + transfer_linksBuild.l(linksModel,years,transfer_techs,vintage)
        $transfer_availTech(linksModel,years,transfer_techs,vintage)
    - transfer_linksDecom.l(linksModel,years,transfer_techs,vintage)
        $transfer_usedTech(linksModel,years,transfer_techs,vintage);
);

* Set initial state for planned units
transfer_linksTotal.l(linksModel,years,transfer_techs,vintage)
  = transfer_linksPlanned(linksModel,years,transfer_techs,vintage)

* Calculate if planned links expansion is bound by upper and lower limits
set transfer_linkBoundsFixed(linksModel,years,transfer_techs);
transfer_linkBoundsFixed(linksModel,years,transfer_techs)
  $(sum(vintage$transfer_usedTech(linksModel,years,transfer_techs,vintage),
        transfer_linksPlanned(linksModel,years,transfer_techs,vintage))
    = transfer_linksParam(linksModel,years,transfer_techs,"linksUpperLimit")
  and sum(vintage$transfer_usedTech(linksModel,years,transfer_techs,vintage),
        transfer_linksPlanned(linksModel,years,transfer_techs,vintage))
    = transfer_linksParam(linksModel,years,transfer_techs,"linksLowerLimit"))
  = yes;

* Fix linksBuild, linksDecom, linksTotal if levels are predetermined by upper and lower limits
transfer_linksBuild.fx(linksModel,years,transfer_techs,vintage)
  $(transfer_availTech(linksModel,years,transfer_techs,vintage)
    and transfer_linkBoundsFixed(linksModel,years,transfer_techs))
  = transfer_linksBuild.l(linksModel,years,transfer_techs,vintage);
transfer_linksDecom.fx(linksModel,years,transfer_techs,vintage)
  $(transfer_usedTech(linksModel,years,transfer_techs,vintage)
    and transfer_linkBoundsFixed(linksModel,years,transfer_techs))
  = transfer_linksDecom.l(linksModel,years,transfer_techs,vintage);
transfer_linksTotal.fx(linksModel,years,transfer_techs,vintage)
  $(transfer_usedTech(linksModel,years,transfer_techs,vintage)
    and transfer_linkBoundsFixed(linksModel,years,transfer_techs))
  = transfer_linksTotal.l(linksModel,years,transfer_techs,vintage);

transfer_linksTotal_MIP.up(linksModel,years,transfer_techs,vintage)
    $(transfer_usedTech(linksModel,years,transfer_techs,vintage)
      and transfer_techParam(transfer_techs,vintage,"miplinks") = 1)
    = transfer_linksParam(linksModel,years,transfer_techs,"linksUpperLimit");

* Add parameter for fixing capacities during myopic runs
parameter transfer_linksDelta_upper(linksModel,years,transfer_techs);
parameter transfer_linksDelta_lower(linksModel,years,transfer_techs);


* ==== declaration of equations ====

equations
  Eq_transfer_linksBalance(linksModel,years,transfer_techs,vintage
    ) "Ensures the links balance over the planning period."
  Eq_transfer_linksFixedDecom(linksModel,years,transfer_techs,vintage
    ) "Restricts the fixed decommissioning of links over the planning period."
  Eq_transfer_linksFreeDecom(linksModel,years,transfer_techs,vintage
    ) "Restricts the free decommissioning of links over the planning period."
  Eq_transfer_linksUpperLimit(linksModel,years,transfer_techs
    ) "Upper bound for the total number of links."
  Eq_transfer_linksLowerLimit(linksModel,years,transfer_techs
    ) "Lower bound for the total number of links."
  Eq_transfer_linksTotalMIP(linksModel,years,transfer_techs,vintage
    ) "Fixes the total number of links to the corresponding integer variable."

  Eq_transfer_flowAlongUpperLimit(timeModel,linksModel,years,transfer_techs,vintage
    ) "Upper bound for the flow along the transfer links."
  Eq_transfer_flowAgainstUpperLimit(timeModel,linksModel,years,transfer_techs,vintage
    ) "Upper bound for the flow against the transfer links."
  Eq_transfer_dcopf_angleFlows(timeModel,linksModel,years,gridSegments
    )
  Eq_transfer_dcopf_cycleFlows(timeModel,years,cycles,gridSegments
    )
  ;


* ==== equation definition ====
* // ## Equations
* // ### Transfer Links Balance
* // Ensures that the transfer between nodes is balanced.
* // {Eq_transfer_linksBalance}
Eq_transfer_linksBalance(linksModelToCalc,yearsSel,transfer_techs,vintage)
    $((transfer_usedTech(linksModelToCalc,yearsSel,transfer_techs,vintage)
          or sum(years$sameas(years,yearsSel), transfer_usedTech(linksModelToCalc,years-1,transfer_techs,vintage)))
        and not transfer_linkBoundsFixed(linksModelToCalc,yearsSel,transfer_techs))
    ..
    transfer_linksTotal(linksModelToCalc,yearsSel,transfer_techs,vintage)
    =e=
    sum(yearsToCalc$(ord(yearsToCalc) = 1 and sameas(yearsToCalc, yearsSel)),
      sum(years$sameas(years, yearsToCalc),
        transfer_linksTotal(linksModelToCalc,years-1,transfer_techs,vintage)
          $transfer_usedTech(linksModelToCalc,years-1,transfer_techs,vintage)))
    + sum((yearsToCalc)$(ord(yearsToCalc) > 1 and sameas(yearsToCalc, yearsSel)),
      transfer_linksTotal(linksModelToCalc,yearsToCalc-1,transfer_techs,vintage)
        $transfer_usedTech(linksModelToCalc,yearsToCalc-1,transfer_techs,vintage))
    + transfer_linksBuild(linksModelToCalc,yearsSel,transfer_techs,vintage)
        $transfer_availTech(linksModelToCalc,yearsSel,transfer_techs,vintage)
    - transfer_linksDecom(linksModelToCalc,yearsSel,transfer_techs,vintage)
        $transfer_usedTech(linksModelToCalc,yearsSel,transfer_techs,vintage);

* // ### Transfer Links Fixed Decommissioning
* // Balances fixed link decommissioning.
* // {Eq_transfer_linksFixedDecom}
Eq_transfer_linksFixedDecom(linksModelToCalc,yearsSel,transfer_techs,vintage)
    $(transfer_decomTech(linksModelToCalc,yearsSel,transfer_techs,vintage)
        and not transfer_techParam(transfer_techs,vintage,"freeDecom") = 1
        and not transfer_linkBoundsFixed(linksModelToCalc,yearsSel,transfer_techs))
    ..
    transfer_linksDecom(linksModelToCalc,yearsSel,transfer_techs,vintage)
    =e=
    sum(years$
        (transfer_availTech(linksModelToCalc,years,transfer_techs,vintage)
            and years.val > sum(yearsToCalc$sameas(yearsToCalc+1, yearsSel), yearsToCalc.val) - transfer_techParam(transfer_techs,vintage,'lifeTime')
            and years.val <= yearsSel.val - transfer_techParam(transfer_techs,vintage,'lifeTime')),
        transfer_linksBuild(linksModelToCalc,years,transfer_techs,vintage));

* // ### Transfer Links Free Decommissioning
* // Balances free link decommissioning.
* // {Eq_transfer_linksFreeDecom}
Eq_transfer_linksFreeDecom(linksModelToCalc,yearsSel,transfer_techs,vintage)
    $((transfer_decomTech(linksModelToCalc,yearsSel,transfer_techs,vintage)
        or transfer_usedTech(linksModelToCalc,yearsSel,transfer_techs,vintage))
        and transfer_techParam(transfer_techs,vintage,"freeDecom") = 1)
    ..
    sum(years$
          ((transfer_decomTech(linksModelToCalc,years,transfer_techs,vintage)
            or transfer_usedTech(linksModelToCalc,years,transfer_techs,vintage))
            and years.val < sum(yearsToCalc_a$(ord(yearsToCalc_a)=1), yearsToCalc_a.val)),
        transfer_linksDecom(linksModelToCalc,years,transfer_techs,vintage))
    + sum(yearsToCalc$
          ((transfer_decomTech(linksModelToCalc,yearsToCalc,transfer_techs,vintage)
            or transfer_usedTech(linksModelToCalc,yearsToCalc,transfer_techs,vintage))
            and yearsToCalc.val >= sum(yearsToCalc_a$(ord(yearsToCalc_a)=1), yearsToCalc_a.val)
            and yearsToCalc.val <= yearsSel.val),
        transfer_linksDecom(linksModelToCalc,yearsToCalc,transfer_techs,vintage))
    =g=
    sum(years$
          (transfer_availTech(linksModelToCalc,years,transfer_techs,vintage)
            and years.val < sum(yearsToCalc$(ord(yearsToCalc)=1), yearsToCalc.val) - transfer_techParam(transfer_techs,vintage,'lifeTime')),
        transfer_linksBuild(linksModelToCalc,years,transfer_techs,vintage))
    + sum(yearsToCalc$
          (transfer_availTech(linksModelToCalc,yearsToCalc,transfer_techs,vintage)
            and yearsToCalc.val >= sum(yearsToCalc_a$(ord(yearsToCalc_a)=1), yearsToCalc_a.val)
            and yearsToCalc.val <= yearsSel.val - transfer_techParam(transfer_techs,vintage,'lifeTime')),
        transfer_linksBuild(linksModelToCalc,yearsToCalc,transfer_techs,vintage));

* // ### Transfer Links Lower Limit
* // Ensures link capacity lower limits.
* // {Eq_transfer_linksLowerLimit}
Eq_transfer_linksLowerLimit(linksModelToCalc,yearsSel,transfer_techs)
    $(transfer_linksParam(linksModelToCalc,yearsSel,transfer_techs,'linksLowerLimit') > 0 )
    ..
    sum(vintage$transfer_usedTech(linksModelToCalc,yearsSel,transfer_techs,vintage),
        transfer_linksTotal(linksModelToCalc,yearsSel,transfer_techs,vintage))
    =g=
    transfer_linksParam(linksModelToCalc,yearsSel,transfer_techs,"linksLowerLimit");

* // ### Transfer Links Upper Limit
* // Ensures link capacity upper limits.
* // {Eq_transfer_linksUpperLimit}
Eq_transfer_linksUpperLimit(linksModelToCalc,yearsSel,transfer_techs)
    $(transfer_linksParam(linksModelToCalc,yearsSel,transfer_techs,'linksUpperLimit') >= 0
        and transfer_linksParam(linksModelToCalc,yearsSel,transfer_techs,'linksUpperLimit') < +inf )
    ..
    sum(vintage$transfer_usedTech(linksModelToCalc,yearsSel,transfer_techs,vintage),
        transfer_linksTotal(linksModelToCalc,yearsSel,transfer_techs,vintage))
    =l=
    transfer_linksParam(linksModelToCalc,yearsSel,transfer_techs,"linksUpperLimit");

* // ### Transfer Links Total MIP
* // Ensures number of MIP links.
* // {Eq_transfer_linksTotalMIP}
Eq_transfer_linksTotalMIP(linksModelToCalc,yearsSel,transfer_techs,vintage)
    $(transfer_usedTech(linksModelToCalc,yearsSel,transfer_techs,vintage)
        and transfer_techParam(transfer_techs,vintage,"mipLinks"))
    ..
    transfer_linksTotal(linksModelToCalc,yearsSel,transfer_techs,vintage)
    =e=
    transfer_linksTotal_MIP(linksModelToCalc,yearsSel,transfer_techs,vintage);

* // ### Transfer Links Flow Along Upper Limit
* // Ensures links flow along upper limit.
* // {Eq_transfer_flowAlongUpperLimit}
Eq_transfer_flowAlongUpperLimit(timeModelSel,linksModelToCalc,yearsSel,transfer_techs,vintage)
    $transfer_usedTech(linksModelToCalc,yearsSel,transfer_techs,vintage)
    ..
    transfer_flowAlong(timeModelSel,linksModelToCalc,yearsSel,transfer_techs,vintage)
    =l=
    transfer_flowProfile(timeModelSel,linksModelToCalc,yearsSel,transfer_techs,vintage,"along")
    * transfer_linksTotal(linksModelToCalc,yearsSel,transfer_techs,vintage)
    * transfer_techParam(transfer_techs,vintage,"flowUpperLimit")
    ;

* // ### Transfer Links Flow Against Upper Limit
* // Ensures links flow against upper limit.
* // {Eq_transfer_flowAgainstUpperLimit}
Eq_transfer_flowAgainstUpperLimit(timeModelSel,linksModelToCalc,yearsSel,transfer_techs,vintage)
    $transfer_usedTech(linksModelToCalc,yearsSel,transfer_techs,vintage)
    ..
    transfer_flowAgainst(timeModelSel,linksModelToCalc,yearsSel,transfer_techs,vintage)
    =l=
    transfer_flowProfile(timeModelSel,linksModelToCalc,yearsSel,transfer_techs,vintage,"against")
    * transfer_linksTotal(linksModelToCalc,yearsSel,transfer_techs,vintage)
    * transfer_techParam(transfer_techs,vintage,"flowUpperLimit")
    ;

set transfer_usedOpf(linksModel,years,gridSegments,transfer_techs,vintage);
transfer_usedOpf(linksModel,yearsToCalc,gridSegments,transfer_techs,vintage)
    $(transfer_usedTech(linksModel,yearsToCalc,transfer_techs,vintage)
        and gridSegments_dcopf(linksModel,transfer_techs,gridSegments))
    = yes;

set transfer_usedOpfLinks(linksModel,years,gridSegments);
option transfer_usedOpfLinks < transfer_usedOpf;

set transfer_usedOpfSegments(years,gridSegments);
option transfer_usedOpfSegments < transfer_usedOpf;

$onVerbatim
$iftheni.opfmethod %opfmethod%==kirchhoff
$offVerbatim
* // ### Transfer DC optimal flow cycle flows
* // Optimal flow cycles.
* // {Eq_transfer_dcopf_cycleFlows}
Eq_transfer_dcopf_cycleFlows(timeModelSel,yearsSel,cycles,gridSegments)
    $transfer_usedOpfSegments(yearsSel,gridSegments)
    ..
    sum ((linksModelToCalc,transfer_techs,vintage)
            $transfer_usedOpf(linksModelToCalc,yearsSel,gridSegments,transfer_techs,vintage),
        transfer_Cmat(cycles,linksModelToCalc,yearsSel,gridSegments)
        * transfer_dcopf_Xtech(linksModelToCalc,yearsSel,transfer_techs,vintage,gridSegments)
        * (transfer_flowAlong(timeModelSel,linksModelToCalc,yearsSel,transfer_techs,vintage)
            - transfer_flowAgainst(timeModelSel,linksModelToCalc,yearsSel,transfer_techs,vintage)) )
    =e=
    0;

$onVerbatim
$elseifi.opfmethod %opfmethod%==angle
$offVerbatim
* // ### Transfer DC optimal flow angle flows
* // Angle flows.
* // {Eq_transfer_dcopf_cycleFlows}
Eq_transfer_dcopf_angleFlows(timeModelSel,linksModelToCalc,yearsSel,gridSegments)
    $transfer_usedOpfLinks(linksModelToCalc,yearsSel,gridSegments)
    ..
    sum ((transfer_techs,vintage)
            $transfer_usedOpf(linksModelToCalc,yearsSel,gridSegments,transfer_techs,vintage),
        ( transfer_flowAlong(timeModelSel,linksModelToCalc,yearsSel,transfer_techs,vintage)
            - transfer_flowAgainst(timeModelSel,linksModelToCalc,yearsSel,transfer_techs,vintage) )
        * transfer_dcopf_Xtech(linksModelToCalc,yearsSel,transfer_techs,vintage,gridSegments) )
    =e=
    sum ( nodesModelSel,
        - transfer_incidenceModel(nodesModelSel,linksModelToCalc)
        * transfer_dcopf_voltageAngle(timeModelSel,nodesModelSel,yearsSel,gridSegments));
$onVerbatim
$endif.opfmethod
$offVerbatim


* ==== model definition ====

Model M_transfer
/
  Eq_transfer_linksBalance
  Eq_transfer_linksFixedDecom
  Eq_transfer_linksFreeDecom
  Eq_transfer_linksLowerLimit
  Eq_transfer_linksUpperLimit
  Eq_transfer_linksTotalMIP
  Eq_transfer_flowAlongUpperLimit
  Eq_transfer_flowAgainstUpperLimit
$onVerbatim
$iftheni.opfmethod %opfmethod%==kirchhoff
$offVerbatim
  Eq_transfer_dcopf_cycleFlows
$onVerbatim
$elseifi.opfmethod %opfmethod%==angle
$offVerbatim
  Eq_transfer_dcopf_angleFlows
$onVerbatim
$endif.opfmethod
$offVerbatim
/;
