* SPDX-FileCopyrightText: Copyright (c) 2023-2024 German Aerospace Center (DLR)
* SPDX-License-Identifier: BSD-3-Clause

* // # accounting_input
* // The parameters in this file describe the accounting of indicators in the model.

* // ## Variables
* // {special_table_input_variables}
$offListing
$onEmpty

* The table headers in the parameter input sets represent the following:
* "Title | Description| Constraints | Type |  Units | Ontology"

* // ## Input Files
** // INPUT: accounting_converterUnits | IAO:0000100:data set
* // ### accounting_converterUnits
* // Title: Accounting Converter Units
* // Description: Comprises coefficients to account for indicators of converter units.
* // {table_accounting_converterUnits}
set pc_accounting_converterUnits
    /
    perUnitBuild        "Indicator per new unit | Indicator coefficient per new unit | | number | {indicator} | OEO_00000350:quantity value"
    perUnitDecom        "Indicator per decommissioned unit | Indicator coefficient per decommissioned unit | | number | {indicator} | OEO_00000350:quantity value"
    perUnitTotal        "Indicator per existing unit | Indicator coefficient per existing unit | | number | {indicator} | OEO_00000350:quantity value"
    useAnnuity          "Use Annuity | Specify if the cost should be calculated as annuities instead of absolute cost | | boolean | {none} | OEO_00000339:program parameter"
    amorTime            "Amortization Time | (Only required for annuity) The time used for the calculation of the annuity factors | minimum:0;needs:['useAnnuity'] | integer | {span} | MISSING_TERM:amortization time"
    interest            "Interest | (Only required for annuity) The interest rate used for the calculation of the annuity factors | needs:['useAnnuity'] | number | {none} | OEO_00240021:discount rate"
    /;
table accounting_converterUnitsIn(indicator,accNodesData,converter_techs,vintage,pc_accounting_converterUnits)
$onDelim
$if exist "%instancedir%/accounting_converterunits.csv" $include "%instancedir%/accounting_converterunits.csv"
$offDelim
$if not exist "%instancedir%/accounting_converterunits.csv" $log "No accounting for converter units included"
;

** // INPUT: accounting_converterActivity | IAO:0000100:data set
* // ### accounting_converterActivity
* // Title: Accounting Converter Activity
* // Description: Comprises coefficients to account for indicators that arise through activities performed by units.
* // {table_accounting_converterActivity}
set pc_accounting_converterActivity
    /
    perActivity         "Indicator per Activity | Indicator coefficient per unit activity | | number | {indicator} | OEO_00000350:quantity value"
    /;
table accounting_converterActivityIn(indicator,accNodesData,converter_techs,vintage,activity,pc_accounting_converterActivity)
$onDelim
$if exist "%instancedir%/accounting_converteractivity.csv" $include "%instancedir%/accounting_converteractivity.csv"
$offDelim
$if not exist "%instancedir%/accounting_converteractivity.csv" $log "No accounting for converter activities included"
;

** // INPUT: accounting_converterStartup | IAO:0000100:data set
* // ### accounting_converterStartup
* // Title: Accounting Converter Startup
* // Description: Comprises coefficients to account for costs arising from converter unit startups and/or ramping.
* // {table_accounting_converterStartup}
set pc_accounting_converterStartup
    /
    perStartup          "Indicator per Startup | Indicator coefficient per unit startup | | number | {indicator} | OEO_00000350:quantity value"
    perRamp             "Indicator per ramping | Indicator coefficient per unit ramping | | number | {indicator} | OEO_00000350:quantity value"
    perRampPos          "Indicator per positive ramping | Indicator coefficient per unit positive ramping | | number | {indicator} | OEO_00000350:quantity value"
    perRampNeg          "Indicator per negative ramping | Indicator coefficient per unit negative ramping | | number | {indicator} | OEO_00000350:quantity value"
    /;
table accounting_converterStartupIn(indicator,accNodesData,converter_techs,vintage,pc_accounting_converterStartup)
$onDelim
$if exist "%instancedir%/accounting_converterstartup.csv" $include "%instancedir%/accounting_converterstartup.csv"
$offDelim
$if not exist "%instancedir%/accounting_converterstartup.csv" $log "No accounting for converter startups included"
;

** // INPUT: accounting_storageUnits | IAO:0000100:data set
* // ### accounting_storageUnits
* // Title: Accounting Storage Units
* // Description: Comprises coefficients to account for costs of storage units.
* // {table_accounting_storageUnits}
set pc_accounting_storageUnits
    /
    perUnitBuild        "Indicator per new unit | Indicator coefficient per new unit. | | number | {indicator} | OEO_00000350:quantity value"
    perUnitDecom        "Indicator per decommissioned unit | Indicator coefficient per decommissioned unit. | | number | {indicator} | OEO_00000350:quantity value"
    perUnitTotal        "Indicator per existing unit | Indicator coefficient per existing unit. | | number | {indicator} | OEO_00000350:quantity value"
    useAnnuity          "Use Annuity | Specify if the cost should be calculated as annuities instead of absolute cost | | boolean | {none} | OEO_00000339:program parameter"
    amorTime            "Amortisation Time | (Only required for annuity) The time used for the calculation of the annuity factors | minimum:0;needs:['useAnnuity'] | integer | {span} | MISSING_TERM:amortization time"
    interest            "Interest | (Only required for annuity) The interest rate used for the calculation of the annuity factors | needs:['useAnnuity'] | number | {none} | OEO_00240021:discount rate"
    /;
table accounting_storageUnitsIn(indicator,accNodesData,storage_techs,vintage,pc_accounting_storageUnits)
$onDelim
$if exist "%instancedir%/accounting_storageunits.csv" $include "%instancedir%/accounting_storageunits.csv"
$offDelim
$if not exist "%instancedir%/accounting_storageunits.csv" $log "No accounting for storage units included"
;

** // INPUT: accounting_transferLinks | IAO:0000100:data set
* // ### accounting_transferLinks
* // Title: Accounting Transfer Per Link
* // Description: Comprises coefficients to account for costs arising from building and using transfers per link built.
* // {table_accounting_transferLinks}
set pc_accounting_transferLinks
    /
    perLinkBuild        "Indicator per new link | Indicator coefficient per new link | | number | {indicator} | OEO_00000350:quantity value"
    perLinkDecom        "Indicator per decommissioned link | Indicator coefficient per decommissioned link | | number | {indicator} | OEO_00000350:quantity value"
    perLinkTotal        "Indicator per existing link | Indicator coefficient per existing link | | number | {indicator} | OEO_00000350:quantity value"
    perFlow             "Indicator per flow on a link | Indicator coefficient per flow on a link independent of direction | | number | {indicator} | OEO_00000350:quantity value"
    perFlowAlong        "Indicator per flow along the link | Indicator coefficient per flow along the direction of the link | | number | {indicator} | OEO_00000350:quantity value"
    perFlowAgainst      "Indicator per flow against the link | Indicator coefficient per flow against the direction of the link | | number | {indicator} | OEO_00000350:quantity value"
    useAnnuity          "Use Annuity | Specify if the cost should be calculated as annuities instead of absolute cost | | boolean | {none} | OEO_00000339:program parameter"
    amorTime            "Amortisation Time | (Only required for annuity) The time used for the calculation of the annuity factors | minimum:0;needs:['useAnnuity'] | integer | {span} | MISSING_TERM:amortization time"
    interest            "Interest | (Only required for annuity) The interest rate used for the calculation of the annuity factors | needs:['useAnnuity'] | number | {none} | OEO_00240021:discount rate"
    /;
table accounting_transferLinksIn(indicator,accLinksData,transfer_techs,vintage,pc_accounting_transferLinks)
$onDelim
$if exist "%instancedir%/accounting_transferlinks.csv" $include "%instancedir%/accounting_transferlinks.csv"
$offDelim
$if not exist "%instancedir%/accounting_transferlinks.csv" $log "No accounting for transfers per link included"
;

** // INPUT: accounting_transferPerLength | IAO:0000100:data set
* // ### accounting_transferPerLength
* // Title: Accounting Transfer Per Length
* // Description: Comprises coefficients to account for costs arising from building and using transfers relative to the link lengths.
* // {table_accounting_transferPerLength}
set pc_accounting_transferPerLength
    /
    perLengthBuild    "Indicator per new link and length | Indicator coefficient per new link times length | | number | {indicator} | OEO_00000350:quantity value"
    perLengthDecom    "Indicator per decommissioned link and length | Indicator coefficient per decommissioned link times length | | number | {indicator} | OEO_00000350:quantity value"
    perLengthTotal    "Indicator per existing link and length | Indicator coefficient per existing link times length | | number | {indicator} | OEO_00000350:quantity value"
    perFlow             "Indicator per flow and length on a link | Indicator coefficient per flow times length on a link independent of direction | | number | {indicator} | OEO_00000350:quantity value"
    perFlowAlong        "Indicator per flow and length along the link | Indicator coefficient per flow times length along the direction of the link | | number | {indicator} | OEO_00000350:quantity value"
    perFlowAgainst      "Indicator per flow and length against the link | Indicator coefficient per flow times length against the direction of the link | | number | {indicator} | OEO_00000350:quantity value"
    useAnnuity          "Use Annuity | Specify if the cost should be calculated as annuities instead of absolute cost | | boolean | {none} | OEO_00000339:program parameter"
    amorTime            "Amortization Time | (Only required for annuity) The time used for the calculation of the annuity factors | minimum:0;needs:['useAnnuity'] | integer | {span} | MISSING_TERM:amortization time"
    interest            "Interest | (Only required for annuity) The interest rate used for the calculation of the annuity factors | needs:['useAnnuity'] | number | {none} | OEO_00240021:discount rate"
    /;
table accounting_transferPerLengthIn(indicator,accLinksData,transfer_techs,vintage,link_types,pc_accounting_transferPerLength)
$onDelim
$if exist "%instancedir%/accounting_transferperlength.csv" $include "%instancedir%/accounting_transferperlength.csv"
$offDelim
$if not exist "%instancedir%/accounting_transferperlength.csv" $log "No accounting for transfers per length included"
;

** // INPUT: accounting_sourcesinkFlow | IAO:0000100:data set
* // ### accounting_sourcesinkFlow
* // Title: Accounting Sources and Sinks Flow
* // Description: Comprises coefficients to account for costs arising from sources and sinks per commodity flow or per peak.
* // {table_accounting_sourcesinkFlow}
set pc_accounting_sourcesinkFlow
    /
    perFlow             "Indicator per commodity flow | Indicator coefficient per commodity flow across system boundaries | | number | {indicator} | OEO_00000350:quantity value"
    perPeak             "Indicator per commodity peak | Indicator coefficient per peak (only applied to rescalable profiles) | | number | {indicator} | OEO_00000350:quantity value"
    /;
table accounting_sourcesinkFlowIn(indicator,accNodesData,years,sourcesink_techs,commodity,pc_accounting_sourcesinkFlow)
$onDelim
$if exist "%instancedir%/accounting_sourcesinkflow.csv" $include "%instancedir%/accounting_sourcesinkflow.csv"
$offDelim
$if not exist "%instancedir%/accounting_sourcesinkflow.csv" $log "No accounting for sourcesink flows included"
;

** // INPUT: accounting_perIndicator | IAO:0000100:data set
* // ### accounting_perIndicator
* // Title: Accounting per Indicator
* // Description: Accounting per indicator describe accounting coefficients per indicator used for endpoint indicators.
* // {table_accounting_perIndicator}
set pc_accounting_perIndicator
    /
    perIndicator        "Indicator per indicator | Indicator coefficient per indicator used for endpoint indicators | | number | {indicator} | OEO_00000350:quantity value"
    /;
table accounting_perIndicatorIn(indicator,indicator,accNodesData,accYears,pc_accounting_perIndicator)
$onDelim
$if exist "%instancedir%/accounting_perindicator.csv" $include "%instancedir%/accounting_perindicator.csv"
$offDelim
$if not exist "%instancedir%/accounting_perindicator.csv" $log "No accounting for indicator based indicators included"
;

** // INPUT: accounting_indicatorBounds | IAO:0000100:data set
* // ### accounting_indicatorBounds
* // Title: Accounting Indicator Bounds
* // Description: Accounting indicator bounds describe overall indicator bounds such as the objective value, social discount rate, lower and upper bounds for the indicator etc.
* // {table_accounting_indicatorBounds}
set pc_accounting_indicatorBounds
    /
    obj                 "Indicator Objective | Set the indicator to be the objective value, -1 for minimization, 1 for maximization | enum:[1,0,-1] | integer | {none} | OEO_00240004:objective variable"
    isVariable          "Indicator as variable | Allows to declare an indicator as a free variable instead of an indirect variable | | boolean | | OEO_00000435:variable"
    useFixed            "Use Fixed Value | Use a fixed value for the indicator | | boolean | {none} | OEO_00000339:program parameter"
    fixedValue          "Fixed Value | Value to be fixed for the indicator | | number | {indicator} | OEO_00000104:constraint"
    useLower            "Use Lower Value | Use a lower bound for the indicator | | boolean | {none} | OEO_00000339:program parameter"
    lowerValue          "Lower Value | Value of the lower bound for the indicator | | number | {indicator} | OEO_00000104:constraint"
    useUpper            "Use Upper | Use an upper bound for the indicatorr | | boolean | {none} | OEO_00000339:program parameter"
    upperValue          "Upper Value | Value of the upper bound for the indicator | | number | {indicator} | OEO_00000104:constraint"
    discount            "Discount | The social discount rate for the indicator | | number | {indicator} | MISSING_TERM:social discount rate"
    integral            "Integral Indicator | If enabled, the indicator will also be accounted for the years between the year slices | | boolean | {none} | OEO_00000339:program parameter"
    endyear             "End Year for Integral | (Only required for integral indicators) The number of years used for the last year slice or perpetuity | | integer | {none} | OEO_00000339:program parameter"
    mga                 "Modeling to Generate Alternatives (MGA) | Set the indicator to a MGA distance to be maximized | | boolean | {none} | OEO_00000339:program parameter"
    pareto              "Pareto-front indicator | Set the secondary indicator for the Pareto-front indicator, -1 for minimization, 1 for maximization | enum:[1,0,-1] | integer | {none} | OEO_00000339:program parameter"
    /;
table accounting_indicatorBoundsIn(accNodesData,accYears,indicator,pc_accounting_indicatorBounds)
$onDelim
$if exist "%instancedir%/accounting_indicatorbounds.csv" $include "%instancedir%/accounting_indicatorbounds.csv"
$offDelim
$if not exist "%instancedir%/accounting_indicatorbounds.csv" $log "No bounds for indicator included - no objective function set!"
;

** // INPUT: accounting_indicatorBounds_links | IAO:0000100:data set
* // ### accounting_indicatorBounds_links
* // Title: Accounting Indicator Bounds Links
* // Description: Accounting indicator bounds per links describe fixed, upper and lower values for the indicator for each link.
* // {table_accounting_indicatorBounds_links}
set pc_accounting_indicatorBounds_links
    /
    useFixed            "Use fixed value | Use a fixed value for the indicator | | boolean | {none} | OEO_00000339:program parameter"
    fixedValue          "Fixed value | Value to be fixed for the indicator | | number | {indicator} | OEO_00000104:constraint"
    useLower            "Use lower value | Use a lower bound for the indicator | | boolean | {none} | OEO_00000339:program parameter"
    lowerValue          "Lower value | Value of the lower bound for the indicator | | number | {indicator} | OEO_00000104:constraint"
    useUpper            "Use upper value | Use an upper bound for the indicator | | boolean | {none} | OEO_00000339:program parameter"
    upperValue          "Upper value | Value of the upper bound for the indicator | | number | {indicator} | OEO_00000104:constraint"
    /;
table accounting_indicatorBounds_linksIn(linksData,years,indicator,pc_accounting_indicatorBounds_links)
$onDelim
$if exist "%instancedir%/accounting_indicatorbounds_links.csv" $include "%instancedir%/accounting_indicatorbounds_links.csv"
$offDelim
$if not exist "%instancedir%/accounting_indicatorbounds_links.csv" $log "No bounds for indicator per links included"
;
$offEmpty
$onListing

* Aggregate accounting parameters
$batinclude %aggregateAccountingMean% accounting_converterUnits indicator nodes converter_techs,vintage 1
$batinclude %aggregateAccountingMean% accounting_converterActivity indicator nodes converter_techs,vintage,activity 0
$batinclude %aggregateAccountingMean% accounting_converterStartup indicator nodes converter_techs,vintage 0
$batinclude %aggregateAccountingMean% accounting_storageUnits indicator nodes storage_techs,vintage 1
$batinclude %aggregateAccountingMean% accounting_transferLinks indicator links transfer_techs,vintage 1
$batinclude %aggregateAccountingMean% accounting_transferPerLength indicator links transfer_techs,vintage,link_types 1
$batinclude %aggregateAccountingMean% accounting_sourcesinkFlow indicator nodes years,sourcesink_techs,commodity 0

parameter accounting_indicatorBounds_links(linksModel,years,indicator,pc_accounting_indicatorBounds_links);
accounting_indicatorBounds_links(linksModelToCalc,yearsToCalc,indicator,pc_accounting_indicatorBounds_links)
    = sum((linksData)$sameas(linksModelToCalc,linksData),
            accounting_indicatorBounds_linksIn(linksData,yearsToCalc,indicator,pc_accounting_indicatorBounds_links));

accounting_indicatorBounds_links(linksModelToCalc,yearsToCalc,indicator,pc_accounting_indicatorBounds_links)
    = sum((linksData)$links_aggregate(linksModelToCalc,linksData),
            accounting_indicatorBounds_linksIn(linksData,yearsToCalc,indicator,pc_accounting_indicatorBounds_links));

parameter accounting_indicatorBounds(accNodesModel,accYears,indicator,pc_accounting_indicatorBounds);
accounting_indicatorBounds(accNodesModel,accYears,indicator,pc_accounting_indicatorBounds)
    = sum (accNodesData$sameas(accNodesModel,accNodesData),
        accounting_indicatorBoundsIn(accNodesData,accYears,indicator,pc_accounting_indicatorBounds));

accounting_indicatorBounds(accNodesModel,accYears,indicator,pc_accounting_indicatorBounds)
    $(sum (nodesModel$sameas(accNodesModel,nodesModel), 1) > 0)
    = sum ((accNodesData,nodesData,nodesModel)
            $(aggregateNodesModel(nodesData,nodesModel) and sameas(accNodesModel,nodesModel) and sameas(accNodesData,nodesData)),
        accounting_indicatorBoundsIn(accNodesData,accYears,indicator,pc_accounting_indicatorBounds));

set accounting_perIndicatorNonzero(indicator,indicator_a,accNodesData,accYears);
accounting_perIndicatorNonzero(indicator,indicator_a,accNodesData,accYears)
    $sum(pc_accounting_perIndicator$accounting_perIndicatorIn(indicator,indicator_a,accNodesData,accYears,pc_accounting_perIndicator), 1)
    = yes;

parameter accounting_perIndicatorAgg(indicator,indicator_a,nodesModel,accYears,pc_accounting_perIndicator);
accounting_perIndicatorAgg(indicator,indicator_a,nodesModel,accYears,pc_accounting_perIndicator)
    $sum((nodesData,accnodesData)
            $(aggregatenodesModel(nodesData,nodesModel) and sameas(nodesData,accnodesData)
                and accounting_perIndicatorNonzero(indicator,indicator_a,accnodesData,accYears)), 1)
    = sum((nodesData,accnodesData)
            $(aggregatenodesModel(nodesData,nodesModel) and sameas(nodesData,accnodesData)
                and accounting_perIndicatorNonzero(indicator,indicator_a,accnodesData,accYears)),
        accounting_perIndicatorIn(indicator,indicator_a,accnodesData,accYears,pc_accounting_perIndicator))
    / sum((nodesData,accnodesData)
            $(aggregatenodesModel(nodesData,nodesModel) and sameas(nodesData,accnodesData)
                and accounting_perIndicatorNonzero(indicator,indicator_a,accnodesData,accYears)),
        1);

parameter accounting_perIndicator(indicator,indicator_a,accNodesModel,accYears,pc_accounting_perIndicator);

loop(accNodes,
accounting_perIndicator(indicator,indicator_a,accNodesModel,accYears,pc_accounting_perIndicator)
    $(sum(accNodesData$(sameas(accNodes,accNodesData)
            and accounting_perIndicatorIn(indicator,indicator_a,accNodesData,"horizon",pc_accounting_perIndicator)), 1)
        and map_accNodes(accNodesModel,accNodes)
            )
    = sum(accNodesData$(sameas(accNodes,accNodesData)),
            accounting_perIndicatorIn(indicator,indicator_a,accNodesData,"horizon",pc_accounting_perIndicator));

accounting_perIndicator(indicator,indicator_a,accNodesModel,accYears,pc_accounting_perIndicator)
    $(sum(accNodesData$(sameas(accNodes,accNodesData)
            and accounting_perIndicatorIn(indicator,indicator_a,accNodesData,accYears,pc_accounting_perIndicator)), 1)
        and map_accNodes(accNodesModel,accNodes)
        and not sameas(accYears, "horizon"))
    = sum(accNodesData$(sameas(accNodes,accNodesData)),
            accounting_perIndicatorIn(indicator,indicator_a,accNodesData,accYears,pc_accounting_perIndicator));
);

accounting_perIndicator(indicator,indicator_a,accNodesModel,accYears,pc_accounting_perIndicator)
    $sum(nodesModel$(accounting_perIndicatorAgg(indicator,indicator_a,nodesModel,"horizon",pc_accounting_perIndicator)
        and sameas(accNodesModel,nodesModel)), 1)
    = sum(nodesModel$sameas(accNodesModel,nodesModel),
        accounting_perIndicatorAgg(indicator,indicator_a,nodesModel,"horizon",pc_accounting_perIndicator));

accounting_perIndicator(indicator,indicator_a,accNodesModel,accYears,pc_accounting_perIndicator)
    $sum(nodesModel$(accounting_perIndicatorAgg(indicator,indicator_a,nodesModel,accYears,pc_accounting_perIndicator)
        and sameas(accNodesModel,nodesModel) and not sameas(accYears, "horizon")), 1)
    = sum(nodesModel$sameas(accNodesModel,nodesModel),
        accounting_perIndicatorAgg(indicator,indicator_a,nodesModel,accYears,pc_accounting_perIndicator));


* Make sure all required indicators end up in the model
set activeIndicators(accNodesModel,accYears,indicator);
activeIndicators(accNodesModel,accYears,indicator)
   $((accounting_indicatorBounds(accNodesModel,accYears,indicator,"obj") <> 0
$iftheni.mga %method%==mga
      or accounting_indicatorBounds(accNodesModel,accYears,indicator,"mga") <> 0
$endif.mga
$iftheni.pareto %method%==pareto
      or accounting_indicatorBounds(accNodesModel,accYears,indicator,"pareto") <> 0
$endif.pareto
      or accounting_indicatorBounds(accNodesModel,accYears,indicator,"useFixed") <> 0
      or accounting_indicatorBounds(accNodesModel,accYears,indicator,"useLower") <> 0
      or accounting_indicatorBounds(accNodesModel,accYears,indicator,"useUpper") <> 0 )
        and not accounting_indicatorBounds(accNodesModel,accYears,indicator,"isVariable"))
   = yes;

set activeIndicators_links(linksModel,years,indicator);
activeIndicators_links(linksModelToCalc,yearsToCalc,indicator)
   $(accounting_indicatorBounds_links(linksModelToCalc,yearsToCalc,indicator,"useFixed") <> 0
      or accounting_indicatorBounds_links(linksModelToCalc,yearsToCalc,indicator,"useLower") <> 0
      or accounting_indicatorBounds_links(linksModelToCalc,yearsToCalc,indicator,"useUpper") <> 0 )
   = yes;

set variableIndicators(accNodesModel,accYears,indicator);
variableIndicators(accNodesModel,accYears,indicator)
   $(accounting_indicatorBounds(accNodesModel,accYears,indicator,"obj") = 0
      and accounting_indicatorBounds(accNodesModel,accYears,indicator,"isVariable"))
   = yes;

* Calculate length of years and discount rates per indicator
parameter yearFactor(accNodesModel,accYears,indicator,accYears);
yearFactor(accNodesModel,accYears,indicator,accYears_a)
    $(map_accYears(accYears_a,accYears)
        and sum(yearsToCalc$(sameas(yearsToCalc,accYears_a)), 1))
    = 1;

yearFactor(accNodesModel,accYears,indicator,accYears_a)
    $(map_accYears(accYears_a,accYears)
        and sum(yearsToCalc$(sameas(yearsToCalc,accYears_a)), 1)
        and accounting_indicatorBounds(accNodesModel,accYears,indicator,"integral"))
    = sum(yearsToCalc$sameas(accYears_a,yearsToCalc), yearsLen(yearsToCalc));

yearFactor(accNodesModel,accYears,indicator,accYears_a)
    $(map_accYears(accYears_a,accYears)
        and sum(yearsToCalc$(sameas(yearsToCalc,accYears_a)), 1)
        and accounting_indicatorBounds(accNodesModel,accYears,indicator,"integral")
        and accounting_indicatorBounds(accNodesModel,accYears,indicator,"endyear")
        and yearFactor(accNodesModel,accYears,indicator,accYears_a) = inf)
    = accounting_indicatorBounds(accNodesModel,accYears,indicator,"endyear");

yearFactor(accNodesModel,accYears,indicator,accYears_a)
    $(map_accYears(accYears_a,accYears)
        and sum(yearsToCalc$(sameas(yearsToCalc,accYears_a)), 1)
        and accounting_indicatorBounds(accNodesModel,accYears,indicator,"integral")
        and accounting_indicatorBounds(accNodesModel,accYears,indicator,"discount"))
    = yearFactor(accNodesModel,accYears,indicator,accYears_a)
        * (1 - accounting_indicatorBounds(accNodesModel,accYears,indicator,"discount"))
                ** (accYears_a.val - smin(yearsToCalc, yearsToCalc.val));


* ==== calculate compound indicators for the optimization ====
scalar compIndicators_pre;
scalar compIndicators_post;

parameter compoundIndicators(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a);
compoundIndicators(accNodesModel,accYears,indicator,accNodesModel,accYears,indicator)
    $activeIndicators(accNodesModel,accYears,indicator) = 1;

set compoundIndicators_act(accNodesModel,accYears,indicator);
option compoundIndicators_act < compoundIndicators;

parameter compoundIndicatorsExt(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a,accNodesModel_aa,accYears_aa,indicator_aa);
parameter compoundIndicatorsExt_r(accNodesModel,accYears,indicator,accNodesModel_aa,accYears_aa,indicator_aa,accNodesModel_a,accYears_a,indicator_a);
compoundIndicatorsExt(accNodesModel,accYears,indicator,accNodesModel,accYears,indicator,accNodesModel,accYears,indicator)
    $compoundIndicators_act(accNodesModel,accYears,indicator) = 1;

compIndicators_pre = 0;
compIndicators_post = 1;

option sparseval = 1;
while(compIndicators_pre < compIndicators_post,
    compIndicators_pre = sum((accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a,accNodesModel_aa,accYears_aa,indicator_aa)
        $compoundIndicatorsExt(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a,accNodesModel_aa,accYears_aa,indicator_aa), 1);

    compoundIndicatorsExt(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a,accNodesModel_aa,accYears_aa,indicator_aa)
        $(compoundIndicators_act(accNodesModel,accYears,indicator)
            and compoundIndicators(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a)
            and not accounting_indicatorBounds(accNodesModel_a,accYears_a,indicator_a,"isVariable")
            and map_accNodes(accNodesModel_aa,accNodesModel_a)
            and map_accYears(accYears_aa,accYears_a))
        = compoundIndicators(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a)
            * accounting_perIndicator(indicator_a,indicator_aa,accNodesModel_aa,accYears_aa,"perIndicator");

    option compoundIndicatorsExt_r < compoundIndicatorsExt;

    compoundIndicators(accNodesModel,accYears,indicator,accNodesModel_aa,accYears_aa,indicator_aa)
        $(compoundIndicators_act(accNodesModel,accYears,indicator)
            and sum((accNodesModel_a,accYears_a,indicator_a)
                    $compoundIndicatorsExt_r(accNodesModel,accYears,indicator,accNodesModel_aa,accYears_aa,indicator_aa,accNodesModel_a,accYears_a,indicator_a), 1) = 1)
        = sum((accNodesModel_a,accYears_a,indicator_a)
                $(compoundIndicatorsExt_r(accNodesModel,accYears,indicator,accNodesModel_aa,accYears_aa,indicator_aa,accNodesModel_a,accYears_a,indicator_a)),
            compoundIndicatorsExt_r(accNodesModel,accYears,indicator,accNodesModel_aa,accYears_aa,indicator_aa,accNodesModel_a,accYears_a,indicator_a));

    compoundIndicators(accNodesModel,accYears,indicator,accNodesModel_aa,accYears_aa,indicator_aa)
        $(compoundIndicators_act(accNodesModel,accYears,indicator)
            and sum((accNodesModel_a,accYears_a,indicator_a)
                    $compoundIndicatorsExt_r(accNodesModel,accYears,indicator,accNodesModel_aa,accYears_aa,indicator_aa,accNodesModel_a,accYears_a,indicator_a), 1) > 1)
        = sum((accNodesModel_a,accYears_a,indicator_a)
                $(compoundIndicatorsExt_r(accNodesModel,accYears,indicator,accNodesModel_aa,accYears_aa,indicator_aa,accNodesModel_a,accYears_a,indicator_a)
                    and not (sameas(accNodesModel, accNodesModel_a) and sameas(accYears, accYears_a) and sameas(indicator, indicator_a))),
            compoundIndicatorsExt_r(accNodesModel,accYears,indicator,accNodesModel_aa,accYears_aa,indicator_aa,accNodesModel_a,accYears_a,indicator_a))
        / sum((accNodesModel_a,accYears_a,indicator_a)
                $(compoundIndicatorsExt_r(accNodesModel,accYears,indicator,accNodesModel_aa,accYears_aa,indicator_aa,accNodesModel_a,accYears_a,indicator_a)
                    and not (sameas(accNodesModel, accNodesModel_a) and sameas(accYears, accYears_a) and sameas(indicator, indicator_a))), 1);

    option compoundIndicators_act < compoundIndicators;
    compIndicators_post = sum((accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a,accNodesModel_aa,accYears_aa,indicator_aa)
        $compoundIndicatorsExt(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a,accNodesModel_aa,accYears_aa,indicator_aa), 1);
);
option sparseval = 0;

* map to sub-years and sub-nodes
compoundIndicatorsExt(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a,accNodesModel_aa,accYears_aa,indicator_a)
    $(compoundIndicators_act(accNodesModel,accYears,indicator)
        and compoundIndicators(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a)
        and map_accYears(accYears_aa,accYears_a)
        and map_accNodes(accNodesModel_aa,accNodesModel_a))
    = compoundIndicators(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a)
        * yearFactor(accNodesModel,accYears,indicator,accYears_aa);

* Remove all accounting regions and accounting years
compoundIndicatorsExt(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a,accNodesModel_aa,accYears_aa,indicator_aa)
    $(compoundIndicatorsExt(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a,accNodesModel_aa,accYears_aa,indicator_aa)
        and not accounting_indicatorBounds(accNodesModel_aa,accYears_aa,indicator_aa,"isVariable")
        and sum(accNodes$(sameas(accNodes,accNodesModel_aa)), 1)
            or sameas("horizon",accYears_aa))
    = 0;

* Map from extended parameter to final parameter
option compoundIndicatorsExt_r < compoundIndicatorsExt;
compoundIndicators(accNodesModel,accYears,indicator,accNodesModel_aa,accYears_aa,indicator_aa)
    $(compoundIndicators_act(accNodesModel,accYears,indicator)
        and sum((accNodesModel_a,accYears_a,indicator_a)
                $compoundIndicatorsExt_r(accNodesModel,accYears,indicator,accNodesModel_aa,accYears_aa,indicator_aa,accNodesModel_a,accYears_a,indicator_a), 1) = 1)
    = sum((accNodesModel_a,accYears_a,indicator_a)
            $(compoundIndicatorsExt_r(accNodesModel,accYears,indicator,accNodesModel_aa,accYears_aa,indicator_aa,accNodesModel_a,accYears_a,indicator_a)),
        compoundIndicatorsExt_r(accNodesModel,accYears,indicator,accNodesModel_aa,accYears_aa,indicator_aa,accNodesModel_a,accYears_a,indicator_a));

compoundIndicators(accNodesModel,accYears,indicator,accNodesModel_aa,accYears_aa,indicator_aa)
    $(compoundIndicators_act(accNodesModel,accYears,indicator)
        and sum((accNodesModel_a,accYears_a,indicator_a)
                $compoundIndicatorsExt_r(accNodesModel,accYears,indicator,accNodesModel_aa,accYears_aa,indicator_aa,accNodesModel_a,accYears_a,indicator_a), 1) > 1)
    = sum((accNodesModel_a,accYears_a,indicator_a)
            $(compoundIndicatorsExt_r(accNodesModel,accYears,indicator,accNodesModel_aa,accYears_aa,indicator_aa,accNodesModel_a,accYears_a,indicator_a)
                and not (sameas(accNodesModel, accNodesModel_a) and sameas(accYears, accYears_a) and sameas(indicator, indicator_a))),
        compoundIndicatorsExt_r(accNodesModel,accYears,indicator,accNodesModel_aa,accYears_aa,indicator_aa,accNodesModel_a,accYears_a,indicator_a))
    / sum((accNodesModel_a,accYears_a,indicator_a)
            $(compoundIndicatorsExt_r(accNodesModel,accYears,indicator,accNodesModel_aa,accYears_aa,indicator_aa,accNodesModel_a,accYears_a,indicator_a)
                and not (sameas(accNodesModel, accNodesModel_a) and sameas(accYears, accYears_a) and sameas(indicator, indicator_a))), 1);

* Remove all accounting regions and accounting years
compoundIndicators(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a)
    $(compoundIndicators(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a)
        and not accounting_indicatorBounds(accNodesModel_a,accYears_a,indicator_a,"isVariable")
        and (sum(accNodes$(sameas(accNodes,accNodesModel_a)), 1)
            or sameas("horizon",accYears_a)))
    = 0;

* Remove all slack indicators except the ones declared in indicatorBounds
compoundIndicators(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a)
    $(compoundIndicators(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a)
        and sum((accNodesModel_aa,accYears_aa)$accounting_indicatorBounds(accNodesModel_aa,accYears_aa,indicator_a,"isVariable"), 1)
        and not accounting_indicatorBounds(accNodesModel_a,accYears_a,indicator_a,"isVariable"))
    = 0;

* ==== compound indicators for the post calculation ====
parameter compoundIndicatorsFull(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a);
compoundIndicatorsFull(accNodesModel,accYears,indicator,accNodesModel,accYears,indicator)
    $(sum(nodesModelToCalc$map_accNodesToCalc(accNodesModel,nodesModelToCalc), 1)
        and sum(yearsToCalc$map_accYearsToCalc(accYears,yearsToCalc), 1)) = 1;

set compoundIndicatorsFull_act(accNodesModel,accYears,indicator);
option compoundIndicatorsFull_act < compoundIndicatorsFull;

parameter compoundIndicatorsFullExt(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a,accNodesModel_aa,accYears_aa,indicator_aa);
parameter compoundIndicatorsFullExt_r(accNodesModel,accYears,indicator,accNodesModel_aa,accYears_aa,indicator_aa,accNodesModel_a,accYears_a,indicator_a);
compoundIndicatorsFullExt(accNodesModel,accYears,indicator,accNodesModel,accYears,indicator,accNodesModel,accYears,indicator) = 1;

compIndicators_pre = 0;
compIndicators_post = 1;

option sparseval = 1;
while(compIndicators_pre < compIndicators_post,
    compIndicators_pre = sum((accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a,accNodesModel_aa,accYears_aa,indicator_aa)
        $compoundIndicatorsFullExt(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a,accNodesModel_aa,accYears_aa,indicator_aa), 1);

    compoundIndicatorsFullExt(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a,accNodesModel_aa,accYears_aa,indicator_aa)
        $(compoundIndicatorsFull_act(accNodesModel,accYears,indicator)
            and compoundIndicatorsFull(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a)
            and not accounting_indicatorBounds(accNodesModel_a,accYears_a,indicator_a,"isVariable")
            and map_accNodes(accNodesModel_aa,accNodesModel_a)
            and map_accYears(accYears_aa,accYears_a))
        = compoundIndicatorsFull(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a)
            * accounting_perIndicator(indicator_a,indicator_aa,accNodesModel_aa,accYears_aa,"perIndicator");

    option compoundIndicatorsFullExt_r < compoundIndicatorsFullExt;

    compoundIndicatorsFull(accNodesModel,accYears,indicator,accNodesModel_aa,accYears_aa,indicator_aa)
        $(compoundIndicatorsFull_act(accNodesModel,accYears,indicator)
            and sum((accNodesModel_a,accYears_a,indicator_a)
                    $compoundIndicatorsFullExt_r(accNodesModel,accYears,indicator,accNodesModel_aa,accYears_aa,indicator_aa,accNodesModel_a,accYears_a,indicator_a), 1) = 1)
        = sum((accNodesModel_a,accYears_a,indicator_a)
                $(compoundIndicatorsFullExt_r(accNodesModel,accYears,indicator,accNodesModel_aa,accYears_aa,indicator_aa,accNodesModel_a,accYears_a,indicator_a)),
            compoundIndicatorsFullExt_r(accNodesModel,accYears,indicator,accNodesModel_aa,accYears_aa,indicator_aa,accNodesModel_a,accYears_a,indicator_a));

    compoundIndicatorsFull(accNodesModel,accYears,indicator,accNodesModel_aa,accYears_aa,indicator_aa)
        $(compoundIndicatorsFull_act(accNodesModel,accYears,indicator)
            and sum((accNodesModel_a,accYears_a,indicator_a)
                    $compoundIndicatorsFullExt_r(accNodesModel,accYears,indicator,accNodesModel_aa,accYears_aa,indicator_aa,accNodesModel_a,accYears_a,indicator_a), 1) > 1)
        = sum((accNodesModel_a,accYears_a,indicator_a)
                $(compoundIndicatorsFullExt_r(accNodesModel,accYears,indicator,accNodesModel_aa,accYears_aa,indicator_aa,accNodesModel_a,accYears_a,indicator_a)
                    and not (sameas(accNodesModel, accNodesModel_a) and sameas(accYears, accYears_a) and sameas(indicator, indicator_a))),
            compoundIndicatorsFullExt_r(accNodesModel,accYears,indicator,accNodesModel_aa,accYears_aa,indicator_aa,accNodesModel_a,accYears_a,indicator_a))
        / sum((accNodesModel_a,accYears_a,indicator_a)
                $(compoundIndicatorsFullExt_r(accNodesModel,accYears,indicator,accNodesModel_aa,accYears_aa,indicator_aa,accNodesModel_a,accYears_a,indicator_a)
                    and not (sameas(accNodesModel, accNodesModel_a) and sameas(accYears, accYears_a) and sameas(indicator, indicator_a))), 1);

    option compoundIndicatorsFull_act < compoundIndicatorsFull;
    compIndicators_post = sum((accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a,accNodesModel_aa,accYears_aa,indicator_aa)
        $compoundIndicatorsFullExt(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a,accNodesModel_aa,accYears_aa,indicator_aa), 1);
);
option sparseval = 0;

* map to sub-years and sub-nodes
compoundIndicatorsFullExt(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a,accNodesModel_aa,accYears_aa,indicator_a)
    $(compoundIndicatorsFull_act(accNodesModel,accYears,indicator)
        and compoundIndicatorsFull(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a)
        and map_accYears(accYears_aa,accYears_a)
        and map_accNodes(accNodesModel_aa,accNodesModel_a))
    = compoundIndicatorsFull(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a)
        * yearFactor(accNodesModel,accYears,indicator,accYears_aa);

* Remove all accounting regions and accounting years
compoundIndicatorsFullExt(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a,accNodesModel_aa,accYears_aa,indicator_aa)
    $(compoundIndicatorsFullExt(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a,accNodesModel_aa,accYears_aa,indicator_aa)
        and not accounting_indicatorBounds(accNodesModel_aa,accYears_aa,indicator_aa,"isVariable")
        and sum(accNodes$(sameas(accNodes,accNodesModel_aa)), 1)
            or sameas("horizon",accYears_aa))
    = 0;

* Map from extended parameter to final parameter
option compoundIndicatorsFullExt_r < compoundIndicatorsFullExt;
compoundIndicatorsFull(accNodesModel,accYears,indicator,accNodesModel_aa,accYears_aa,indicator_aa)
    $(compoundIndicatorsFull_act(accNodesModel,accYears,indicator)
        and sum((accNodesModel_a,accYears_a,indicator_a)
                $compoundIndicatorsFullExt_r(accNodesModel,accYears,indicator,accNodesModel_aa,accYears_aa,indicator_aa,accNodesModel_a,accYears_a,indicator_a), 1) = 1)
    = sum((accNodesModel_a,accYears_a,indicator_a)
            $(compoundIndicatorsFullExt_r(accNodesModel,accYears,indicator,accNodesModel_aa,accYears_aa,indicator_aa,accNodesModel_a,accYears_a,indicator_a)),
        compoundIndicatorsFullExt_r(accNodesModel,accYears,indicator,accNodesModel_aa,accYears_aa,indicator_aa,accNodesModel_a,accYears_a,indicator_a));

compoundIndicatorsFull(accNodesModel,accYears,indicator,accNodesModel_aa,accYears_aa,indicator_aa)
    $(compoundIndicatorsFull_act(accNodesModel,accYears,indicator)
        and sum((accNodesModel_a,accYears_a,indicator_a)
                $compoundIndicatorsFullExt_r(accNodesModel,accYears,indicator,accNodesModel_aa,accYears_aa,indicator_aa,accNodesModel_a,accYears_a,indicator_a), 1) > 1)
    = sum((accNodesModel_a,accYears_a,indicator_a)
            $(compoundIndicatorsFullExt_r(accNodesModel,accYears,indicator,accNodesModel_aa,accYears_aa,indicator_aa,accNodesModel_a,accYears_a,indicator_a)
                and not (sameas(accNodesModel, accNodesModel_a) and sameas(accYears, accYears_a) and sameas(indicator, indicator_a))),
        compoundIndicatorsFullExt_r(accNodesModel,accYears,indicator,accNodesModel_aa,accYears_aa,indicator_aa,accNodesModel_a,accYears_a,indicator_a))
    / sum((accNodesModel_a,accYears_a,indicator_a)
            $(compoundIndicatorsFullExt_r(accNodesModel,accYears,indicator,accNodesModel_aa,accYears_aa,indicator_aa,accNodesModel_a,accYears_a,indicator_a)
                and not (sameas(accNodesModel, accNodesModel_a) and sameas(accYears, accYears_a) and sameas(indicator, indicator_a))), 1);

* Remove all accounting regions and accounting years
compoundIndicatorsFull(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a)
    $(compoundIndicatorsFull(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a)
        and not accounting_indicatorBounds(accNodesModel_a,accYears_a,indicator_a,"isVariable")
        and (sum(accNodes$(sameas(accNodes,accNodesModel_a)), 1)
            or sameas("horizon",accYears_a)))
    = 0;

* Remove all slack indicators except the ones declared in indicatorBounds
compoundIndicatorsFull(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a)
    $(compoundIndicatorsFull(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a)
        and sum((accNodesModel_aa,accYears_aa)$accounting_indicatorBounds(accNodesModel_aa,accYears_aa,indicator_a,"isVariable"), 1)
        and not accounting_indicatorBounds(accNodesModel_a,accYears_a,indicator_a,"isVariable"))
    = 0;

* Check if the objective value and optimization sense is set correctly
scalar opti_values, opti_sense;
opti_values = sum((accNodesModel,accYears,indicator)$(accounting_indicatorBounds(accNodesModel,accYears,indicator,"obj") <> 0), 1);
opti_sense = sum((accNodesModel,accYears,indicator), accounting_indicatorBounds(accNodesModel,accYears,indicator,"obj"));
abort$(opti_values < 1)
    "Accounting: No indicator specified as objective value"
abort$(opti_values > 1)
    "Accounting: Too many indicators specified as objective value"
abort$(opti_sense <> -1 and opti_sense <> 1)
    "Accounting: Optimization sense has to be either -1 or 1"
