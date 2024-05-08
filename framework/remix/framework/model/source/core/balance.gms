* SPDX-FileCopyrightText: Copyright (c) 2023-2024 German Aerospace Center (DLR)
* SPDX-License-Identifier: BSD-3-Clause

* // # balance
* // The equations in this file describe the commodity balancing in the model.

Equations
  Eq_balance_commodities(timeModel,nodesModel,years,commodity
    ) "Balance for each commodity used in each time step, region, and year";


* ==== calculation of mappings ====

set balance_techComm(techs,commodity);
balance_techComm(converter_techs(techs),commodity)
    $sum((vintage,activity)$converter_coefficient(converter_techs,vintage,activity,commodity,"coefficient"), 1) = yes;
balance_techComm(storage_techs(techs),commodity)
    $sum((vintage)$storage_sizeParam(storage_techs,vintage,commodity,"size"), 1) = yes;
balance_techComm(transfer_techs(techs),commodity)
    $sum((vintage)$transfer_coefficient(transfer_techs,vintage,commodity,"coefficient"), 1) = yes;
balance_techComm(sourcesink_techs(techs),commodity)
    $sum((nodesModelToCalc,yearsToCalc)$sourcesink_enabled(nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity), 1) = yes;

set balance_usedConverter(nodesModel,years,commodity);
set balance_usedStorage(nodesModel,years,commodity);
set balance_usedTransfer(nodesModel,years,commodity);
set balance_usedSourceSink(nodesModel,years,commodity);
set balance_usedBalance(nodesModel,years,commodity);

balance_usedConverter(nodesModel,years,commodity)
    = sum ((converter_techs,vintage,activity)
            $( converter_coefficient(converter_techs,vintage,activity,commodity,"coefficient") <> 0
                and converter_usedTech(nodesModel,years,converter_techs,vintage) ), 1);

balance_usedStorage(nodesModel,years,commodity)
    = sum ((storage_techs,vintage,activity)
            $(storage_sizeParam(storage_techs,vintage,commodity,"size") <> 0
                and storage_usedTech(nodesModel,years,storage_techs,vintage) ), 1);

balance_usedTransfer(nodesModel,years,commodity)
    = sum ((linksModel,transfer_techs,vintage)
            $(transfer_coefficient(transfer_techs,vintage,commodity,"coefficient") <> 0
                and transfer_incidenceModel(nodesModel,linksModel) <> 0), 1);

option balance_usedSourceSink < sourcesink_enabled;

balance_usedBalance(nodesModel,years,commodity)
    $(balance_usedConverter(nodesModel,years,commodity)
        or balance_usedStorage(nodesModel,years,commodity)
        or balance_usedTransfer(nodesModel,years,commodity)
        or balance_usedSourceSink(nodesModel,years,commodity) )
    = yes;


* ==== equation definition ====

* // ### Balance commodities
* // Balancing of commodities for all model regions, time steps, and years
* // {Eq_balance_commodities}
Eq_balance_commodities(timeModelSel(timeModelToCalc),nodesModelSel,yearsSel,commodity)
    $balance_usedBalance(nodesModelSel,yearsSel,commodity)
    ..
* converter
    sum((converter_techs,vintage,activity)
            $( converter_coefficientProfile(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,activity,commodity) <> 0
                AND converter_usedTech(nodesModelSel,yearsSel,converter_techs,vintage) ),
        converter_activity(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,activity)
            * converter_coefficientProfile(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,activity,commodity)
            * timeLength(timeModelSel)
        + converter_unitsUsingActivity_MIP(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,activity)
            * converter_coefficient(converter_techs,vintage,activity,commodity,"constant")
            * timeLength(timeModelSel))

* storages
    + sum((storage_techs,vintage)
            $( storage_sizeParam(storage_techs,vintage,commodity,"size") <> 0
                and storage_usedTech(nodesModelSel,yearsSel,storage_techs,vintage) ),
        storage_level(timeModelToCalc--1,nodesModelSel,yearsSel,storage_techs,vintage,commodity)
        - storage_level(timeModelSel,nodesModelSel,yearsSel,storage_techs,vintage,commodity)
$iftheni.pips %method%==pips
        + (storage_level(timeModelToCalc--1,nodesModelSel,yearsSel,storage_techs,vintage,commodity)
            * storage_sizeParam(storage_techs,vintage,commodity,"selfdischarge")
        - storage_sizeParam(storage_techs,vintage,commodity,"selfdischargeAbs"))
        * timeLength(timeModelToCalc)
$else.pips
        - storage_losses(timeModelSel,nodesModelSel,yearsSel,storage_techs,vintage,commodity)
        * timeLength(timeModelSel)
$endif.pips
        )

* transfer
    + sum((linksModel,transfer_techs,vintage)
            $(transfer_usedTech(linksModel,yearsSel,transfer_techs,vintage)
                and linksModelToCalc(linksModel)),
        (   transfer_flowAlong(timeModelSel,linksModel,yearsSel,transfer_techs,vintage)
                $(transfer_incidenceModel(nodesModelSel,linksModel) > 0)
          + transfer_flowAgainst(timeModelSel,linksModel,yearsSel,transfer_techs,vintage)
                $(transfer_incidenceModel(nodesModelSel,linksModel) < 0) )
        * transfer_coefficient(transfer_techs,vintage,commodity,"coefficient")
        * timeLength(timeModelSel))

    - sum((linksModel,transfer_techs,vintage)
            $(transfer_usedTech(linksModel,yearsSel,transfer_techs,vintage)
                and linksModelToCalc(linksModel)),
        (   transfer_flowAlong(timeModelSel,linksModel,yearsSel,transfer_techs,vintage)
                $(transfer_incidenceModel(nodesModelSel,linksModel) < 0)
          + transfer_flowAgainst(timeModelSel,linksModel,yearsSel,transfer_techs,vintage)
                $(transfer_incidenceModel(nodesModelSel,linksModel) > 0) )
        * transfer_coefficient(transfer_techs,vintage,commodity,"coefficient")
        * timeLength(timeModelSel))

    + 0.5 * sum((linksModel,transfer_techs,vintage)
            $(transfer_usedTech(linksModel,yearsSel,transfer_techs,vintage)
                AND linksModelToCalc(linksModel)),
        (   transfer_flowAlong(timeModelSel,linksModel,yearsSel,transfer_techs,vintage)
                $(transfer_incidenceModel(nodesModelSel,linksModel) <> 0)
          + transfer_flowAgainst(timeModelSel,linksModel,yearsSel,transfer_techs,vintage)
                $(transfer_incidenceModel(nodesModelSel,linksModel) <> 0) )
        * timeLength(timeModelSel)
        * ( transfer_coefPerFlow(transfer_techs,vintage,commodity,"coefPerFlow")
            + sum(link_types,
                transfer_coefPerLength(transfer_techs,vintage,commodity,link_types,"coefPerLength")
                * transfer_lengthParam(linksModel,link_types,"length"))))

* sourcesink
    + sum((sourcesink_techs)
            $sourcesink_enabled(nodesModelSel,yearsSel,sourcesink_techs,commodity),
        sourcesink_flow(timeModelSel,nodesModelSel,yearsSel,sourcesink_techs,commodity)
        * timeLength(timeModelSel))
    =e=
    0;


* ==== model definition ====

Model M_balance
/
  Eq_balance_commodities
/;
