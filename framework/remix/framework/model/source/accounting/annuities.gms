* SPDX-FileCopyrightText: Copyright (c) 2023-2024 German Aerospace Center (DLR)
* SPDX-License-Identifier: BSD-3-Clause

* Calculation of annuities

abort$(sum((indicator,nodesModelToCalc,converter_techs,vintage)
        $(accounting_converterUnits(indicator,nodesModelToCalc,converter_techs,vintage,"useAnnuity") = 1
        and accounting_converterUnits(indicator,nodesModelToCalc,converter_techs,vintage,"amorTime") < 1), 1) > 0 )
    "Accounting: Some converter technologies use annuities but have no amortization time"

parameter accounting_annuityFactor_converter(indicator,nodesModel,converter_techs,vintage);
accounting_annuityFactor_converter(indicator,nodesModelToCalc,converter_techs,vintage)
    $accounting_converterUnits(indicator,nodesModelToCalc,converter_techs,vintage,"useAnnuity")
    = 
    accounting_converterUnits(indicator,nodesModelToCalc,converter_techs,vintage,"interest")
        * (1 + accounting_converterUnits(indicator,nodesModelToCalc,converter_techs,vintage,"interest"))
        ** accounting_converterUnits(indicator,nodesModelToCalc,converter_techs,vintage,"amorTime")
    / ((1 + accounting_converterUnits(indicator,nodesModelToCalc,converter_techs,vintage,"interest"))
        ** accounting_converterUnits(indicator,nodesModelToCalc,converter_techs,vintage,"amorTime") - 1);


abort$(sum((indicator,nodesModelToCalc,storage_techs,vintage)
        $(accounting_storageUnits(indicator,nodesModelToCalc,storage_techs,vintage,"useAnnuity") = 1
        and accounting_storageUnits(indicator,nodesModelToCalc,storage_techs,vintage,"amorTime") < 1), 1) > 0 )
    "Accounting: Some storage technologies use annuities but have no amortization time"

parameter accounting_annuityFactor_storage(indicator,nodesModel,storage_techs,vintage);
accounting_annuityFactor_storage(indicator,nodesModelToCalc,storage_techs,vintage)
    $accounting_storageUnits(indicator,nodesModelToCalc,storage_techs,vintage,"useAnnuity")
    = 
    accounting_storageUnits(indicator,nodesModelToCalc,storage_techs,vintage,"interest")
        * (1 + accounting_storageUnits(indicator,nodesModelToCalc,storage_techs,vintage,"interest"))
        ** accounting_storageUnits(indicator,nodesModelToCalc,storage_techs,vintage,"amorTime")
    / ((1 + accounting_storageUnits(indicator,nodesModelToCalc,storage_techs,vintage,"interest"))
        ** accounting_storageUnits(indicator,nodesModelToCalc,storage_techs,vintage,"amorTime") - 1);


abort$(sum((indicator,linksModelToCalc,transfer_techs,vintage)
        $(accounting_transferLinks(indicator,linksModelToCalc,transfer_techs,vintage,"useAnnuity") = 1
        and accounting_transferLinks(indicator,linksModelToCalc,transfer_techs,vintage,"amorTime") < 1), 1) > 0 )
    "Accounting: Some transfer technologies use annuities but have no amortization time"

parameter accounting_annuityFactor_transferLink(indicator,linksModel,transfer_techs,vintage);
accounting_annuityFactor_transferLink(indicator,linksModelToCalc,transfer_techs,vintage)
    $accounting_transferLinks(indicator,linksModelToCalc,transfer_techs,vintage,"useAnnuity")
    = 
    accounting_transferLinks(indicator,linksModelToCalc,transfer_techs,vintage,"interest")
        * (1 + accounting_transferLinks(indicator,linksModelToCalc,transfer_techs,vintage,"interest"))
        ** accounting_transferLinks(indicator,linksModelToCalc,transfer_techs,vintage,"amorTime")
    / ((1 + accounting_transferLinks(indicator,linksModelToCalc,transfer_techs,vintage,"interest"))
        ** accounting_transferLinks(indicator,linksModelToCalc,transfer_techs,vintage,"amorTime") - 1);



abort$(sum((indicator,linksModelToCalc,transfer_techs,vintage,link_types)
        $(accounting_transferPerLength(indicator,linksModelToCalc,transfer_techs,vintage,link_types,"useAnnuity") = 1
        and accounting_transferPerLength(indicator,linksModelToCalc,transfer_techs,vintage,link_types,"amorTime") < 1), 1) > 0 )
    "Accounting: Some transfer-per-length technologies use annuities but have no amortization time"

parameter accounting_annuityFactor_transferPerLength(indicator,linksModel,transfer_techs,vintage,link_types);
accounting_annuityFactor_transferPerLength(indicator,linksModelToCalc,transfer_techs,vintage,link_types) = 1;
    
accounting_annuityFactor_transferPerLength(indicator,linksModelToCalc,transfer_techs,vintage,link_types)
    $accounting_transferPerLength(indicator,linksModelToCalc,transfer_techs,vintage,link_types,"useAnnuity")
    = 
    accounting_transferPerLength(indicator,linksModelToCalc,transfer_techs,vintage,link_types,"interest")
        * (1 + accounting_transferPerLength(indicator,linksModelToCalc,transfer_techs,vintage,link_types,"interest"))
        ** accounting_transferPerLength(indicator,linksModelToCalc,transfer_techs,vintage,link_types,"amorTime")
    / ((1 + accounting_transferPerLength(indicator,linksModelToCalc,transfer_techs,vintage,link_types,"interest"))
        ** accounting_transferPerLength(indicator,linksModelToCalc,transfer_techs,vintage,link_types,"amorTime") - 1);

$ifthene.roundcoefs %roundcoefs%=1
accounting_annuityFactor_converter(indicator,nodesModelToCalc,converter_techs,vintage)
    = round(accounting_annuityFactor_converter(indicator,nodesModelToCalc,converter_techs,vintage), 5);
accounting_annuityFactor_storage(indicator,nodesModelToCalc,storage_techs,vintage)
    = round(accounting_annuityFactor_storage(indicator,nodesModelToCalc,storage_techs,vintage), 5);
accounting_annuityFactor_transferLink(indicator,linksModelToCalc,transfer_techs,vintage)
    = round(accounting_annuityFactor_transferLink(indicator,linksModelToCalc,transfer_techs,vintage), 5);
accounting_annuityFactor_transferPerLength(indicator,linksModelToCalc,transfer_techs,vintage,link_types)
    = round(accounting_annuityFactor_transferPerLength(indicator,linksModelToCalc,transfer_techs,vintage,link_types), 5);
$endif.roundcoefs
