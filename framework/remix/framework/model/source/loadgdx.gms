* SPDX-FileCopyrightText: Copyright (c) 2023-2024 German Aerospace Center (DLR)
* SPDX-License-Identifier: BSD-3-Clause

* ==== global settings ====
$if not set capsfromgdx     $setglobal capsfromgdx       None

$iftheni.cfgdx not %capsfromgdx%==None
$log "Loading capacities from file %capsfromgdx%"

$gdxin %capsfromgdx%
$load converter_units, storage_units, transfer_links

converter_units(%selscen%accNodesModel,accYears,converter_techs,vintage,capType)
    = ceil(converter_units(%selscen%accNodesModel,accYears,converter_techs,vintage,capType) * 1e3) / 1e3;
storage_units(%selscen%accNodesModel,accYears,storage_techs,vintage,capType)
    = ceil(storage_units(%selscen%accNodesModel,accYears,storage_techs,vintage,capType) * 1e3) / 1e3;
transfer_links(%selscen%linksModel,years,transfer_techs,vintage,capType)
    = ceil(transfer_links(%selscen%linksModel,years,transfer_techs,vintage,capType) * 1e3) / 1e3;
$endif.cfgdx