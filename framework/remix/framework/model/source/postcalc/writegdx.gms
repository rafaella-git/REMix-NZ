* SPDX-FileCopyrightText: Copyright (c) 2023-2024 German Aerospace Center (DLR)
* SPDX-License-Identifier: BSD-3-Clause

$onVerbatim
$ifthene.run_postcalc %run_postcalc%

execute_unload "%resultdir%%system.dirsep%%resultfile%.gdx"
$iftheni.gdx_metadata %metadata% == 0
$elseIf.gdx_metadata a == a
    metadata
$endif.gdx_metadata
    timeModel
    timeModelToCalc
    nodesModel
    linksModel
    indicator
    commodity
    techs
    accNodesModel
    accYears

    map_nodesModel
    map_linksModel
    map_nodesAccounting

    indicator_accounting
    indicator_accounting_ref
    indicator_accounting_comp
    indicator_accounting_detailed
    indicator_accounting_links

    converter_caps
    converter_units
$iftheni.gdx_mipconverter %gdx_mipconverter%==1
    converter_unitsOnline
    converter_unitsUsingActivity_MIP
    converter_activity
    converter_unitStartups
$endif.gdx_mipconverter

    transfer_caps
    transfer_links
    transfer_flows
    transfer_flows_annual
    transfer_losses
    transfer_losses_annual

    storage_caps
    storage_units
    storage_flows
    storage_level_out
    storage_flows_annual
    storage_losses_out
    storage_losses_annual
$iftheni.gdx_mipstorage %gdx_mipstorage%==1
    storage_unitsStateTracker
    storage_unitsStateTrackerDecom
    storage_levelPerAge
    storage_chargePerAge
    storage_unitsSoC
    storage_unitsBuild
    storage_unitsDecom
$endif.gdx_mipstorage

    commodity_balance
    commodity_balance_annual

    marginals_balance
    marginals_sourcesink_profile
    marginals_sourcesink_sum
    marginals_indicator_bounds

$iftheni.r2a %gdx_r2a%==1
    r2a_annuity_cost_converter
    r2a_spec_cost_converter
    r2a_spec_cost_fuel
    r2a_spec_cost_indicator
    r2a_converter_efficiencies
    r2a_converter_avail_factor
    r2a_converter_avail_profile
    r2a_storage_selfdischarge
    r2a_storage_e2p
$endif.r2a
    diagnostics
    ;
$endif.run_postcalc
$offVerbatim