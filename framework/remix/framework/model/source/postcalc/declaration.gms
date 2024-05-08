* SPDX-FileCopyrightText: Copyright (c) 2023-2024 German Aerospace Center (DLR)
* SPDX-License-Identifier: BSD-3-Clause
* // # Output data
* // This is the reference data model of the REMix output.
* //
* // ## Standard output
$onVerbatim
$ifthene.run_postcalc %run_postcalc%

$if not set resultdir        $setglobal resultdir               %gams.workdir%%system.dirsep%results
$if not set resultfile       $setglobal resultfile              remix
$if not set gdx_roundts      $setglobal gdx_roundts             1
$if not set gdx_mipconverter $setglobal gdx_mipconverter        0
$if not set gdx_mipstorage   $setglobal gdx_mipstorage          0
$if not set gdx_r2a          $setglobal gdx_r2a                 0

$if not dexist "%resultdir%" put_utility 'exec' / 'mkdir -p %resultdir%'

set capType / "build", "decom", "total", "lowerLimit", "upperLimit", "total_degraded" /;
set balanceType / "netto", "brutto", "positive", "negative", "flh" /;
set profileType / "upper", "fixed", "lower" /;
set r2a_has_converter_cost(indicator,nodesModel,years,techs,vintage,commodity);

** // OUTPUT: indicator_accounting | OEO_00000350:quantity value
* // ### indicator_accounting
* // Title: Accounting indicators Post-calculation
parameter indicator_accounting(%scenidx%accNodesModel,accYears,indicator) "Post-calculation aggregated indicator accounting" ;
* //
** // OUTPUT: indicator_accounting_ref | OEO_00000350:quantity value
* // ### indicator_accounting_ref
* // Title: Accounting indicators reference
parameter indicator_accounting_ref(%scenidx%accNodesModel,accYears,indicator) "Indicator accounting optimization levels" ;
* //
** // OUTPUT: indicator_accounting_comp | OEO_00000350:quantity value
* // ### indicator_accounting_comp
* // Title: Accounting indicators composition
parameter indicator_accounting_comp(%scenidx%accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a) "Indicator accounting with other indicator composition" ;
* //
** // OUTPUT: indicator_accounting_detailed | OEO_00000350:quantity value
* // ### indicator_accounting_detailed
* // Title: Accounting indicators detail
parameter indicator_accounting_detailed(%scenidx%indicator,nodesModel,years,techs) "Indicator accounting with technological composition" ;
* //
** // OUTPUT: indicator_accounting_links | OEO_00000350:quantity value
* // ### indicator_accounting_links
* // Title: Accounting indicators per transfer link
parameter indicator_accounting_links(%scenidx%indicator,nodesModel,nodesModel_a,linksModel,years,techs) "Indicator accounting of transfer technologies" ;
* //
parameter converter_ratedInput(techs,vintage,commodity);
parameter converter_ratedOutput(techs,vintage,commodity);
parameter converter_ratedOutput_min(nodesModel,years,techs,commodity);
parameter converter_ratedOutput_max(nodesModel,years,techs,commodity);

** // OUTPUT: converter_caps | OEO_00010257:power capacity
* // ### converter_caps
* // Title: Coverter capacities post-calculation
parameter converter_caps(%scenidx%accNodesModel,accYears,techs,commodity,capType) "Converter capacities, aggregated total" ;
* //
parameter converter_caps_ext(%scenidx%accNodesModel,nodesModel,accYears,years,techs,vintage,commodity,capType);

** // OUTPUT: converter_units | OEO_00000350:quantity value
* // ### converter_units
* // Title: Coverter units post-calculation
parameter converter_units(%scenidx%accNodesModel,accYears,techs,vintage,capType) "Converter units" ;
* //
parameter converter_units_ext(%scenidx%accNodesModel,nodesModel,accYears,years,techs,vintage,capType);

set transfer_usedStartEnd(nodesModel,nodesModel,linksModel,years,techs);
parameter max_transfer_coefficient(techs,commodity);
parameter min_transfer_coefficient(techs,commodity);
** // OUTPUT: transfer_links | OEO_00000350:quantity value
* // ### transfer_links
* // Title: Transfer post-calculation
parameter transfer_links(%scenidx%linksModel,years,techs,vintage,capType) "Transfer links" ;
* //
** // OUTPUT: transfer_caps | OEO_00010257:power capacity
* // ### transfer_caps
* // Title: Transfer capacities post-calculation
parameter transfer_caps(%scenidx%nodesModel_start,nodesModel_end,linksModel,years,techs,commodity,capType) "Transfer capacities, aggregated total" ;
* //
** // OUTPUT: transfer_flows | OEO_00050019:energy amount value
* // ### transfer_flows
* // Title: Transfer hourly flows
parameter transfer_flows(%scenidx%timeModel,nodesModel_start,nodesModel_end,linksModel,years,techs,commodity) "Hourly link flows" ;
* //
** // OUTPUT: transfer_flows_annual | OEO_00050019:energy amount value
* // ### transfer_flows_annual
* // Title: Transfer annual flows
parameter transfer_flows_annual(%scenidx%nodesModel_start,nodesModel_end,linksModel,years,techs,commodity,balanceType) "Annual link flows" ;
* //
** // OUTPUT: transfer_losses | OEO_00050019:energy amount value
* // ### transfer_losses
* // Title: Transfer hourly losses
parameter transfer_losses(%scenidx%timeModel,nodesModel_start,nodesModel_end,linksModel,years,techs,commodity) "Hourly link losses" ;
* //
** // OUTPUT: transfer_losses_annual | OEO_00050019:energy amount value
* // ### transfer_losses_annual
* // Title: Transfer annual losses
parameter transfer_losses_annual(%scenidx%nodesModel_start,nodesModel_end,linksModel,years,techs,commodity,balanceType) "Annual link losses" ;
* //
parameter storage_size_max(nodesModel,years,techs,commodity);
parameter storage_size_min(nodesModel,years,techs,commodity);
** // OUTPUT: storage_units |  OEO_00000350:quantity value
* // ### storage_units
* // Title: Storage units post-calculation
parameter storage_units(%scenidx%accNodesModel,accYears,techs,vintage,capType) "Storage units" ;
* //
parameter storage_units_ext(%scenidx%accNodesModel,nodesModel,accYears,years,techs,vintage,capType);
** // OUTPUT: storage_caps | OEO_00230000:energy storage capacity
* // ### storage_caps
* // Title: Storage capacities post-calculation
parameter storage_caps(%scenidx%accNodesModel,accYears,techs,commodity,capType) "Storage capacities, aggregated total" ;
* //
parameter storage_caps_ext(%scenidx%accNodesModel,nodesModel,accYears,years,techs,vintage,commodity,capType);
** // OUTPUT: storage_level_out | OEO_00330012:energy storage content
* // ### storage_level_out
* // Title: Storage level
parameter storage_level_out(%scenidx%timeModel,accNodesModel,accYears,techs,commodity) "Storage level per time step" ;
* //
parameter storage_level_out_ext(%scenidx%timeModel,accNodesModel,nodesModel,accYears,years,techs,vintage,commodity);
** // OUTPUT: storage_flows | OEO_00050019:energy amount value
* // ### storage_flows
* // Title: Storage hourly flows
parameter storage_flows(%scenidx%timeModel,accNodesModel,accYears,techs,commodity) "Storage flows per time step" ;
* //
parameter storage_flows_ext(%scenidx%timeModel,accNodesModel,nodesModel,accYears,years,techs,vintage,commodity);
** // OUTPUT: storage_flows_annual | OEO_00050019:energy amount value
* // ### storage_flows_annual
* // Title: Storage annual flows
parameter storage_flows_annual(%scenidx%accNodesModel,accYears,techs,commodity,balanceType) "Aggregated storage flows" ;
* //
** // OUTPUT: storage_losses_out | OEO_00050019:energy amount value
* // ### storage_losses_out
* // Title: Storage hourly losses
parameter storage_losses_out(%scenidx%timeModel,accNodesModel,accYears,techs,commodity) "Storage losses per time step" ;
* //
** // OUTPUT: storage_losses_annual | OEO_00050019:energy amount value
* // ### storage_losses_annual
* // Title: Storage annual losses
parameter storage_losses_annual(%scenidx%accNodesModel,accYears,techs,commodity,balanceType) "Aggregated storage losses" ;
* //
** // OUTPUT: commodity_balance | OEO_00000350:quantity value
* // ### commodity_balance
* // Title: Hourly commodity balance
parameter commodity_balance(%scenidx%timeModel,accNodesModel,accYears,techs,commodity) "Commodity balances per model hour" ;
parameter commodity_balance_ext(%scenidx%timeModel,accNodesModel,nodesModel,accYears,years,techs,vintage,commodity);
* //
** // OUTPUT: commodity_balance_annual | OEO_00000350:quantity value
* // ### commodity_balance_annual
* // Title: Annual commodity balance
parameter commodity_balance_annual(%scenidx%accNodesModel,accYears,techs,commodity,balanceType) "Aggregated commodity balances" ;
* //
** // OUTPUT: marginals_sourcesink_profile | OEO_00040008:marginal cost
* // ### marginals_sourcesink_profile
* // Title: Source-sink marginals
parameter marginals_sourcesink_profile(%scenidx%timeModel,nodesModel,years,techs,commodity) "Source and sink flow marginal values" ;
* //
** // OUTPUT: marginals_balance | OEO_00040008:marginal cost
* // ### marginals_balance
* // Title: Nodel balance marginals
parameter marginals_balance(%scenidx%timeModel,nodesModel,years,commodity) "Nodal balance marginal values" ;
* //
** // OUTPUT: marginals_sourcesink_sum | OEO_00040008:marginal cost
* // ### marginals_sourcesink_sum
* // Title: Source-sink annual marginals
parameter marginals_sourcesink_sum(%scenidx%nodesModel,years,techs,commodity) "Annual source and sink flow marginal values" ;
* //
** // OUTPUT: marginals_indicator_bounds | OEO_00040008:marginal cost
* // ### marginals_indicator_bounds
* // Title: Indicator bounds marginals
parameter marginals_indicator_bounds(%scenidx%accNodesModel,accYears,indicator) "Accounting indicator marginals" ;
* //
* // ## REMix-AMIRIS interface output
* //
** // OUTPUT: r2a_annuity_cost_converter | OEO_00040009:cost
* // ### r2a_annuity_cost_converter
* // Title: Converter annuity cost
parameter r2a_annuity_cost_converter(indicator,accNodesModel,accYears,techs,vintage,commodity) "Annuity cost of converters" ;
* //
** // OUTPUT: r2a_spec_cost_converter | OEO_00040009:cost
* // ### r2a_spec_cost_converter
* // Title: Converter specific cost
parameter r2a_spec_cost_converter(indicator,accNodesModel,accYears,techs,vintage,commodity) "Specific cost of converters" ;
* //
** // OUTPUT: r2a_spec_cost_fuel | OEO_00040009:cost
* // ### r2a_spec_cost_fuel
* // Title: Commodity specific cost
parameter r2a_spec_cost_fuel(indicator,accNodesModel,accYears,techs,commodity) "Specific cost of commodities, used particularly for fuels" ;
* //
** // OUTPUT: r2a_spec_cost_indicator | OEO_00040009:cost
* // ### r2a_spec_cost_indicator
* // Title: Indicator specific cost
parameter r2a_spec_cost_indicator(indicator,indicator_a,accNodesModel,accYears) "Specific indicator cost" ;
* //
** // OUTPUT: r2a_converter_efficiencies | OEO_00140050:efficiency value
* // ### r2a_converter_efficiencies
* // Title: Indicator specific cost
parameter r2a_converter_efficiencies(techs,vintage,activity,commodity,commodity_a) "Converter-specific efficiency" ;
* //
** // OUTPUT: r2a_converter_avail_factor | OEO_00000350:quantity value
* // ### r2a_converter_avail_factor
* // Title: Converter availability factor
parameter r2a_converter_avail_factor(accNodesModel,accYears,techs,vintage) "Availability factor of converter technologies" ;
parameter r2a_converter_avail_profile(timeModel,accNodesModel,accYears,techs,commodity,profileType) "Availability of converter technologies" ;
parameter r2a_converter_avail_profile_ext(timeModel,accNodesModel,nodesModel,accYears,years,techs,vintage,commodity,profileType);
* //
** // OUTPUT: r2a_storage_e2p | OEO_00000350:quantity value
* // ### r2a_storage_e2p
* // Title: E2P Storage
parameter r2a_storage_e2p(%scenidx%accNodesModel,accYears,techs,vintage,commodity,capType) "Energy-to-power ratio" ;
* //
** // OUTPUT: r2a_storage_selfdischarge | OEO_00000350:quantity value
* // ### r2a_storage_selfdischarge
* // Title: Storage self discharge rates
parameter r2a_storage_selfdischarge(accNodesModel,accYears,techs,vintage,commodity);

$endif.run_postcalc
$offVerbatim