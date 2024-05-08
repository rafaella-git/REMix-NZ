* SPDX-FileCopyrightText: Copyright (c) 2023-2024 German Aerospace Center (DLR)
* SPDX-License-Identifier: BSD-3-Clause

$ifthene.run_postcalc %run_postcalc%

* ==== indicator accounting ====

indicator_accounting_detailed(%selscen%indicator,nodesModelToCalc,yearsToCalc,techs)
    =
* == converters ==
    sum ((converter_techs(techs),vintage)
                $(converter_availTech(nodesModelToCalc,yearsToCalc,converter_techs,vintage)
                    and accounting_converterUnits(indicator,nodesModelToCalc,converter_techs,vintage,"useAnnuity") = 0),
        converter_unitsBuild.l(nodesModelToCalc,yearsToCalc,converter_techs,vintage)
        * accounting_converterUnits(indicator,nodesModelToCalc,converter_techs,vintage,"perUnitBuild"))

    + sum ((years_a,converter_techs(techs),vintage)
                $(converter_availTech(nodesModelToCalc,years_a,converter_techs,vintage)
                    and accounting_converterUnits(indicator,nodesModelToCalc,converter_techs,vintage,"useAnnuity") = 1
                    and years_a.val + accounting_converterUnits(indicator,nodesModelToCalc,converter_techs,vintage,"amorTime") > yearsToCalc.val
                    and years_a.val <= yearsToCalc.val ),
        converter_unitsBuild.l(nodesModelToCalc,years_a,converter_techs,vintage)
        * accounting_converterUnits(indicator,nodesModelToCalc,converter_techs,vintage,"perUnitBuild")
        * accounting_annuityFactor_converter(indicator,nodesModelToCalc,converter_techs,vintage) )

    + sum ((converter_techs(techs),vintage)
                $converter_usedTech(nodesModelToCalc,yearsToCalc,converter_techs,vintage),
        converter_unitsDecom.l(nodesModelToCalc,yearsToCalc,converter_techs,vintage)
        * accounting_converterUnits(indicator,nodesModelToCalc,converter_techs,vintage,"perUnitDecom")

        + converter_unitsTotal.l(nodesModelToCalc,yearsToCalc,converter_techs,vintage)
        * accounting_converterUnits(indicator,nodesModelToCalc,converter_techs,vintage,"perUnitTotal") )

    + sum ((timeModelToCalc,converter_techs(techs),vintage,activity)
                $converter_usedTechAct(nodesModelToCalc,yearsToCalc,converter_techs,vintage,activity),
        converter_activity.l(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage,activity)
        * timeLength(timeModelToCalc)
        * accounting_converterActivity(indicator,nodesModelToCalc,converter_techs,vintage,activity,"perActivity") )

    + sum ((timeModelToCalc,converter_techs(techs),vintage)
                $converter_usedTech(nodesModelToCalc,yearsToCalc,converter_techs,vintage),
        converter_unitStartups.l(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage)
        * accounting_converterStartup(indicator,nodesModelToCalc,converter_techs,vintage,"perStartup") )

    + sum ((timeModelToCalc,converter_techs(techs),vintage)
                $converter_usedTech(nodesModelToCalc,yearsToCalc,converter_techs,vintage),
        converter_rampPos.l(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage)
        * (accounting_converterStartup(indicator,nodesModelToCalc,converter_techs,vintage,"perRamp")
            + accounting_converterStartup(indicator,nodesModelToCalc,converter_techs,vintage,"perRampPos"))

        + converter_rampNeg.l(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage)
        * (accounting_converterStartup(indicator,nodesModelToCalc,converter_techs,vintage,"perRamp")
            + accounting_converterStartup(indicator,nodesModelToCalc,converter_techs,vintage,"perRampNeg")))

* == storage ==
    + sum ((storage_techs(techs),vintage)
                $(storage_availTech(nodesModelToCalc,yearsToCalc,storage_techs,vintage)
                    and accounting_storageUnits(indicator,nodesModelToCalc,storage_techs,vintage,"useAnnuity") = 0),
        storage_unitsBuild.l(nodesModelToCalc,yearsToCalc,storage_techs,vintage)
        * accounting_storageUnits(indicator,nodesModelToCalc,storage_techs,vintage,"perUnitBuild") )

    + sum ((years_a,storage_techs(techs),vintage)
                $(storage_availTech(nodesModelToCalc,yearsToCalc,storage_techs,vintage)
                    and accounting_storageUnits(indicator,nodesModelToCalc,storage_techs,vintage,"useAnnuity") = 1
                    and years_a.val + accounting_storageUnits(indicator,nodesModelToCalc,storage_techs,vintage,"amorTime") > yearsToCalc.val
                    and years_a.val <= yearsToCalc.val ),
        storage_unitsBuild.l(nodesModelToCalc,years_a,storage_techs,vintage)
        * accounting_storageUnits(indicator,nodesModelToCalc,storage_techs,vintage,"perUnitBuild")
        * accounting_annuityFactor_storage(indicator,nodesModelToCalc,storage_techs,vintage) )

    + sum ((storage_techs(techs),vintage)
                $storage_usedTech(nodesModelToCalc,yearsToCalc,storage_techs,vintage),
        storage_unitsDecom.l(nodesModelToCalc,yearsToCalc,storage_techs,vintage)
        * accounting_storageUnits(indicator,nodesModelToCalc,storage_techs,vintage,"perUnitDecom")

        + storage_unitsTotal.l(nodesModelToCalc,yearsToCalc,storage_techs,vintage)
        * accounting_storageUnits(indicator,nodesModelToCalc,storage_techs,vintage,"perUnitTotal") )

* == transfer ==
    + sum ((linksModel,transfer_techs(techs),vintage)
                $(transfer_availTech(linksModel,yearsToCalc,transfer_techs,vintage)
                    and linksModelToCalc(linksModel)
                    and transfer_incidenceModel(nodesModelToCalc,linksModel) <> 0
                    and accounting_transferLinks(indicator,linksModel,transfer_techs,vintage,"useAnnuity") = 0),
        0.5
        * transfer_linksBuild.l(linksModel,yearsToCalc,transfer_techs,vintage)
        * accounting_transferLinks(indicator,linksModel,transfer_techs,vintage,"perLinkBuild") )

    + sum ((linksModel,years_a,transfer_techs(techs),vintage)
                $(transfer_availTech(linksModel,yearsToCalc,transfer_techs,vintage)
                    and linksModelToCalc(linksModel)
                    and transfer_incidenceModel(nodesModelToCalc,linksModel) <> 0
                    and accounting_transferLinks(indicator,linksModel,transfer_techs,vintage,"useAnnuity") = 1
                    and years_a.val + accounting_transferLinks(indicator,linksModel,transfer_techs,vintage,"amorTime") > yearsToCalc.val
                    and years_a.val <= yearsToCalc.val ),
        0.5
        * transfer_linksBuild.l(linksModel,years_a,transfer_techs,vintage)
        * accounting_transferLinks(indicator,linksModel,transfer_techs,vintage,"perLinkBuild")
        * accounting_annuityFactor_transferLink(indicator,linksModel,transfer_techs,vintage) )

    + sum ((linksModel,transfer_techs(techs),vintage,link_types)
                $(transfer_availTech(linksModel,yearsToCalc,transfer_techs,vintage)
                    and linksModelToCalc(linksModel)
                    and transfer_incidenceModel(nodesModelToCalc,linksModel) <> 0
                    and accounting_transferPerLength(indicator,linksModel,transfer_techs,vintage,link_types,"useAnnuity") = 0 ),
        0.5
        * transfer_linksBuild.l(linksModel,yearsToCalc,transfer_techs,vintage)
        * transfer_lengthParam(linksModel,link_types,"length")
        * accounting_transferPerLength(indicator,linksModel,transfer_techs,vintage,link_types,"perLengthBuild") )

    + sum ((linksModel,years_a,transfer_techs(techs),vintage,link_types)
                $(transfer_availTech(linksModel,yearsToCalc,transfer_techs,vintage)
                    and linksModelToCalc(linksModel)
                    and transfer_incidenceModel(nodesModelToCalc,linksModel) <> 0
                    and accounting_transferPerLength(indicator,linksModel,transfer_techs,vintage,link_types,"useAnnuity") = 1
                    and years_a.val + accounting_transferPerLength(indicator,linksModel,transfer_techs,vintage,link_types,"amorTime") > yearsToCalc.val
                    and years_a.val <= yearsToCalc.val ),
        0.5
        * transfer_linksBuild.l(linksModel,years_a,transfer_techs,vintage)
        * transfer_lengthParam(linksModel,link_types,"length")
        * accounting_transferPerLength(indicator,linksModel,transfer_techs,vintage,link_types,"perLengthBuild")
        * accounting_annuityFactor_transferPerLength(indicator,linksModel,transfer_techs,vintage,link_types) )

    + sum ((linksModel,transfer_techs(techs),vintage)
                $(transfer_usedTech(linksModel,yearsToCalc,transfer_techs,vintage)
                    and linksModelToCalc(linksModel)
                    and transfer_incidenceModel(nodesModelToCalc,linksModel) <> 0 ),
        0.5
        * transfer_linksDecom.l(linksModel,yearsToCalc,transfer_techs,vintage)
        * accounting_transferLinks(indicator,linksModel,transfer_techs,vintage,"perLinkDecom")

        + 0.5
        * transfer_linksTotal.l(linksModel,yearsToCalc,transfer_techs,vintage)
        * accounting_transferLinks(indicator,linksModel,transfer_techs,vintage,"perLinkTotal")

        + 0.5
        * sum (link_types,
            transfer_linksDecom.l(linksModel,yearsToCalc,transfer_techs,vintage)
            * transfer_lengthParam(linksModel,link_types,"length")
            * accounting_transferPerLength(indicator,linksModel,transfer_techs,vintage,link_types,"perLengthDecom")

            + transfer_linksTotal.l(linksModel,yearsToCalc,transfer_techs,vintage)
            * transfer_lengthParam(linksModel,link_types,"length")
            * accounting_transferPerLength(indicator,linksModel,transfer_techs,vintage,link_types,"perLengthTotal"))

        + 0.5
        * sum (timeModelToCalc,
            transfer_flowAlong.l(timeModelToCalc,linksModel,yearsToCalc,transfer_techs,vintage)
            * timeLength(timeModelToCalc)
            * ( accounting_transferLinks(indicator,linksModel,transfer_techs,vintage,"perFlow")
                + accounting_transferLinks(indicator,linksModel,transfer_techs,vintage,"perFlowAlong"))

            + transfer_flowAgainst.l(timeModelToCalc,linksModel,yearsToCalc,transfer_techs,vintage)
            * timeLength(timeModelToCalc)
            * ( accounting_transferLinks(indicator,linksModel,transfer_techs,vintage,"perFlow")
                + accounting_transferLinks(indicator,linksModel,transfer_techs,vintage,"perFlowAgainst")))

        + 0.5
        * sum ((timeModelToCalc, link_types),
            transfer_flowAlong.l(timeModelToCalc,linksModel,yearsToCalc,transfer_techs,vintage)
            * timeLength(timeModelToCalc)
            * transfer_lengthParam(linksModel,link_types,"length")
            * (accounting_transferPerLength(indicator,linksModel,transfer_techs,vintage,link_types,"perFlow")
                + accounting_transferPerLength(indicator,linksModel,transfer_techs,vintage,link_types,"perFlowAlong"))

            + transfer_flowAgainst.l(timeModelToCalc,linksModel,yearsToCalc,transfer_techs,vintage)
            * timeLength(timeModelToCalc)
            * transfer_lengthParam(linksModel,link_types,"length")
            * (accounting_transferPerLength(indicator,linksModel,transfer_techs,vintage,link_types,"perFlow")
                + accounting_transferPerLength(indicator,linksModel,transfer_techs,vintage,link_types,"perFlowAgainst"))))


* == sources / sinks ==
    + sum ((timeModelToCalc,sourcesink_techs(techs),commodity)
            $sourcesink_enabled(nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity),
        sourcesink_flow.l(timeModelToCalc,nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity)
        * timeLength(timeModelToCalc)
        * accounting_sourcesinkFlow(indicator,nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity,"perFlow"))
    ;


* ==== full set of main indicators ====

indicator_accounting(%selscen%accNodesModel,accYears,indicator)
    $sum((accNodesModel_a,accYears_a,indicator_a)
            $compoundIndicatorsFull(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a), 1)
    =
    + sum((accNodesModel_a,accYears_a,indicator_a)
            $(compoundIndicatorsFull(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a)
                and variableIndicators(accNodesModel_a,accYears_a,indicator_a)),
        compoundIndicatorsFull(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a)
        * accounting_indicator.l(accNodesModel_a,accYears_a,indicator_a))

    + sum((accNodesModel_a,accYears_a,indicator_a)
            $(compoundIndicatorsFull(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a)),
        compoundIndicatorsFull(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a)
        * (sum((nodesModelToCalc,yearsToCalc,techs)
                $(sameas(accNodesModel_a, nodesModelToCalc)
                    and sameas(accYears_a, yearsToCalc)),
                indicator_accounting_detailed(%selscen%indicator_a,nodesModelToCalc,yearsToCalc,techs))))
    ;

indicator_accounting_comp(%selscen%accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a)
    $compoundIndicatorsFull(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a)
    =
        compoundIndicatorsFull(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a)
        * accounting_indicator.l(accNodesModel_a,accYears_a,indicator_a)
            $variableIndicators(accNodesModel_a,accYears_a,indicator_a)

        + compoundIndicatorsFull(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a)
        * (sum((nodesModelToCalc,yearsToCalc,techs)
                $(sameas(accNodesModel_a, nodesModelToCalc)
                    and sameas(accYears_a, yearsToCalc)),
                indicator_accounting_detailed(%selscen%indicator_a,nodesModelToCalc,yearsToCalc,techs)))
    ;

indicator_accounting_ref(%selscen%accNodesModel,accYears,indicator)
    $accounting_indicator.l(accNodesModel,accYears,indicator)
    = accounting_indicator.l(accNodesModel,accYears,indicator);

* ==== link-specific indicator accounting ====

indicator_accounting_links(%selscen%indicator,nodesModelToCalc,nodesModelToCalc_a,linksModelToCalc,yearsToCalc,techs)
    $(transfer_incidenceModel(nodesModelToCalc,linksModelToCalc) < 0
        and transfer_incidenceModel(nodesModelToCalc_a,linksModelToCalc) > 0)
    =
    sum ((transfer_techs(techs),vintage)
                $(transfer_availTech(linksModelToCalc,yearsToCalc,transfer_techs,vintage)
                    and accounting_transferLinks(indicator,linksModelToCalc,transfer_techs,vintage,"useAnnuity") = 0),
        transfer_linksBuild.l(linksModelToCalc,yearsToCalc,transfer_techs,vintage)
        * accounting_transferLinks(indicator,linksModelToCalc,transfer_techs,vintage,"perLinkBuild") )

    + sum ((years_a,transfer_techs(techs),vintage)
                $(transfer_availTech(linksModelToCalc,yearsToCalc,transfer_techs,vintage)
                    and accounting_transferLinks(indicator,linksModelToCalc,transfer_techs,vintage,"useAnnuity") = 1
                    and years_a.val + accounting_transferLinks(indicator,linksModelToCalc,transfer_techs,vintage,"amorTime") > yearsToCalc.val
                    and years_a.val <= yearsToCalc.val ),
        transfer_linksBuild.l(linksModelToCalc,years_a,transfer_techs,vintage)
        * accounting_transferLinks(indicator,linksModelToCalc,transfer_techs,vintage,"perLinkBuild")
        * accounting_annuityFactor_transferLink(indicator,linksModelToCalc,transfer_techs,vintage) )

    + sum ((transfer_techs(techs),vintage,link_types)
                $(transfer_availTech(linksModelToCalc,yearsToCalc,transfer_techs,vintage)
                    and accounting_transferPerLength(indicator,linksModelToCalc,transfer_techs,vintage,link_types,"useAnnuity") = 0 ),
        transfer_linksBuild.l(linksModelToCalc,yearsToCalc,transfer_techs,vintage)
        * transfer_lengthParam(linksModelToCalc,link_types,"length")
        * accounting_transferPerLength(indicator,linksModelToCalc,transfer_techs,vintage,link_types,"perLengthBuild") )

    + sum ((years_a,transfer_techs(techs),vintage,link_types)
                $(transfer_availTech(linksModelToCalc,yearsToCalc,transfer_techs,vintage)
                    and accounting_transferPerLength(indicator,linksModelToCalc,transfer_techs,vintage,link_types,"useAnnuity") = 1
                    and years_a.val + accounting_transferPerLength(indicator,linksModelToCalc,transfer_techs,vintage,link_types,"amorTime") > yearsToCalc.val
                    and years_a.val <= yearsToCalc.val ),
        transfer_linksBuild.l(linksModelToCalc,years_a,transfer_techs,vintage)
        * transfer_lengthParam(linksModelToCalc,link_types,"length")
        * accounting_transferPerLength(indicator,linksModelToCalc,transfer_techs,vintage,link_types,"perLengthBuild")
        * accounting_annuityFactor_transferPerLength(indicator,linksModelToCalc,transfer_techs,vintage,link_types) )

    + sum ((transfer_techs(techs),vintage)
                $(transfer_usedTech(linksModelToCalc,yearsToCalc,transfer_techs,vintage)),
        transfer_linksDecom.l(linksModelToCalc,yearsToCalc,transfer_techs,vintage)
        * accounting_transferLinks(indicator,linksModelToCalc,transfer_techs,vintage,"perLinkDecom")

        + transfer_linksTotal.l(linksModelToCalc,yearsToCalc,transfer_techs,vintage)
        * accounting_transferLinks(indicator,linksModelToCalc,transfer_techs,vintage,"perLinkTotal")

        + sum (link_types,
            transfer_linksDecom.l(linksModelToCalc,yearsToCalc,transfer_techs,vintage)
            * transfer_lengthParam(linksModelToCalc,link_types,"length")
            * accounting_transferPerLength(indicator,linksModelToCalc,transfer_techs,vintage,link_types,"perLengthDecom")

            + transfer_linksTotal.l(linksModelToCalc,yearsToCalc,transfer_techs,vintage)
            * transfer_lengthParam(linksModelToCalc,link_types,"length")
            * accounting_transferPerLength(indicator,linksModelToCalc,transfer_techs,vintage,link_types,"perLengthTotal"))

        + sum (timeModelToCalc,
            transfer_flowAlong.l(timeModelToCalc,linksModelToCalc,yearsToCalc,transfer_techs,vintage)
            * timeLength(timeModelToCalc)
            * ( accounting_transferLinks(indicator,linksModelToCalc,transfer_techs,vintage,"perFlow")
                + accounting_transferLinks(indicator,linksModelToCalc,transfer_techs,vintage,"perFlowAlong"))

            + transfer_flowAgainst.l(timeModelToCalc,linksModelToCalc,yearsToCalc,transfer_techs,vintage)
            * timeLength(timeModelToCalc)
            * ( accounting_transferLinks(indicator,linksModelToCalc,transfer_techs,vintage,"perFlow")
                + accounting_transferLinks(indicator,linksModelToCalc,transfer_techs,vintage,"perFlowAgainst")))

        + sum ((timeModelToCalc, link_types),
            transfer_flowAlong.l(timeModelToCalc,linksModelToCalc,yearsToCalc,transfer_techs,vintage)
            * timeLength(timeModelToCalc)
            * transfer_lengthParam(linksModelToCalc,link_types,"length")
            * (accounting_transferPerLength(indicator,linksModelToCalc,transfer_techs,vintage,link_types,"perFlow")
                + accounting_transferPerLength(indicator,linksModelToCalc,transfer_techs,vintage,link_types,"perFlowAlong"))

            + transfer_flowAgainst.l(timeModelToCalc,linksModelToCalc,yearsToCalc,transfer_techs,vintage)
            * timeLength(timeModelToCalc)
            * transfer_lengthParam(linksModelToCalc,link_types,"length")
            * (accounting_transferPerLength(indicator,linksModelToCalc,transfer_techs,vintage,link_types,"perFlow")
                + accounting_transferPerLength(indicator,linksModelToCalc,transfer_techs,vintage,link_types,"perFlowAgainst")))
    );


* == converter capacities ==

converter_ratedInput(converter_techs,vintage,commodity)
    $(converter_techParam(converter_techs,vintage,"lifeTime") > 0)
    = smin(activity$converter_usedAct(converter_techs,vintage,activity),
            converter_coefficient(converter_techs,vintage,activity,commodity,"coefficient"));
converter_ratedInput(converter_techs,vintage,commodity)
    $(converter_ratedInput(converter_techs,vintage,commodity) > 0)
    = 0;

converter_ratedOutput(converter_techs,vintage,commodity)
    $(converter_techParam(converter_techs,vintage,"lifeTime") > 0)
    = smax(activity$converter_usedAct(converter_techs,vintage,activity),
            converter_coefficient(converter_techs,vintage,activity,commodity,"coefficient"));
converter_ratedOutput(converter_techs,vintage,commodity)
    $(converter_ratedOutput(converter_techs,vintage,commodity) < 0)
    = 0;

converter_ratedOutput_min(nodesModelToCalc,years,converter_techs(techs),commodity)
    $sum(vintage$(converter_usedTech(nodesModelToCalc,years,converter_techs,vintage)
        and converter_ratedOutput(converter_techs,vintage,commodity)), 1)
    = smin(vintage$(converter_usedTech(nodesModelToCalc,years,converter_techs,vintage)
            and converter_ratedOutput(converter_techs,vintage,commodity)),
        converter_ratedOutput(converter_techs,vintage,commodity));

converter_ratedOutput_max(nodesModelToCalc,years,converter_techs(techs),commodity)
    $sum(vintage$(converter_usedTech(nodesModelToCalc,years,converter_techs,vintage)
        and converter_ratedOutput(converter_techs,vintage,commodity)), 1)
    = smax(vintage$(converter_usedTech(nodesModelToCalc,years,converter_techs,vintage)
            and converter_ratedOutput(converter_techs,vintage,commodity)),
        converter_ratedOutput(converter_techs,vintage,commodity));

converter_units_ext(%selscen%map_accNodesPostCalc(accNodesModel,nodesModelToCalc),map_accYearsPostCalc(accYears,years),converter_techs(techs),vintage,"build")
    = converter_unitsBuild.l(nodesModelToCalc,years,converter_techs,vintage);

converter_units_ext(%selscen%map_accNodesPostCalc(accNodesModel,nodesModelToCalc),map_accYearsPostCalc(accYears,years),converter_techs(techs),vintage,"decom")
    = converter_unitsDecom.l(nodesModelToCalc,years,converter_techs,vintage);

converter_units_ext(%selscen%map_accNodesPostCalc(accNodesModel,nodesModelToCalc),map_accYearsPostCalc(accYears,years),converter_techs(techs),vintage,"total")
    = converter_unitsTotal.l(nodesModelToCalc,years,converter_techs,vintage);

converter_caps_ext(%selscen%map_accNodesPostCalc(accNodesModel,nodesModelToCalc),map_accYearsPostCalc(accYears,years),converter_techs(techs),vintage,commodity,capType)
    $(converter_ratedOutput(converter_techs,vintage,commodity))
    = converter_units_ext(%selscen%accNodesModel,nodesModelToCalc,accYears,years,techs,vintage,capType)
        * converter_ratedOutput(converter_techs,vintage,commodity);

converter_caps_ext(%selscen%map_accNodesPostCalc(accNodesModel,nodesModelToCalc),map_accYearsPostCalc(accYears,years),converter_techs(techs),vintage,commodity,"lowerLimit")
    $(converter_usedTech(nodesModelToCalc,years,converter_techs,vintage)
        and converter_capacityParam(nodesModelToCalc,years,converter_techs,"unitsLowerLimit") > 0
        and converter_ratedOutput_min(nodesModelToCalc,years,converter_techs,commodity))
    = converter_capacityParam(nodesModelToCalc,years,converter_techs,"unitsLowerLimit")
        * converter_ratedOutput_min(nodesModelToCalc,years,converter_techs,commodity);

converter_caps_ext(%selscen%map_accNodesPostCalc(accNodesModel,nodesModelToCalc),map_accYearsPostCalc(accYears,years),converter_techs(techs),vintage,commodity,"upperLimit")
    $(converter_usedTech(nodesModelToCalc,years,converter_techs,vintage)
        and converter_capacityParam(nodesModelToCalc,years,converter_techs,"unitsUpperLimit") < inf
        and converter_ratedOutput_min(nodesModelToCalc,years,converter_techs,commodity))
    = converter_capacityParam(nodesModelToCalc,years,converter_techs,"unitsUpperLimit")
        * converter_ratedOutput_max(nodesModelToCalc,years,converter_techs,commodity);

converter_caps(%selscen%accNodesModel,accYears,converter_techs(techs),commodity,capType)
    = sum((nodesModelToCalc,years,vintage)
            $(map_accNodesPostCalc(accNodesModel,nodesModelToCalc)
                and map_accYearsPostCalc(accYears,years)),
        converter_caps_ext(%selscen%accNodesModel,nodesModelToCalc,accYears,years,converter_techs,vintage,commodity,capType));
option clear = converter_caps_ext;

converter_units(%selscen%accNodesModel,accYears,converter_techs(techs),vintage,capType)
    = sum((nodesModelToCalc,years)
            $(map_accNodesPostCalc(accNodesModel,nodesModelToCalc)
                and map_accYearsPostCalc(accYears,years)),
        converter_units_ext(%selscen%accNodesModel,nodesModelToCalc,accYears,years,converter_techs,vintage,capType));
option clear = converter_units_ext;


* ==== transfer capacities ====

max_transfer_coefficient(transfer_techs(techs),commodity)
	$(smax(vintage, transfer_coefficient(transfer_techs,vintage,commodity,"coefficient")) > 0)
	= smax(vintage, transfer_coefficient(transfer_techs,vintage,commodity,"coefficient"));

min_transfer_coefficient(transfer_techs(techs),commodity)
	$(smin(vintage, transfer_coefficient(transfer_techs,vintage,commodity,"coefficient")) > 0)
	= smin(vintage, transfer_coefficient(transfer_techs,vintage,commodity,"coefficient"));

transfer_usedStartEnd(nodesModelToCalc_start,nodesModelToCalc_end,linksModel,years,transfer_techs)
    $(transfer_incidenceModel(nodesModelToCalc_start,linksModel) < 0
        and transfer_incidenceModel(nodesModelToCalc_end,linksModel) > 0
        and sum(vintage$transfer_usedTech(linksModel,years,transfer_techs,vintage), 1))
    = yes;

transfer_links(%selscen%linksModel,years,transfer_techs,vintage,"build")
    = transfer_linksBuild.l(linksModel,years,transfer_techs,vintage);

transfer_links(%selscen%linksModel,years,transfer_techs,vintage,"decom")
    = transfer_linksDecom.l(linksModel,years,transfer_techs,vintage);

transfer_links(%selscen%linksModel,years,transfer_techs,vintage,"total")
    = transfer_linksTotal.l(linksModel,years,transfer_techs,vintage);

transfer_caps(%selscen%nodesModelToCalc_start,nodesModelToCalc_end,linksModel,years,transfer_techs(techs),commodity,capType)
    $(transfer_usedStartEnd(nodesModelToCalc_start,nodesModelToCalc_end,linksModel,years,transfer_techs)
        and sum(vintage$(transfer_usedTech(linksModel,years,transfer_techs,vintage)
                            and transfer_coefficient(transfer_techs,vintage,commodity,"coefficient")), 1))
    = sum(vintage,
        transfer_links(%selscen%linksModel,years,transfer_techs,vintage,capType)
        * transfer_coefficient(transfer_techs,vintage,commodity,"coefficient"));

transfer_caps(%selscen%nodesModelToCalc_start,nodesModelToCalc_end,linksModel,years,transfer_techs(techs),commodity,"lowerLimit")
    $(transfer_usedStartEnd(nodesModelToCalc_start,nodesModelToCalc_end,linksModel,years,transfer_techs)
        and sum(vintage$(transfer_usedTech(linksModel,years,transfer_techs,vintage)
                            and transfer_coefficient(transfer_techs,vintage,commodity,"coefficient")), 1))
    = transfer_linksParam(linksModel,years,transfer_techs,'linksLowerLimit')
        * min_transfer_coefficient(transfer_techs,commodity);

transfer_caps(%selscen%nodesModelToCalc_start,nodesModelToCalc_end,linksModel,years,transfer_techs(techs),commodity,"upperLimit")
    $(transfer_usedStartEnd(nodesModelToCalc_start,nodesModelToCalc_end,linksModel,years,transfer_techs)
        and sum(vintage$(transfer_usedTech(linksModel,years,transfer_techs,vintage)
                            and transfer_coefficient(transfer_techs,vintage,commodity,"coefficient")), 1))
    = transfer_linksParam(linksModel,years,transfer_techs,'linksUpperLimit')
        * max_transfer_coefficient(transfer_techs,commodity);


* ==== transfer flows ====

transfer_flows(%selscen%timeModelToCalc,nodesModelToCalc_start,nodesModelToCalc_end,linksModel,yearsToCalc,transfer_techs(techs),commodity)
    $transfer_usedStartEnd(nodesModelToCalc_start,nodesModelToCalc_end,linksModel,yearsToCalc,transfer_techs)
    = sum(vintage,
        ( transfer_flowAlong.l(timeModelToCalc,linksModel,yearsToCalc,transfer_techs,vintage)
            - transfer_flowAgainst.l(timeModelToCalc,linksModel,yearsToCalc,transfer_techs,vintage) )
        * timeLength(timeModelToCalc)
        * transfer_coefficient(transfer_techs,vintage,commodity,"coefficient"));

transfer_flows_annual(%selscen%nodesModelToCalc_start,nodesModelToCalc_end,linksModel,yearsToCalc,transfer_techs(techs),commodity,"netto")
    $transfer_usedStartEnd(nodesModelToCalc_start,nodesModelToCalc_end,linksModel,yearsToCalc,transfer_techs)
    = sum(timeModelToCalc,
        transfer_flows(%selscen%timeModelToCalc,nodesModelToCalc_start,nodesModelToCalc_end,linksModel,yearsToCalc,transfer_techs,commodity));

transfer_flows_annual(%selscen%nodesModelToCalc_start,nodesModelToCalc_end,linksModel,yearsToCalc,transfer_techs(techs),commodity,"positive")
    $transfer_usedStartEnd(nodesModelToCalc_start,nodesModelToCalc_end,linksModel,yearsToCalc,transfer_techs)
    = sum(timeModelToCalc
            $(transfer_flows(%selscen%timeModelToCalc,nodesModelToCalc_start,nodesModelToCalc_end,linksModel,yearsToCalc,transfer_techs,commodity) > 0),
        transfer_flows(%selscen%timeModelToCalc,nodesModelToCalc_start,nodesModelToCalc_end,linksModel,yearsToCalc,transfer_techs,commodity));

transfer_flows_annual(%selscen%nodesModelToCalc_start,nodesModelToCalc_end,linksModel,yearsToCalc,transfer_techs(techs),commodity,"negative")
    $transfer_usedStartEnd(nodesModelToCalc_start,nodesModelToCalc_end,linksModel,yearsToCalc,transfer_techs)
    = sum(timeModelToCalc
            $(transfer_flows(%selscen%timeModelToCalc,nodesModelToCalc_start,nodesModelToCalc_end,linksModel,yearsToCalc,transfer_techs,commodity) < 0),
        transfer_flows(%selscen%timeModelToCalc,nodesModelToCalc_start,nodesModelToCalc_end,linksModel,yearsToCalc,transfer_techs,commodity));

transfer_flows_annual(%selscen%nodesModelToCalc_start,nodesModelToCalc_end,linksModel,yearsToCalc,transfer_techs(techs),commodity,"brutto")
    $transfer_usedStartEnd(nodesModelToCalc_start,nodesModelToCalc_end,linksModel,yearsToCalc,transfer_techs)
    = transfer_flows_annual(%selscen%nodesModelToCalc_start,nodesModelToCalc_end,linksModel,yearsToCalc,transfer_techs,commodity,"positive")
        - transfer_flows_annual(%selscen%nodesModelToCalc_start,nodesModelToCalc_end,linksModel,yearsToCalc,transfer_techs,commodity,"negative");

transfer_flows_annual(%selscen%nodesModelToCalc_start,nodesModelToCalc_end,linksModel,yearsToCalc,transfer_techs(techs),commodity,"flh")
    $(transfer_usedStartEnd(nodesModelToCalc_start,nodesModelToCalc_end,linksModel,yearsToCalc,transfer_techs)
        and transfer_flows_annual(%selscen%nodesModelToCalc_start,nodesModelToCalc_end,linksModel,yearsToCalc,transfer_techs,commodity,"brutto") > 0
        and transfer_caps(%selscen%nodesModelToCalc_start,nodesModelToCalc_end,linksModel,yearsToCalc,transfer_techs,commodity,"total") > 0)
    = transfer_flows_annual(%selscen%nodesModelToCalc_start,nodesModelToCalc_end,linksModel,yearsToCalc,transfer_techs,commodity,"brutto")
        / transfer_caps(%selscen%nodesModelToCalc_start,nodesModelToCalc_end,linksModel,yearsToCalc,transfer_techs,commodity,"total");

transfer_losses(%selscen%timeModelToCalc,nodesModelToCalc_start,nodesModelToCalc_end,linksModel,yearsToCalc,transfer_techs(techs),commodity)
    $(transfer_usedStartEnd(nodesModelToCalc_start,nodesModelToCalc_end,linksModel,yearsToCalc,transfer_techs)
        and (sum(vintage$transfer_coefPerFlow(transfer_techs,vintage,commodity,"coefPerFlow"), 1)
            or sum((vintage, link_types)$transfer_coefPerLength(transfer_techs,vintage,commodity,link_types,"coefPerLength"), 1)))
    = -1 * abs(sum(vintage,
        ( transfer_flowAlong.l(timeModelToCalc,linksModel,yearsToCalc,transfer_techs,vintage)
            - transfer_flowAgainst.l(timeModelToCalc,linksModel,yearsToCalc,transfer_techs,vintage) )
        * ( transfer_coefPerFlow(transfer_techs,vintage,commodity,"coefPerFlow")
            + sum(link_types,
                transfer_coefPerLength(transfer_techs,vintage,commodity,link_types,"coefPerLength")
                * transfer_lengthParam(linksModel,link_types,"length")))))
        * timeLength(timeModelToCalc);

transfer_losses_annual(%selscen%nodesModelToCalc_start,nodesModelToCalc_end,linksModel,yearsToCalc,transfer_techs(techs),commodity,"netto")
    $(transfer_usedStartEnd(nodesModelToCalc_start,nodesModelToCalc_end,linksModel,yearsToCalc,transfer_techs)
        and (sum(vintage$transfer_coefPerFlow(transfer_techs,vintage,commodity,"coefPerFlow"), 1)
            or sum((vintage, link_types)$transfer_coefPerLength(transfer_techs,vintage,commodity,link_types,"coefPerLength"), 1)))
    = sum(timeModelToCalc,
        transfer_losses(%selscen%timeModelToCalc,nodesModelToCalc_start,nodesModelToCalc_end,linksModel,yearsToCalc,transfer_techs,commodity));


* ==== storage capacities ====

storage_size_max(nodesModelToCalc,years,storage_techs(techs),commodity)
    $sum(vintage$(storage_usedTech(nodesModelToCalc,years,storage_techs,vintage)
            and storage_usedCom(storage_techs,vintage,commodity)
            and storage_sizeParam(storage_techs,vintage,commodity,"size")), 1)
    = smax(vintage$(storage_usedTech(nodesModelToCalc,years,storage_techs,vintage)
            and storage_usedCom(storage_techs,vintage,commodity)
            and storage_sizeParam(storage_techs,vintage,commodity,"size")),
        storage_sizeParam(storage_techs,vintage,commodity,"size"));

storage_size_min(nodesModelToCalc,years,storage_techs(techs),commodity)
    $sum(vintage$(storage_usedTech(nodesModelToCalc,years,storage_techs,vintage)
            and storage_usedCom(storage_techs,vintage,commodity)
            and storage_sizeParam(storage_techs,vintage,commodity,"size")), 1)
    = smin(vintage$(storage_usedTech(nodesModelToCalc,years,storage_techs,vintage)
            and storage_usedCom(storage_techs,vintage,commodity)
            and storage_sizeParam(storage_techs,vintage,commodity,"size")),
        storage_sizeParam(storage_techs,vintage,commodity,"size"));

storage_units_ext(%selscen%map_accNodesPostCalc(accNodesModel,nodesModelToCalc),map_accYearsPostCalc(accYears,years),storage_techs(techs),vintage,"build")
    $storage_usedTech(nodesModelToCalc,years,storage_techs,vintage)
    = storage_unitsBuild.l(nodesModelToCalc,years,storage_techs,vintage);

storage_units_ext(%selscen%map_accNodesPostCalc(accNodesModel,nodesModelToCalc),map_accYearsPostCalc(accYears,years),storage_techs(techs),vintage,"decom")
    $storage_usedTech(nodesModelToCalc,years,storage_techs,vintage)
    = storage_unitsDecom.l(nodesModelToCalc,years,storage_techs,vintage);

storage_units_ext(%selscen%map_accNodesPostCalc(accNodesModel,nodesModelToCalc),map_accYearsPostCalc(accYears,years),storage_techs(techs),vintage,"total")
    $storage_usedTech(nodesModelToCalc,years,storage_techs,vintage)
    = storage_unitsTotal.l(nodesModelToCalc,years,storage_techs,vintage);

storage_caps_ext(%selscen%map_accNodesPostCalc(accNodesModel,nodesModelToCalc),map_accYearsPostCalc(accYears,years),storage_techs(techs),vintage,commodity,capType)
    $storage_sizeParam(storage_techs,vintage,commodity,"size")
    = storage_units_ext(%selscen%accNodesModel,nodesModelToCalc,accYears,years,storage_techs,vintage,capType)
        * storage_sizeParam(storage_techs,vintage,commodity,"size");

storage_caps_ext(%selscen%map_accNodesPostCalc(accNodesModel,nodesModelToCalc),map_accYearsPostCalc(accYears,years),storage_techs(techs),vintage,commodity,"lowerLimit")
    $(storage_usedTech(nodesModelToCalc,years,storage_techs,vintage)
        and storage_reservoirParam(nodesModelToCalc,years,storage_techs,"unitsLowerLimit") > 0
        and storage_size_min(nodesModelToCalc,years,storage_techs,commodity))
    = storage_reservoirParam(nodesModelToCalc,years,storage_techs,"unitsLowerLimit")
        * storage_size_min(nodesModelToCalc,years,storage_techs,commodity);

storage_caps_ext(%selscen%map_accNodesPostCalc(accNodesModel,nodesModelToCalc),map_accYearsPostCalc(accYears,years),storage_techs(techs),vintage,commodity,"upperLimit")
    $(storage_usedTech(nodesModelToCalc,years,storage_techs,vintage)
        and storage_reservoirParam(nodesModelToCalc,years,storage_techs,"unitsUpperLimit") < inf
        and storage_size_max(nodesModelToCalc,years,storage_techs,commodity))
    = storage_reservoirParam(nodesModelToCalc,years,storage_techs,"unitsUpperLimit")
        * storage_size_max(nodesModelToCalc,years,storage_techs,commodity);

storage_caps_ext(%selscen%map_accNodesPostCalc(accNodesModel,nodesModelToCalc),map_accYearsPostCalc(accYears,years),storage_techs(techs),vintage,commodity,"total_degraded")
    $(storage_usedCom(storage_techs,vintage,commodity)
        and (storage_techParam(storage_techs,vintage,"annualDegradation") > 0
             or storage_techParam(storage_techs,vintage,"usageDegradation")))
    = sum((degradation_states,yearsCom), (storage_degradationParam(storage_techs,vintage,degradation_states,"remainingCapacity")
                                    - (years.val - yearsCom.val) * storage_techParam(storage_techs,vintage,"annualDegradation"))
                                    * storage_unitsStateTracker.l(nodesModelToCalc,years,yearsCom,storage_techs,vintage,degradation_states)
                                    * storage_sizeParam(storage_techs,vintage,commodity,"size"));

storage_caps(%selscen%accNodesModel,accYears,storage_techs(techs),commodity,capType)
    = sum((nodesModelToCalc,years,vintage)
            $(map_accNodesPostCalc(accNodesModel,nodesModelToCalc)
                and map_accYearsPostCalc(accYears,years)),
        storage_caps_ext(%selscen%accNodesModel,nodesModelToCalc,accYears,years,storage_techs,vintage,commodity,capType));
option clear = storage_caps_ext;

storage_units(%selscen%accNodesModel,accYears,storage_techs(techs),vintage,capType)
    = sum((nodesModelToCalc,years)
            $(map_accNodesPostCalc(accNodesModel,nodesModelToCalc)
                and map_accYearsPostCalc(accYears,years)),
        storage_units_ext(%selscen%accNodesModel,nodesModelToCalc,accYears,years,storage_techs,vintage,capType));
option clear = storage_units_ext;


* ==== storage levels and flows ====

storage_flows_ext(%selscen%timeModelToCalc,map_accNodesPostCalc(accNodesModel,nodesModelToCalc),map_accYearsPostCalc(accYears,yearsToCalc),storage_techs(techs),vintage,commodity)
    $(storage_usedTech(nodesModelToCalc,yearsToCalc,storage_techs,vintage) and balance_techComm(storage_techs,commodity))
    = storage_level.l(timeModelToCalc,nodesModelToCalc,yearsToCalc,storage_techs,vintage,commodity)
        - storage_level.l(timeModelToCalc--1,nodesModelToCalc,yearsToCalc,storage_techs,vintage,commodity)
$iftheni.pips %method%==solpoint
        - (storage_level.l(timeModelToCalc--1,nodesModelToCalc,yearsToCalc,storage_techs,vintage,commodity)
            * storage_sizeParam(storage_techs,vintage,commodity,"selfdischarge")
        + storage_sizeParam(storage_techs,vintage,commodity,"selfdischargeAbs"))
        * timeLength(timeModelToCalc)
$else.pips
    + storage_losses.l(timeModelToCalc,nodesModelToCalc,yearsToCalc,storage_techs,vintage,commodity)
        * timeLength(timeModelToCalc)
$endif.pips
    ;

storage_flows(%selscen%timeModelToCalc,accNodesModel,accYears,storage_techs(techs),commodity)
    = sum((nodesModelToCalc,yearsToCalc,vintage)
            $(storage_usedTech(nodesModelToCalc,yearsToCalc,storage_techs,vintage) and balance_techComm(storage_techs,commodity)),
        storage_flows_ext(%selscen%timeModelToCalc,accNodesModel,nodesModelToCalc,accYears,yearsToCalc,storage_techs,vintage,commodity));
option clear = storage_flows_ext;

storage_level_out_ext(%selscen%timeModelToCalc,map_accNodesPostCalc(accNodesModel,nodesModelToCalc),map_accYearsPostCalc(accYears,yearsToCalc),storage_techs(techs),vintage,commodity)
    $(storage_usedTech(nodesModelToCalc,yearsToCalc,storage_techs,vintage) and balance_techComm(storage_techs,commodity))
    = storage_level.l(timeModelToCalc,nodesModelToCalc,yearsToCalc,storage_techs,vintage,commodity);

storage_level_out(%selscen%timeModelToCalc,accNodesModel,accYears,storage_techs(techs),commodity)
    = sum((nodesModelToCalc,yearsToCalc,vintage)
            $(storage_usedTech(nodesModelToCalc,yearsToCalc,storage_techs,vintage) and balance_techComm(storage_techs,commodity)),
        storage_level_out_ext(%selscen%timeModelToCalc,accNodesModel,nodesModelToCalc,accYears,yearsToCalc,storage_techs,vintage,commodity));
option clear = storage_level_out_ext;

storage_flows_annual(%selscen%accNodesModel,accYears,storage_techs(techs),commodity,"netto")
    = sum(timeModelToCalc,
        storage_flows(%selscen%timeModelToCalc,accNodesModel,accYears,storage_techs,commodity));

storage_flows_annual(%selscen%accNodesModel,accYears,storage_techs(techs),commodity,"positive")
    = sum(timeModelToCalc
            $(storage_flows(%selscen%timeModelToCalc,accNodesModel,accYears,storage_techs,commodity) > 0),
        storage_flows(%selscen%timeModelToCalc,accNodesModel,accYears,storage_techs,commodity));

storage_flows_annual(%selscen%accNodesModel,accYears,storage_techs(techs),commodity,"negative")
    = sum(timeModelToCalc
            $(storage_flows(%selscen%timeModelToCalc,accNodesModel,accYears,storage_techs,commodity) < 0),
        storage_flows(%selscen%timeModelToCalc,accNodesModel,accYears,storage_techs,commodity));

storage_flows_annual(%selscen%accNodesModel,accYears,storage_techs(techs),commodity,"brutto")
    = storage_flows_annual(%selscen%accNodesModel,accYears,storage_techs,commodity,"positive")
        + storage_flows_annual(%selscen%accNodesModel,accYears,storage_techs,commodity,"negative");

storage_flows_annual(%selscen%accNodesModel,accYears,storage_techs(techs),commodity,"flh")
    $(storage_flows_annual(%selscen%accNodesModel,accYears,storage_techs,commodity,"brutto") > 0
        and storage_caps(%selscen%accNodesModel,accYears,storage_techs,commodity,"total") > 0)
    = storage_flows_annual(%selscen%accNodesModel,accYears,storage_techs,commodity,"brutto")
        / storage_caps(%selscen%accNodesModel,accYears,storage_techs,commodity,"total");

storage_losses_out(%selscen%timeModelToCalc,accNodesModel,accYears,storage_techs(techs),commodity)
    = sum((nodesModelToCalc,yearsToCalc,vintage)
            $(storage_losses.l(timeModelToCalc,nodesModelToCalc,yearsToCalc,storage_techs,vintage,commodity)
                and map_accNodesPostCalc(accNodesModel,nodesModelToCalc)
                and map_accYearsPostCalc(accYears,yearsToCalc)),
        storage_losses.l(timeModelToCalc,nodesModelToCalc,yearsToCalc,storage_techs,vintage,commodity)
        * timeLength(timeModelToCalc));

storage_losses_annual(%selscen%accNodesModel,accYears,storage_techs(techs),commodity,"netto")
    = sum(timeModelToCalc,
        storage_losses_out(%selscen%timeModelToCalc,accNodesModel,accYears,storage_techs,commodity));


* == commodity balance ==

commodity_balance_ext(%selscen%timeModelToCalc,map_accNodesPostCalc(accNodesModel,nodesModelToCalc),map_accYearsPostCalc(accYears,yearsToCalc),converter_techs(techs),vintage,commodity)
    $(converter_usedTech(nodesModelToCalc,yearsToCalc,converter_techs,vintage)
        and balance_usedConverter(nodesModelToCalc,yearsToCalc,commodity))
    = sum((activity)
            $(converter_coefficientProfile(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage,activity,commodity)),
        converter_activity.l(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage,activity)
            * timeLength(timeModelToCalc)
            * converter_coefficientProfile(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage,activity,commodity)
        + converter_unitsUsingActivity_MIP.l(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage,activity)
            * timeLength(timeModelToCalc)
            * converter_coefficient(converter_techs,vintage,activity,commodity,"constant"));

commodity_balance_ext(%selscen%timeModelToCalc,map_accNodesPostCalc(accNodesModel,nodesModelToCalc),map_accYearsPostCalc(accYears,yearsToCalc),storage_techs(techs),vintage,commodity)
    $(storage_usedTech(nodesModelToCalc,yearsToCalc,storage_techs,vintage)
        and storage_sizeParam(storage_techs,vintage,commodity,"size"))
    = commodity_balance_ext(%selscen%timeModelToCalc,accNodesModel,nodesModelToCalc,accYears,yearsToCalc,storage_techs,vintage,commodity)
    + storage_level.l(timeModelToCalc--1,nodesModelToCalc,yearsToCalc,storage_techs,vintage,commodity)
        - storage_level.l(timeModelToCalc,nodesModelToCalc,yearsToCalc,storage_techs,vintage,commodity)
    - storage_losses.l(timeModelToCalc,nodesModelToCalc,yearsToCalc,storage_techs,vintage,commodity)
    * timeLength(timeModelToCalc);

commodity_balance_ext(%selscen%timeModelToCalc,map_accNodesPostCalc(accNodesModel,nodesModelToCalc),map_accYearsPostCalc(accYears,yearsToCalc),sourcesink_techs(techs),vintage,commodity)
    $(sourcesink_enabled(nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity)
        and sameas(yearsToCalc,vintage))
    = commodity_balance_ext(%selscen%timeModelToCalc,accNodesModel,nodesModelToCalc,accYears,yearsToCalc,sourcesink_techs,vintage,commodity)
    + sourcesink_flow.l(timeModelToCalc,nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity)
    * timeLength(timeModelToCalc);

commodity_balance_ext(%selscen%timeModelToCalc,map_accNodesPostCalc(accNodesModel,nodesModelToCalc),map_accYearsPostCalc(accYears,yearsToCalc),transfer_techs(techs),vintage,commodity)
    $sum((linksModel)
        $(transfer_incidenceModel(nodesModelToCalc,linksModel)
                and transfer_usedTech(linksModel,yearsToCalc,transfer_techs,vintage)
                and transfer_coefficient(transfer_techs,vintage,commodity,"coefficient")), 1)
    = commodity_balance_ext(%selscen%timeModelToCalc,accNodesModel,nodesModelToCalc,accYears,yearsToCalc,transfer_techs,vintage,commodity)
    + sum((linksModel)
            $(transfer_incidenceModel(nodesModelToCalc,linksModel)
                and transfer_usedTech(linksModel,yearsToCalc,transfer_techs,vintage)
                and transfer_coefficient(transfer_techs,vintage,commodity,"coefficient")),

            (transfer_flowAlong.l(timeModelToCalc,linksModel,yearsToCalc,transfer_techs,vintage)
                    $(transfer_incidenceModel(nodesModelToCalc,linksModel) > 0)
              + transfer_flowAgainst.l(timeModelToCalc,linksModel,yearsToCalc,transfer_techs,vintage)
                    $(transfer_incidenceModel(nodesModelToCalc,linksModel) < 0))
            * timeLength(timeModelToCalc)
            * transfer_coefficient(transfer_techs,vintage,commodity,"coefficient")

            - (transfer_flowAlong.l(timeModelToCalc,linksModel,yearsToCalc,transfer_techs,vintage)
                    $(transfer_incidenceModel(nodesModelToCalc,linksModel) < 0)
              + transfer_flowAgainst.l(timeModelToCalc,linksModel,yearsToCalc,transfer_techs,vintage)
                    $(transfer_incidenceModel(nodesModelToCalc,linksModel) > 0))
            * timeLength(timeModelToCalc)
            * transfer_coefficient(transfer_techs,vintage,commodity,"coefficient")

            + 0.5 * ((transfer_flowAlong.l(timeModelToCalc,linksModel,yearsToCalc,transfer_techs,vintage)
                        $(transfer_incidenceModel(nodesModelToCalc,linksModel) <> 0)
                    + transfer_flowAgainst.l(timeModelToCalc,linksModel,yearsToCalc,transfer_techs,vintage)
                        $(transfer_incidenceModel(nodesModelToCalc,linksModel) <> 0))
                    * timeLength(timeModelToCalc)
                    * ( transfer_coefPerFlow(transfer_techs,vintage,commodity,"coefPerFlow")
                        + sum(link_types,
                            transfer_coefPerLength(transfer_techs,vintage,commodity,link_types,"coefPerLength")
                            * transfer_lengthParam(linksModel,link_types,"length")))))
    ;

commodity_balance(%selscen%timeModelToCalc,accNodesModel,accYears,balance_techComm(techs,commodity))
    = sum((nodesModelToCalc,yearsToCalc,vintage)
            $(commodity_balance_ext(%selscen%timeModelToCalc,accNodesModel,nodesModelToCalc,accYears,yearsToCalc,techs,vintage,commodity)
                and map_accNodesPostCalc(accNodesModel,nodesModelToCalc)
                and map_accYearsPostCalc(accYears,yearsToCalc)),
        commodity_balance_ext(%selscen%timeModelToCalc,accNodesModel,nodesModelToCalc,accYears,yearsToCalc,techs,vintage,commodity));
option clear = commodity_balance_ext;


* ==== annual commodity sums ====

commodity_balance_annual(%selscen%accNodesModel,accYears,techs,commodity,"netto")
    = sum(timeModelToCalc,
        commodity_balance(%selscen%timeModelToCalc,accNodesModel,accYears,techs,commodity));

commodity_balance_annual(%selscen%accNodesModel,accYears,techs,commodity,"positive")
    = sum(timeModelToCalc
            $(commodity_balance(%selscen%timeModelToCalc,accNodesModel,accYears,techs,commodity) > 0),
        commodity_balance(%selscen%timeModelToCalc,accNodesModel,accYears,techs,commodity));

commodity_balance_annual(%selscen%accNodesModel,accYears,techs,commodity,"negative")
    = sum(timeModelToCalc
            $(commodity_balance(%selscen%timeModelToCalc,accNodesModel,accYears,techs,commodity) < 0),
        commodity_balance(%selscen%timeModelToCalc,accNodesModel,accYears,techs,commodity));

commodity_balance_annual(%selscen%accNodesModel,accYears,techs,commodity,"brutto")
    = commodity_balance_annual(%selscen%accNodesModel,accYears,techs,commodity,"positive")
        - commodity_balance_annual(%selscen%accNodesModel,accYears,techs,commodity,"negative");

commodity_balance_annual(%selscen%accNodesModel,accYears,techs,commodity,"flh")
    $(commodity_balance_annual(%selscen%accNodesModel,accYears,techs,commodity,"brutto") > 0
        and converter_caps(%selscen%accNodesModel,accYears,techs,commodity,"total") > 0 )
    = commodity_balance_annual(%selscen%accNodesModel,accYears,techs,commodity,"brutto")
        / converter_caps(%selscen%accNodesModel,accYears,techs,commodity,"total");


* ==== marginal information ====

marginals_balance(%selscen%timeModelToCalc,nodesModelToCalc,yearsToCalc,commodity)
  $balance_usedBalance(nodesModelToCalc,yearsToCalc,commodity)
  =
  Eq_balance_commodities.m(timeModelToCalc,nodesModelToCalc,yearsToCalc,commodity);

marginals_sourcesink_profile(%selscen%timeModelToCalc,nodesModelToCalc,yearsToCalc,sourcesink_techs(techs),commodity)
    $( sourcesink_flow.m(timeModelToCalc,nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity) <> 0 )
    =
    sourcesink_flow.m(timeModelToCalc,nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity);

marginals_sourcesink_sum(%selscen%nodesModelToCalc,yearsToCalc,sourcesink_techs(techs),commodity)
    $( ( Eq_sourcesink_useLowerSum.m(nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity)
        + Eq_sourcesink_useUpperSum.m(nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity)
        + Eq_sourcesink_useFixedSum.m(nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity) ) <> eps )
    =
    ( Eq_sourcesink_useLowerSum.m(nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity)
    + Eq_sourcesink_useUpperSum.m(nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity)
    + Eq_sourcesink_useFixedSum.m(nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity) );

marginals_indicator_bounds(%selscen%accNodesModel,accYears,indicator)
    $( accounting_indicator.m(accNodesModel,accYears,indicator) <> 0 )
    =
    accounting_indicator.m(accNodesModel,accYears,indicator);


* ==== R2A postcalc information ====

$ifthene.r2a %gdx_r2a%
r2a_has_converter_cost(indicator,nodesModelToCalc,yearsToCalc,converter_techs(techs),vintage,commodity)
    = (sum(activity$(accounting_converterUnits(indicator,nodesModelToCalc,converter_techs,vintage,"perUnitBuild")
                or accounting_converterUnits(indicator,nodesModelToCalc,converter_techs,vintage,"perUnitDecom")
                or accounting_converterUnits(indicator,nodesModelToCalc,converter_techs,vintage,"perUnitTotal")
                or accounting_converterActivity(indicator,nodesModelToCalc,converter_techs,vintage,activity,"perActivity")), 1)
        and converter_ratedOutput(converter_techs,vintage,commodity) <> 0);

r2a_annuity_cost_converter(indicator,accNodesModel,accYears,converter_techs(techs),vintage,commodity)
    $(sum((accNodesModel_a,nodesModelToCalc,yearsToCalc)$(map_accNodes(accNodesModel_a,accNodesModel)
        and sameas(accNodesModel_a,nodesModelToCalc) and sameas(accYears,yearsToCalc)
        and r2a_has_converter_cost(indicator,nodesModelToCalc,yearsToCalc,techs,vintage,commodity)), 1)
    and sum((nodesModelToCalc,accNodesModel_a,yearsToCalc,activity)$(map_accNodes(accNodesModel_a,accNodesModel)
        and sameas(accNodesModel_a,nodesModelToCalc) and sameas(yearsToCalc,accYears)
        and converter_usedTechAct(nodesModelToCalc,yearsToCalc,converter_techs,vintage,activity)), 1) > 0)
    =
    ( sum((nodesModelToCalc,accNodesModel_a,yearsToCalc,activity)$(map_accNodes(accNodesModel_a,accNodesModel)
            and sameas(accNodesModel_a,nodesModelToCalc) and sameas(yearsToCalc,accYears)),
        accounting_converterUnits(indicator,nodesModelToCalc,converter_techs,vintage,"perUnitBuild")
            * accounting_annuityFactor_converter(indicator,nodesModelToCalc,converter_techs,vintage)
        + accounting_converterUnits(indicator,nodesModelToCalc,converter_techs,vintage,"perUnitDecom")
        + accounting_converterUnits(indicator,nodesModelToCalc,converter_techs,vintage,"perUnitTotal"))
        / sum((nodesModelToCalc,accNodesModel_a,yearsToCalc,activity)$(map_accNodes(accNodesModel_a,accNodesModel)
            and sameas(accNodesModel_a,nodesModelToCalc) and sameas(yearsToCalc,accYears)), 1)

    + ( sum((nodesModelToCalc,accNodesModel_a,yearsToCalc,activity)$(map_accNodes(accNodesModel_a,accNodesModel)
            and sameas(accNodesModel_a,nodesModelToCalc) and sameas(yearsToCalc,accYears)
            and converter_usedTechAct(nodesModelToCalc,yearsToCalc,converter_techs,vintage,activity)),
        accounting_converterActivity(indicator,nodesModelToCalc,converter_techs,vintage,activity,"perActivity"))
        / sum((nodesModelToCalc,accNodesModel_a,yearsToCalc,activity)$(map_accNodes(accNodesModel_a,accNodesModel)
            and sameas(accNodesModel_a,nodesModelToCalc) and sameas(yearsToCalc,accYears)
            and converter_usedTechAct(nodesModelToCalc,yearsToCalc,converter_techs,vintage,activity)), 1))

    ) / converter_ratedOutput(converter_techs,vintage,commodity)
    ;

r2a_spec_cost_converter(indicator,accNodesModel,accYears,converter_techs(techs),vintage,commodity)
    $(sum((accNodesModel_a,nodesModelToCalc,yearsToCalc)$(map_accNodes(accNodesModel_a,accNodesModel)
        and sameas(accNodesModel_a,nodesModelToCalc) and sameas(accYears,yearsToCalc)
        and r2a_has_converter_cost(indicator,nodesModelToCalc,yearsToCalc,techs,vintage,commodity)), 1)
    and sum((nodesModelToCalc,accNodesModel_a,yearsToCalc,activity)$(map_accNodes(accNodesModel_a,accNodesModel)
        and sameas(accNodesModel_a,nodesModelToCalc) and sameas(yearsToCalc,accYears)
        and converter_usedTechAct(nodesModelToCalc,yearsToCalc,converter_techs,vintage,activity)), 1) > 0)
    =
    ( sum((nodesModelToCalc,accNodesModel_a,yearsToCalc,activity)$(map_accNodes(accNodesModel_a,accNodesModel)
            and sameas(accNodesModel_a,nodesModelToCalc) and sameas(yearsToCalc,accYears)),
        accounting_converterUnits(indicator,nodesModelToCalc,converter_techs,vintage,"perUnitBuild")
        + accounting_converterUnits(indicator,nodesModelToCalc,converter_techs,vintage,"perUnitDecom")
        + accounting_converterUnits(indicator,nodesModelToCalc,converter_techs,vintage,"perUnitTotal"))
        / sum((nodesModelToCalc,accNodesModel_a,yearsToCalc,activity)$(map_accNodes(accNodesModel_a,accNodesModel)
            and sameas(accNodesModel_a,nodesModelToCalc) and sameas(yearsToCalc,accYears)), 1)

    + ( sum((nodesModelToCalc,accNodesModel_a,yearsToCalc,activity)$(map_accNodes(accNodesModel_a,accNodesModel)
            and sameas(accNodesModel_a,nodesModelToCalc) and sameas(yearsToCalc,accYears)
            and converter_usedTechAct(nodesModelToCalc,yearsToCalc,converter_techs,vintage,activity)),
        accounting_converterActivity(indicator,nodesModelToCalc,converter_techs,vintage,activity,"perActivity"))
        / sum((nodesModelToCalc,accNodesModel_a,yearsToCalc,activity)$(map_accNodes(accNodesModel_a,accNodesModel)
            and sameas(accNodesModel_a,nodesModelToCalc) and sameas(yearsToCalc,accYears)
            and converter_usedTechAct(nodesModelToCalc,yearsToCalc,converter_techs,vintage,activity)), 1))

    ) / converter_ratedOutput(converter_techs,vintage,commodity)
    ;

r2a_spec_cost_fuel(indicator,accNodesModel,accYears,sourcesink_techs(techs),commodity)
    $(sum((accNodesModel_a, nodesModelToCalc, yearsToCalc)$(map_accNodes(accNodesModel_a,accNodesModel)
        and sameas(accNodesModel_a,nodesModelToCalc) and sameas(accYears,yearsToCalc)
        and sourcesink_enabled(nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity)), 1)
    and sum(yearsToCalc$sameas(yearsToCalc,accYears), 1) > 0)
    =
    sum((accNodesModel_a, nodesModelToCalc, yearsToCalc)$(map_accNodes(accNodesModel_a,accNodesModel)
        and sameas(accNodesModel_a,nodesModelToCalc) and sameas(accYears,yearsToCalc)),
        accounting_sourcesinkFlow(indicator,nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity,"perFlow"))
    / sum((accNodesModel_a, nodesModelToCalc, yearsToCalc)$(map_accNodes(accNodesModel_a,accNodesModel)
        and sameas(accNodesModel_a,nodesModelToCalc) and sameas(accYears,yearsToCalc)), 1)
    ;

r2a_spec_cost_indicator(indicator,indicator_a,accNodesModel,accYears)
    $(accounting_perIndicator(indicator,indicator_a,accNodesModel,accYears,"perIndicator") <> 0)
    =
    accounting_perIndicator(indicator,indicator_a,accNodesModel,accYears,"perIndicator");

r2a_converter_efficiencies(converter_techs(techs),vintage,activity,commodity,commodity_a)
    $(converter_coefficient(converter_techs,vintage,activity,commodity,"coefficient") < 0
        and converter_coefficient(converter_techs,vintage,activity,commodity_a,"coefficient") > 0)
    =
    - converter_coefficient(converter_techs,vintage,activity,commodity_a,"coefficient")
    / converter_coefficient(converter_techs,vintage,activity,commodity,"coefficient")
    ;

r2a_converter_avail_factor(accNodesModel,accYears,converter_techs(techs),vintage)
    $(sum((accNodesModel_a, nodesModelToCalc, yearsToCalc)$(map_accNodes(accNodesModel_a,accNodesModel)
        and sameas(accNodesModel_a,nodesModelToCalc) and sameas(accYears,yearsToCalc)
        and converter_usedTech(nodesModelToCalc,yearsToCalc,converter_techs,vintage)
        and not converter_activity_hasProfile(nodesModelToCalc,yearsToCalc,converter_techs,"upper")
        and not converter_activity_hasProfile(nodesModelToCalc,yearsToCalc,converter_techs,"fixed")), 1))
    =
    converter_techParam(converter_techs,vintage,"activityUpperLimit");

r2a_converter_avail_profile_ext(timeModelToCalc,map_accNodesPostCalc(accNodesModel,nodesModelToCalc),map_accYearsPostCalc(accYears,yearsToCalc),converter_usedCom(converter_techs,vintage,commodity),"lower")
    $(converter_usedTech(nodesModelToCalc,yearsToCalc,converter_techs,vintage)
        and (converter_activity_hasProfile(nodesModelToCalc,yearsToCalc,converter_techs,"lower")
                or converter_activity_hasProfile(nodesModelToCalc,yearsToCalc,converter_techs,"fixed")))
    =
    converter_activityProfile(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage,"lower")
    * converter_unitsTotal.l(nodesModelToCalc,yearsToCalc,converter_techs,vintage)
    * converter_ratedOutput(converter_techs,vintage,commodity);

r2a_converter_avail_profile_ext(timeModelToCalc,map_accNodesPostCalc(accNodesModel,nodesModelToCalc),map_accYearsPostCalc(accYears,yearsToCalc),converter_usedCom(converter_techs,vintage,commodity),"upper")
    $(converter_usedTech(nodesModelToCalc,yearsToCalc,converter_techs,vintage)
        and (converter_activity_hasProfile(nodesModelToCalc,yearsToCalc,converter_techs,"upper")
                or converter_activity_hasProfile(nodesModelToCalc,yearsToCalc,converter_techs,"fixed")))
    =
    converter_activityProfile(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage,"upper")
    * converter_unitsTotal.l(nodesModelToCalc,yearsToCalc,converter_techs,vintage)
    * converter_ratedOutput(converter_techs,vintage,commodity);

r2a_converter_avail_profile(timeModelToCalc,accNodesModel,accYears,converter_techs(techs),commodity,profileType)
    = sum((nodesModelToCalc,yearsToCalc,vintage)
            $(converter_usedTech(nodesModelToCalc,yearsToCalc,converter_techs,vintage)
                and map_accNodesPostCalc(accNodesModel,nodesModelToCalc)
                and map_accYearsPostCalc(accYears,yearsToCalc)),
            r2a_converter_avail_profile_ext(timeModelToCalc,accNodesModel,nodesModelToCalc,accYears,yearsToCalc,converter_techs,vintage,commodity,profileType));
option clear = r2a_converter_avail_profile_ext;

r2a_storage_selfdischarge(accNodesModel,accYears,storage_techs(techs),vintage,commodity)
    $(sum((accNodesModel_a, nodesModelToCalc, yearsToCalc)$(map_accNodes(accNodesModel_a,accNodesModel)
        and sameas(accNodesModel_a,nodesModelToCalc) and sameas(accYears,yearsToCalc)
        and storage_usedTech(nodesModelToCalc,yearsToCalc,storage_techs,vintage)), 1))
    =
    storage_sizeParam(storage_techs,vintage,commodity,"selfdischarge");

r2a_storage_e2p(%selscen%accNodesModel,accYears,storage_techs(techs),vintage,commodity,"total")
    $(storage_caps(%selscen%accNodesModel,accYears,storage_techs,commodity,"total") > 0
    and sum(converter_techs, converter_caps(%selscen%accNodesModel,accYears,converter_techs,commodity,"total")) > 0)
    =
    storage_caps(%selscen%accNodesModel,accYears,storage_techs,commodity,"total")
    / (sum(converter_techs, converter_caps(%selscen%accNodesModel,accYears,converter_techs,commodity,"total"))
        / sum(converter_techs$converter_caps(%selscen%accNodesModel,accYears,converter_techs,commodity,"total"), 1));

$endif.r2a


* ==== round profiles to reduce size of gdx ====

$ifthene.roundts %gdx_roundts%
commodity_balance(%selscen%timeModelToCalc,accNodesModel,accYears,techs,commodity)
    $commodity_balance(%selscen%timeModelToCalc,accNodesModel,accYears,techs,commodity)
    = round(commodity_balance(%selscen%timeModelToCalc,accNodesModel,accYears,techs,commodity), 6);

transfer_flows(%selscen%timeModelToCalc,nodesModelToCalc,nodesModelToCalc_a,linksModel,yearsToCalc,transfer_techs,commodity)
    $transfer_flows(%selscen%timeModelToCalc,nodesModelToCalc,nodesModelToCalc_a,linksModel,yearsToCalc,transfer_techs,commodity)
    = round(transfer_flows(%selscen%timeModelToCalc,nodesModelToCalc,nodesModelToCalc_a,linksModel,yearsToCalc,transfer_techs,commodity), 6);

transfer_losses(%selscen%timeModelToCalc,nodesModelToCalc,nodesModelToCalc_a,linksModel,yearsToCalc,transfer_techs,commodity)
    $transfer_losses(%selscen%timeModelToCalc,nodesModelToCalc,nodesModelToCalc_a,linksModel,yearsToCalc,transfer_techs,commodity)
    = round(transfer_losses(%selscen%timeModelToCalc,nodesModelToCalc,nodesModelToCalc_a,linksModel,yearsToCalc,transfer_techs,commodity), 6);

storage_flows(%selscen%timeModelToCalc,accNodesModel,accYears,storage_techs,commodity)
    $storage_flows(%selscen%timeModelToCalc,accNodesModel,accYears,storage_techs,commodity)
    = round(storage_flows(%selscen%timeModelToCalc,accNodesModel,accYears,storage_techs,commodity), 6);

storage_level_out(%selscen%timeModelToCalc,accNodesModel,accYears,storage_techs,commodity)
    $storage_level_out(%selscen%timeModelToCalc,accNodesModel,accYears,storage_techs,commodity)
    = round(storage_level_out(%selscen%timeModelToCalc,accNodesModel,accYears,storage_techs,commodity), 6);

marginals_balance(%selscen%timeModelToCalc,nodesModelToCalc,yearsToCalc,commodity)
    $marginals_balance(%selscen%timeModelToCalc,nodesModelToCalc,yearsToCalc,commodity)
    = round(marginals_balance(%selscen%timeModelToCalc,nodesModelToCalc,yearsToCalc,commodity), 6);

marginals_sourcesink_profile(%selscen%timeModelToCalc,nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity)
    $marginals_sourcesink_profile(%selscen%timeModelToCalc,nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity)
    = round(marginals_sourcesink_profile(%selscen%timeModelToCalc,nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity), 6);
$endif.roundts

$endif.run_postcalc
