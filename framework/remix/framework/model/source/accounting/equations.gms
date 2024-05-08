* SPDX-FileCopyrightText: Copyright (c) 2023-2024 German Aerospace Center (DLR)
* SPDX-License-Identifier: BSD-3-Clause

* // # accounting_equations

* ==== declaration of variables ====

variables
accounting_indicator(accNodesModel,accYears,indicator)
accounting_indicator_links(linksModel,years,indicator)
accounting_objective
  ;


* ==== definition of variables ====

* set the variable levels to be fixed for years before the optimization years
accounting_indicator.l(accNodesModel,accYears,indicator)
    $activeIndicators(accNodesModel,accYears,indicator)
    = 0;

accounting_indicator.lo(accNodesModel,accYears,indicator)
    $(accounting_indicatorBounds(accNodesModel,accYears,indicator,"useLower") <> 0 )
    = accounting_indicatorBounds(accNodesModel,accYears,indicator,"lowerValue");

accounting_indicator.up(accNodesModel,accYears,indicator)
    $(accounting_indicatorBounds(accNodesModel,accYears,indicator,"useUpper") <> 0 )
    = accounting_indicatorBounds(accNodesModel,accYears,indicator,"upperValue");

accounting_indicator.fx(accNodesModel,accYears,indicator)
    $(accounting_indicatorBounds(accNodesModel,accYears,indicator,"useFixed") <> 0 )
    = accounting_indicatorBounds(accNodesModel,accYears,indicator,"fixedValue");


accounting_indicator_links.lo(linksModelToCalc,yearsToCalc,indicator)
    $(accounting_indicatorBounds_links(linksModelToCalc,yearsToCalc,indicator,"useLower") <> 0 )
    = accounting_indicatorBounds_links(linksModelToCalc,yearsToCalc,indicator,"lowerValue");

accounting_indicator_links.up(linksModelToCalc,yearsToCalc,indicator)
    $(accounting_indicatorBounds_links(linksModelToCalc,yearsToCalc,indicator,"useUpper") <> 0 )
    = accounting_indicatorBounds_links(linksModelToCalc,yearsToCalc,indicator,"upperValue");

accounting_indicator_links.fx(linksModelToCalc,yearsToCalc,indicator)
    $(accounting_indicatorBounds_links(linksModelToCalc,yearsToCalc,indicator,"useFixed") <> 0 )
    = accounting_indicatorBounds_links(linksModelToCalc,yearsToCalc,indicator,"fixedValue");


* ==== declaration of equations ====

equations
Eq_accounting_indicatorCalc(accNodesModel,accYears,indicator
    ) "Calculates the level of an indicator per accounting region"
Eq_accounting_indicatorCalc_links(linksModel,years,indicator
    ) "Calculates the level of an indicator per model link"
Eq_accounting_objective "Calculates the objective value based on the specified indicator"
  ;


* ==== equations definition ====
* // ## Equations
* // ### Accounting Indicator Calculation
* // Calculates the indicators for each model node for converters, sources and sinks, transfer, storage and variable indicators.
* // {Eq_accounting_indicatorCalc}
Eq_accounting_indicatorCalc(accNodesModel,accYearsSel(accYears),indicator)
    $activeIndicators(accNodesModel,accYears,indicator)
    ..
    accounting_indicator(accNodesModel,accYears,indicator)
    =e=

* == variable indicators ==
    sum((accNodesModel_a,accYears_a,indicator_a)
        $(compoundIndicators(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a)
            and variableIndicators(accNodesModel_a,accYears_a,indicator_a)),
        compoundIndicators(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a)
        * accounting_indicator(accNodesModel_a,accYears_a,indicator_a))

* == converters ==
    + sum ((accNodesModel_a,nodesModelSel,accYears_a,yearsSel,indicator_a)
            $( compoundIndicators(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a)
                and sameas(nodesModelSel,accNodesModel_a) and sameas(yearsSel,accYears_a)),
        compoundIndicators(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a)
        *
        ( sum ((converter_techs,vintage)
                    $(converter_availTech(nodesModelSel,yearsSel,converter_techs,vintage)
                        and accounting_converterUnits(indicator_a,nodesModelSel,converter_techs,vintage,"useAnnuity") = 0),
            converter_unitsBuild(nodesModelSel,yearsSel,converter_techs,vintage)
            * accounting_converterUnits(indicator_a,nodesModelSel,converter_techs,vintage,"perUnitBuild"))

        + sum ((years_a,converter_techs,vintage)
                    $(converter_availTech(nodesModelSel,years_a,converter_techs,vintage)
                        and years_a.val < sum(yearsToCalc_a$(ord(yearsToCalc_a)=1), yearsToCalc_a.val)
                        and accounting_converterUnits(indicator_a,nodesModelSel,converter_techs,vintage,"useAnnuity") = 1
                        and years_a.val + accounting_converterUnits(indicator_a,nodesModelSel,converter_techs,vintage,"amorTime") > yearsSel.val
                        and years_a.val <= yearsSel.val ),
            converter_unitsBuild(nodesModelSel,years_a,converter_techs,vintage)
            * accounting_converterUnits(indicator_a,nodesModelSel,converter_techs,vintage,"perUnitBuild")
            * accounting_annuityFactor_converter(indicator_a,nodesModelSel,converter_techs,vintage) )

        + sum ((yearsToCalc,converter_techs,vintage)
                    $(converter_availTech(nodesModelSel,yearsToCalc,converter_techs,vintage)
                        and yearsToCalc.val >= sum(yearsToCalc_a$(ord(yearsToCalc_a)=1), yearsToCalc_a.val)
                        and accounting_converterUnits(indicator_a,nodesModelSel,converter_techs,vintage,"useAnnuity") = 1
                        and yearsToCalc.val + accounting_converterUnits(indicator_a,nodesModelSel,converter_techs,vintage,"amorTime") > yearsSel.val
                        and yearsToCalc.val <= yearsSel.val ),
            converter_unitsBuild(nodesModelSel,yearsToCalc,converter_techs,vintage)
            * accounting_converterUnits(indicator_a,nodesModelSel,converter_techs,vintage,"perUnitBuild")
            * accounting_annuityFactor_converter(indicator_a,nodesModelSel,converter_techs,vintage) )

        + sum ((converter_techs,vintage)
                    $converter_usedTech(nodesModelSel,yearsSel,converter_techs,vintage),
            converter_unitsDecom(nodesModelSel,yearsSel,converter_techs,vintage)
            * accounting_converterUnits(indicator_a,nodesModelSel,converter_techs,vintage,"perUnitDecom")

            + converter_unitsTotal(nodesModelSel,yearsSel,converter_techs,vintage)
            * accounting_converterUnits(indicator_a,nodesModelSel,converter_techs,vintage,"perUnitTotal") )

        + sum ((timeModelSel,converter_techs,vintage,activity)
                    $converter_usedTechAct(nodesModelSel,yearsSel,converter_techs,vintage,activity),
            converter_activity(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage,activity)
            * timeLength(timeModelSel)
            * accounting_converterActivity(indicator_a,nodesModelSel,converter_techs,vintage,activity,"perActivity") )

        + sum ((timeModelSel,converter_techs,vintage)
                    $converter_usedTech(nodesModelSel,yearsSel,converter_techs,vintage),
            converter_unitStartups(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage)
            * accounting_converterStartup(indicator_a,nodesModelSel,converter_techs,vintage,"perStartup") )

        + sum ((timeModelSel,converter_techs,vintage)
                    $converter_usedTech(nodesModelSel,yearsSel,converter_techs,vintage),
            converter_rampPos(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage)
            * (accounting_converterStartup(indicator_a,nodesModelSel,converter_techs,vintage,"perRamp")
               + accounting_converterStartup(indicator_a,nodesModelSel,converter_techs,vintage,"perRampPos"))

            + converter_rampNeg(timeModelSel,nodesModelSel,yearsSel,converter_techs,vintage)
            * (accounting_converterStartup(indicator_a,nodesModelSel,converter_techs,vintage,"perRamp")
               + accounting_converterStartup(indicator_a,nodesModelSel,converter_techs,vintage,"perRampNeg")))
        )
    )

* == storage ==
    + sum ((accNodesModel_a,nodesModelSel,accYears_a,yearsSel,indicator_a)
            $( compoundIndicators(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a)
                and sameas(nodesModelSel,accNodesModel_a) and sameas(yearsSel,accYears_a)),
        compoundIndicators(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a)
        *
        ( sum ((storage_techs,vintage)
                    $(storage_availTech(nodesModelSel,yearsSel,storage_techs,vintage)
                        and accounting_storageUnits(indicator_a,nodesModelSel,storage_techs,vintage,"useAnnuity") = 0),
            storage_unitsBuild(nodesModelSel,yearsSel,storage_techs,vintage)
            * accounting_storageUnits(indicator_a,nodesModelSel,storage_techs,vintage,"perUnitBuild") )

        + sum ((years_a,storage_techs,vintage)
                    $(storage_availTech(nodesModelSel,years_a,storage_techs,vintage)
                        and years_a.val < sum(yearsToCalc_a$(ord(yearsToCalc_a)=1), yearsToCalc_a.val)
                        and accounting_storageUnits(indicator_a,nodesModelSel,storage_techs,vintage,"useAnnuity") = 1
                        and years_a.val + accounting_storageUnits(indicator_a,nodesModelSel,storage_techs,vintage,"amorTime") > yearsSel.val
                        and years_a.val <= yearsSel.val ),
            storage_unitsBuild(nodesModelSel,years_a,storage_techs,vintage)
            * accounting_storageUnits(indicator_a,nodesModelSel,storage_techs,vintage,"perUnitBuild")
            * accounting_annuityFactor_storage(indicator_a,nodesModelSel,storage_techs,vintage) )

        + sum ((yearsToCalc,storage_techs,vintage)
                    $(storage_availTech(nodesModelSel,yearsToCalc,storage_techs,vintage)
                        and yearsToCalc.val >= sum(yearsToCalc_a$(ord(yearsToCalc_a)=1), yearsToCalc_a.val)
                        and accounting_storageUnits(indicator_a,nodesModelSel,storage_techs,vintage,"useAnnuity") = 1
                        and yearsToCalc.val + accounting_storageUnits(indicator_a,nodesModelSel,storage_techs,vintage,"amorTime") > yearsSel.val
                        and yearsToCalc.val <= yearsSel.val ),
            storage_unitsBuild(nodesModelSel,yearsToCalc,storage_techs,vintage)
            * accounting_storageUnits(indicator_a,nodesModelSel,storage_techs,vintage,"perUnitBuild")
            * accounting_annuityFactor_storage(indicator_a,nodesModelSel,storage_techs,vintage) )

        + sum ((storage_techs,vintage)
                    $storage_usedTech(nodesModelSel,yearsSel,storage_techs,vintage),
            storage_unitsDecom(nodesModelSel,yearsSel,storage_techs,vintage)
            * accounting_storageUnits(indicator_a,nodesModelSel,storage_techs,vintage,"perUnitDecom")

            + storage_unitsTotal(nodesModelSel,yearsSel,storage_techs,vintage)
            * accounting_storageUnits(indicator_a,nodesModelSel,storage_techs,vintage,"perUnitTotal") )
        )
    )


* == transfer ==
    + sum ((accNodesModel_a,nodesModelSel,accYears_a,yearsSel,indicator_a)
            $( compoundIndicators(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a)
                and sameas(nodesModelSel,accNodesModel_a) and sameas(yearsSel,accYears_a)),
        compoundIndicators(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a)
        *
        ( sum ((linksModelToCalc,transfer_techs,vintage)
                    $(transfer_availTech(linksModelToCalc,yearsSel,transfer_techs,vintage)
                        and transfer_incidenceModel(nodesModelSel,linksModelToCalc) <> 0
                        and accounting_transferLinks(indicator_a,linksModelToCalc,transfer_techs,vintage,"useAnnuity") = 0),
            0.5
            * transfer_linksBuild(linksModelToCalc,yearsSel,transfer_techs,vintage)
            * accounting_transferLinks(indicator_a,linksModelToCalc,transfer_techs,vintage,"perLinkBuild") )

        + sum ((linksModelToCalc,years_a,transfer_techs,vintage)
                    $(transfer_availTech(linksModelToCalc,years_a,transfer_techs,vintage)
                        and years_a.val < sum(yearsToCalc_a$(ord(yearsToCalc_a)=1), yearsToCalc_a.val)
                        and transfer_incidenceModel(nodesModelSel,linksModelToCalc) <> 0
                        and accounting_transferLinks(indicator_a,linksModelToCalc,transfer_techs,vintage,"useAnnuity") = 1
                        and years_a.val + accounting_transferLinks(indicator_a,linksModelToCalc,transfer_techs,vintage,"amorTime") > yearsSel.val
                        and years_a.val <= yearsSel.val ),
            0.5
            * transfer_linksBuild(linksModelToCalc,years_a,transfer_techs,vintage)
            * accounting_transferLinks(indicator_a,linksModelToCalc,transfer_techs,vintage,"perLinkBuild")
            * accounting_annuityFactor_transferLink(indicator_a,linksModelToCalc,transfer_techs,vintage) )

        + sum ((linksModelToCalc,yearsToCalc,transfer_techs,vintage)
                    $(transfer_availTech(linksModelToCalc,yearsToCalc,transfer_techs,vintage)
                        and yearsToCalc.val >= sum(yearsToCalc_a$(ord(yearsToCalc_a)=1), yearsToCalc_a.val)
                        and transfer_incidenceModel(nodesModelSel,linksModelToCalc) <> 0
                        and accounting_transferLinks(indicator_a,linksModelToCalc,transfer_techs,vintage,"useAnnuity") = 1
                        and yearsToCalc.val + accounting_transferLinks(indicator_a,linksModelToCalc,transfer_techs,vintage,"amorTime") > yearsSel.val
                        and yearsToCalc.val <= yearsSel.val ),
            0.5
            * transfer_linksBuild(linksModelToCalc,yearsToCalc,transfer_techs,vintage)
            * accounting_transferLinks(indicator_a,linksModelToCalc,transfer_techs,vintage,"perLinkBuild")
            * accounting_annuityFactor_transferLink(indicator_a,linksModelToCalc,transfer_techs,vintage) )

        + sum ((linksModelToCalc,transfer_techs,vintage,link_types)
                    $(transfer_availTech(linksModelToCalc,yearsSel,transfer_techs,vintage)
                        and transfer_incidenceModel(nodesModelSel,linksModelToCalc) <> 0
                        and accounting_transferPerLength(indicator_a,linksModelToCalc,transfer_techs,vintage,link_types,"useAnnuity") = 0 ),
            0.5
            * transfer_linksBuild(linksModelToCalc,yearsSel,transfer_techs,vintage)
            * transfer_lengthParam(linksModelToCalc,link_types,"length")
            * accounting_transferPerLength(indicator_a,linksModelToCalc,transfer_techs,vintage,link_types,"perLengthBuild") )

        + sum ((linksModelToCalc,years_a,transfer_techs,vintage,link_types)
                    $(transfer_availTech(linksModelToCalc,years_a,transfer_techs,vintage)
                        and years_a.val < sum(yearsToCalc_a$(ord(yearsToCalc_a)=1), yearsToCalc_a.val)
                        and transfer_incidenceModel(nodesModelSel,linksModelToCalc) <> 0
                        and accounting_transferPerLength(indicator_a,linksModelToCalc,transfer_techs,vintage,link_types,"useAnnuity") = 1
                        and years_a.val + accounting_transferPerLength(indicator_a,linksModelToCalc,transfer_techs,vintage,link_types,"amorTime") > yearsSel.val
                        and years_a.val <= yearsSel.val ),
            0.5
            * transfer_linksBuild(linksModelToCalc,years_a,transfer_techs,vintage)
            * transfer_lengthParam(linksModelToCalc,link_types,"length")
            * accounting_transferPerLength(indicator_a,linksModelToCalc,transfer_techs,vintage,link_types,"perLengthBuild")
            * accounting_annuityFactor_transferPerLength(indicator_a,linksModelToCalc,transfer_techs,vintage,link_types) )

        + sum ((linksModelToCalc,yearsToCalc,transfer_techs,vintage,link_types)
                    $(transfer_availTech(linksModelToCalc,yearsToCalc,transfer_techs,vintage)
                        and yearsToCalc.val >= sum(yearsToCalc_a$(ord(yearsToCalc_a)=1), yearsToCalc_a.val)
                        and transfer_incidenceModel(nodesModelSel,linksModelToCalc) <> 0
                        and accounting_transferPerLength(indicator_a,linksModelToCalc,transfer_techs,vintage,link_types,"useAnnuity") = 1
                        and yearsToCalc.val + accounting_transferPerLength(indicator_a,linksModelToCalc,transfer_techs,vintage,link_types,"amorTime") > yearsSel.val
                        and yearsToCalc.val <= yearsSel.val ),
            0.5
            * transfer_linksBuild(linksModelToCalc,yearsToCalc,transfer_techs,vintage)
            * transfer_lengthParam(linksModelToCalc,link_types,"length")
            * accounting_transferPerLength(indicator_a,linksModelToCalc,transfer_techs,vintage,link_types,"perLengthBuild")
            * accounting_annuityFactor_transferPerLength(indicator_a,linksModelToCalc,transfer_techs,vintage,link_types) )

        + sum ((linksModelToCalc,transfer_techs,vintage)
                    $(transfer_usedTech(linksModelToCalc,yearsSel,transfer_techs,vintage)
                        and transfer_incidenceModel(nodesModelSel,linksModelToCalc) <> 0 ),
            0.5
            * transfer_linksDecom(linksModelToCalc,yearsSel,transfer_techs,vintage)
            * accounting_transferLinks(indicator_a,linksModelToCalc,transfer_techs,vintage,"perLinkDecom")

            + 0.5
            * transfer_linksTotal(linksModelToCalc,yearsSel,transfer_techs,vintage)
            * accounting_transferLinks(indicator_a,linksModelToCalc,transfer_techs,vintage,"perLinkTotal")

            + 0.5
            * sum (link_types,
                transfer_linksDecom(linksModelToCalc,yearsSel,transfer_techs,vintage)
                * transfer_lengthParam(linksModelToCalc,link_types,"length")
                * accounting_transferPerLength(indicator_a,linksModelToCalc,transfer_techs,vintage,link_types,"perLengthDecom")

                + transfer_linksTotal(linksModelToCalc,yearsSel,transfer_techs,vintage)
                * transfer_lengthParam(linksModelToCalc,link_types,"length")
                * accounting_transferPerLength(indicator_a,linksModelToCalc,transfer_techs,vintage,link_types,"perLengthTotal"))

            + 0.5
            * sum (timeModelSel,
                transfer_flowAlong(timeModelSel,linksModelToCalc,yearsSel,transfer_techs,vintage)
                * timeLength(timeModelSel)
                * ( accounting_transferLinks(indicator_a,linksModelToCalc,transfer_techs,vintage,"perFlow")
                    + accounting_transferLinks(indicator_a,linksModelToCalc,transfer_techs,vintage,"perFlowAlong"))

                + transfer_flowAgainst(timeModelSel,linksModelToCalc,yearsSel,transfer_techs,vintage)
                * timeLength(timeModelSel)
                * ( accounting_transferLinks(indicator_a,linksModelToCalc,transfer_techs,vintage,"perFlow")
                    + accounting_transferLinks(indicator_a,linksModelToCalc,transfer_techs,vintage,"perFlowAgainst")))

            + 0.5
            * sum ((timeModelSel, link_types),
                transfer_flowAlong(timeModelSel,linksModelToCalc,yearsSel,transfer_techs,vintage)
                * timeLength(timeModelSel)
                * transfer_lengthParam(linksModelToCalc,link_types,"length")
                * (accounting_transferPerLength(indicator_a,linksModelToCalc,transfer_techs,vintage,link_types,"perFlow")
                    + accounting_transferPerLength(indicator_a,linksModelToCalc,transfer_techs,vintage,link_types,"perFlowAlong"))

                + transfer_flowAgainst(timeModelSel,linksModelToCalc,yearsSel,transfer_techs,vintage)
                * timeLength(timeModelSel)
                * transfer_lengthParam(linksModelToCalc,link_types,"length")
                * (accounting_transferPerLength(indicator_a,linksModelToCalc,transfer_techs,vintage,link_types,"perFlow")
                    + accounting_transferPerLength(indicator_a,linksModelToCalc,transfer_techs,vintage,link_types,"perFlowAgainst")))
            )
        )
    )


* == sources / sinks ==
    + sum ((accNodesModel_a,nodesModelSel,accYears_a,yearsSel,indicator_a)
            $( compoundIndicators(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a)
                and sameas(nodesModelSel,accNodesModel_a) and sameas(yearsSel,accYears_a)),
        compoundIndicators(accNodesModel,accYears,indicator,accNodesModel_a,accYears_a,indicator_a)
        *
        sum ((timeModelSel,sourcesink_techs,commodity)
                $sourcesink_enabled(nodesModelSel,yearsSel,sourcesink_techs,commodity),
            sourcesink_flow(timeModelSel,nodesModelSel,yearsSel,sourcesink_techs,commodity)
            * timeLength(timeModelSel)
            * accounting_sourcesinkFlow(indicator_a,nodesModelSel,yearsSel,sourcesink_techs,commodity,"perFlow") )
    );

* // ### Accounting Indicator Calculation Links
* // Calculates the indicators for each transfer for converters, sources and sinks, transfer, storage and variable indicators.
* // {Eq_accounting_indicatorCalc_links}
Eq_accounting_indicatorCalc_links(linksModelToCalc,yearsSel,indicator)
    $activeIndicators_links(linksModelToCalc,yearsSel,indicator)
    ..
    accounting_indicator_links(linksModelToCalc,yearsSel,indicator)
    =e=
    sum ((transfer_techs,vintage)
                $(transfer_availTech(linksModelToCalc,yearsSel,transfer_techs,vintage)
                    and accounting_transferLinks(indicator,linksModelToCalc,transfer_techs,vintage,"useAnnuity") = 0),
        transfer_linksBuild(linksModelToCalc,yearsSel,transfer_techs,vintage)
        * accounting_transferLinks(indicator,linksModelToCalc,transfer_techs,vintage,"perLinkBuild") )

    + sum ((years_a,transfer_techs,vintage)
                $(transfer_availTech(linksModelToCalc,yearsSel,transfer_techs,vintage)
                    and accounting_transferLinks(indicator,linksModelToCalc,transfer_techs,vintage,"useAnnuity") = 1
                    and years_a.val + accounting_transferLinks(indicator,linksModelToCalc,transfer_techs,vintage,"amorTime") > yearsSel.val
                    and years_a.val <= yearsSel.val ),
        transfer_linksBuild(linksModelToCalc,years_a,transfer_techs,vintage)
        * accounting_transferLinks(indicator,linksModelToCalc,transfer_techs,vintage,"perLinkBuild")
        * accounting_annuityFactor_transferLink(indicator,linksModelToCalc,transfer_techs,vintage) )

    + sum ((transfer_techs,vintage,link_types)
                $(transfer_availTech(linksModelToCalc,yearsSel,transfer_techs,vintage)
                    and accounting_transferPerLength(indicator,linksModelToCalc,transfer_techs,vintage,link_types,"useAnnuity") = 0 ),
        transfer_linksBuild(linksModelToCalc,yearsSel,transfer_techs,vintage)
        * transfer_lengthParam(linksModelToCalc,link_types,"length")
        * accounting_transferPerLength(indicator,linksModelToCalc,transfer_techs,vintage,link_types,"perLengthBuild") )

    + sum ((years_a,transfer_techs,vintage,link_types)
                $(transfer_availTech(linksModelToCalc,yearsSel,transfer_techs,vintage)
                    and accounting_transferPerLength(indicator,linksModelToCalc,transfer_techs,vintage,link_types,"useAnnuity") = 1
                    and years_a.val + accounting_transferPerLength(indicator,linksModelToCalc,transfer_techs,vintage,link_types,"amorTime") > yearsSel.val
                    and years_a.val <= yearsSel.val ),
        transfer_linksBuild(linksModelToCalc,years_a,transfer_techs,vintage)
        * transfer_lengthParam(linksModelToCalc,link_types,"length")
        * accounting_transferPerLength(indicator,linksModelToCalc,transfer_techs,vintage,link_types,"perLengthBuild")
        * accounting_annuityFactor_transferPerLength(indicator,linksModelToCalc,transfer_techs,vintage,link_types) )

    + sum ((transfer_techs,vintage)
                $(transfer_usedTech(linksModelToCalc,yearsSel,transfer_techs,vintage)),
        transfer_linksDecom(linksModelToCalc,yearsSel,transfer_techs,vintage)
        * accounting_transferLinks(indicator,linksModelToCalc,transfer_techs,vintage,"perLinkDecom")

        + transfer_linksTotal(linksModelToCalc,yearsSel,transfer_techs,vintage)
        * accounting_transferLinks(indicator,linksModelToCalc,transfer_techs,vintage,"perLinkTotal")

        + sum (link_types,
            transfer_linksDecom(linksModelToCalc,yearsSel,transfer_techs,vintage)
            * transfer_lengthParam(linksModelToCalc,link_types,"length")
            * accounting_transferPerLength(indicator,linksModelToCalc,transfer_techs,vintage,link_types,"perLengthDecom")

            + transfer_linksTotal(linksModelToCalc,yearsSel,transfer_techs,vintage)
            * transfer_lengthParam(linksModelToCalc,link_types,"length")
            * accounting_transferPerLength(indicator,linksModelToCalc,transfer_techs,vintage,link_types,"perLengthTotal"))

        + sum (timeModelSel,
            transfer_flowAlong(timeModelSel,linksModelToCalc,yearsSel,transfer_techs,vintage)
            * timeLength(timeModelSel)
            * ( accounting_transferLinks(indicator,linksModelToCalc,transfer_techs,vintage,"perFlow")
                + accounting_transferLinks(indicator,linksModelToCalc,transfer_techs,vintage,"perFlowAlong"))

            + transfer_flowAgainst(timeModelSel,linksModelToCalc,yearsSel,transfer_techs,vintage)
            * timeLength(timeModelSel)
            * ( accounting_transferLinks(indicator,linksModelToCalc,transfer_techs,vintage,"perFlow")
                + accounting_transferLinks(indicator,linksModelToCalc,transfer_techs,vintage,"perFlowAgainst")))

        + sum ((timeModelSel, link_types),
            transfer_flowAlong(timeModelSel,linksModelToCalc,yearsSel,transfer_techs,vintage)
            * timeLength(timeModelSel)
            * transfer_lengthParam(linksModelToCalc,link_types,"length")
            * (accounting_transferPerLength(indicator,linksModelToCalc,transfer_techs,vintage,link_types,"perFlow")
                + accounting_transferPerLength(indicator,linksModelToCalc,transfer_techs,vintage,link_types,"perFlowAlong"))

            + transfer_flowAgainst(timeModelSel,linksModelToCalc,yearsSel,transfer_techs,vintage)
            * timeLength(timeModelSel)
            * transfer_lengthParam(linksModelToCalc,link_types,"length")
            * (accounting_transferPerLength(indicator,linksModelToCalc,transfer_techs,vintage,link_types,"perFlow")
                + accounting_transferPerLength(indicator,linksModelToCalc,transfer_techs,vintage,link_types,"perFlowAgainst")))
        );

* // ### Accounting Objective
* // Calculates the indicators for the objective.
* // {Eq_accounting_objective}
Eq_accounting_objective
    ..
    accounting_objective
    =e=
    sum ((accNodesModel,accYears,indicator)
            $(accounting_indicatorBounds(accNodesModel,accYears,indicator,"obj") <> 0 ),
        accounting_indicator(accNodesModel,accYears,indicator) )


* ==== model definition ====

Model M_accounting
/
Eq_accounting_indicatorCalc_links
Eq_accounting_indicatorCalc
Eq_accounting_objective
/;
