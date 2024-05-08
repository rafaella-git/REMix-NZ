* SPDX-FileCopyrightText: Copyright (c) 2023-2024 German Aerospace Center (DLR)
* SPDX-License-Identifier: BSD-3-Clause

$onVerbatim
$iftheni.method %method%==iternodes

$ifi not "%gams.optdir%"=="" $if not dexist "%gams.optdir%" put_utility 'exec' / 'mkdir -p %gams.optdir%'

* ==== global options ====
$if not set solver     $set solver            cplex
$if not set optfile    $set optfile           1
$if not set crossover  $set crossover         0
$if not set threads    $set threads           4
$if not set submits    $set submits           8
$if not set accuracy   $set accuracy          1e-6
$if not set names      $set names             0
$if not set rerun      $set rerun             0
$if not set pathopt    $set pathopt           foresight
$if not set epgap      $set epgap             1e-3

* ==== setup optimization ====
if ((sum(nodesModelToCalc, 1)<40 and sum(timeModelToCalc, 1)<25),
   option limrow=100000, limcol=100000;
else
   option limrow=0, limcol=0, solprint=off;
);

$setenv GDXCOMPRESS 1

option LP = %solver%;
option MIP = %solver%;
option reslim = 259200;
option optcr = %epgap%;
remix.optfile = %optfile%;


* ==== configure optionfiles ====

$iftheni.solver %solver%==cplex
file opt / "%gams.optdir%cplex.opt" /;
put opt;
$ife %debug%=0 $ife %names%=0 put "names no" /;
$ife %rerun%=0 put "rerun no" /;
put "threads %threads%" /;
put "epgap %epgap%" /;
put "lpmethod 4" /;
put "barepcomp %ACCURACY%" /;
$if set cpumask put "cpumask %cpumask%" /;
putclose;

$else.solver
$abort "No valid solver specified. Please use the CPLEX solver for the iternodes method."

$endif.solver
$offVerbatim

* ==== solve the problem ====

* mapping from optimization frame to years
parameter h(nodesModel);
h(nodesModel) = 0;
parameter sol(nodesModel);
sol(nodesModel) = 0;
scalar submit;
submit = 0;

remix.solveLink = 3;

loop ( optiframeToCalc,
    yearsSel(years) = no;
    yearsSel(years)$map_optiframe(optiframeToCalc,years) = yes;
    yearsToFix(years) = no;
    yearsToFix(years)$(years.val < smin(years_a$yearsSel(years_a), years_a.val)) = yes;
    accYearsSel(accYears) = no;
    accYearsSel("horizon") = yes;
    accYearsSel(accYears)$(sum(yearsSel$sameas(accYears,yearsSel), 1)) = yes;
    accYearsToFix(accYears) = no;
    accYearsToFix(accYears)$(sum(years$(sameas(years,accYears) and years.val < smin(years_a$yearsSel(years_a), years_a.val)), 1) > 0) = yes;
    timeModelSel(timeModel) = no;
    timeModelSel(timeModel)$timeModelToCalc(timeModel) = yes;

* Fix decision for years previously optimized in case of myopic or foresight
    converter_unitsBuild.fx(nodesModelToCalc,yearsToFix,converter_techs,vintage)
        = converter_unitsBuild.l(nodesModelToCalc,yearsToFix,converter_techs,vintage);
    converter_unitsDecom.fx(nodesModelToCalc,yearsToFix,converter_techs,vintage)
        = converter_unitsDecom.l(nodesModelToCalc,yearsToFix,converter_techs,vintage);

    storage_unitsBuild.fx(nodesModelToCalc,yearsToFix,storage_techs,vintage)
        = storage_unitsBuild.l(nodesModelToCalc,yearsToFix,storage_techs,vintage);
    storage_unitsDecom.fx(nodesModelToCalc,yearsToFix,storage_techs,vintage)
        = storage_unitsDecom.l(nodesModelToCalc,yearsToFix,storage_techs,vintage);

    transfer_linksBuild.fx(linksModelToCalc,yearsToFix,transfer_techs,vintage)
        = transfer_linksBuild.l(linksModelToCalc,yearsToFix,transfer_techs,vintage);
    transfer_linksDecom.fx(linksModelToCalc,yearsToFix,transfer_techs,vintage)
        = transfer_linksDecom.l(linksModelToCalc,yearsToFix,transfer_techs,vintage);

    accounting_indicator.fx(accNodesModel,accYearsToFix,indicator)
        = accounting_indicator.l(accNodesModel,accYearsToFix,indicator);

    repeat
        submit = 0;
        loop(nodesModel$(nodesModelToCalc(nodesModel) and h(nodesModel) = 0 and sol(nodesModel) = 0),
            if ((sum(nodesModel_a$h(nodesModel_a), 1) < %submits% and submit < 2),
                nodesModelSel(nodesModel)$nodesModelToCalc(nodesModel) = yes;

                if (opti_sense < 0,
                solve remix minimizing accounting_objective using mip;
                else
                solve remix maximizing accounting_objective using mip;
                );

                nodesModelSel(nodesModel) = no;
                h(nodesModel) = remix.handle;
                submit = submit + 1;
            );
        );

        loop(nodesModel$handleCollect(h(nodesModel)),
            display$handleDelete(h(nodesModel)) 'trouble deleting handles' ;
            h(nodesModel) = 0;
            sol(nodesModel) = 1;
        );
    until sum(nodesModel_a$sol(nodesModel_a), 1) = card(nodesModelToCalc);
);

$onVerbatim
$endif.method
$offVerbatim
