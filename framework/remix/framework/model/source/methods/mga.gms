* SPDX-FileCopyrightText: Copyright (c) 2023-2024 German Aerospace Center (DLR)
* SPDX-License-Identifier: BSD-3-Clause

$onVerbatim
$iftheni.method %method%==mga

$ifi not "%gams.optdir%"=="" $if not dexist "%gams.optdir%" put_utility 'exec' / 'mkdir -p %gams.optdir%'

* ==== global options ====
$if not set solver        $set solver            cplex
$if not set holdfixed     $set holdfixed         1
$if not set optfile       $set optfile           1
$if not set crossover     $set crossover         0
$if not set threads       $set threads           8
$if not set accuracy      $set accuracy          1e-6
$if not set names         $set names             0
$if not set rerun         $set rerun             0
$if not set pathopt       $set pathopt           foresight


if ((sum(nodesModelToCalc, 1)<40 and sum(timeModelToCalc, 1)<25),
   option limrow=100000, limcol=100000;
else
   option limrow=0, limcol=0, solprint=off;
);

$setenv GDXCOMPRESS 1

option MIP = %solver%;
option reslim = 259200;
option optcr = 1e-2;
remix.optfile = %optfile%;


* ==== configure optionfiles ====

$iftheni.solver %solver%==cplex
file opt / "%gams.optdir%cplex.opt" /;
put opt;
$ife %debug%=0 $ife %names%=0 put "names no" /;
put "lpmethod 4" /;
put "rerun no" /;
$ife %crossover%=0 put "solutiontype 2" /;
put "threads %THREADS%" /;
put "barepcomp %ACCURACY%" /;
$if set cpumask put "cpumask %cpumask%" /;
putclose;

file opt_mga / "%gams.optdir%cplex.o99" /;
put opt_mga;
$ife %debug%=0 $ife %names%=0 put "names no" /;
put "advind 0" /;
put "lpmethod 4" /;
put "startalg 4" /;
put "rerun no" /;
$ife %crossover%=0 put "solutiontype 2" /;
put "threads %THREADS%" /;
put "barepcomp 1e-3" /;
$if set cpumask put "cpumask %cpumask%" /;
putclose;

$elseifi.solver %solver%==gurobi
file opt / "%gams.optdir%gurobi.opt" /;
put opt;
$ife %debug%=0 $ife %names%=0 put "names no" /;
$ife %rerun%=1 put "rerun 1" /;
put "method 2" /;
put "mipgap 1e-3" /;
put "nonConvex 2" /;
put "threads %THREADS%" /;
put "barconvtol %ACCURACY%" /;
putclose;

file opt_mga / "%gams.optdir%gurobi.o99" /;
put opt_mga;
$ife %debug%=0 $ife %names%=0 put "names no" /;
$ife %rerun%=1 put "rerun 1" /;
put "method 2" /;
put "threads %THREADS%" /;
put "barconvtol 1e-3" /;
putclose;


$else.solver
$abort "No valid solver specified. Used solvers are CPLEX, Gurobi, XPRESS, or Convert"

$endif.solver


* ==== initial solution ====

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
    nodesModelSel(nodesModel) = no;
    nodesModelSel(nodesModel)$nodesModelToCalc(nodesModel) = yes;

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

* Optimize and log values
    remix.holdFixed = %holdfixed%;

put_utility 'log' / 'Running base optimization ';

    if (opti_sense < 0,
    solve remix minimizing accounting_objective using MIP;
    else
    solve remix maximizing accounting_objective using MIP;
    );

    put_utility 'log' / 'Model status ' remix.modelstat:0:0;
    put_utility 'log' / 'Objectiv value ' accounting_objective.l:0:3;

);

$include "%sourcedir%/postcalc/definition.gms"
mga_comp(mga)$mga_act(mga) = yes;
mga_act(mga) = no;


* ==== modify the model and solve for alternatives ====
variable mga_objective;

positive variable mga_dist(mga,accNodesModel,accYears,indicator);
positive variable mga_dist_pos(mga,accNodesModel,accYears,indicator);
positive variable mga_dist_neg(mga,accNodesModel,accYears,indicator);

parameter mga_indicatorPoints(mga,accNodesModel,accYears,indicator);
mga_indicatorPoints(mga,accNodesModel,accYears,indicator)
    $(accounting_indicatorBounds(accNodesModel,accYears,indicator,"mga")
        and ord(mga) = 1)
    = accounting_indicator.l(accNodesModel,accYears,indicator);

equation Eq_mga_limitObjective(accNodesModel,accYears,indicator);
Eq_mga_limitObjective(accNodesModel,accYears,indicator)
    $(accounting_indicatorBounds(accNodesModel,accYears,indicator,"obj") <> 0 )
    ..
    accounting_indicator(accNodesModel,accYears,indicator)
    =l=
    accounting_objective.l * %mgafactor%;


$iftheni.mgamethod %mgamethod%==linear
* LINEAR FORMULATION, PUSHES TOWARDS HIGHER INDICATORS

equation Eq_mga_obj;
Eq_mga_obj
    ..
    mga_objective
    =e=
    sum((mga_comp,accNodesModel,accYears,indicator)
            $accounting_indicatorBounds(accNodesModel,accYears,indicator,"mga"),
        mga_dist(mga_comp,accNodesModel,accYears,indicator))
    / sum((mga_comp,accNodesModel,accYears,indicator)
            $accounting_indicatorBounds(accNodesModel,accYears,indicator,"mga"),
        1);

equation Eq_mga_distance(mga,accNodesModel,accYears,indicator);
Eq_mga_distance(mga_comp,accNodesModel,accYears,indicator)
    $accounting_indicatorBounds(accNodesModel,accYears,indicator,"mga")
    ..
    mga_dist(mga_comp,accNodesModel,accYears,indicator)
    =e=
    ( accounting_indicator(accNodesModel,accYears,indicator))
    / mga_indicatorPoints(mga_comp,accNodesModel,accYears,indicator)

model remix_mga
    /
    remix
    - Eq_accounting_objective
    + Eq_mga_limitObjective
    + Eq_mga_distance
    + Eq_mga_obj
    /;

$set modeltype LP


$elseifi.mgamethod %mgamethod%==binary
* BINARY FORMULATION, DISTANCE AS ABSOLUTE VALUE
binary variable mga_dist_bool(mga,accNodesModel,accYears,indicator);

equation Eq_mga_obj;
Eq_mga_obj
    ..
    mga_objective
    =e=
    sum((mga_comp,accNodesModel,accYears,indicator)
            $accounting_indicatorBounds(accNodesModel,accYears,indicator,"mga"),
        mga_dist_pos(mga_comp,accNodesModel,accYears,indicator)
        + mga_dist_neg(mga_comp,accNodesModel,accYears,indicator))
    / sum((mga_comp,accNodesModel,accYears,indicator)
            $accounting_indicatorBounds(accNodesModel,accYears,indicator,"mga"),
        1);

equation Eq_mga_bool_pos(mga,accNodesModel,accYears,indicator);
Eq_mga_bool_pos(mga_comp,accNodesModel,accYears,indicator)
    $accounting_indicatorBounds(accNodesModel,accYears,indicator,"mga")
    ..
    mga_dist_pos(mga_comp,accNodesModel,accYears,indicator)
    =l=
    10 * mga_dist_bool(mga_comp,accNodesModel,accYears,indicator)

equation Eq_mga_bool_neg(mga,accNodesModel,accYears,indicator);
Eq_mga_bool_neg(mga_comp,accNodesModel,accYears,indicator)
    $accounting_indicatorBounds(accNodesModel,accYears,indicator,"mga")
    ..
    mga_dist_neg(mga_comp,accNodesModel,accYears,indicator)
    =l=
    10 * (1 - mga_dist_bool(mga_comp,accNodesModel,accYears,indicator))

equation Eq_mga_distance(mga,accNodesModel,accYears,indicator);
Eq_mga_distance(mga_comp,accNodesModel,accYears,indicator)
    $accounting_indicatorBounds(accNodesModel,accYears,indicator,"mga")
    ..
    mga_dist_pos(mga_comp,accNodesModel,accYears,indicator)
    - mga_dist_neg(mga_comp,accNodesModel,accYears,indicator)
    =e=
    ( accounting_indicator(accNodesModel,accYears,indicator))
    / mga_indicatorPoints(mga_comp,accNodesModel,accYears,indicator)

model remix_mga
    /
    remix
    - Eq_accounting_objective
    + Eq_mga_limitObjective
    + Eq_mga_bool_pos
    + Eq_mga_bool_neg
    + Eq_mga_distance
    + Eq_mga_obj
    /;

$set modeltype MIP


$elseifi.mgamethod %mgamethod%==quadratic
* QUADRATIC FORMULATION, DISTANCE AS SQUARE VALUE

mga_indicatorPoints(mga,accNodesModel,accYears,indicator)
    $(accounting_indicatorBounds(accNodesModel,accYears,indicator,"mga")
        and ord(mga) = 1)
    = accounting_indicator.l(accNodesModel,accYears,indicator);

equation Eq_mga_obj;
Eq_mga_obj
    ..
    mga_objective
    =e=
    sum((mga_comp,accNodesModel,accYears,indicator)
            $accounting_indicatorBounds(accNodesModel,accYears,indicator,"mga"),
        mga_dist_pos(mga_comp,accNodesModel,accYears,indicator)
        + mga_dist_neg(mga_comp,accNodesModel,accYears,indicator))
    / sum((mga_comp,indicator)
            $accounting_indicatorBounds(accNodesModel,accYears,indicator,"mga"),
        1);

equation Eq_mga_bool_pos(mga,accNodesModel,accYears,indicator);
Eq_mga_bool_pos(mga_comp,accNodesModel,accYears,indicator)
    $accounting_indicatorBounds(accNodesModel,accYears,indicator,"mga")
    ..
    mga_dist_pos(mga_comp,accNodesModel,accYears,indicator)
    =l=
    10 * mga_dist_bool(mga_comp,accNodesModel,accYears,indicator)

equation Eq_mga_bool_neg(mga,accNodesModel,accYears,indicator);
Eq_mga_bool_neg(mga_comp,accNodesModel,accYears,indicator)
    $accounting_indicatorBounds(accNodesModel,accYears,indicator,"mga")
    ..
    mga_dist_neg(mga_comp,accNodesModel,accYears,indicator)
    =l=
    10 * (1 - mga_dist_bool(mga_comp,accNodesModel,accYears,indicator))

equation Eq_mga_distance(mga,accNodesModel,accYears,indicator);
Eq_mga_distance(mga_comp,accNodesModel,accYears,indicator)
    $accounting_indicatorBounds(accNodesModel,accYears,indicator,"mga")
    ..
    mga_dist_pos(mga_comp,accNodesModel,accYears,indicator)
    - mga_dist_neg(mga_comp,accNodesModel,accYears,indicator)
    =e=
    ( accounting_indicator(accNodesModel,accYears,indicator))
    / mga_indicatorPoints(mga_comp,accNodesModel,accYears,indicator)

model remix_mga
    /
    remix
    - Eq_accounting_objective
    + Eq_mga_limitObjective
    + Eq_mga_bool_pos
    + Eq_mga_bool_neg
    + Eq_mga_distance
    + Eq_mga_obj
    /;

$set modeltype QCP


$elseifi.mgamethod %mgamethod%==hypersphere

set mga_indicators(indicator);
mga_indicators(indicator)
    $sum((accNodesModel,accYears)$accounting_indicatorBounds(accNodesModel,accYears,indicator,"mga"), 1)
 = yes;

parameter mga_weights(mga, indicator);

embeddedCode Python:

    try:
        from remix.framework.tools.mga import uniform_hypersphere
    except ImportError:
        from sys import path
        model_dir = Path(r'%sourcedir% '.strip()).parents[3].as_posix()
        if model_dir not in path:
            path.append(model_dir)
        from remix.framework.tools.mga import uniform_hypersphere


    points = list(gams.get("mga"))
    indicators = list(gams.get("mga_indicators"))

    data = [[round(i,4) for i in j] for j in uniform_hypersphere(len(indicators), len(points))]
    weights = [((p,m), data[i][j]) for i, p in enumerate(points) for j, m in enumerate(indicators)]

    gams.set("mga_weights", weights)
endEmbeddedCode mga_weights

equation Eq_mga_hyperray(accNodesModel,accYears,indicator);
Eq_mga_hyperray(accNodesModel,accYears,indicator)
    $accounting_indicatorBounds(accNodesModel,accYears,indicator,"mga")
    ..
    accounting_indicator(accNodesModel,accYears,indicator)
    =e=
    mga_indicatorPoints("mga0",accNodesModel,accYears,indicator)
    + sum(mga_act, mga_objective * mga_weights(mga_act,indicator))
    ;

model remix_mga
    /
    remix
    - Eq_accounting_objective
    + Eq_mga_limitObjective
    + Eq_mga_hyperray
    /;

$set modeltype MIP

$endif.mgamethod


* run the loop for the configured MGA method

option %modeltype% = %solver%;

remix_mga.holdFixed = %holdfixed%;
remix_mga.optfile = 99;

loop(mga$(ord(mga) > 1),
    mga_act(mga) = yes;

put_utility 'log' / 'Running MGA point ' (ord(mga)-1):0:0 ;

    solve remix_mga maximizing mga_objective using %modeltype%;

    mga_indicatorPoints(mga,accNodesModel,accYears,indicator)
    $(accounting_indicatorBounds(accNodesModel,accYears,indicator,"mga"))
    = accounting_indicator.l(accNodesModel,accYears,indicator);

$include "%sourcedir%/postcalc/definition.gms"
    mga_comp(mga) = yes;
    mga_act(mga) = no;
);

$include "%sourcedir%/postcalc/writegdx.gms"
$setglobal run_postcalc 0

$onVerbatim
$endif.method
$offVerbatim
