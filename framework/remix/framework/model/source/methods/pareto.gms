* SPDX-FileCopyrightText: Copyright (c) 2023-2024 German Aerospace Center (DLR)
* SPDX-License-Identifier: BSD-3-Clause

$onVerbatim
$iftheni.method %method%==pareto

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
put "threads %threads%" /;
put "barepcomp %accuracy%" /;
$if set cpumask put "cpumask %cpumask%" /;
putclose;

file opt_pareto / "%gams.optdir%cplex.o99" /;
put opt_pareto;
$ife %debug%=0 $ife %names%=0 put "names no" /;
* put "advind 0" /;
put "lpmethod 4" /;
put "rerun no" /;
$ife %crossover%=0 put "solutiontype 2" /;
put "threads %threads%" /;
put "barepcomp %accuracy%" /;
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
put "threads %threads%" /;
put "barconvtol %accuracy%" /;
putclose;

file opt_pareto / "%gams.optdir%gurobi.o99" /;
put opt_pareto;
$ife %debug%=0 $ife %names%=0 put "names no" /;
$ife %rerun%=1 put "rerun 1" /;
put "method 2" /;
put "threads %threads%" /;
put "barconvtol %accuracy%" /;
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


* ==== modify the model and solve for pareto points ====

* After writing the initial solution to pareto0 reset the active set and run the pareto loop
pareto_act(pareto) = no;

variable pareto_objective;

parameter pareto_points(pareto);
pareto_points(pareto) = (%paretofactor% - 1) / %paretopoints% * (ord(pareto) - 1);

equation Eq_pareto_limitObjective(accNodesModel,accYears,indicator);
Eq_pareto_limitObjective(accNodesModel,accYears,indicator)
    $(accounting_indicatorBounds(accNodesModel,accYears,indicator,"obj") <> 0 )
    ..
    accounting_indicator(accNodesModel,accYears,indicator)
    =l=
    accounting_objective.l * sum(pareto_act, 1 + pareto_points(pareto_act));

equation Eq_pareto_obj;
Eq_pareto_obj
    ..
    pareto_objective
    =e=
    sum((accNodesModel,accYears,indicator)
            $accounting_indicatorBounds(accNodesModel,accYears,indicator,"pareto"),
        accounting_indicatorBounds(accNodesModel,accYears,indicator,"pareto")
        * accounting_indicator(accNodesModel,accYears,indicator));

model remix_pareto
    /
    remix
    - Eq_accounting_objective
    + Eq_pareto_limitObjective
    + Eq_pareto_obj
    /;

$set modeltype MIP

option %modeltype% = %solver%;

remix_pareto.holdFixed = %holdfixed%;
remix_pareto.optfile = 99;

loop(pareto$(ord(pareto) > 1),
    pareto_act(pareto) = yes;

put_utility 'log' / 'Running pareto point ' (ord(pareto)-1):0:0 ;

    solve remix_pareto maximizing pareto_objective using %modeltype%;

$include "%sourcedir%/postcalc/definition.gms"
    pareto_act(pareto) = no;
);

$include "%sourcedir%/postcalc/writegdx.gms"
$setglobal run_postcalc 0

$onVerbatim
$endif.method
$offVerbatim
