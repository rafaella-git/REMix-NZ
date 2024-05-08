* SPDX-FileCopyrightText: Copyright (c) 2023-2024 German Aerospace Center (DLR)
* SPDX-License-Identifier: BSD-3-Clause

$onVerbatim
$iftheni.method %method%==solve

$ifi not "%gams.optdir%"=="" $if not dexist "%gams.optdir%" put_utility 'exec' / 'mkdir -p %gams.optdir%'

* ==== global options ====
$if not set solver          $set solver          cplex
$if not set equlist         $set equlist         0
$if not set solvelink       $set solvelink       0
$if not set optfile         $set optfile         1
$if not set holdfixed       $set holdfixed       1
$if not set pathopt         $set pathopt         foresight

* ==== general solver options ====
$if not set barrier         $set barrier         1
$if not set lpmethod        $set lpmethod        6
$if not set crossover       $set crossover       0
$if not set threads         $set threads         4
$if not set accuracy        $set accuracy        1e-6
$if not set names           $set names           0
$if not set iis             $set iis             0
$if not set epgap           $set epgap           1e-3
$if not set barorder        $set barorder        -2

* ==== cplex options ====
$if not set preind          $set preind          1
$if not set scaind          $set scaind          0
$if not set predual         $set predual         -1
$if not set parallel        $set parallel        1
$if not set baralg          $set baralg          0
$if not set barstartalg     $set barstartalg     1
$if not set barcolnz        $set barcolnz        0
$if not set rerun           $set rerun           0

* ==== gurobi options ====


* ==== copt options ====
$if not set dualize         $set dualize         -1
$if not set barhomogeneous  $set barhomogeneous  -1
$if not set presolve        $set presolve        -1

* ==== debug options ====
$ife %debug%<>0             $set barrier         0
$ife %debug%<>0             $set names           1


* ==== setup optimization ====
if ((sum(nodesModelToCalc, 1)>40 or sum(timeModelToCalc, 1)>50) and not %equlist%,
   option limRow=0, limCol=0, solPrint=off;
else
   option limRow=100000, limCol=100000, solPrint=on;
);

$setenv GDXCOMPRESS 1

option mip = %solver%;
option reslim = 1209600;
option optcr = %epgap%;
remix.threads = %threads%;
remix.optFile = %optfile%;
remix.solveLink = %solvelink%;
remix.holdFixed = %holdfixed%;


* ==== configure option files ====

$iftheni.solver %solver%==cplex
$ife %barorder%=-2 $set barorder 3

file opt / "%gams.optdir%cplex.opt" /;
put opt;
$ife %names%=0 put "names no" /;
$ife %rerun%=0 put "rerun no" /;
$ife %iis%>0 put "iis %iis%" /;
$ife %barrier%=1 put "lpmethod 4" /;
$ife %barrier%=0 put "lpmethod %lpmethod%" /;
put "barorder %barorder%" /;
put "preind %preind%" /;
put "scaind %scaind%" /;
put "predual %predual%" /;
put "baralg %baralg%" /;
put "barstartalg %barstartalg%" /;
put "barepcomp %accuracy%" /;
$ife %crossover%=0 put "solutiontype 2" /;
$if set datacheck put "datacheck %datacheck%" /;
put "startalg 4" /;
put "epgap %epgap%" /;
put "quality 1" /;
put "barcolnz %barcolnz%" /;
put "threads %threads%" /;
put "parallelmode %parallel%" /;
$if set randomseed put "randomseed %randomseed%" /;
$if set cpumask put "cpumask %cpumask%" /;
putclose;

$elseifi.solver %solver%==gurobi
$ife %barorder%=-2 $set barorder 1

file opt / "%gams.optdir%gurobi.opt" /;
put opt;
$ife %names%=0 put "names no" /;
$ife %rerun%=1 put "rerun 1" /;
$ife %barrier%=1 put "method 2" /;
$ife %iis%>0 put "iis %iis%" /;
put "barorder %barorder%" /;
put "presolve %presolve%" /;
put "mipgap %epgap%" /;
$ife %crossover%=0 put "crossover 0" /;
put "threads %threads%" /;
put "barconvtol %accuracy%" /;
putclose;

$elseifi.solver %solver%==copt
$ife %barorder%=-2 $set barorder 1

file opt / "%gams.optdir%copt.opt" /;
put opt;
$ife %barrier%=1 put 'lpmethod 2' /;
$ife %crossover%=0 put "crossover 0" /;
put "absgap 1e-4" /;
put "relgap %accuracy%" /;
put "barorder %barorder%" /;
put "presolve %presolve%" /;
put "dualize %dualize%" /;
put "barhomogeneous %barhomogeneous%" /;
put "threads %threads%" /;
putclose;

$elseifi.solver %solver%==xpress
file opt / "%gams.optdir%xpress.opt" /;
put opt;
$ife %rerun%=1 put "rerun 1" /;
$ife %barrier%=1 put "algorithm barrier" /;
$ife %crossover%=0 put "crossover 0" /;
put "threads %threads%" /;
put "barGapStop %accuracy%" /;
putclose;

$elseifi.solver %solver%==highs
file opt / "%gams.optdir%highs.opt" /;
put opt;
$ife %barrier%=1 put 'solver = ipm' /;
put 'ipm_optimality_tolerance = %accuracy%' /;
$ife %crossover%=0 put 'run_crossover = false' /;
put "threads = %threads%" /;
put "parallel = on" /;
putclose;

$elseifi.solver %solver%==convert
file opt / "%gams.optdir%convert.opt" /;
put opt;
put "CplexLP %gams.optdir%remix.lp" /;
put "CplexMPS %gams.optdir%remix.mps" /;
put "Dict %gams.optdir%dict.txt" /;
putclose;

$elseifi.solver %solver%==scip
file opt / "%gams.optdir%scip.opt" /;
put opt;
put 'gams/interactive = "write prob remix.cip quit"' /;
putclose;

$else.solver
$abort "No valid solver specified. Available solvers are CPLEX, Gurobi, XPRESS, COPT, Convert, HiGHS, or SCIP."

$endif.solver


* ==== solve the problem ====

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
    converter_unitsDelta_upper(nodesModelToCalc,yearsToFix,converter_techs)
        $(sum(yearsToCalc$sameas(yearsToFix, yearsToCalc), 1))
        = sum(vintage, converter_unitsTotal.l(nodesModelToCalc,yearsToFix,converter_techs,vintage))
            - converter_capacityParam(nodesModelToCalc,yearsToFix,converter_techs,"unitsUpperLimit");
    converter_unitsDelta_upper(nodesModelToCalc,yearsToFix,converter_techs)
        $(converter_unitsDelta_upper(nodesModelToCalc,yearsToFix,converter_techs) < 0) = 0;

    converter_unitsDelta_lower(nodesModelToCalc,yearsToCalc,converter_techs,vintage)
        = converter_unitsDecom.lo(nodesModelToCalc,yearsToCalc,converter_techs,vintage)
            - converter_unitsTotal.l(nodesModelToCalc,yearsToCalc-1,converter_techs,vintage);
    converter_unitsDelta_lower(nodesModelToCalc,yearsToCalc,converter_techs,vintage)
        $(converter_unitsDelta_lower(nodesModelToCalc,yearsToCalc,converter_techs,vintage) < 0) = 0;

    converter_unitsBuild.l(nodesModelToCalc,yearsToFix,converter_techs,vintage)
        $converter_availTech(nodesModelToCalc,yearsToFix,converter_techs,vintage)
        = converter_unitsBuild.l(nodesModelToCalc,yearsToFix,converter_techs,vintage)
            - converter_unitsDelta_upper(nodesModelToCalc,yearsToFix,converter_techs);

    converter_unitsDecom.lo(nodesModelToCalc,yearsToCalc,converter_techs,vintage)
        $converter_usedTech(nodesModelToCalc,yearsToCalc,converter_techs,vintage)
        = converter_unitsDecom.lo(nodesModelToCalc,yearsToCalc,converter_techs,vintage)
            - converter_unitsDelta_lower(nodesModelToCalc,yearsToCalc,converter_techs,vintage);

    converter_unitsBuild.l(nodesModelToCalc,yearsToFix,converter_techs,vintage)
        $(converter_unitsBuild.l(nodesModelToCalc,yearsToFix,converter_techs,vintage) < 0) = 0;
    converter_unitsBuild.fx(nodesModelToCalc,yearsToFix,converter_techs,vintage)
        = converter_unitsBuild.l(nodesModelToCalc,yearsToFix,converter_techs,vintage);
    converter_unitsDecom.l(nodesModelToCalc,yearsToFix,converter_techs,vintage)
        $(converter_unitsDecom.l(nodesModelToCalc,yearsToFix,converter_techs,vintage) < 0) = 0;
    converter_unitsDecom.fx(nodesModelToCalc,yearsToFix,converter_techs,vintage)
        = converter_unitsDecom.l(nodesModelToCalc,yearsToFix,converter_techs,vintage);
    converter_unitsTotal.fx(nodesModelToCalc,yearsToFix,converter_techs,vintage)
        = converter_unitsTotal.l(nodesModelToCalc,yearsToFix,converter_techs,vintage);


    storage_unitsDelta_upper(nodesModelToCalc,yearsToFix,storage_techs)
        $(sum(yearsToCalc$sameas(yearsToFix, yearsToCalc), 1))
        = sum(vintage, storage_unitsTotal.l(nodesModelToCalc,yearsToFix,storage_techs,vintage))
            - storage_reservoirParam(nodesModelToCalc,yearsToFix,storage_techs,"unitsUpperLimit");
    storage_unitsDelta_upper(nodesModelToCalc,yearsToFix,storage_techs)
        $(storage_unitsDelta_upper(nodesModelToCalc,yearsToFix,storage_techs) < 0) = 0;

    storage_unitsDelta_lower(nodesModelToCalc,yearsToFix,storage_techs)
        $(sum(yearsToCalc$sameas(yearsToFix, yearsToCalc), 1))
        = storage_reservoirParam(nodesModelToCalc,yearsToFix,storage_techs,"unitsLowerLimit")
            - sum(vintage, storage_unitsTotal.l(nodesModelToCalc,yearsToFix,storage_techs,vintage));
    storage_unitsDelta_lower(nodesModelToCalc,yearsToFix,storage_techs)
        $(storage_unitsDelta_lower(nodesModelToCalc,yearsToFix,storage_techs) < 0) = 0;

    storage_unitsBuild.l(nodesModelToCalc,yearsToFix,storage_techs,vintage)
        $storage_availTech(nodesModelToCalc,yearsToFix,storage_techs,vintage)
        = storage_unitsBuild.l(nodesModelToCalc,yearsToFix,storage_techs,vintage)
            - storage_unitsDelta_upper(nodesModelToCalc,yearsToFix,storage_techs);

    storage_unitsDecom.l(nodesModelToCalc,yearsToFix,storage_techs,vintage)
        $storage_usedTech(nodesModelToCalc,yearsToFix,storage_techs,vintage)
        = storage_unitsDecom.l(nodesModelToCalc,yearsToFix,storage_techs,vintage)
            - storage_unitsDelta_lower(nodesModelToCalc,yearsToFix,storage_techs);

    storage_unitsBuild.l(nodesModelToCalc,yearsToFix,storage_techs,vintage)
        $(storage_unitsBuild.l(nodesModelToCalc,yearsToFix,storage_techs,vintage) < 0) = 0;
    storage_unitsBuild.fx(nodesModelToCalc,yearsToFix,storage_techs,vintage)
        = storage_unitsBuild.l(nodesModelToCalc,yearsToFix,storage_techs,vintage);
    storage_unitsDecom.l(nodesModelToCalc,yearsToFix,storage_techs,vintage)
        $(storage_unitsDecom.l(nodesModelToCalc,yearsToFix,storage_techs,vintage) < 0) = 0;
    storage_unitsDecom.fx(nodesModelToCalc,yearsToFix,storage_techs,vintage)
        = storage_unitsDecom.l(nodesModelToCalc,yearsToFix,storage_techs,vintage);
    storage_unitsTotal.fx(nodesModelToCalc,yearsToFix,storage_techs,vintage)
        = storage_unitsTotal.l(nodesModelToCalc,yearsToFix,storage_techs,vintage);


    transfer_linksDelta_upper(linksModelToCalc,yearsToFix,transfer_techs)
        $(sum(yearsToCalc$sameas(yearsToFix, yearsToCalc), 1))
        = sum(vintage, transfer_linksTotal.l(linksModelToCalc,yearsToFix,transfer_techs,vintage))
            - transfer_linksParam(linksModelToCalc,yearsToFix,transfer_techs,"linksUpperLimit");
    transfer_linksDelta_upper(linksModelToCalc,yearsToFix,transfer_techs)
        $(transfer_linksDelta_upper(linksModelToCalc,yearsToFix,transfer_techs) < 0) = 0;

    transfer_linksDelta_lower(linksModelToCalc,yearsToFix,transfer_techs)
        $(sum(yearsToCalc$sameas(yearsToFix, yearsToCalc), 1))
        = transfer_linksParam(linksModelToCalc,yearsToFix,transfer_techs,"linksLowerLimit")
            - sum(vintage, transfer_linksTotal.l(linksModelToCalc,yearsToFix,transfer_techs,vintage));
    transfer_linksDelta_lower(linksModelToCalc,yearsToFix,transfer_techs)
        $(transfer_linksDelta_lower(linksModelToCalc,yearsToFix,transfer_techs) < 0) = 0;

    transfer_linksBuild.l(linksModelToCalc,yearsToFix,transfer_techs,vintage)
        $transfer_availTech(linksModelToCalc,yearsToFix,transfer_techs,vintage)
        = transfer_linksBuild.l(linksModelToCalc,yearsToFix,transfer_techs,vintage)
            - transfer_linksDelta_upper(linksModelToCalc,yearsToFix,transfer_techs);

    transfer_linksDecom.l(linksModelToCalc,yearsToFix,transfer_techs,vintage)
        $transfer_usedTech(linksModelToCalc,yearsToFix,transfer_techs,vintage)
        = transfer_linksDecom.l(linksModelToCalc,yearsToFix,transfer_techs,vintage)
            - transfer_linksDelta_lower(linksModelToCalc,yearsToFix,transfer_techs);

    transfer_linksBuild.l(linksModelToCalc,yearsToFix,transfer_techs,vintage)
        $(transfer_linksBuild.l(linksModelToCalc,yearsToFix,transfer_techs,vintage) < 0) = 0;
    transfer_linksBuild.fx(linksModelToCalc,yearsToFix,transfer_techs,vintage)
        = transfer_linksBuild.l(linksModelToCalc,yearsToFix,transfer_techs,vintage);
    transfer_linksDecom.l(linksModelToCalc,yearsToFix,transfer_techs,vintage)
        $(transfer_linksDecom.l(linksModelToCalc,yearsToFix,transfer_techs,vintage) < 0) = 0;
    transfer_linksDecom.fx(linksModelToCalc,yearsToFix,transfer_techs,vintage)
        = transfer_linksDecom.l(linksModelToCalc,yearsToFix,transfer_techs,vintage);
    transfer_linksTotal.fx(linksModelToCalc,yearsToFix,transfer_techs,vintage)
        = transfer_linksTotal.l(linksModelToCalc,yearsToFix,transfer_techs,vintage);

    accounting_indicator.fx(accNodesModel,accYearsToFix,indicator)
        = accounting_indicator.l(accNodesModel,accYearsToFix,indicator);

* Optimize and log values
    if (opti_sense < 0,
    solve remix minimizing accounting_objective using mip;
    else
    solve remix maximizing accounting_objective using mip;
    );

    put_utility 'log' / 'Model status ' remix.modelstat:0:0;
    put_utility 'log' / 'Objective value ' accounting_objective.l:0:3;

);

$onVerbatim
$endif.method
$offVerbatim
