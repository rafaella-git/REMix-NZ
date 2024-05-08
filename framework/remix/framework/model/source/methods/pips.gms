* SPDX-FileCopyrightText: Copyright (c) 2023-2024 German Aerospace Center (DLR)
* SPDX-License-Identifier: BSD-3-Clause

$onVerbatim
$iftheni.method %method%==pips

* ==== global options ====
$if not set holdfixed     $set holdfixed     1
$if not set equlist       $set equlist       0
$if not set solvelink     $set solvelink     0
$if not set timeblocks    $set timeblocks    2190
$if not set splitgdx      $set splitgdx      0
$if not set gdxname       $set gdxname       remix
$if not set blockdir      $set blockdir      %gams.workdir%%system.dirsep%blocks
$if not set writelp       $set writelp       0
$if not set writemps      $set writemps      0
$if not set gdxnames      $set gdxnames      0
$if not set gdxuels       $set gdxuels       0
$if not set gdxdict       $set gdxdict       1
$if not set optfile       $set optfile       1

$ife %checkanno%==1       $set gdxuels       1

$setEnv GDXCONVERT v7
$setEnv GDXCOMPRESS 1
$offVerbatim

* ==== block mappings ====

$eval numberYearBlocksEval card(years)
$eval yearsBlockSize ceil(card(years)/%numberYearBlocksEval%)
set yearsHelper / yH1*yH%yearsBlockSize% /;
set yearsBlocks / yB1*yB%numberYearBlocksEval% /;
set yearsBlockMappingHelper(yearsBlocks, yearsHelper, years) / ( #yearsBlocks . #yearsHelper ) : #years /;
set yearsBlockMapping(yearsBlocks, years);
option yearsBlockMapping<yearsBlockMappingHelper;

$eval numberTimeBlocksEval ceil(%timeBlocks%)
$eval timeModelBlockSize card(timeModel)/%timeBlocks%
set timeModelBlocks / tmB1*tmB%numberTimeBlocksEval% /;
set timeBlockMapping(timeModelBlocks, timeModel);
timeBlockMapping(timeModelBlocks, timeModel)$(ord(timeModelBlocks)
    = ceil(ord(timeModel) * %timeBlocks% / card(timeModel))) = yes;

$eval lastBlock %numberTimeBlocksEval%*%numberYearBlocksEval%
set blocks / bl1*bl%lastBlock% /;
set blocksToCalc(blocks) / bl1*bl%lastBlock% /;
set blockSelected(blocks);

set blockMapping(blocks,yearsBlocks,timeModelBlocks) / #blocks : ( #yearsBlocks . #timeModelBlocks ) /;
set blockAssignmentsHelper(blocks,yearsBlocks,timeModelBlocks,years,timeModel);
blockAssignmentsHelper(blocks,yearsBlocks,timeModelBlocks,years,timeModel)
    $( blockMapping(blocks,yearsBlocks,timeModelBlocks)
        and yearsBlockMapping(yearsBlocks, years)
        and timeBlockMapping(timeModelBlocks, timeModel) )
    = yes;

set blockAssignments(blocks,years,timeModel);
option blockAssignments < blockAssignmentsHelper;

blockAssignments(blocks,years,timeModel)
    $(not yearsToCalc(years) or not timeModelToCalc(timeModel))
    = no;
option blocksToCalc < blockAssignments;

set blockAssignments_time(blocks,timeModel);
option blockAssignments_time < blockAssignments;

set blockAssignments_years(blocks,years);
option blockAssignments_years < blockAssignments;

$offorder
scalar annot_blockLast;
annot_blockLast = card(blocksToCalc) + 2;

parameter annot_blockStage(timeModel,years);
parameter annot_blockStage_x(years,timeModel);
annot_blockStage_x(yearsToCalc,timeModelToCalc)
    = sum(blockAssignments(blocksToCalc,yearsToCalc,timeModelToCalc), ord(blocksToCalc)) + 1;
option annot_blockStage < annot_blockStage_x;
$onorder

* ==== time model linked (tml) mappings ====
set tml_ext(blocks,timeModel);
set tml(timeModel);
tml_ext(blocks,timeModelToCalc)
    $(blockAssignments_time(blocks,timeModelToCalc) xor blockAssignments_time(blocks,timeModelToCalc--1)) = yes;
option tml<tml_ext;

* minUptime-1 due to consideration of current timestep
set tml_uptime_ext(blocks,timeModel,converter_techs,vintage);
set tml_uptime(timeModel,converter_techs,vintage);
tml_uptime_ext(blocks,timeModelToCalc,converter_hasMinDowntime(converter_techs,vintage))
    $(blockAssignments_time(blocks,timeModelToCalc) xor blockAssignments_time(blocks,timeModelToCalc--(converter_techParam(converter_techs,vintage,"minUptime")-1)))
    = yes;
option tml_uptime < tml_uptime_ext;

* minDowntime-1 due to consideration of current timestep
set tml_downtime_ext(blocks,timeModel,converter_techs,vintage);
set tml_downtime(timeModel,converter_techs,vintage);
tml_downtime_ext(blocks,timeModelToCalc,converter_hasMinUptime(converter_techs,vintage))
    $(blockAssignments_time(blocks,timeModelToCalc) xor blockAssignments_time(blocks,timeModelToCalc--(converter_techParam(converter_techs,vintage,"minDowntime")-1)))
    = yes;
option tml_downtime < tml_downtime_ext;


* ==== annotation ====

converter_unitsOnline.stage(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage)
        $converter_usedTech(nodesModelToCalc,yearsToCalc,converter_techs,vintage)
    = annot_blockStage(timeModelToCalc,yearsToCalc);
converter_activity.stage(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage,activity)
        $converter_usedTechAct(nodesModelToCalc,yearsToCalc,converter_techs,vintage,activity)
    = annot_blockStage(timeModelToCalc,yearsToCalc);

converter_unitsOnline_MIP.stage(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage)
        $(converter_usedTech(nodesModelToCalc,yearsToCalc,converter_techs,vintage)
            and converter_techParam(converter_techs,vintage,"mipDispatch"))
    = annot_blockStage(timeModelToCalc,yearsToCalc);
converter_unitsUsingActivity_MIP.stage(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage,activity)
        $(converter_usedTechAct(nodesModelToCalc,yearsToCalc,converter_techs,vintage,activity)
            and converter_techParam(converter_techs,vintage,"mipDispatch"))
    = annot_blockStage(timeModelToCalc,yearsToCalc);

converter_unitStartups.stage(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage)
        $(converter_usedTech(nodesModelToCalc,yearsToCalc,converter_techs,vintage)
            and converter_techParam(converter_techs,vintage,"mipDispatch"))
    = annot_blockStage(timeModelToCalc,yearsToCalc);
converter_unitShutdowns.stage(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage)
        $(converter_usedTech(nodesModelToCalc,yearsToCalc,converter_techs,vintage)
            and converter_techParam(converter_techs,vintage,"mipDispatch"))
    = annot_blockStage(timeModelToCalc,yearsToCalc);

Eq_converter_unitsOnline.stage(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage)
        $converter_usedTech(nodesModelToCalc,yearsToCalc,converter_techs,vintage)
    = annot_blockStage(timeModelToCalc,yearsToCalc);
Eq_converter_activityUpperLimit.stage(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage)
        $converter_usedTech(nodesModelToCalc,yearsToCalc,converter_techs,vintage)
    = annot_blockStage(timeModelToCalc,yearsToCalc);
Eq_converter_activityLowerLimit.stage(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage)
        $converter_usedTech(nodesModelToCalc,yearsToCalc,converter_techs,vintage)
    = annot_blockStage(timeModelToCalc,yearsToCalc);
Eq_converter_activityFixedLimit.stage(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage)
        $converter_usedTech(nodesModelToCalc,yearsToCalc,converter_techs,vintage)
    = annot_blockStage(timeModelToCalc,yearsToCalc);

Eq_converter_unitsOnlineMIP.stage(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage)
        $(converter_usedTech(nodesModelToCalc,yearsToCalc,converter_techs,vintage)
            and converter_techParam(converter_techs,vintage,"mipDispatch"))
    = annot_blockStage(timeModelToCalc,yearsToCalc);
Eq_converter_unitsOnlineUC.stage(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage)
        $(converter_usedTech(nodesModelToCalc,yearsToCalc,converter_techs,vintage)
            and converter_techParam(converter_techs,vintage,"mipDispatch"))
    = annot_blockStage(timeModelToCalc,yearsToCalc);
Eq_converter_noOnlineIdle.stage(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage)
        $(converter_usedTech(nodesModelToCalc,yearsToCalc,converter_techs,vintage)
            and converter_techParam(converter_techs,vintage,"mipDispatch"))
    = annot_blockStage(timeModelToCalc,yearsToCalc);

Eq_converter_minUptime.stage(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage)
        $(converter_usedTech(nodesModelToCalc,yearsToCalc,converter_techs,vintage)
            and converter_techParam(converter_techs,vintage,"mipDispatch")
            and converter_techParam(converter_techs,vintage,"minUptime"))
    = annot_blockStage(timeModelToCalc,yearsToCalc);
Eq_converter_minUptime.stage(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage)
        $(converter_usedTech(nodesModelToCalc,yearsToCalc,converter_techs,vintage)
            and converter_techParam(converter_techs,vintage,"mipDispatch")
            and converter_techParam(converter_techs,vintage,"minUptime")
            and tml_uptime(timeModelToCalc,converter_techs,vintage))
    = annot_blockLast;
Eq_converter_minDowntime.stage(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage)
        $(converter_usedTech(nodesModelToCalc,yearsToCalc,converter_techs,vintage)
            and converter_techParam(converter_techs,vintage,"mipDispatch")
            and converter_techParam(converter_techs,vintage,"minDowntime"))
    = annot_blockStage(timeModelToCalc,yearsToCalc);
Eq_converter_minDowntime.stage(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage)
        $(converter_usedTech(nodesModelToCalc,yearsToCalc,converter_techs,vintage)
            and converter_techParam(converter_techs,vintage,"mipDispatch")
            and converter_techParam(converter_techs,vintage,"minDowntime")
            and tml_downtime(timeModelToCalc,converter_techs,vintage))
    = annot_blockLast;

Eq_converter_activityUpperLimitPartLoad.stage(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage,activity)
        $(converter_usedTechAct(nodesModelToCalc,yearsToCalc,converter_techs,vintage,activity)
            and converter_techParam(converter_techs,vintage,"mipDispatch"))
    = annot_blockStage(timeModelToCalc,yearsToCalc);
Eq_converter_activityLowerLimitPartLoad.stage(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage,activity)
        $(converter_usedTechAct(nodesModelToCalc,yearsToCalc,converter_techs,vintage,activity)
            and converter_techParam(converter_techs,vintage,"mipDispatch"))
    = annot_blockStage(timeModelToCalc,yearsToCalc);

Eq_converter_activityStartups.stage(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage)
        $converter_usedTech(nodesModelToCalc,yearsToCalc,converter_techs,vintage)
    = annot_blockStage(timeModelToCalc,yearsToCalc);
Eq_converter_activityStartups.stage(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage)
        $(converter_usedTech(nodesModelToCalc,yearsToCalc,converter_techs,vintage)
            and tml(timeModelToCalc))
    = annot_blockLast;
Eq_converter_activityShutdowns.stage(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage)
        $converter_usedTech(nodesModelToCalc,yearsToCalc,converter_techs,vintage)
    = annot_blockStage(timeModelToCalc,yearsToCalc);
Eq_converter_activityShutdowns.stage(timeModelToCalc,nodesModelToCalc,yearsToCalc,converter_techs,vintage)
        $(converter_usedTech(nodesModelToCalc,yearsToCalc,converter_techs,vintage)
            and tml(timeModelToCalc))
    = annot_blockLast;


sourcesink_flow.stage(timeModelToCalc,nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity)
        $sourcesink_enabled(nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity)
    = annot_blockStage(timeModelToCalc,yearsToCalc);

Eq_sourcesink_useFixedSum.stage(nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity)
        $sourcesink_enabled(nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity)
    = annot_blockLast;
Eq_sourcesink_useLowerSum.stage(nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity)
        $sourcesink_enabled(nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity)
    = annot_blockLast;
Eq_sourcesink_useUpperSum.stage(nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity)
        $sourcesink_enabled(nodesModelToCalc,yearsToCalc,sourcesink_techs,commodity)
    = annot_blockLast;


storage_level.stage(timeModelToCalc,nodesModelToCalc,yearsToCalc,storage_techs,vintage,commodity)
        $storage_usedTechCom(nodesModelToCalc,yearsToCalc,storage_techs,vintage,commodity)
    = annot_blockStage(timeModelToCalc,yearsToCalc);
storage_losses.stage(timeModelToCalc,nodesModelToCalc,yearsToCalc,storage_techs,vintage,commodity)
        $storage_usedTechCom(nodesModelToCalc,yearsToCalc,storage_techs,vintage,commodity)
    = annot_blockStage(timeModelToCalc,yearsToCalc);

Eq_storage_levelLowerLimit.stage(timeModelToCalc,nodesModelToCalc,yearsToCalc,storage_techs,vintage,commodity)
        $storage_usedTechCom(nodesModelToCalc,yearsToCalc,storage_techs,vintage,commodity)
    = annot_blockStage(timeModelToCalc,yearsToCalc);
Eq_storage_levelUpperLimit.stage(timeModelToCalc,nodesModelToCalc,yearsToCalc,storage_techs,vintage,commodity)
        $storage_usedTechCom(nodesModelToCalc,yearsToCalc,storage_techs,vintage,commodity)
    = annot_blockStage(timeModelToCalc,yearsToCalc);
Eq_storage_losses.stage(timeModelToCalc,nodesModelToCalc,yearsToCalc,storage_techs,vintage,commodity)
        $storage_usedTechCom(nodesModelToCalc,yearsToCalc,storage_techs,vintage,commodity)
    = annot_blockStage(timeModelToCalc,yearsToCalc);


transfer_flowAlong.stage(timeModelToCalc,linksModelToCalc,yearsToCalc,transfer_techs,vintage)
    $transfer_usedTech(linksModelToCalc,yearsToCalc,transfer_techs,vintage)
    = annot_blockStage(timeModelToCalc,yearsToCalc);
transfer_flowAgainst.stage(timeModelToCalc,linksModelToCalc,yearsToCalc,transfer_techs,vintage)
    $transfer_usedTech(linksModelToCalc,yearsToCalc,transfer_techs,vintage)
    = annot_blockStage(timeModelToCalc,yearsToCalc);

Eq_transfer_flowAlongUpperLimit.stage(timeModelToCalc,linksModelToCalc,yearsToCalc,transfer_techs,vintage)
    $transfer_usedTech(linksModelToCalc,yearsToCalc,transfer_techs,vintage)
    = annot_blockStage(timeModelToCalc,yearsToCalc);
Eq_transfer_flowAgainstUpperLimit.stage(timeModelToCalc,linksModelToCalc,yearsToCalc,transfer_techs,vintage)
    $transfer_usedTech(linksModelToCalc,yearsToCalc,transfer_techs,vintage)
    = annot_blockStage(timeModelToCalc,yearsToCalc);

$onVerbatim
$iftheni.opfmethod %opfmethod%==kirchhoff
$offVerbatim
set grid_segments_sub(years,gridSegments);
grid_segments_sub(yearsToCalc,gridSegments)
    = sum((linksModelToCalc,transfer_techs,vintage)
            $(transfer_usedTech(linksModelToCalc,yearsToCalc,transfer_techs,vintage)
                and gridSegments_dcopf(linksModelToCalc,transfer_techs,gridSegments)), 1);

Eq_transfer_dcopf_cycleFlows.stage(timeModelToCalc,yearsToCalc,cycles,gridSegments)
   $grid_segments_sub(yearsToCalc,gridSegments)
    = annot_blockStage(timeModelToCalc,yearsToCalc);

$onVerbatim
$elseifi.opfmethod %opfmethod%==angle
$offVerbatim
transfer_dcopf_voltageAngle.stage(timeModelToCalc,nodesModelToCalc,yearsToCalc,gridSegments)
    $sum((linksModelToCalc)$transfer_incidenceSegments(nodesModelToCalc,linksModelToCalc,yearsToCalc,gridSegments), 1)
    = annot_blockStage(timeModelToCalc,yearsToCalc);

Eq_transfer_dcopf_angleFlows.stage(timeModelToCalc,linksModelToCalc,yearsToCalc,gridSegments)
    $sum((transfer_techs,vintage)$(transfer_usedTech(linksModelToCalc,yearsToCalc,transfer_techs,vintage)
                            and gridSegments_dcopf(linksModelToCalc,transfer_techs,gridSegments)), 1)
    = annot_blockStage(timeModelToCalc,yearsToCalc);

$onVerbatim
$endif.opfmethod
$offVerbatim

Eq_balance_commodities.stage(timeModelToCalc,nodesModelToCalc,yearsToCalc,commodity)
    = annot_blockStage(timeModelToCalc,yearsToCalc);
Eq_balance_commodities.stage(timeModelToCalc,nodesModelToCalc,yearsToCalc,commodity)
    $( balance_usedStorage(nodesModelToCalc,yearsToCalc,commodity)
        and tml(timeModelToCalc) )
    = annot_blockLast;

Eq_accounting_indicatorCalc.stage(accNodesModel,accYears,indicator)
    $activeIndicators(accNodesModel,accYears,indicator)
    = annot_blockLast;


* ==== write jacobian ====
yearsSel(years)$yearsToCalc(years) = yes;
timeModelSel(timeModel)$timeModelToCalc(timeModel) = yes;
nodesModelSel(nodesModel)$nodesModelToCalc(nodesModel) = yes;

if ((sum(nodesModelToCalc, 1)>40 or sum(timeModelToCalc, 1)>50) and not %equlist%,
   option limRow=0, limCol=0, solPrint=off;
else
   option limRow=100000, limCol=100000, solPrint=on;
);

$setenv GDXCOMPRESS 1

option MIP = Convert;

$onVerbatim
remix.optFile = %optfile%;
remix.solveLink = %solvelink%;
remix.holdFixed = %holdfixed%;
remix.priorOpt = 1;

$if not dexist "%blockdir%" put_utility 'exec' / 'mkdir -p %blockdir%'
file opt / "%gams.optdir%convert.opt" /;
put opt;
put "DumpGDX %blockdir%/%gdxname%.gdx" /;
$ife %WRITELP%==1 put "cplexLP %blockdir%/%gdxname%.lp" /;
$ife %WRITEMPS%==1 put "cplexMPS %blockdir%/%gdxname%.mps" /;
$ife %gdxnames%=0 put "GDXNames 0" /;
$ife %gdxuels%=0 put "GDXUELs 0" /;
put "GDXHessian 0" /;
put "GDXQuadratic 0" /;
putclose;
$offVerbatim

* Always calculate all yearsToCalc
yearsSel(years)$yearsToCalc(years) = yes;
yearsToFix(years)$(years.val < smin(years_a$yearsSel(years_a), years_a.val)) = yes;
accYearsSel("horizon") = yes;
accYearsSel(accYears)$(sum(yearsSel$sameas(accYears,yearsSel), 1)) = yes;
accYearsToFix(accYears)$(sum(years$(sameas(years,accYears) and years.val < smin(years_a$yearsSel(years_a), years_a.val)), 1) > 0) = yes;
timeModelSel(timeModel)$timeModelToCalc(timeModel) = yes;
nodesModelSel(nodesModel)$nodesModelToCalc(nodesModel) = yes;

* Fix variables before yearsToCalc
converter_unitsBuild.fx(nodesModelToCalc,yearsToFix,converter_techs,vintage)
    = converter_unitsBuild.l(nodesModelToCalc,yearsToFix,converter_techs,vintage);
converter_unitsDecom.fx(nodesModelToCalc,yearsToFix,converter_techs,vintage)
    = converter_unitsDecom.l(nodesModelToCalc,yearsToFix,converter_techs,vintage);
converter_unitsTotal.fx(nodesModelToCalc,yearsToFix,converter_techs,vintage)
    = converter_unitsTotal.l(nodesModelToCalc,yearsToFix,converter_techs,vintage);

storage_unitsBuild.fx(nodesModelToCalc,yearsToFix,storage_techs,vintage)
    = storage_unitsBuild.l(nodesModelToCalc,yearsToFix,storage_techs,vintage);
storage_unitsDecom.fx(nodesModelToCalc,yearsToFix,storage_techs,vintage)
    = storage_unitsDecom.l(nodesModelToCalc,yearsToFix,storage_techs,vintage);
storage_unitsTotal.fx(nodesModelToCalc,yearsToFix,storage_techs,vintage)
    = storage_unitsTotal.l(nodesModelToCalc,yearsToFix,storage_techs,vintage);

transfer_linksBuild.fx(linksModelToCalc,yearsToFix,transfer_techs,vintage)
    = transfer_linksBuild.l(linksModelToCalc,yearsToFix,transfer_techs,vintage);
transfer_linksDecom.fx(linksModelToCalc,yearsToFix,transfer_techs,vintage)
    = transfer_linksDecom.l(linksModelToCalc,yearsToFix,transfer_techs,vintage);
transfer_linksTotal.fx(linksModelToCalc,yearsToFix,transfer_techs,vintage)
    = transfer_linksTotal.l(linksModelToCalc,yearsToFix,transfer_techs,vintage);

accounting_indicator.fx(accNodesModel,accYearsToFix,indicator)
    = accounting_indicator.l(accNodesModel,accYearsToFix,indicator);

if (opti_sense<0,
solve remix minimizing accounting_objective using MIP;
else
solve remix maximizing accounting_objective using MIP;
);


* ==== check and split jacobian ====

$onVerbatim
$ifthene.checkanno %checkanno%==1
scalar maxStages; maxStages = sum(blocksToCalc, 1) + 2;
put_utility 'shell' / 'gams %gams.scrDir%checkanno.gms profile=%gams.profile% --jacFileName=%blockdir%/%gdxname%.gdx lo=4 --maxStages=' maxStages:0:0 ' --skipFix=1';
$endif.checkanno

$ifthene.gdxdict %gdxdict%==1
execute "mv -f  %gams.scrDir%/gamsdict.dat %blockdir%/%gdxname%_dict.gdx"
$endif.gdxdict

$ifthene.splitgdx %splitgdx%==1
scalar nbBlocks; nbBlocks = sum(blocksToCalc, 1) + 1;
scalar nbFiles; nbFiles = nbBlocks + 100;
put_utility 'shell' / 'ulimit -S -n ' nbFiles:0:0 '; gmschk -T -X -g %gams.sysdir% ' nbBlocks:0:0 ' %blockdir%/%gdxname%.gdx';
put_utility 'exec' / 'rm %blockdir%/%gdxname%.gdx'
put_utility 'exec' / 'rm %blockdir%/%gdxname%.map'
$endif.splitgdx

$endif.method
$offVerbatim
