* SPDX-FileCopyrightText: Copyright (c) 2023-2024 German Aerospace Center (DLR)
* SPDX-License-Identifier: BSD-3-Clause

* ==== global settings ====
$onVerbatim
$if not set method          $setglobal method            solve
$if not set pathopt         $setglobal pathopt           foresight
$if not set roundcoefs      $setglobal roundcoefs        0
$if not set debug           $setglobal debug             0
$if not set scenidx         $setglobal scenidx
$if not set selscen         $setglobal selscen

$setglobal run_postcalc                                  1
$if set postcalc            $setglobal run_postcalc      %postcalc%

$offVerbatim

$setglobal aggregateNodes          %sourcedir%/battools/aggregateNodes.gms
$setglobal aggregateAccountingMean %sourcedir%/battools/aggregateAccountingMean.gms

* ==== include modules ====
$include "%sourcedir%/sets.gms"
$include "%sourcedir%/accounting/input.gms"
$include "%sourcedir%/accounting/annuities.gms"

$include "%sourcedir%/methods/mga_pre.gms"
$include "%sourcedir%/methods/pareto_pre.gms"
$include "%sourcedir%/postcalc/declaration.gms"
$include "%sourcedir%/loadgdx.gms"

$include "%sourcedir%/core/converter.gms"
$include "%sourcedir%/core/storage.gms"
$include "%sourcedir%/core/transfer.gms"
$include "%sourcedir%/core/sourcesink.gms"
$include "%sourcedir%/core/balance.gms"

$include "%sourcedir%/accounting/equations.gms"
$include "%sourcedir%/optiframe.gms"

Model remix /
  M_converter
  M_storage
  M_transfer
  M_sourcesink
  M_balance
  M_accounting
/;


* ==== include methods  ====

* check if the method is valid
$onVerbatim
$iftheni.method %method%==solve
$log "Solving the REMix model"
$elseifi.method %method%==pareto
$log "Solving the REMix model and generating a pareto front"
$elseifi.method %method%==mga
$log "Solving the REMix model and generating alternatives (MGA)"
$elseifi.method %method%==iternodes
$log "Solving the REMix model by iteration of individual nodes"
$elseifi.method %method%==pips
$log "Building annotated gdx file for PIPS-IPM++"
$setGlobal run_postcalc 0
$else.method
$abort "REMix: No valid method chosen"
$endif.method
$offVerbatim


* if method is pips write the checkanno tool
$include "%sourcedir%/methods/checkanno.gms"

* if method is pips generate the annotated gdx file
$include "%sourcedir%/methods/pips.gms"

* if method is solve run the model with commercial solvers
$include "%sourcedir%/methods/solve.gms"

* if method is pareto solve the model once and afterwards run multiple points along the pareto front
$include "%sourcedir%/methods/pareto.gms"

* if method is mga solve the model once and afterwards maximize the length metric of indicators
$include "%sourcedir%/methods/mga.gms"

* if method is iternodes run the model by iterating through all nodesToCalc
$include "%sourcedir%/methods/iternodes.gms"

* if postcalc is one write the results gdx
$include "%sourcedir%/postcalc/definition.gms"
$include "%sourcedir%/postcalc/writegdx.gms"
