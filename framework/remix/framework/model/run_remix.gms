* SPDX-FileCopyrightText: Copyright (c) 2023-2024 German Aerospace Center (DLR)
* SPDX-License-Identifier: BSD-3-Clause

* ==== REMix version number ====
$setglobal remixversion 0.9.7

* ==== global settings ====
$if not set sourcedir           $setglobal sourcedir                 %gams.workdir%%system.dirsep%source
$if not set datadir             $setglobal datadir                   %gams.workdir%
$if not set scendir             $setglobal scendir                   .
$if not set instancedir         $setglobal instancedir               %gams.scrdir%

$if not dexist "%datadir%/%scendir%" $abort "Error: Data directory %datadir%/%scendir% not found!"

* ==== write metadata and inherit dataset ====
$include "%sourcedir%/metadata.gms"
$include "%sourcedir%/data_inheritance.gms"

* ==== run remix ====
$onVerbatim
$iftheni.test_inputdata %test%==inputdata
$log "Testing input data, model will not be solved"
$include "%sourcedir%/methods/test_inputdata.gms"
$else.test_inputdata
$offVerbatim
$include "%sourcedir%/remix.gms"
$onVerbatim
$endif.test_inputdata
$offVerbatim
