* SPDX-FileCopyrightText: Copyright (c) 2023-2024 German Aerospace Center (DLR)
* SPDX-License-Identifier: BSD-3-Clause

$setargs paramOut paramIn operator weightingFactor

$if "%operator%"==""        $setlocal operator           sum
$if "%weightingFactor%"=="" $setlocal weightingFactor    1

%paramOut%$nodesModelToCalc(nodesModel)
    = %operator%(nodesData$aggregateNodesModel(nodesData,nodesModel), %paramIn% * %weightingFactor%) ;
