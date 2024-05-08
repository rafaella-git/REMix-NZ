* SPDX-FileCopyrightText: Copyright (c) 2023-2024 German Aerospace Center (DLR)
* SPDX-License-Identifier: BSD-3-Clause

$onVerbatim
$iftheni.method %method%==mga

$if not set alternatives  $setglobal alternatives      5
$if not set mgafactor     $setglobal mgafactor         1.03
$if not set mgamethod     $setglobal mgamethod         hypersphere

$setglobal scenidx mga,
$setglobal selscen mga_act,

set mga / mga0 * mga%alternatives% /;
set mga_act(mga);
set mga_comp(mga);

mga_act(mga)$(ord(mga) = 1) = yes;
$endif.method