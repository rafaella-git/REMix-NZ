* SPDX-FileCopyrightText: Copyright (c) 2023-2024 German Aerospace Center (DLR)
* SPDX-License-Identifier: BSD-3-Clause

$onVerbatim
$iftheni.method %method%==pareto

$if not set paretopoints  $setglobal paretopoints      5
$if not set paretofactor  $setglobal paretofactor      1.05

$setglobal scenidx pareto,
$setglobal selscen pareto_act,

set pareto / pareto0 * pareto%paretopoints% /;
set pareto_act(pareto);
pareto_act(pareto)$(ord(pareto) = 1) = yes;
$endif.method