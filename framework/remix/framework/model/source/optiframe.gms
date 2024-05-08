* SPDX-FileCopyrightText: Copyright (c) 2023-2024 German Aerospace Center (DLR)
* SPDX-License-Identifier: BSD-3-Clause

* ==== optimization frames ====
$onVerbatim

* mapping from optimization frame to years
$eval numYearsToCalc card(yearsToCalc)
set optiframe / of1 * of%numYearsToCalc% /;
set optiframeToCalc(optiframe);

$iftheni.pathopt %pathopt%==foresight
set map_optiframe(optiframe,years) / of1 . #yearsToCalc /;
$elseifi.pathopt %pathopt%==myopic
set map_optiframe(optiframe,years) / #optiframe : #yearsToCalc /;
$elseifi.pathopt %pathopt%==target
set map_optiframe(optiframe,years);
map_optiframe(optiframe, years)$(ord(optiframe) = 1 and years.val = smax(years_a$yearsToCalc(years_a), years_a.val)) = yes;
$else.pathopt
$abort "No valid pathopt mode specified. Use either foresight, myopic or target."
$endif.pathopt

$offVerbatim

option optiframeToCalc < map_optiframe;