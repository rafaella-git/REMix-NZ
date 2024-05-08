* SPDX-FileCopyrightText: Copyright (c) 2023-2024 German Aerospace Center (DLR)
* SPDX-License-Identifier: BSD-3-Clause

$onVerbatim
$if not set resultdir        $set resultdir               %gams.workdir%/results
$if not set resultfile       $set resultfile              remix
$offVerbatim

$offListing
$onEmpty

table largenumbers_ts(*,*,*)
$onDelim
$if exist "%gams.scrDir%/largenumbers_profile.csv" $include "%gams.scrDir%/largenumbers_profile.csv"
$offDelim
;
display largenumbers_ts;

table movingdecimals_ts(*,*,*)
$onDelim
$if exist "%gams.scrDir%/movingdecimals_profile.csv" $include "%gams.scrDir%/movingdecimals_profile.csv"
$offDelim
;
display movingdecimals_ts;

table normal_parameters(*,*,*)
$onDelim
$if exist "%gams.scrDir%/normal_parameters.csv" $include "%gams.scrDir%/normal_parameters.csv"
$offDelim
;
display normal_parameters;

table onlydecimals_ts(*,*,*)
$onDelim
$if exist "%gams.scrDir%/onlydecimals_profile.csv" $include "%gams.scrDir%/onlydecimals_profile.csv"
$offDelim
;
display onlydecimals_ts;

table onlyintegers_ts(*,*,*)
$onDelim
$if exist "%gams.scrDir%/onlyintegers_profile.csv" $include "%gams.scrDir%/onlyintegers_profile.csv"
$offDelim
;
display onlyintegers_ts;

table zerotails_ts(*,*,*)
$onDelim
$if exist "%gams.scrDir%/zerotails_profile.csv" $include "%gams.scrDir%/zerotails_profile.csv"
$offDelim
;
display zerotails_ts;

set testset(*) /
$onDelim
$if exist "%gams.scrDir%/set_testset.csv" $include "%gams.scrDir%/set_testset.csv"
$offDelim
/;
display testset;

$offEmpty
$onListing

$onVerbatim
execute_unload "%resultdir%/%resultfile%.gdx"
    largenumbers_ts,
    movingdecimals_ts,
    normal_parameters,
    onlydecimals_ts,
    onlyintegers_ts,
    zerotails_ts,
    testset;
$offVerbatim