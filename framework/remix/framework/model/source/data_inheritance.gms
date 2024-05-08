* SPDX-FileCopyrightText: Copyright (c) 2023-2024 German Aerospace Center (DLR)
* SPDX-License-Identifier: BSD-3-Clause

$onVerbatim
$iftheni.doInheritance not "%system.FE%"==".dmp"

$if not set roundts  $setglobal roundts   0
$if not set cache    $setglobal cache     0
$if not set format    $setglobal format   "dat"
$offVerbatim

* pandas based data inheritance reading in the dat files, merging them and writing csv files in the scrdir
$onEmbeddedCode  Python:
    from pathlib import Path

    try:
        from remix.framework.tools.inheritance import main_inheritance
    except ImportError:
        from sys import path
        model_dir = Path(r'%sourcedir% '.strip()).parents[3].as_posix()
        if model_dir not in path:
            path.append(model_dir)
        from remix.framework.tools.inheritance import main_inheritance

    scr_dir = Path(r'%gams.scrDir% '.strip()).as_posix()
    instance_dir = Path(r'%instancedir% '.strip()).as_posix()
    base_path = Path(r'%datadir% '.strip()).as_posix()
    scen_dir = Path(r'%scendir% '.strip()).as_posix()

    if instance_dir == scr_dir:
        main_inheritance(instance_dir, base_path, scen_dir, roundts=%roundts%, _logger=gams.printLog)

$offEmbeddedCode
;

$onVerbatim
$endif.doInheritance
$offVerbatim
