default_test_remix_kwargs = {
    "specific": "",
    "mode": "solve",
    "junitxml": "report.xml",
    "keep": 0,
    "timelimit": 7200,
}

default_test_remix_help = {
    "specific":
        "Filter for test case selection.",

    "timelimit":
        "Filter TCs based on expected runtime.",

    "mode":
        "Decide if test instances are run (solve) or metrics are checked (metrics).",

    "junitxml":
        "File name of the final testing report.",

    "keep":
        "Instruct GAMS to keep the scratch directory (225a) after the run is "
        "finished.",
}
