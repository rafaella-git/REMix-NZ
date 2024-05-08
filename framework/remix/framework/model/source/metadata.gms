$log ### Metadata ###
$log Using REMix version %remixversion%
$log Using GAMS version %system.GamsReleaseMaint%
$log Using data directory %datadir% and scenario directory %scendir%
$log Using framework path %sourcedir%
$log Using project path %datadir%%scendir%
$log Using instance directory %instancedir%

$call 'printf "remix_version \"%remixversion%\"\n" > %gams.scrdir%%system.dirsep%metadata';
$call 'printf "gams_version \"%system.GamsReleaseMaint%\"\n" >> %gams.scrdir%%system.dirsep%metadata';
$call 'printf "timestamp \"%system.DATE% %system.TIME%\"\n" >> %gams.scrdir%%system.dirsep%metadata';

$call 'printf "framework_path " >> %gams.scrdir%%system.dirsep%metadata';
$iftheni.os %system.FileSys%==msnt
$call 'cd %sourcedir% && cd| sed "s/^/\"/;s/$/\"/" >> %gams.scrdir%%system.dirsep%metadata';
$else.os
$call 'cd %sourcedir% && pwd| sed "s/^/\"/;s/$/\"/" >> %gams.scrdir%%system.dirsep%metadata';
$endif.os
$call 'printf "\n" >> %gams.scrdir%%system.dirsep%metadata';
$call 'printf "framework_hash " >> %gams.scrdir%%system.dirsep%metadata';
$iftheni.os %system.FileSys%==msnt
$call 'cd %sourcedir% && git rev-parse HEAD >nul 2>&1 && git rev-parse HEAD >> %gams.scrdir%%system.dirsep%metadata';
$else.os
$call 'cd %sourcedir% && git rev-parse HEAD > /dev/null 2>&1 && git rev-parse HEAD >> %gams.scrdir%%system.dirsep%metadata';
$endif.os
$call 'printf "\n" >> %gams.scrdir%%system.dirsep%metadata';
$call 'printf "framework_branch " >> %gams.scrdir%%system.dirsep%metadata';
$iftheni.os %system.FileSys%==msnt
$call 'cd %sourcedir% && git rev-parse --abbrev-ref HEAD >nul 2>&1 && git rev-parse --abbrev-ref HEAD | sed "s/^/\"/;s/$/\"/" >> %gams.scrdir%%system.dirsep%metadata';
$else.os
$call 'cd %sourcedir% && git rev-parse --abbrev-ref HEAD > /dev/null 2>&1 && git rev-parse --abbrev-ref HEAD | sed "s/^/\"/;s/$/\"/" >> %gams.scrdir%%system.dirsep%metadata';
$endif.os
$call 'printf "\n" >> %gams.scrdir%%system.dirsep%metadata';

$call 'printf "project_path " >> %gams.scrdir%%system.dirsep%metadata';
$iftheni.os %system.FileSys%==msnt
$call 'cd %datadir%%system.dirsep%%scendir% && cd| sed "s/^/\"/;s/$/\"/" >> %gams.scrdir%%system.dirsep%metadata';
$else.os
$call 'cd %datadir%%system.dirsep%%scendir% && pwd| sed "s/^/\"/;s/$/\"/" >> %gams.scrdir%%system.dirsep%metadata';
$endif.os
$call 'printf "\n" >> %gams.scrdir%%system.dirsep%metadata';
$call 'printf "project_hash " >> %gams.scrdir%%system.dirsep%metadata';
$iftheni.os %system.FileSys%==msnt
$call 'cd %datadir%%system.dirsep%%scendir% && git rev-parse HEAD >nul 2>&1 && git rev-parse HEAD | sed "s/^/\"/;s/$/\"/" >> %gams.scrdir%%system.dirsep%metadata';
$else.os
$call 'cd %datadir%%system.dirsep%%scendir% && git rev-parse HEAD > /dev/null 2>&1 && git rev-parse HEAD | sed "s/^/\"/;s/$/\"/" >> %gams.scrdir%%system.dirsep%metadata';
$endif.os
$call 'printf "\n" >> %gams.scrdir%%system.dirsep%metadata';
$call 'printf "project_branch " >> %gams.scrdir%%system.dirsep%metadata';
$iftheni.os %system.FileSys%==msnt
$call 'cd %datadir%%system.dirsep%%scendir% && git rev-parse --abbrev-ref HEAD >nul 2>&1 && git rev-parse --abbrev-ref HEAD | sed "s/^/\"/;s/$/\"/" >> %gams.scrdir%%system.dirsep%metadata';
$else.os
$call 'cd %datadir%%system.dirsep%%scendir% && git rev-parse --abbrev-ref HEAD > /dev/null 2>&1 && git rev-parse --abbrev-ref HEAD | sed "s/^/\"/;s/$/\"/" >> %gams.scrdir%%system.dirsep%metadata';
$endif.os
$call 'printf "\n" >> %gams.scrdir%%system.dirsep%metadata';

set metadata(*) /
$include "%gams.scrdir%%system.dirsep%metadata"
/;