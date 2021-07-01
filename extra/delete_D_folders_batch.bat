REM Deletes all .D folders in all subdirectories of current directory
REM Place the .bat file in the parent folder of the folders containing .D files
ECHO OFF
setlocal EnableDelayedExpansion


:PROMPT
SET /P AREYOUSURE=This will delete all .D folders in all subdirectories, this cannot be reversed, are you sure? ([Y]/[N])?
IF /I "%AREYOUSURE%" NEQ "Y" (
    echo ------- Aborted Process -------
    GOTO END
)

for /D %%d in (.\*) do (
@echo ---------------------------Processing folder: %%~nd----------------------------

    for /D %%D in ("%%d"\*) do (
        set "folder_extension=%%~xD"
        if !folder_extension!==.D RD /S /Q "%%D"
        
    )

)
:END

pause