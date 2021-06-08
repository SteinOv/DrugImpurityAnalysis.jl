REM Converts all .D folders in all subdirectories of current directory to mzxml
REM Additionally, retrieves sample_info.xml files from .D folders

@ECHO OFF
setlocal EnableDelayedExpansion

for /D %%d in (.\*) do (
@echo ---------------------------Processing folder: %%~nd----------------------------

    "P:\ProteoWizard\ProteoWizard 3.0.21063.9e153e63f\MSconvert.exe" "%%d\*.D" "--mzXML" "--32" "-o" "%%~nd"
    @echo ----------------------------Retrieving sample info----------------------------
    for /D %%D in ("%%d"\*) do (
        copy "%%D\AcqData\sample_info.xml" "%%d\%%~nD_sample_info.xml"
    )

)
pause