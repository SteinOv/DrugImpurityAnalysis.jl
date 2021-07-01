REM Converts all .D folders in all subdirectories of current directory to mzxml
REM Additionally, retrieves sample_info.xml files from .D folders
REM Replace {MS_convert.exe location} with location of MS_convert.exe,
REM example location: "C:\Program Files\ProteoWizard\ProteoWizard 3.0.21180.d45de83ec\msconvert.exe"
REM Then place the .bat file in the parent folder of the folders containing .D files
@ECHO OFF
setlocal EnableDelayedExpansion

for /D %%d in (.\*) do (
@echo ---------------------------Processing folder: %%~nd----------------------------

    "{MS_convert.exe location}" "%%d\*.D" "--mzXML" "--32" "-o" "%%~nd"
    @echo ----------------------------Retrieving sample info----------------------------
    for /D %%D in ("%%d"\*) do (
        copy "%%D\AcqData\sample_info.xml" "%%d\%%~nD_sample_info.xml"
    )

)
pause