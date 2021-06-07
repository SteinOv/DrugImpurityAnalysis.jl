REM Converts all .D folders in current directory to mzxml
REM Additionally, retrieves sample_info.xml files from .D folders

@ECHO OFF
setlocal EnableDelayedExpansion

"P:\ProteoWizard\ProteoWizard 3.0.21063.9e153e63f\MSconvert.exe" "*.D" "--mzXML" "--32"
@echo ----------------------------Retrieving sample info----------------------------
for /D %%D in (.\*) do (
    @REM set "dir_name=%%~nD"
    echo "%%D\AcqData\sample_info.xml"
    copy "%%D\AcqData\sample_info.xml" "%%~nD_sample_info.xml"
)


pause