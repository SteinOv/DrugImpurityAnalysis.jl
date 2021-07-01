REM Converts all .D folders in current directory to mzxml
REM Additionally, retrieves sample_info.xml files from .D folders
REM Replace {MS_convert.exe location} with location of MS_convert.exe
REM example location: "C:\Program Files\ProteoWizard\ProteoWizard 3.0.21180.d45de83ec\msconvert.exe"
REM Then place the .bat file in the same folder as the .D files
@ECHO OFF
setlocal EnableDelayedExpansion

"{MS_convert.exe location}" "*.D" "--mzXML" "--32"
@echo ----------------------------Retrieving sample info----------------------------
for /D %%D in (.\*) do (
    echo "%%D\AcqData\sample_info.xml"
    copy "%%D\AcqData\sample_info.xml" "%%~nD_sample_info.xml"
)


pause