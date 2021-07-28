# DrugImpurityAnalysis.jl

## -------- README in process --------

<br/>
This program can process gc-ms files of drugs, with the focus on cocaine. The program can detect adulterants and impurities in the drug sample, and determine the fraction compared to the major drug. This can be used to find similarities between, and link different samples.
<br/><br/>

### Requirements and Usage

This program has only been tested on Windows and on Julia version 1.6.1.

#### Conversion to .mzxml

The gc-ms files must have the .mzXML extension, gc-ms files can be converted to this extension using ProteoWizard: `http://proteowizard.sourceforge.net/download.html`\
Download ProteoWizard for your operation system, install it and note down the installation directory.
In the "extra" folder, .bat files are present to automatically convert agilent .D folders to .mzXML after installation of ProteoWizard, but the installation location of ProteoWizard must be set to the right folder in the .bat files.\
In the "convert_to_mzxml" bat files, replace `{MS_convert.exe location}` with the location of the MS_convert.exe file located in the installation directory, keep the qoutes. For example line 9 in `convert_to_mzxml.bat` could be: 
`"C:\Program Files\ProteoWizard\ProteoWizard 3.0.21193.ccb3e0136\" "*.D" "--mzXML" "--32"` 

The program supports importing sample metadata, the metadata is not present in the .mzxml files. So the .bat file automatically extracts the file containing the metadata from the .D files. The program can read this file. Optionally, if this file is not present, the program will try to find the metadata in the .D files, if present.

#### compounds.csv
compounds.csv should contain information about the major drug, impurities and adulterants.
- RT: The expected retention time of the compound.
- mz: Mass over charge values that the program should take into account, separated by a semicolon. Alternatively, the mz values can be put in brackets and separated by commas, the program will then take the sum of those mz values, instead of processing them separately.
- overlap: Manually define RT of peak that overlaps with this compound.
<br/><br/>


#### settings.json
In settings.json the following settings can me modified:
main_compound: Name of the main compound (must match a compound in compounds.csv)
internal_standard: Name of the internal standard (must match a compound in compounds.csv)
main/IS min ratio: Minimum ratio between integrals of main compound and internal standard. Samples with a ratio below this value will be ignored and not be present in the output.

output_metadata: Which metadata should be shown in the output.
type: Type of metadata, currently only "sample_info_xml" is supported.
header: Header under which the metadata is output.
element_name: Name of the element in the sample_info.xml file.

exclude_from_impurity_profile: Compounds listed here are not shown in the output
compound_groups: Custom grouping of compounds, can be used in analysis.


### Functions
