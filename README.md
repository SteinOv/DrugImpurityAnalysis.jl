# DrugImpurityAnalysis.jl

## -------- README in process --------

<br/>
This program can process gc-ms files of drugs, with the focus on cocaine. The program can detect adulterants and impurities in the drug sample, and determine the fraction compared to the major drug. This can be used to find similarities between, and link different samples.
<br/><br/>

### Requirements and Usage

This program has only been tested on Windows.

The gc-ms files must have the .mzXML extension, gc-ms files can be converted to this extension using ProteoWizard: `http://proteowizard.sourceforge.net/download.html`\
In the "extra" folder, .bat files are present to automatically convert agilent .D folders to .mzXML after installation of ProteoWizard, but the installation location of ProteoWizard must be set to the right folder in the .bat files.\
The program supports importing sample names of .D files, to use this feature, the .mzxml and .D files should both be in the same directory or the sample_info extracted using the .bat file. Open the .bat files as text for more detailed instructions\
Function create_impurity_profile creates impurity profiles from mzXML files in pathin folder. create_impurity_profiles_batch creates impurity profiles from all mzXML files in subfolders.
<br/><br/>

compounds.csv should contain information about the major drug, impurities and adulterants.
- RT: The expected retention time of the compound.
- mz: Mass over charge values that the program should take into account, separated by a semicolon. Alternatively, the mz values can be put in brackets and separated by commas, the program will then take the sum of those mz values, instead of processing them separately.
- overlap: Manually define RT of peak that overlaps with this compound.
<br/><br/>

In settings.json the following settings can me modified:

