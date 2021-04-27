# Bachelor Project

© 2021 "All Rights Reserved"

This program can process gc-ms files of drugs, with the focus on cocaine. The program can detect adulterants and impurities in the drug sample, and determine the fraction compared to the major drug. This can be used to find similarities between, and link different samples.


### Requirements

This program has only been tested on Windows.
The package MS_Import.jl is needed to use this program, which can be found here: https://bitbucket.org/SSamanipour/ms_import.jl/src/master/

Run the following to install the MS_Import.jl package:
'''
using Pkg
Pkg.add(PackageSpec(url="https://bitbucket.org/SSamanipour/ms_import.jl/src/master/"))
'''

The gc-ms files must have the .mzxml extension, gc-ms files can be converted to this extension using ProteoWizard: http://proteowizard.sourceforge.net/download.html
The program supports importing sample names of .D files, to use this feature the .mzxml and .D files should both be in the same directory.
The gc-ms files should be under data/[subfolder] and the subfolder name specified in main.jl.

compounds.csv should contain information about the major drug, impurities and adulterants.
- RT: The expected retention time of the compound
- mz: Mass over charge values that the program should take into account, separated by a semicolon. Alternatively, the mz values can be put in brackets and separated by commas, the program will then take the sum of those mz values, instead of processing them separate.
- type: Whether the compound is the major compound (major), internal standard (internal_standard), impurity or an adulterant.
- type_int: The type displayed as an integer; major: -1, internal_standard: -10, impurity: 0, adulterant: 1.
- ratio: The normalized ratio's of the mz values to the mz value with the highest intensity. Not used at the moment, except for the internal standard; the ratio defined here is the minimum integral ratio {highest mz of IS / highest mz major compound} to process the sample.
- overlap: Manually define RT of peak that overlaps with this compound.


### Program in detail


In short, the program consists of five steps
- Read files
  For each spectrum:
    - Process metadata
    - Process major compound 
    - Process impurities and adulterants
- Write results to .csv file


#### Read files
1) The gc-ms files are imported using the MS_Import.jl package. Additionally, the filename, folder name and if .D folder is present; sample name are retrieved.
2) compounds.csv is read.

#### Process metadata
1) The filename, folder name and sample name are stored in a dataframe.

#### Process major compound
1) The actual retention time of the major compound is determined by searching for the highest intensity around the predicted RT in the ion-extracted chromatogram (XIC) of the highest mz value.
2) The retention time modifier is determined {actual RT / predicted RT}.
3) The highest mz value for the major compound is integrated.
4) The highest mz value for the internal standard is integrated.
5) The {IS / major compound} ratio is determined.
6.1) If the {IS / major compound} ratio is higher than the ratio defined in compounds.csv, the major integral is set to 0 and the sample is not processed further.
6.2) Else; The sample is processed further, all mz values of the major compound are integrated and the total intensity determined.

#### Process impurities and adulterants
For each impurity/adulterant:
  1) The actual RT is determined by {predicted RT * RT modifier}
  For each mz value:
    2) An XIC spectrum is created of the specified mz value (or the sum of several mz values if defined in compounds.csv).
    3) The peak is integrated
  4) The total intensity is determined
  5) The ratio to the major compound is calculated: {total intensity impurity/adulterant / total intensity major compound * factor}
  6) This ratio is added to the dataframe containing the metadata


#### Write results to .csv file
1) The dataframe containing the metadata and ratio's of all spectra is written to a .csv file.


### Peak integration in detail



