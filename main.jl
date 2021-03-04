using MS_Import



folder = "P:/Git/bachelor_project"
mzxml_path = "data/mzxml"
pathin = joinpath(folder, mzxml_path)
# filenames=["B039695.mzXML", "B039699.mzXML"]
filenames=["B039695.mzXML"]

# mz_values,mz_int,t0,t_end,file_name,pathin,
# msModel,msIonisation,msManufacturer,polarity,Rt=import_files_MS1(pathin,filenames) # Only MS1

chrom=import_files(pathin,filenames)


