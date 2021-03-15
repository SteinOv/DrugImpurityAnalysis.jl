using MS_Import
using Plots
function main()

	# Specify path
	folder = "P:/Git/bachelor_project"
	mzxml_path = "data/mzxml"
	pathin = joinpath(folder, mzxml_path)

	# Import spectra
	spectra = batch_import(pathin)

	spectrum = spectra[10]["MS1"]
	mass = [81.9, 82.1]
	mass = [182, 183] # 182.10000610351562 (5 @[1930, 137]), 182.10000610351562 (10 @[1930, 138])

	spectrum_filtered = filter_mass(spectrum, mass)
	println(spectrum["Rt"][findmax(spectrum_filtered)[2]])
	plot(spectrum["Rt"], spectrum_filtered, xminorticks=0.1)
	spectrum["Mz_values"][1930, :]
	plot(spectrum["Mz_values"][1930, :], spectrum["Mz_intensity"][1930, :])
	findmax(spectrum["Mz_intensity"][1930, :])
	findmax(spectrum["Mz_values"][1930, 138])

	# print Rt and intensity of max for each spectrum
	for i=1:length(spectra)
		spectrum = spectra[i]["MS1"]
		spectrum_filtered = filter_mass(spectrum, mass)
		max_Rt = spectrum["Rt"][findmax(spectrum_filtered)[2]]
		max_int = findmax(spectrum_filtered)[1]
		println("$i, $max_Rt, $max_int")
	end
end


function filter_mass(spectrum, mass)
	"""Returns filtered spectrum based on mass range"""
	
	# Store which indices are between mass range
	filter = (spectrum["Mz_values"] .>= mass[1]) .& (spectrum["Mz_values"] .<= mass[2])

	# For every scan (Rt), sum mz intensities at chosen mz values
	spectrum_filtered = zeros(Int, size(spectrum["Mz_values"][:, 1]))
	for i=1:length(spectrum_filtered)
		for j=1:length(filter[i, :])
			if filter[i, j] == true
				spectrum_filtered[i] += spectrum["Mz_intensity"][i, j]
			end
		end
	end

	return spectrum_filtered

end


function batch_import(pathin)
	"""Returns list containing all imported spectra from folder"""
	mz_thresh = [0, 0]
	Int_thresh = 0
	files = readdir(pathin)
	files_supported = zeros(Bool, length(files),1)

	# Determine supported files
	for i=1:length(files)
		if lowercase(files[i][end-4:end]) == "mzxml" || lowercase(files[i][end-3:end]) == "cdf"
			files_supported[i] = true
		end
	end

	# Save indices of supported files
	supported_indices=findall(x -> x == true, files_supported)

	# Create array for storing spectra
	spectra = Array{Any}(undef,size(supported_indices,1))

	# Load files into spectra array
	for i in supported_indices
		filenames = [files[i]]
		spectra[i] = import_files(pathin,filenames,mz_thresh,Int_thresh)
	end

	return spectra
end

main()

