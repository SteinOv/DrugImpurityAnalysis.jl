using MS_Import
using Plots
function main()

	# Specify path
	folder = "P:/Git/bachelor_project"
	mzxml_path = "data/mzxml"
	pathin = joinpath(folder, mzxml_path)

	# Import spectra
	spectra = batch_import(pathin)
	instr = spectrum["MS_Instrument"]
	spectrum = spectra[1]["MS1"]


	spectrum_filtered = filter_mass(spectrum, mass)
	println(spectrum["Rt"][findmax(spectrum_filtered)[2]])
	plot(spectrum["Rt"], spectrum_filtered, xminorticks=0.1)
	findmax(spectrum_filtered[1940:end])

	spectrum_filtered[2490]

	println(spectrum["Mz_intensity"][2490, :])
	spectrum["Mz_values"][2490, 110]

	plot(spectrum["Mz_values"][2490, :], spectrum["Mz_intensity"][2490, :])
	findmax(spectrum["Mz_intensity"][2490, :])
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


function cocaine_intensity(spectrum)
	"""Returns intensity of cocaine peak if present, else 0""" 

	# Cocaine approximate peak locations
	RT = 6.66 # minutes
	max_RT_deviation = 0.1 # minutes
	mass_vals = [82.07, 182.12, 83.07, 94.07, 77.04, 105.03, 96.08, 42.03, 303.15]
	max_mass_deviation = 0.05 

	RT_range = [RT - max_RT_deviation, RT + max_RT_deviation]
	RT_range_index = RT_indices(spectrum, RT_range)

	for mass in mass_vals
		mass_range = [mass - max_mass_deviation, mass + max_mass_deviation]
		spectrum_XIC = filter_mass(spectrum, mass_range)
		max = findmax(spectrum_XIC[RT_range_index[1]:RT_range_index[2]]) # max[1] = int, [2] = index
		println("$mass: $(max[1])")
	end

end

function RT_indices(spectrum, RT)
	"""Returns index range of given retention time range"""

	# Only possible if Rt is sorted
	if !issorted(spectrum["Rt"])
		error("Error: Retention time array is not in order")
	end
	
	return [findfirst(i -> i >= RT[1], spectrum["Rt"]), findlast(i -> i <= RT[2], spectrum["Rt"])]
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
		println("Read spectra $i")
	end

	return spectra
end

main()

