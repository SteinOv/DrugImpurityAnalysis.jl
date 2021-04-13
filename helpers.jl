"""
Helper functions

RT_indices: Converts retention time to indexes in spectrum
filter_XIC: Returns spectrum filtered at chosen mz values
batch_import: Imports spectra present in given folder
retrieve_sample_name: Retrieves sample names for the spectra from the agilent .D folder


"""

function RT_indices(spectrum, RT)
	"""Returns index range of given retention time range"""

	# Only possible if Rt is sorted
	if !issorted(spectrum["Rt"])
		error("Error: Retention time array is not in order")
	end
	
	return [findfirst(i -> i >= RT[1], spectrum["Rt"]), findlast(i -> i <= RT[2], spectrum["Rt"])]
end

function filter_XIC(spectrum, mass_values)
	"""Returns filtered spectrum based on given mass (mz) values"""

	# Return complete spectrum
	if mass_values == 0
		return reduce(+, spectrum["Mz_intensity"], dims=2)
	# Convert to tuple
	elseif typeof(mass_values) == Int
		mass_values = (mass_values,)
	end

	indices = Array{CartesianIndex{2}, 1}()

	# Add indices that correspond to given mz values
	for mass in mass_values
		mass_range = [mass - MAX_MASS_DEVIATION, mass + MAX_MASS_DEVIATION]

		# Store which indices are between mass range
		append!(indices, findall(mz -> (mz .>= mass_range[1]) .& (mz .<= mass_range[2]), spectrum["Mz_values"]))
	end

	# Create XIC spectrum
	spectrum_XIC = @MVector zeros(Int, size(spectrum["Mz_values"], 1))

	# Filter spectrum using stored indices
	for I in indices
		spectrum_XIC[I[1]] += spectrum["Mz_intensity"][I]
	end

	return spectrum_XIC
end


function batch_import(pathin)
	"""Returns list containing all imported spectra from folder"""
	mz_thresh = [0, 0]
	Int_thresh = 0
	files = readdir(pathin, sort=false)
	files_supported = zeros(Bool, length(files),1)

	# Determine supported files
	for i=1:length(files)
		if lowercase(files[i][end-4:end]) == "mzxml" || lowercase(files[i][end-3:end]) == "cdf"
			files_supported[i] = true
		end
	end

	# Save indices of supported files
	supported_indices=findall(x -> x == true, files_supported)
	num_of_spectra = length(supported_indices)

	# Create array for storing spectra
	spectra = Array{Any}(undef,size(supported_indices,1))

	# Load files into spectra array
	for (i, file_index) in enumerate(supported_indices)
		filename = files[file_index]
		spectra[i] = import_files(pathin,[filename],mz_thresh,Int_thresh)
		spectra[i]["MS1"]["Filename"] = [filename]
		spectra[i]["MS1"]["Sample Name"] = [retrieve_sample_name(pathin, filename)]
		

		println("Read spectra $(i[1]) of $num_of_spectra")
	end

	return spectra
end

function retrieve_sample_name(pathin, filename)
	# Path to xml file which contains Sample Name
	folder_name = filename[begin:end-5] * "D"
	folder_loc = joinpath(pathin, folder_name)
	sample_info_file = joinpath(folder_loc, "AcqData\\sample_info.xml")

	# File does not exist
	if !isfile(sample_info_file)
		return ""
	end

	# Load xml file into array
	xdoc = parse_file(sample_info_file)
	xroot = root(xdoc)
	xarray = xroot["Field"]

	# Find index where Sample Name is stored and retrieve its value
	element_index = findfirst(i -> content(find_element(i, "Name")) == "Sample Name", xarray)
	sample_name = content(find_element(xarray[element_index], "Value"))

	# Free memory
	free(xdoc)

	return sample_name
end
