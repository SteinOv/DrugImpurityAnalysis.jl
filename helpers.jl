"""
Helper functions

RT_to_scans: Converts retention time to indexes in spectrum
filter_XIC: Returns spectrum filtered at chosen mz values
batch_import: Imports spectra present in given folder
retrieve_sample_name: Retrieves sample names for the spectra from the agilent .D folder


"""

"""Returns index range of given retention time range"""
function RT_to_scans(spectrum, RT_range)
	# Only possible if Rt is sorted
	if !issorted(spectrum["Rt"])
		error("Error: Retention time array is not in order")
	end

	# Convert Integer to tuple
	if RT_range isa(Number)
		RT_range = (RT_range, RT_range)
	end

	return searchsortedfirst(spectrum["Rt"], RT_range[1]), searchsortedlast(spectrum["Rt"], RT_range[2])
end

"""Returns filtered spectrum based on given mass (mz) values"""
function filter_XIC(spectrum, mass_values)
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
		mass_range = [mass - MAX_MZ_DEVIATION, mass + MAX_MZ_DEVIATION]

		# Store which indices are between mass range
		append!(indices, findall(mz -> (mz .>= mass_range[1]) .& (mz .<= mass_range[2]), spectrum["Mz_values"]))
	end

	# Create XIC spectrum
	spectrum_XIC = zeros(Int, size(spectrum["Mz_values"], 1))

	# Filter spectrum using stored indices
	for I in indices
		spectrum_XIC[I[1]] += spectrum["Mz_intensity"][I]
	end

	return spectrum_XIC
end

"""
	parse_data(row)
Parses mz values and RT overlap values from string to array
return mz_values, RT_overlap_values
"""
function parse_data(row)
	# Parse mz values
	mz_values_string = split(string(row.mz), ";")
	mz_values = Array{Any,1}(undef, length(mz_values_string))
	for (i, mz) in enumerate(mz_values_string)

		# Tuple of mz values
		if startswith(mz, "(")
			mz_tuple_string = split(mz[2:end - 1], ",")
			mz_values[i] = Tuple(([parse(Float32, sub_mass) for sub_mass in mz_tuple_string]))

		# Separate mz value
		else
			mz_values[i] = parse(Float32, mz)
		end
	end

	# Parse overlap values
	RT_overlap_values_string = split(string(row.overlap), ",")
	RT_overlap_values = zeros(Float16, 2)
	for RT_overlap in RT_overlap_values_string
		RT_overlap = RT_overlap == "missing" ? 0 : parse(Float16, RT_overlap)
		index = RT_overlap < row.RT ? 1 : 2
		RT_overlap_values[index] = RT_overlap
	end

	return mz_values, RT_overlap_values

end

"""
	batch_import(pathin)
Returns list containing all imported spectra from folder and metadata headers
"""
function batch_import(pathin, settings_json)
	
	files = readdir(pathin, sort=false)
	files_supported = zeros(Bool, length(files),1)

	# Determine supported files
	for i=1:length(files)
		if length(files[i]) > length(".mzxml") && lowercase(files[i][end-4:end]) == "mzxml"
			files_supported[i] = true
		end
	end

	# Save indices of supported files
	supported_indices=findall(x -> x == true, files_supported)
	num_of_spectra = length(supported_indices)

	# Retrieve which metadata to save
	info_containing = settings_json[Symbol("output_metadata")]
	sample_info_metadata = Dict(i["element_name"] => i["header"] for i in info_containing if i["type"] == "sample_info_xml")
	metadata_headers = append!([], values(sample_info_metadata), ["Filename", "Folder Name"])
	# Retrieve filters
	filters = settings_json[Symbol("filters")]
	filters_must_contain = Dict(i["element_name"] => i["value"] for i in filters if i["filter"] == "must_contain")



	spectra = []
	# Load files into spectra array
	for (i, file_index) in enumerate(supported_indices)
		filename = files[file_index]

		# Retrieve sample info and check filter requirements
		element_names = append!([], keys(sample_info_metadata), keys(filters_must_contain))
		sample_info = retrieve_sample_info(pathin, filename, element_names)
		filters_accepted = true
		for (element_name, value) in filters_must_contain
			if !contains(sample_info[element_name], value)
				filters_accepted = false
				break
			end
		end

		if filters_accepted == false
			println("Spectrum $(i[1]) of $num_of_spectra did not meet filter requirements")
			continue
		end

		# Read spectrum
		spectrum = import_files_light(pathin,[filename])
		push!(spectra, spectrum)

		# Add metadata
		spectrum["MS1"]["Filename"] = filename
		merge!(spectrum["MS1"], Dict(header => sample_info[element_name] for (element_name, header) in sample_info_metadata))
		spectrum["MS1"]["Folder Name"] = split(pathin, "\\")[end]
		

		println("Read spectra $(i[1]) of $num_of_spectra")
	end

	return spectra, metadata_headers
end

function retrieve_sample_info(pathin, filename, element_names)
	# Path to xml file which contains Sample Name
	folder_name = filename[begin:end-5] * "D"
	folder_loc = joinpath(pathin, folder_name)
	sample_info_file = joinpath(folder_loc, "AcqData\\sample_info.xml")

	# File does not exist
	if !isfile(sample_info_file)
		@warn "$(sample_info_file) does not exist, metadata will not be saved"
		return ""
	end

	# Load xml file into array
	xdoc = parse_file(sample_info_file)
	xroot = root(xdoc)
	xarray = xroot["Field"]

	
	# Find index where element is stored and retrieve its value
	elements = Dict{String, String}()
	for element_name in element_names
		element_index = findfirst(i -> content(find_element(i, "Name")) == element_name, xarray)
		element_value = content(find_element(xarray[element_index], "Value"))
		elements[element_name] = element_value
	end

	# Free memory
	free(xdoc)

	return elements
end
