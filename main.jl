using MS_Import
using Plots
using Printf
using CSV
using DataFrames
using StaticArrays
using LightXML

using BenchmarkTools

# Minimum intensity of cocaine to process sample
const MIN_INTENSITY = 10^5
const NORM_CONSTANT = 10000
const MAX_MASS_DEVIATION = 0.5
const MAX_RT_DEVIATION = 0.08
const NOISE_CUTOFF = 1000 # TODO Maybe determined dynamically



"""
TODO
- Major compound (cocaine) must be the first compound in compounds.csv,
  could implement sorting of DataFrame to allow different order in compounds.csv
"""
function main()

	# Specify path
	data_folder = "P:/Git/bachelor_project/data"
	subfolder = "200006"
	pathin = joinpath(data_folder, subfolder)
	csvout = joinpath(pathin, "impurity_profile.csv")

	# Import spectra
	load_time = @elapsed spectra = batch_import(pathin)
	
	# @btime begin #  3.168 s (7339762 allocations: 277.15 MiB)
		# Import RT and mz info of valid compounds into DataFrame
		compounds = CSV.read("compounds.csv", DataFrame)
		filter!(row -> !(any(ismissing, (row.RT, row.mz)) || any((row.RT, row.mz) .== 0)), compounds)


		# Create DataFrame for storing impurity profile (output)
		imp_profile = DataFrame()
		insertcols!(imp_profile, :item_number => String[])
		insertcols!(imp_profile, :filename => String[])
		insertcols!(imp_profile, :folder => String[])
		for name in compounds[!, "compound"]
			insertcols!(imp_profile, Symbol(name) => Int32[])
		end

		# Analyse all spectra
		sample_profile = zeros(Int32, size(compounds, 1))
		for i=1:length(spectra)
			sample_metadata = Array{Any,1}(undef, 3)

			println("Analysing spectrum $i...")
			spectrum = spectra[i]["MS1"]

			# Save metadata
			sample_metadata[1] = spectrum["Sample Name"][1]
			sample_metadata[2] = spectrum["Filename"][1]
			sample_metadata[3] = subfolder

			# Determine intensity of major compound and RT shift (modifier)
			major_intensity, RT_modifier, major_compound_name = process_major(compounds, spectrum)

			# Check if major compounds sufficiently present
			if major_intensity > MIN_INTENSITY
				sample_profile[1] = NORM_CONSTANT
			else
				sample_profile[:] .= 0
				push!(imp_profile, append!(sample_metadata, sample_profile))
				println("Done (major compound not present)")
				continue
			end

			# Integrate peaks of all compounds
			for (j, row) in enumerate(eachrow(filter(row -> row.compound != major_compound_name, compounds)))
				compound = row.compound
				RT = row.RT * RT_modifier

				mass_vals = split(row.mz, ";")
				mass_vals = [parse(Float32, mass) for mass in mass_vals]

				# Integrate all mz values
				mass_integral = integrate_peaks(spectrum, RT, mass_vals)

				# Determine final intensity for use in RT_deviation
				intensity = determine_intensity(mass_integral)
				
				sample_profile[j + 1] = round(Int, intensity/major_intensity * NORM_CONSTANT)
			end

			push!(imp_profile, append!(sample_metadata, sample_profile))
			println("Done")
		end
	# end
	CSV.write(csvout, imp_profile)


	# spectrum = spectra[4]["MS1"]




end

function process_major(compounds, spectrum)
	major_compound = filter(row -> row.type_bool == -1, compounds)
	major_compound_name = major_compound.compound[1]
	RT = major_compound.RT[1]

	mass_vals = split(major_compound.mz[1], ";")
	mass_vals = [parse(Float32, mass) for mass in mass_vals]

	# Find peak using highest mz value within +/- 1 min of known RT
	mass = maximum(mass_vals)
	RT_range = [RT - 1, RT + 1]
	RT_range_index = RT_indices(spectrum, RT_range)
	maximum_RT = findmax(filter_XIC(spectrum, mass)[RT_range_index[1]:RT_range_index[2]])
	real_RT = spectrum["Rt"][maximum_RT[2] + RT_range_index[1] - 1]

	# Peak below noise
	if maximum_RT[1] < NOISE_CUTOFF
		return (0, 0, major_compound_name)
	end

	# Compute RT modifier
	RT_modifier = round(real_RT / RT, digits = 3)
	RT = real_RT

	# Integrate all mz values
	mass_integral = integrate_peaks(spectrum, RT, mass_vals)

	# Determine final intensity for use in RT_deviation
	intensity = determine_intensity(mass_integral)

	return (intensity, RT_modifier, major_compound_name)
end


function determine_intensity(mass_integral)
	"""
	Determines final intensity for use in ratios
	Just sums all integrals for now
	"""
	return sum(mass_integral[:, 2])
end


# TEMPORARY
function mass_ratios(compound, spectra=spectra, compounds=compounds)
"""Print intensities of mz peaks of specific compound for each spectrum"""


	compound_row = findfirst(comp -> comp == compound, compounds.compound)

	normalized_list = zeros(Int16, length(spectra))

	row = compounds[compound_row, :]
	mass_vals = split(row[3], ";")
	mass_vals = [parse(Float32, mass) for mass in mass_vals]
	RT = row[2]

	for i=1:length(spectra)
		spectrum = spectra[i]["MS1"]
		RT_modifier = process_major(compounds, spectrum)[2]

		if RT_modifier == 0
			continue
		end
		mass_intensity = integrate_peaks(spectrum, RT * RT_modifier, mass_vals)
		norm = mass_intensity[1, 2]

		if norm == 0
			continue
		end

		println("Spectrum $i: \t Mass \t   Intensity     Normalised Intensity")

		


		for (mass, intensity) in zip(mass_intensity[:, 1], mass_intensity[:, 2])
			@printf("\t\t%3.2f\t| %10.3E  |  %4i\n", mass, intensity, intensity/norm * 1000)

			if mass == mass_vals[1]
				normalized_list[i] = round(Int16, intensity/norm * 1000)
			end

		end

		println("----------------------------")
	end

	println("Normalized list second mass:")
	for i in normalized_list
		@printf("%4i\n", i)
	end

end


# TEMPORARY
function cocaine_intensity(spectrum)
	"""Returns intensity of cocaine peak if present, else 0"""

	# Cocaine approximate peak locations
	RT = 6.66 # minutes
	mass_vals = [82.07, 182.12, 83.07, 94.07, 77.04, 105.03, 96.08, 42.03, 303.15]
	
	return integrate_peaks(spectrum, RT, mass_vals)

end


"""
TODO
- Improve overlapping peak integration by predicting actual integral after peak overlap
- Noise cut off determined dynamically
- Peak shape should be the same for each mass, only determine peak range for one mass

"""
function integrate_peaks(spectrum, RT, mass_vals)
	"""
	Integrates peaks of specific compound

	Returns matrix mass_integral, where row 1 are the mz values
	and row 2 the corresponding intensities
	"""

	RT_range = [RT - MAX_RT_DEVIATION, RT + MAX_RT_DEVIATION]

	RT_range_index = RT_indices(spectrum, RT_range)

	# Determine peak range based on first mass
	mass = mass_vals[1]
	spectrum_XIC = filter_XIC(spectrum, mass)

	# Determine maximum within RT and mass range
	(max_intensity, max_index) = findmax(spectrum_XIC[RT_range_index[1]:RT_range_index[2]])
	max_index += RT_range_index[1] - 1

	mass_integral = @MMatrix zeros(Float64, length(mass_vals), 2)

	# Below noise cutoff, set intensity to zero
	if max_intensity < NOISE_CUTOFF
		mass_integral[:, 1] .= mass_vals
		return mass_integral
	end

	peak_range = @MVector zeros(Int16, 2)

	for (j, direction) in enumerate([-1, 1])
		peak_end = find_end_of_peak(spectrum_XIC, max_intensity, max_index, direction)
		peak_range[j] = peak_end[1]
		if peak_end[2]
			println("WARNING: Peak overlap at file: $(spectrum["Filename"][1])")
			println("RT: $(spectrum["Rt"][peak_end[1]]), mass: $mass, index: $(peak_end[1])")
		end
	end


	# Integrate peak for all mz values
	for (i, mass) in enumerate(mass_vals)
		mass_integral[i, 1] = mass
		spectrum_XIC = filter_XIC(spectrum, mass)
		
		integral = sum(spectrum_XIC[peak_range[1]:peak_range[2]])
		mass_integral[i,2] = integral
			
	end



	return mass_integral
end

function find_end_of_peak(spectrum_XIC, max_intensity, max_index, direction)

	last_intensity = max_intensity
	current_index = max_index
	below_noise_cutoff = false
	reached_noise = false

	peak_overlap = false

	# For storing last intensity changes
	history_size = 5
	last_intensity_changes = @MVector zeros(Int64, history_size)

	i = 0
	while current_index > 1 && current_index < (length(spectrum_XIC) - 1)
		i = i % history_size + 1
		current_index += direction

		current_intensity = spectrum_XIC[current_index]
		mean_intensity_change = (last_intensity_changes[1] + last_intensity_changes[end]) / 2

		# Large intensity increase, peak overlap
		if mean_intensity_change > 20*NOISE_CUTOFF
			peak_overlap = true
			current_index -= (history_size - 1) * direction
			break
		# Below noise cut off and intensity decreasing
		elseif current_intensity < NOISE_CUTOFF && current_intensity < last_intensity
			below_noise_cutoff = true
		# Below noise cut off and intensity increasing or reached zero
		elseif (below_noise_cutoff && current_intensity > last_intensity) || current_intensity <= 0
			reached_noise = true
			current_index -= direction
			break
		end
		
		last_intensity_changes[i] = current_intensity - last_intensity
		last_intensity = current_intensity
		
	end

	return current_index, peak_overlap
end


function plt(mass=0, RT=0, RT_range_size=0.5)
	"""Plot spectrum of given mass
	plt((mass), (RT), (RT_deviation))

	mass: XIC at specific mass
	mass = 0: Spectrum not filtered

	RT given: Zooms in around given RT
	RT = 0: Shows complete RT spectrum
	RT_range_size: Deviation around RT (default 0.5)
	"""

	# TODO gives error at start and end

	RT_range = [RT - RT_range_size, RT + RT_range_size]
	RT_range_index = RT_indices(spectrum, RT_range)

	if mass != 0
		spectrum_XIC = filter_XIC(spectrum, mass)
	else
		spectrum_XIC = filter_XIC(spectrum, 0)
	end

	if RT > 0
		add_left = 100
		add_right = 100
		plot(spectrum["Rt"][RT_range_index[1]:RT_range_index[2]], spectrum_XIC[(RT_range_index[1]):RT_range_index[2]])
	else
		plot(spectrum["Rt"], spectrum_XIC)
	end
end

function plt!(mass=0, RT=0, RT_range_size=0.5)
	"""Plot spectrum of given mass
	plt((mass), (RT), (RT_deviation))

	mass: XIC at specific mass
	mass = 0: Spectrum not filtered

	RT given: Zooms in around given RT
	RT = 0: Shows complete RT spectrum
	RT_range_size: Deviation around RT (default 0.5)
	"""

	# TODO gives error at start and end

	RT_range = [RT - RT_range_size, RT + RT_range_size]
	RT_range_index = RT_indices(spectrum, RT_range)

	if mass != 0
		spectrum_XIC = filter_XIC(spectrum, mass)
	else
		spectrum_XIC = filter_XIC(spectrum, 0)
	end

	if RT > 0
		add_left = 100
		add_right = 100
		plot!(spectrum["Rt"][RT_range_index[1]:RT_range_index[2]], spectrum_XIC[(RT_range_index[1]):RT_range_index[2]])
	else
		plot!(spectrum["Rt"], spectrum_XIC)
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


# main()
