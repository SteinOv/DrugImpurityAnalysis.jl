using MS_Import
using Plots
using Printf
using CSV
using DataFrames
using StaticArrays

using BenchmarkTools

# Minimum intensity of cocaine to process sample
MIN_INTENSITY = 10^5
NORM_CONSTANT = 10000



"""
TODO
- Major compound (cocaine) must be the first compound in compounds.csv,
  could implement sorting of DataFrame to allow different order in compounds.csv
"""
function main()

	# Specify path
	data_folder = "P:/Git/bachelor_project/data"
	subfolder = "mzxml"
	pathin = joinpath(data_folder, subfolder)
	csvout = joinpath(pathin, "impurity_profile.csv")

	# Import spectra
	spectra = batch_import(pathin)

	# Import RT and mz info of valid compounds into DataFrame
	compounds = CSV.read("compounds.csv", DataFrame)
	filter!(row -> !(any(ismissing, (row.RT, row.mz)) || any((row.RT, row.mz) .== 0)), compounds)


	# Create DataFrame for storing impurity profile (output)
	imp_profile = DataFrame()
	insertcols!(imp_profile, :sample => String[])
	for name in compounds[!, "compound"]
		insertcols!(imp_profile, Symbol(name) => Int16[])
	end

	# Analyse all spectra
	sample_profile = Vector(undef, size(imp_profile, 2))
	for i=1:length(spectra)

		println("Analysing spectrum $i...")
		spectrum = spectra[i]["MS1"]
		sample_profile[1] = spectrum["Filename"][1]
		# sample_profile = zeros(Int16, size(compounds, 1))
		major_intensity = 0

		# Integrate peaks of all compounds
		for (j, row) in enumerate(eachrow(compounds))
			compound = row[1]
			RT = row[2]

			mass_vals = split(row[3], ";")
			mass_vals = [parse(Float32, mass) for mass in mass_vals]

			# Integrate all mz values
			mass_integral = integrate_peaks(spectrum, RT, mass_vals)


			# Determine final intensity for use in RT_deviation
			intensity = determine_intensity(mass_integral)

			# For major compound (cocaine)
			if row[5] == -1 && intensity > MIN_INTENSITY
				major_intensity = intensity
				sample_profile[2] = NORM_CONSTANT
				continue
			# Does not contain major compound (cocaine)
			elseif row[5] == -1
				sample_profile[2:end] .= 0
				break
			end
			
			sample_profile[j + 1] = round(Int, intensity/major_intensity * NORM_CONSTANT)
		end

		push!(imp_profile, sample_profile)
		println("Done")
	end

	CSV.write(csvout, imp_profile)


	# spectrum = spectra[4]["MS1"]




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
		norm = mass_intensity[1, 2]
		
		if norm == 0
			continue
		end

		println("Spectrum $i: \t Mass \t   Intensity     Normalised Intensity")

		mass_intensity = integrate_peaks(spectrum, RT, mass_vals)


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

	# Defines noise intensity
	noise_cutoff = 1000 # TODO Maybe determined dynamically

	max_RT_deviation = 0.08 # minutes
	max_mass_deviation = 0.5 

	RT_range = [RT - max_RT_deviation, RT + max_RT_deviation]

	RT_range_index = RT_indices(spectrum, RT_range)

	# Determine peak range based on first mass
	mass = mass_vals[1]
	mass_range = [mass - max_mass_deviation, mass + max_mass_deviation]
	spectrum_XIC = filter_XIC(spectrum, mass_range)

	# Determine maximum within RT and mass range
	(max_intensity, max_index) = findmax(spectrum_XIC[RT_range_index[1]:RT_range_index[2]])
	max_index += RT_range_index[1] - 1

	mass_integral = @MMatrix zeros(Float64, length(mass_vals), 2)

	# Below noise cutoff, set intensity to zero
	if max_intensity < noise_cutoff
		mass_integral[:, 1] .= mass_vals
		return mass_integral
	end

	peak_range = @MVector zeros(Int16, 2)

	for (j, direction) in enumerate([-1, 1])
		peak_end = find_end_of_peak(spectrum_XIC, max_intensity, max_index, direction, noise_cutoff)
		peak_range[j] = peak_end[1]
		if peak_end[2]
			println("WARNING: Peak overlap at file: $(spectrum["Filename"][1])")
			println("RT: $(spectrum["Rt"][peak_end[1]]), mass: $mass, index: $(peak_end[1])")
		end
	end


	# Integrate peak for all mz values
	for (i, mass) in enumerate(mass_vals)
		mass_integral[i, 1] = mass
		mass_range = [mass - max_mass_deviation, mass + max_mass_deviation]
		spectrum_XIC = filter_XIC(spectrum, mass_range)
		
		integral = sum(spectrum_XIC[peak_range[1]:peak_range[2]])
		mass_integral[i,2] = integral
			
	end



	return mass_integral
end

function find_end_of_peak(spectrum_XIC, max_intensity, max_index, direction, noise_cutoff)

	last_intensity = max_intensity
	current_index = max_index
	below_noise_cutoff = false
	reached_noise = false

	peak_overlap = false

	# For storing last intensity change
	history_size = 5
	last_intensity_changes = @MVector zeros(Int64, history_size)

	i = 0
	while current_index > 1 && current_index < (length(spectrum_XIC) - 1)
		i = i % history_size + 1
		current_index += direction

		current_intensity = spectrum_XIC[current_index]
		mean_intensity_change = (last_intensity_changes[1] + last_intensity_changes[end]) / 2

		# Large intensity increase, peak overlap
		if mean_intensity_change > 20*noise_cutoff
			peak_overlap = true
			current_index -= (history_size - 1) * direction
			break
		# Below noise cut off and intensity decreasing
		elseif current_intensity < noise_cutoff && current_intensity < last_intensity
			below_noise_cutoff = true
		# Below noise cut off and intensity increasing
		elseif below_noise_cutoff && current_intensity > last_intensity
			reached_noise = true
			current_index -= direction
			break
		end
		
		last_intensity_changes[i] = current_intensity - last_intensity
		last_intensity = current_intensity
		
	end

	return current_index, peak_overlap
end

	


function plt(mass=0, mass2=0, RT=0, RT_deviation=0.1)
	"""Plot spectrum of given mass
	plt((mass), (mass2), (RT), (RT_deviation))

	only mass given: XIC at specific mass
	mass and mass2: XIC between range
	mass = 0: Spectrum not filtered

	RT given: Zooms in around given RT
	RT = 0: Shows complete RT spectrum
	RT_deviation: Deviation around RT (default 0.1)
	"""

	# TODO gives error at start and end

	max_mass_deviation = 0.5 
	RT_range = [RT - RT_deviation, RT + RT_deviation]
	RT_range_index = RT_indices(spectrum, RT_range)

	if mass > 0 && mass2 <= 0
		mass_range = [mass - max_mass_deviation, mass + max_mass_deviation]
	elseif mass > 0
		mass_range = [mass, mass2]
	else
		mass_range = [0, Inf]
	end

	spectrum_XIC = filter_XIC(spectrum, mass_range)

	if RT > 0
		add_left = 100
		add_right = 100
		plot(spectrum["Rt"][RT_range_index[1] - add_left:RT_range_index[2] + add_right], spectrum_XIC[(RT_range_index[1] - add_left):(RT_range_index[2] + add_right)])
	else
		plot(spectrum["Rt"], spectrum_XIC)
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

function filter_XIC(spectrum, mass)
	"""Returns filtered spectrum based on mass range"""

	
	# Store which indices are between mass range
	filter = findall(mz -> (mz .>= mass[1]) .& (mz .<= mass[2]), spectrum["Mz_values"])
	filter = getindex.(filter, [1 2])

	# Create XIC spectrum
	spectrum_XIC = @MVector zeros(Int, size(spectrum["Mz_values"], 1))
	for (row, column) in zip(filter[:, 1], filter[:, 2])
		spectrum_XIC[row] += spectrum["Mz_intensity"][row, column]
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
		filename = [files[file_index]]
		spectra[i] = import_files(pathin,filename,mz_thresh,Int_thresh)
		spectra[i]["MS1"]["Filename"] = filename
		println("Read spectra $(i[1]) of $num_of_spectra")
	end

	return spectra
end

main()

