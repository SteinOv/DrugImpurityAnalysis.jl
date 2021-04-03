using MS_Import
using Plots
using Printf
using CSV
using DataFrames
# using StaticArrays

# Minimum intensity of cocaine to process sample
MIN_INTENSITY = 10^5



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
				sample_profile[2] = 1000
				continue
			# Does not contain major compound (cocaine)
			elseif row[5] == -1
				sample_profile[2:end] .= 0
				break
			end
			
			sample_profile[j + 1] = round(Int, intensity/major_intensity * 1000)
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
function cocaine_normalized_list(spectra)
"""Print intensities of mz peaks of cocaine for each spectrum"""

	normalized_list = zeros(Int16, length(spectra))

	for i=1:length(spectra)
		spectrum = spectra[i]["MS1"]
		println("Spectrum $i: \t Mass \t   Intensity     Normalised Intensity")

		mass_intensity = cocaine_intensity(spectrum)
		norm = mass_intensity[1, 2]
		
		for (mass, intensity) in zip(mass_intensity[:, 1], mass_intensity[:, 2])
			# println("$mass: $intensity")
			@printf("\t\t%3.2f\t| %10.3E  |  %4i\n", mass, intensity, intensity/norm)

			if mass == 182.12
				normalized_list[i] = round(Int16, intensity/norm * 1000)
			end

		end

		println("----------------------------")
	end

	println("Normalized list (182.12):")
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

	# Minimum peak height, otherwise intensity is set to 0
	noise_cutoff = 1500 # TODO Maybe determined dynamically

	max_RT_deviation = 0.1 # minutes
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

	mass_integral = zeros(Float64, length(mass_vals), 2)

	# Below noise cutoff, set intensity to zero
	if max_intensity < noise_cutoff
		mass_integral[:, 1] .= mass_vals
		return mass_integral
	end

	peak_range = zeros(Int16, 2)

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
		integral = sum(spectrum_XIC[peak_range[1]:peak_range[2]])
		mass_integral[i,:] = [mass, integral]
		
	end

	return mass_integral
end

function find_end_of_peak(spectrum_XIC, max_intensity, max_index, direction, noise_cutoff)

	current_intensity = max_intensity
	last_intensity = 0
	current_index = max_index
	recurring_intensity_incr = 0
	reached_noise = false
	above_noise_cutoff = true

	while recurring_intensity_incr < 5 && !reached_noise
		current_index += direction

		if current_index == 0 || current_index == length(spectrum_XIC)
			reached_noise = true
			current_index -= direction 
			break
		end
		current_intensity = spectrum_XIC[current_index]

		# Intensity decreasing and above noise cut off
		if current_intensity < last_intensity && above_noise_cutoff && current_intensity > noise_cutoff
			recurring_intensity_incr = 0
		# Intensity increasing and above noise cut off
		elseif current_intensity > last_intensity && above_noise_cutoff
			recurring_intensity_incr += 1
		# Below noise cut off for the first time
		elseif current_intensity < noise_cutoff && above_noise_cutoff
			above_noise_cutoff = false
		# Below noise cut off and intensity increasing (reached noise)
		elseif !above_noise_cutoff && current_intensity > last_intensity
			reached_noise = true
		end
		last_intensity = current_intensity
	end

	# Peak overlap
	if !reached_noise
		current_index -= recurring_intensity_incr * direction
		peak_overlap = true
	else
		peak_overlap = false
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

