using MS_Import
using Plots
using Printf
using CSV
using DataFrames
using StaticArrays
using LightXML

using BenchmarkTools

include("manual_visualisation.jl")
include("helpers.jl")

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
		compounds_in_prof = filter(row -> !(row.type_int == -10), compounds)
		imp_profile = DataFrame()
		insertcols!(imp_profile, :item_number => String[])
		insertcols!(imp_profile, :filename => String[])
		insertcols!(imp_profile, :folder => String[])
		for name in compounds_in_prof.compound
			insertcols!(imp_profile, Symbol(name) => Int32[])
		end

		# Analyse all spectra
		sample_profile = zeros(Int32, size(compounds_in_prof, 1))
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
			if major_intensity > 0
				sample_profile[1] = NORM_CONSTANT
			else
				sample_profile[:] .= 0
				push!(imp_profile, append!(sample_metadata, sample_profile))
				println("Done (major compound not present)")
				continue
			end

			# Integrate peaks of all compounds
			for (j, row) in enumerate(eachrow(filter(row -> !(row.type_int in [-1, -10]), compounds)))
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


	# spectrum = spectra[3]["MS1"]
	# plt(303)



end

function process_major(compounds, spectrum)
	"""Returns intensity of major compound, RT modifier (for RT shift correction) and name of major compound"""
	major_compound = filter(row -> row.type_int == -1, compounds)
	internal_standard = filter(row -> row.type_int == -10, compounds)

	major_compound_name = major_compound.compound[1]
	RT = major_compound.RT[1]
	IS_RT = internal_standard.RT[1]
	IS_major_ratio = parse(Float32, internal_standard.ratio[1])[1]

	mass_vals = split(major_compound.mz[1], ";")
	mass_vals = [parse(Float32, mass) for mass in mass_vals]
	IS_mass_vals = split(internal_standard.mz[1], ";")
	IS_mass_vals = [parse(Float32, mass) for mass in IS_mass_vals]

	# Find peak using highest mz value within +/- 1 min of known RT
	highest_mass = maximum(mass_vals)
	RT_range = [RT - 1, RT + 1]
	RT_range_index = RT_indices(spectrum, RT_range)
	maximum_RT = findmax(filter_XIC(spectrum, highest_mass)[RT_range_index[1]:RT_range_index[2]])
	real_RT = spectrum["Rt"][maximum_RT[2] + RT_range_index[1] - 1]

	# Peak below noise
	if maximum_RT[1] < NOISE_CUTOFF
		return (0, 0, major_compound_name)
	end

	# Compute RT modifier
	RT_modifier = round(real_RT / RT, digits = 3)
	RT = real_RT

	IS_RT = IS_RT * RT_modifier
	IS_integral = integrate_peaks(spectrum, IS_RT, IS_mass_vals[1])[1, 2]

	# Integrate all mz values
	mass_integral = integrate_peaks(spectrum, RT, mass_vals)

	# Required ratio between first mz IS and biggest mass major compound is defined in compounds.csv
	highest_mass_i = findfirst(x -> x == highest_mass, mass_integral[:, 1])
	hi_mz_intensity = mass_integral[highest_mass_i, 2]

	# Ratio too high, no significant amount of major compound present
	if IS_integral / hi_mz_intensity > IS_major_ratio
		return (0, 0, major_compound_name)
	end

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




function determine_noise(spectrum_part)
	"""Returns median of noise or -1 if peak detected"""
	

	# More than 80% zeros, noise level of 0
	if count(x -> x == 0, spectrum_part) / length(spectrum_part) > 0.8 # TODO maybe global?
		return 0
	end

	spectr_median = median(spectrum_part)

	# Count number of median crossings
	crossings_count = 0
	for i in 2:length(spectrum_part)
		crosses_median = (spectrum_part[i] - spectr_median) * (spectrum_part[i - 1] - spectr_median) <= 0
		if crosses_median
			crossings_count += 1
		end
	end

	# Fraction of crossings is smaller than 20%, peak instead of noise
	if crossings_count / length(spectrum_part) < 0.2 # TODO maybe global?
		return -1
	else
		return spectr_median


# main()

