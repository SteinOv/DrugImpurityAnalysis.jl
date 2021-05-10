using MS_Import
using Plots
using Printf
using CSV
using DataFrames
using StaticArrays
using LightXML

using BenchmarkTools

include("helpers.jl")

const NORM_CONSTANT = 1000000
const MAX_MASS_DEVIATION = 0.5
const MAX_RT_DEVIATION = 0.05

# 40 scans is 0.1 minute
const MAX_PEAK_SCANS_DETERMINATION = 20 # Used for differentiating peak from noise
const MAX_PEAK_SCANS_INTEGRATION_LEFT = 40 # Bounds for peak integration (left)
const MAX_PEAK_SCANS_INTEGRATION_RIGHT = 80 # Bounds for peak integration (right)


const NOISE_INTERVAL_SIZE = 20 # Size of interval to search for noise
const NOISE_CROSSINGS_FRACTION = 0.33
const NOISE_ZEROS_FRACTION = 0.8
const OVERLAP_BELOW_MEDIAN = 3



"""
TODO
- Major compound (cocaine) must be the first compound in compounds.csv,
  could implement sorting of DataFrame to allow different order in compounds.csv
- Bring mass and mz under one term
"""
function main()

	# Specify path
	data_folder = joinpath(@__DIR__, "data")
	subfolder = "coca_caf_cal"
	pathin = joinpath(data_folder, subfolder)
	csvout = joinpath(pathin, "impurity_profile.csv")

	# Import spectra
	spectra = batch_import(pathin)
	

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

		# Check if major compound sufficiently present
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

			# Parse mz values
			mass_str = split(row.mz, ";")
			mass_values = Array{Any,1}(undef, length(mass_str))
			for (i, mass) in enumerate(mass_str)
				if startswith(mass, "(")
					mass = split(mass[2:end - 1], ",")
					mass_values[i] = Tuple(([parse(Float32, sub_mass) for sub_mass in mass]))
				else
					mass_values[i] = parse(Float32, mass)
				end	
			end

			overlap_RT = ismissing(row.overlap) ? 0 : row.overlap*RT_modifier # TODO support 2 overlap RT's or return error

			# Integrate all mz values
			mass_integral = integrate_mz_values(spectrum, RT, mass_values, overlap_RT)

			# Determine final intensity
			intensity = determine_intensity(mass_integral)
			
			sample_profile[j + 1] = round(Int, intensity/major_intensity * NORM_CONSTANT)
		end

		push!(imp_profile, append!(sample_metadata, sample_profile))
		println("Done")
	end

	CSV.write(csvout, imp_profile)

end

function process_major(compounds, spectrum)
	"""Returns intensity of major compound, RT modifier (for RT shift correction) and name of major compound"""
	major_compound = filter(row -> row.type_int == -1, compounds)
	internal_standard = filter(row -> row.type_int == -10, compounds)

	major_compound_name = major_compound.compound[1]
	RT = major_compound.RT[1]
	IS_RT = internal_standard.RT[1]
	major_IS_ratio = parse(Float32, internal_standard.ratio[1])[1]

	mass_vals = split(major_compound.mz[1], ";")
	mass_vals = [parse(Float32, mass) for mass in mass_vals]
	IS_mass_vals = split(internal_standard.mz[1], ";")
	IS_mass_vals = [parse(Float32, mass) for mass in IS_mass_vals]

	# Find peak using highest mz value within +/- 1 min of predicted RT
	highest_mass = maximum(mass_vals)
	RT_range = [RT - 1, RT + 1]
	RT_range_index = RT_indices(spectrum, RT_range)
	maximum_RT = findmax(filter_XIC(spectrum, highest_mass)[RT_range_index[1]:RT_range_index[2]])
	real_RT = spectrum["Rt"][maximum_RT[2] + RT_range_index[1] - 1]

	# Compute RT modifier
	RT_modifier = round(real_RT / RT, digits = 3)
	RT = real_RT

	# Determine cocaine/IS ratio
	highest_mz_integral = integrate_mz_values(spectrum, RT, highest_mass, 0)[1, 2]
	IS_RT = IS_RT * RT_modifier
	IS_highest_mz_integral = integrate_mz_values(spectrum, IS_RT, IS_mass_vals[1], 0)[1, 2]


	# Ratio too low, no significant amount of major compound present 
	if highest_mz_integral / IS_highest_mz_integral < major_IS_ratio || highest_mz_integral == 0
		return (0, RT_modifier, major_compound_name)
	elseif IS_highest_mz_integral == 0 && highest_mz_integral != 0
		@warn "Major compound detected, but internal standard not detected"
		return (0, RT_modifier, major_compound_name)
	end



	# Integrate all mz values
	mass_integral = integrate_mz_values(spectrum, RT, mass_vals, 0)

	# Required ratio between first mz IS and biggest mass major compound is defined in compounds.csv
	highest_mass_i = findfirst(x -> x == highest_mass, mass_integral[:, 1])
	hi_mz_intensity = mass_integral[highest_mass_i, 2]

	


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


function integrate_mz_values(spectrum, RT, mz_vals, overlap_RT=0)
	"""
	Integrates all mz values of specific compound

	Returns matrix mass_integral, where row 1 are the mz values
	and row 2 the corresponding intensities
	"""
	
	RT_range = [RT - MAX_RT_DEVIATION, RT + MAX_RT_DEVIATION]
	RT_range_index = RT_indices(spectrum, RT_range)

	mass_integral = Array{Any,2}(undef, length(mz_vals), 2)

	for (i, mz) in enumerate(mz_vals)
		mass_integral[i, 1] = mz
		spectrum_XIC = filter_XIC(spectrum, mz)

		left_index, right_index, noise_median = determine_peak_info(spectrum_XIC, RT_range_index, overlap_RT)

		if left_index == -1
			# No peak
			mass_integral[i, 2] = 0
		else
			# Integrate peak
			spectrum_part = spectrum_XIC[left_index:right_index]
			replace!(x -> x < noise_median ? 0 : x - noise_median, spectrum_part)
			mass_integral[i, 2] = sum(spectrum_part)
		end

	end

	return mass_integral
end

function determine_peak_info(spectrum_XIC, RT_range_index, overlap_RT)
	"""Determines peak range and noise median"""

	peak, left_index, right_index, max_index = find_peak(RT_range_index, spectrum_XIC)
	
	# No peak
	if peak == -1
		return -1, -1, -1
	# Peak without noise
	elseif peak == 0
		return left_index, right_index, 0
	end

	# Broader bounds for integration
	left_index = max_index - round(Int, MAX_PEAK_SCANS_INTEGRATION_LEFT)
	right_index = max_index + round(Int, MAX_PEAK_SCANS_INTEGRATION_RIGHT)
	
	# Check for overlap and adjust bounds if overlap with other peak found
	if overlap_RT > 0
		overlap_RT_range = (overlap_RT - MAX_RT_DEVIATION, overlap_RT + MAX_RT_DEVIATION)
		overlap_RT_range_index = RT_indices(spectrum, overlap_RT_range)
	else
		overlap_RT_range_index = 0
	end
	peak_end = determine_overlap(spectrum_XIC, left_index, right_index, max_index, overlap_RT_range_index)
	left_index = peak_end[1] == -1 ? left_index : peak_end[1]
	right_index = peak_end[2] == -1 ? right_index : peak_end[2]

	# Search for closest noise, left if no overlap defined on left
	direction = overlap_RT > 0 && overlap_RT < RT ? 1 : -1
	noise_end = direction == -1 ? left_index : right_index + NOISE_INTERVAL_SIZE - 1
	noise_start = noise_end - (NOISE_INTERVAL_SIZE - 1)
	noise_median = Int
	# rep_count = 0 # TEMP
	while noise_start > 0 && noise_end < length(spectrum_XIC)
		noise_median = determine_noise(spectrum_XIC[noise_start:noise_end], true)

		# Noise found
		if noise_median >= 0
			break
		end
		
		# rep_count += 1 # TEMP
		noise_start += NOISE_INTERVAL_SIZE * direction
		noise_end += NOISE_INTERVAL_SIZE * direction
	end
	# println("rep_count: $rep_count, mz: $mz, max_index: $max_index") # TEMP

	if noise_median == -1
		error("No noise found (mz: $mz)")
	end

	# Substract noise median around peak max
	spectrum_part = spectrum_XIC[left_index:right_index]
	replace!(x -> x < noise_median ? 0 : x - noise_median, spectrum_part)
	index_shift = left_index - 1

	# Find minimum left and right closest to max
	max_index = findmax(spectrum_part)[2]
	min_left = minimum(spectrum_part[1:max_index])
	min_right = minimum(spectrum_part[max_index:end])
	right_index = peak_end[2] == -1 ? findfirst(x -> x == min_right, spectrum_part[max_index:end]) + max_index - 1 : peak_end[2] - left_index + 1
	left_index = peak_end[1] == -1 ? findlast(x -> x == min_left, spectrum_part[1:max_index]) : peak_end[1] - left_index + 1

	return left_index + index_shift, right_index + index_shift, noise_median
end


function find_peak(RT_range_index, spectrum)
	"""
	Finds highest peak within given RT range in spectrum
	Returns -1 if no peak is present, 0 if peak is present without noise and 
	1 if peak is present with noise, together with left, right and max index.
	"""

	# Determine maximum within RT range
	(max_intensity, max_index) = findmax(spectrum[RT_range_index[1]:RT_range_index[2]])
	max_index += RT_range_index[1] - 1

	# Define bounds around maximum
	left_index = max_index - round(Int, MAX_PEAK_SCANS_DETERMINATION / 2)
	right_index = max_index + round(Int, MAX_PEAK_SCANS_DETERMINATION / 2)

	# A lot of intensities of zero with a peak in between
	if determine_noise(spectrum[left_index:right_index]) == 0 && max_intensity > 0
		left_index = findlast(x -> x == 0, spectrum[left_index:max_index]) + left_index - 1
		right_index = findfirst(x -> x == 0, spectrum[max_index:right_index]) + max_index - 1
		peak = 0 # Peak without noise
	# No peak
	elseif max_intensity == 0
		peak = -1 # No peak
		return (peak, left_index, right_index, max_index)
	# Peak present, with noise
	else
		min_left = minimum(spectrum[left_index:max_index])
		min_right = minimum(spectrum[max_index:right_index])
		left_index = findlast(x -> x == min_left, spectrum[left_index:max_index]) + left_index - 1
		right_index = findfirst(x -> x == min_right, spectrum[max_index:right_index]) + max_index - 1
		peak = 1 # Peak with noise
	end

	# Noise instead of peak
	if determine_noise(spectrum[left_index:right_index]) >= 0
		peak = -1 # No peak
	end

	return (peak, left_index, right_index, max_index)
end



function determine_overlap(spectrum_XIC, left_index, right_index, max_index, overlap_RT_range_index)
	"""
	Determines whether peak has overlap.
	Returns index of minimum between peak and overlapping peak if overlap present.
	"""
	# peak_end[1] is overlap on left, peak_end[2] is overlap on right
	peak_end = Int32[-1, -1]

	# Overlap defined in compounds.csv # TODO, for now only supports one RT
	if overlap_RT_range_index != 0
		overlap_peak, overlap_left_index, overlap_right_index, overlap_max_index = find_peak(overlap_RT_range_index, spectrum_XIC)

		# Peak exists
		if overlap_peak >= 0
			overlap_on_right = max_index < overlap_max_index
			between_peaks = overlap_on_right ? spectrum_XIC[max_index:overlap_max_index] : spectrum_XIC[overlap_max_index:max_index]

			# Index of minimum closest to peak
			between_min = overlap_on_right ? 
						findfirst(x -> x == minimum(between_peaks), between_peaks) + max_index - 1 : 
						findlast(x -> x == minimum(between_peaks), between_peaks) + overlap_max_index - 1
			
			if overlap_on_right
				peak_end[2] = between_min
			else
				peak_end[1] = between_min
			end
		end
	end

	# Check for overlap based on median
	under_median_count = 0
	for (j, index_range) in enumerate(((max_index, left_index), (max_index, right_index)))

		# Peak_end already defined
		if peak_end != -1
			continue
		end
		range_order = index_range[1] < index_range[2] ? (index_range[1], index_range[2]) : (index_range[2], index_range[1])
		range_median = median(spectrum_XIC[range_order[1]:range_order[2]])

		for k in index_range[1]:index_range[2]
			# consecutively below median and subsequently above median
			if spectrum_XIC[k] > range_median && under_median_count >= OVERLAP_BELOW_MEDIAN
				spectrum_part = k < max_index ? spectrum_XIC[k:max_index] : spectrum_XIC[max_index:k]
				peak_end[j] = k < max_index ? findmin(spectrum_part)[2] + k - 1 : findmin(spectrum_part)[2] + max_index - 1
				break
			# Below median
			elseif spectrum_XIC[k] < range_median
				under_median_count += 1
			# Above median, has not been below median consecutively
			else
				under_median_count = 0
			end
		end
	end

	return peak_end
end

function determine_noise(spectrum_part, use_median=false)
	"""
	Returns mean or median of noise or -1 if peak detected
	Uses mean by default, uses median if use_median=true
	"""
	

	# Large amount of intensities of zero
	if count(x -> x == 0, spectrum_part) / length(spectrum_part) > NOISE_ZEROS_FRACTION
		return 0
	end

	if !use_median
		spectrum_part_avg = mean(spectrum_part)
	else
		spectrum_part_avg = median(spectrum_part)
	end

	# Count number of mean/median crossings
	crossings_count = 0
	for i in 2:length(spectrum_part)
		crosses_mean = (spectrum_part[i] - spectrum_part_avg) * (spectrum_part[i - 1] - spectrum_part_avg) <= 0
		if crosses_mean
			crossings_count += 1
		end
	end

	# Small fraction of mean/median crossings, peak instead of noise
	if crossings_count / length(spectrum_part) < NOISE_CROSSINGS_FRACTION
		return -1
	else
		return round(Int, spectrum_part_avg)
	end
end

# main()


# using LinearAlgebra, SparseArrays

# function AsLS(y, lambda, p)
# 	y = vec(y)

# 	# Estimate baseline with asymmetric least squares
# 	m = length(y)
# 	d = zeros(m,m)
# 	d[diagind(d)] .= 1
# 	D = diff(diff(d, dims = 1),dims =1)
# 	w = ones(m, 1)
# 	z = []

# 	for it = 1:10
# 		W = spdiagm(m, m, 0 => vec(w))
# 	    C = cholesky(W + lambda * D' * D)
# 	    z = C.U \ (C.U' \ (w .* y))
# 	    w = p * (y .> z) + (1 - p) * (y .< z)
# 	end

# 	return z
# end

# # Cinnamoylcocaine
# spectrum = spectra[33]["MS1"]
# spectrum_XIC = filter_XIC(spectrum, (82, 182, 83, 96, 94))[2075:2270]
# x = spectrum["Rt"][2075:2270]

# baseline = AsLS(spectrum_XIC, 20000, 0.01)
# plot(x, spectrum_XIC)
# plot!(x, baseline)
# ylims!(0, 10000)
# xlims!(7, 7.5)

# # Norcocaine
# spectrum = spectra[33]["MS1"]
# spectrum_XIC = filter_XIC(spectrum, 168)[1798:1890]
# x = spectrum["Rt"][1798:1890]

# baseline = AsLS(spectrum_XIC, 10000, 0.01)
# plot(x, spectrum_XIC)
# plot!(x, baseline)
# # ylims!(0, 10000)
# # xlims!(7, 7.5)



# plt((168,68,136), 6.5)
# plt!((168), 6.5)
# plt!((68,136), 6.5)

# # High norcocaine
# spectrum = spectra[17]["MS1"]
# ylims!(0, 100000)
# xlims!(6.25, 6.7)

# # Low norcocaine
# spectrum = spectra[29]["MS1"]
# xlims!(6.35, 6.45)
# ylims!(0, 5000)


# include("main.jl")

# spectrum = spectra[3]["MS1"]
# plt(303)
# plt()
# spectrum = spectra[33]["MS1"]
# plt(77)
# plt()