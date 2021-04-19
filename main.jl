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

const NORM_CONSTANT = 1000000
const MAX_MASS_DEVIATION = 0.5
const MAX_RT_DEVIATION = 0.08

# 40 scans is 0.1 minute
const MAX_PEAK_SCANS_DETERMINATION = 20 # Used for differentiating peak from noise
const MAX_PEAK_SCANS_INTEGRATION = 60 # Bounds for peak integration

const NOISE_INTERVAL_SIZE = 20 # Size of interval to search for noise
const NOISE_CROSSINGS_FRACTION = 0.33
const NOISE_ZEROS_FRACTION = 0.8
const OVERLAP_BELOW_MEDIAN = 3



"""
TODO
- Major compound (cocaine) must be the first compound in compounds.csv,
  could implement sorting of DataFrame to allow different order in compounds.csv
"""
function main()

	# Specify path
	# data_folder = "P:/Git/bachelor_project/data"
	data_folder = joinpath(@__DIR__, "data")
	subfolder = "200006"
	pathin = joinpath(data_folder, subfolder)
	csvout = joinpath(pathin, "impurity_profile.csv")

	# Import spectra
	load_time = @elapsed spectra = batch_import(pathin)
	
	# @btime begin # 3.984 s (9058147 allocations: 291.28 MiB) (folder: 200006)
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

				# For tuples of mz values
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

				# append!(mass_values, [parse(Float32, mass) for mass in mass_vals])

				# Integrate all mz values
				mass_integral = integrate_peaks(spectrum, RT, mass_values)

				# Determine final intensity for use in RT_deviation
				intensity = determine_intensity(mass_integral)
				
				sample_profile[j + 1] = round(Int, intensity/major_intensity * NORM_CONSTANT)
			end

			push!(imp_profile, append!(sample_metadata, sample_profile))
			println("Done")
		end
	# end # @btime
	CSV.write(csvout, imp_profile)


	# spectrum = spectra[3]["MS1"]
	# plt(303)
	# plt()
	# spectrum = spectra[33]["MS1"]
	# plt(77)
	# plt()



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

	# Compute RT modifier
	RT_modifier = round(real_RT / RT, digits = 3)
	RT = real_RT

	# Determine IS/cocaine ratio
	highest_mz_intensity = integrate_peaks(spectrum, RT, highest_mass)[1, 2]
	IS_RT = IS_RT * RT_modifier
	IS_integral = integrate_peaks(spectrum, IS_RT, IS_mass_vals[1])[1, 2]
	# Ratio too high, no significant amount of major compound present 
	if IS_integral / highest_mz_intensity > IS_major_ratio || highest_mz_intensity == 0
		return (0, 0, major_compound_name)
	elseif IS_integral == 0
		error("Cocaine peak detected, but internal standard not detected")
	end



	# Integrate all mz values
	mass_integral = integrate_peaks(spectrum, RT, mass_vals)

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


function integrate_peaks(spectrum, RT, mz_vals)
	"""
	Integrates peaks of specific compound

	Returns matrix mass_integral, where row 1 are the mz values
	and row 2 the corresponding intensities
	"""
	
	RT_range = [RT - MAX_RT_DEVIATION, RT + MAX_RT_DEVIATION]
	RT_range_index = RT_indices(spectrum, RT_range)

	mass_integral = Array{Any,2}(undef, length(mz_vals), 2)

	for (i, mz) in enumerate(mz_vals)
		mass_integral[i, 1] = mz
		spectrum_XIC = filter_XIC(spectrum, mz)

		# Determine maximum within RT and mass range
		(max_intensity, max_index) = findmax(spectrum_XIC[RT_range_index[1]:RT_range_index[2]])
		max_index += RT_range_index[1] - 1

		# Define bounds around maximum
		left_index = max_index - round(Int, MAX_PEAK_SCANS_DETERMINATION / 2)
		right_index = max_index + round(Int, MAX_PEAK_SCANS_DETERMINATION / 2)

		# A lot of intensities of zero with a peak in between
		if determine_noise(spectrum_XIC[left_index:right_index]) == 0 && max_intensity > 0
			left_index = findlast(x -> x == 0, spectrum_XIC[left_index:max_index]) + left_index - 1
			right_index = findfirst(x -> x == 0, spectrum_XIC[max_index:right_index]) + max_index - 1
			mass_integral[i, 2] = sum(spectrum_XIC[left_index:right_index])
			continue
		# No peak
		elseif max_intensity == 0
			mass_integral[i, 2] = 0
			continue
		end

		# Define peak between minima left and right of the maximum
		min_left = minimum(spectrum_XIC[left_index:max_index])
		min_right = minimum(spectrum_XIC[max_index:right_index])
		left_index = findlast(x -> x == min_left, spectrum_XIC[left_index:max_index]) + left_index - 1
		right_index = findfirst(x -> x == min_right, spectrum_XIC[max_index:right_index]) + max_index - 1

		# Check if peak is actual peak and not noise
		if determine_noise(spectrum_XIC[left_index:right_index]) >= 0
			mass_integral[i, 2] = 0
			continue
		end

		# Broader bounds for integration
		left_index = max_index - round(Int, MAX_PEAK_SCANS_INTEGRATION / 2)
		right_index = max_index + round(Int, MAX_PEAK_SCANS_INTEGRATION / 2)

		# Check if peak has overlap, if so, define peak end between overlapping peaks
		under_median_count = 0
		peak_end = Int32[-1, -1]
		for (j, index_range) in enumerate(((max_index, left_index), (max_index, right_index)))
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

		# Search for closest noise looking left and right of max
		noise_found = false
		noise_median = Int
		direction = -1
		noise_search = [i == -1 ? max_index : i for i in peak_end]
		while !noise_found # TODO stop searching when end of spectrum is reached
			j = direction == -1 ? 1 : 2
			noise_search[j] += NOISE_INTERVAL_SIZE * direction
			noise_interval_start = noise_search[j] + NOISE_INTERVAL_SIZE * direction
			spectrum_part = direction == -1 ? spectrum_XIC[noise_interval_start:noise_search[j]] : spectrum_XIC[noise_search[j]:noise_interval_start]

			if determine_noise(spectrum_part) >= 0
				noise_found = true
				noise_median = median(spectrum_part)
			end

			direction *= -1
		end

		# Substract noise median around peak max
		spectrum_part = spectrum_XIC[left_index:right_index]
		replace!(x -> x < noise_median ? 0 : x - noise_median, spectrum_part)

		# Find minimum left and right closest to max
		max_index = findmax(spectrum_part)[2]
		min_left = minimum(spectrum_part[1:max_index])
		min_right = minimum(spectrum_part[max_index:end])
		right_index = peak_end[2] == -1 ? findfirst(x -> x == min_right, spectrum_part[max_index:end]) + max_index - 1 : peak_end[2] - left_index + 1
		left_index = peak_end[1] == -1 ? findlast(x -> x == min_left, spectrum_part[1:max_index]) : peak_end[1] - left_index + 1

		# Integrate peak
		mass_integral[i, 2] = sum(spectrum_part[left_index:right_index])

	end

	return mass_integral
end


function determine_noise(spectrum_part)
	"""Returns mean of noise or -1 if peak detected"""
	

	# Large amount of intensities of zero
	if count(x -> x == 0, spectrum_part) / length(spectrum_part) > NOISE_ZEROS_FRACTION
		return 0
	end

	spectrum_part_mean = mean(spectrum_part)

	# Count number of mean crossings
	crossings_count = 0
	for i in 2:length(spectrum_part)
		crosses_mean = (spectrum_part[i] - spectrum_part_mean) * (spectrum_part[i - 1] - spectrum_part_mean) <= 0
		if crosses_mean
			crossings_count += 1
		end
	end

	# Small fraction of mean crossings, peak instead of noise
	if crossings_count / length(spectrum_part) < NOISE_CROSSINGS_FRACTION
		return -1
	else
		return round(Int, spectrum_part_mean)
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
