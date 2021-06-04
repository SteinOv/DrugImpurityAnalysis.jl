using MS_Import
using Plots
using Printf
using CSV
using DataFrames
using LightXML

using BenchmarkTools

include("helpers.jl")
include("manual_inspection.jl")

const MAJOR_RT_DEVIATION = 1 # Absolute maximum shift of cocaine peak

const MAX_MZ_DEVIATION = 0.5 # Maximum deviation from given mz value in XIC spectrum
const MAX_RT_SHIFT = 0.05 # Maximum amount that RT can shift from determined RT

# 40 scans is 0.1 minute
const MAX_SCANS_PEAK_SEARCH = 20 # Used for differentiating peak from noise
const MAX_SCANS_PEAK_LEFT = 40 # Bounds for peak integration (left)
const MAX_SCANS_PEAK_RIGHT = 80 # Bounds for peak integration (right)

# If fraction of median/mean crossings is above either value, noise instead of peak
const NOISE_MEDIAN_CROSSINGS_FRACTION = 0.40
const NOISE_MEAN_CROSSINGS_FRACTION = 0.25

# Amount of scans for determining baseline
const BASELINE_SCANS = 20

# Used for determining overlap
const OVERLAP_CONSECUTIVE_BELOW_MEDIAN = 2
const OVERLAP_CONSECUTIVE_ABOVE_MEDIAN = 5



function main()

	# Specify path
	data_folder = joinpath(@__DIR__, "data")
	subfolder = "cocabricks"
	pathin = joinpath(data_folder, subfolder)
	csvout = joinpath(pathin, "impurity_profile.csv")

	# Import spectra
	spectra = batch_import(pathin)


	# Import RT and mz info of valid compounds into DataFrame
	compounds_csv = CSV.read("compounds.csv", DataFrame)
	filter!(row -> !(any(ismissing, (row.RT, row.mz)) || any((row.RT, row.mz) .== 0)), compounds_csv)

	# Create DataFrame for storing impurity profile (output)
	compounds_in_profile = filter(row -> !(row.type_int == -10), compounds_csv).compound
	impurity_profiles = DataFrame()
	insertcols!(impurity_profiles, :item_number => String[])
	insertcols!(impurity_profiles, :filename => String[])
	insertcols!(impurity_profiles, :folder => String[])
	for compound_name in compounds_in_profile
		insertcols!(impurity_profiles, Symbol(compound_name) => Float32[])
	end

	# Analyse all spectra
	for i=1:length(spectra)
		println("Analysing spectrum $i...")
		spectrum = spectra[i]["MS1"]

		compound_mz_integrals = analyse_spectrum(spectrum, compounds_csv)
		# Create impurity profile and add to DataFrame
		sample_metadata, sample_profile = create_impurity_profile(spectrum, compound_mz_integrals, compounds_csv, compounds_in_profile)
		push!(impurity_profiles, append!(sample_metadata, sample_profile))
	end

	CSV.write(csvout, impurity_profiles)
end

"""
	analyse_spectrum(spectrum, compounds)
Analyses all compounds in spectrum and returns dictionary with mz values for each compound
"""
function analyse_spectrum(spectrum, compounds_csv, compounds_to_analyse=compounds_csv)
	RT_modifier = determine_RT_modifier(spectrum, compounds_csv)

	# Create Dict for storing integral for each compound and mz value
	compound_mz_integrals = Dict(Symbol(compound) => Dict{Symbol, Int}() for compound in compounds_to_analyse.compound)

	# No cocaine peak found
	if RT_modifier == -1
		return 0
	end

	# Iterate over compounds in compounds.csv
	for row in eachrow(compounds_to_analyse)
		RT = row.RT * RT_modifier
		mz_vals, overlap_RT_vals = parse_data(row)
		overlap_RT_vals .*= RT_modifier

		# Iterate over each mz value per compound
		for mz in mz_vals
			spectrum_XIC = filter_XIC(spectrum, mz)
			peak_exists, max_scan = search_peak(spectrum_XIC, RT, spectrum)

			if peak_exists == false
				compound_mz_integrals[Symbol(row.compound)][Symbol(mz)] = 0
				continue
			end
			
			# Determine baseline and bounds
			pre_bounds = max_scan - MAX_SCANS_PEAK_LEFT : max_scan + MAX_SCANS_PEAK_RIGHT
			overlap_max_left, overlap_max_right = determine_overlap(spectrum_XIC, max_scan, overlap_RT_vals, spectrum)
			baseline = determine_baseline(spectrum_XIC, max_scan, overlap_max_left, overlap_max_right, pre_bounds)
			peak_bounds = determine_bounds(spectrum_XIC, max_scan, overlap_max_left, overlap_max_right, baseline, pre_bounds)

			# Calculate and store integral
			peak_integral = integrate_peak(spectrum_XIC, peak_bounds, baseline, pre_bounds)
			compound_mz_integrals[Symbol(row.compound)][Symbol(mz)] = peak_integral
		end
	end

	return compound_mz_integrals
end #FUNCTION


"""
	determine_RT_modifier(spectrum, compounds)
Determines retention time shift of cocaine and returns {RT_actual / RT_predicted}

	return RT_modifier"""
function determine_RT_modifier(spectrum, compounds)

	# Read predicted RT and highest mz value
	major_compound = filter(row -> row.type_int == -1, compounds)[1, :]
	RT_predicted = major_compound.RT
	highest_mz_val = maximum(parse_data(major_compound)[1])

	# Search for maximum around predicted RT
	RT_range = (RT_predicted - MAJOR_RT_DEVIATION, RT_predicted + MAJOR_RT_DEVIATION)
	RT_index_range = RT_to_scans(spectrum, RT_range)
	spectrum_XIC_part = filter_XIC(spectrum, highest_mz_val)[RT_index_range[1] : RT_index_range[2]]

	# No peak at all
	if maximum(spectrum_XIC_part) <= 0
		return -1
	end

	max_scan = findmax(spectrum_XIC_part)[2] + RT_index_range[1] - 1
	RT_actual = spectrum["Rt"][max_scan]

	return round(RT_actual / RT_predicted, digits=3)
end

"""
	search_peak(spectrum_XIC, RT, spectrum)
Searches for peak in ion-extracted chromatogram
	
	return peak_exists, max_scan_number"""
function search_peak(spectrum_XIC, RT, spectrum, left_scan=0, right_scan=0)
	# Range to search peak in

	RT_range = (RT - MAX_RT_SHIFT, RT + MAX_RT_SHIFT)
	scan_range = collect(RT_to_scans(spectrum, RT_range))

	# If scans pre-defined
	scan_range[1] = left_scan == 0 ? scan_range[1] : left_scan
	scan_range[2] = right_scan == 0 ? scan_range[2] : right_scan

	max_intensity, max_scan_number = findmax(spectrum_XIC[scan_range[1] : scan_range[2]])
	max_scan_number += scan_range[1] - 1

	# No peak at all
	if max_intensity <= 0
		return false, -1
	end

	# Bounds used for determining whether peak is noise or not
	left_scan = left_scan == 0 ? max_scan_number - round(Int, MAX_SCANS_PEAK_SEARCH / 2) : left_scan
	right_scan = right_scan == 0 ? max_scan_number + round(Int, MAX_SCANS_PEAK_SEARCH / 2) : right_scan

	# Count number of median crosses
	spectrum_part_median = median(spectrum_XIC[left_scan:right_scan])
	median_crossings_count = 0
	for i in (left_scan + 1):right_scan
		if (spectrum_XIC[i] - spectrum_part_median) * (spectrum_XIC[i - 1] - spectrum_part_median) <= 0
			median_crossings_count += 1
		end
	end

	# Count number of mean crosses
	spectrum_part_mean = mean(spectrum_XIC[left_scan:right_scan])
	mean_crossings_count = 0
	for i in (left_scan + 1):right_scan
		if (spectrum_XIC[i] - spectrum_part_mean) * (spectrum_XIC[i - 1] - spectrum_part_mean) <= 0
			mean_crossings_count += 1
		end
	end

	# More crossings than allowed, noise instead of peak
	if median_crossings_count / length(spectrum_XIC[left_scan:right_scan]) > NOISE_MEDIAN_CROSSINGS_FRACTION ||
	   		mean_crossings_count / length(spectrum_XIC[left_scan:right_scan]) > NOISE_MEAN_CROSSINGS_FRACTION
		return false, -1
	else
		return true, max_scan_number
	end
end


"""
Determines scan numbers of maxima of overlapping peaks
For side at which no overlap is defined, takes a simple and 
non-thorough approach for determining whether overlap is present
"""
function determine_overlap(spectrum_XIC, max_scan, RT_overlap_vals, spectrum)
	overlap_max_scans = zeros(Int, 2)
	for (i, RT_overlap) in enumerate(RT_overlap_vals)

		left_scan, right_scan = 0, 0
		# Overlap not defined
		if RT_overlap == 0
			# Determine iteration sequence
			direction  = i == 1 ? -1 : 1
			scan_upper_bound = i == 1 ? -MAX_SCANS_PEAK_LEFT : MAX_SCANS_PEAK_RIGHT
			spectrum_scan_range = max_scan:direction:(scan_upper_bound + max_scan)

			# Determine mean of range and initialize variables
			range_mean = mean(spectrum_XIC[spectrum_scan_range])
			below_consecutive, above_consecutive = 0, 0
			reached_below = false
			overlap_max_scan = 0

			# Iterate through sequence and search for peak separate from main peak
			for scan=spectrum_scan_range
				if spectrum_XIC[scan] < range_mean
					below_consecutive += 1
					above_consecutive = 0
				else
					below_consecutive = 0
					above_consecutive += 1
				end

				# Consecutively below median
				if below_consecutive >= OVERLAP_CONSECUTIVE_BELOW_MEDIAN
					reached_below = true
				# Previously consecutively below median and now consecutively above median
				elseif reached_below && above_consecutive >= OVERLAP_CONSECUTIVE_ABOVE_MEDIAN
					overlap_range = (scan - OVERLAP_CONSECUTIVE_ABOVE_MEDIAN * direction):direction:spectrum_scan_range[end]
					# overlap_max_scan = overlap_range[begin] + direction * (findmax(spectrum_XIC[overlap_range])[2] - 1)#TEMP
					# Bounds for search_peak
					left_scan = direction == -1 ? overlap_range[end] : overlap_range[begin]
					right_scan = direction == -1 ? overlap_range[begin] : overlap_range[end]
					break
				end
			end

			# Set RT_overlap to 0 if not found, else to RT of maximum
			RT_overlap = overlap_max_scan == 0 ? 0 : spectrum["Rt"][overlap_max_scan]
		else
			overlap_scan = RT_to_scans(spectrum, RT_overlap)[1]
			scan_between = round(Int, (overlap_scan + max_scan) / 2)
			left_scan = i == 1 ? 0 : scan_between
			right_scan = i == 1 ? scan_between : 0
		end

		if RT_overlap > 0
			# Final check whether peak exists
			overlap_max_scans[i] = search_peak(spectrum_XIC, RT_overlap, spectrum, left_scan, right_scan)[2]
		end
	end

	return overlap_max_scans
end


"""
Determines baseline based on range left of peak
"""
function determine_baseline(spectrum_XIC, max_scan, overlap_max_left, overlap_max_right, bounds)
	#TODO take section on right side if overlap on left
	if overlap_max_left > 0
		@warn "Overlap on left side" overlap_max_left, overlap_max_right
	end
	# Noise end as minimum closest to peak
	noise_end = max_scan - findmin(spectrum_XIC[max_scan:-1:max_scan - 30])[2] + 1
	spectrum_part = spectrum_XIC[noise_end - BASELINE_SCANS:noise_end]

	# Calculate median and set median as baseline
	spectrum_part_median = round(Int, median(spectrum_part))
	baseline = fill(spectrum_part_median, length(bounds))

	return baseline
end



#=
"""
Determines baseline within bounds
First sets all values above half of maximum to zero,
then sets baseline to median of remaining spectrum part
"""
function determine_baseline(spectrum_XIC, max_scan, bounds)
	# Set all values above half of the maximum to zero
	spectrum_part = spectrum_XIC[round(Int, max_scan - BASELINE_SCANS / 2) : round( Int, max_scan + BASELINE_SCANS / 2)]
	spectrum_part_cutoff = maximum(spectrum_part) / 2
	spectrum_part = map(x -> x > spectrum_part_cutoff ? 0 : x, spectrum_part)

	# Calculate median and set median as baseline
	spectrum_part_median = median(spectrum_part)
	baseline = fill(spectrum_part_median, length(bounds))

	return baseline
end
=#

"""
Determine bounds of peak
Set bounds to baseline crossing closest to peak, or minimum if it does not cross
"""
function determine_bounds(spectrum_XIC, max_scan, overlap_max_left, overlap_max_right, baseline, pre_bounds)
	# Set bound to max of overlapping peak if it exists
	pre_left_bound = overlap_max_left > 0 ? overlap_max_left : pre_bounds[begin]
	pre_right_bound = overlap_max_right > 0 ? overlap_max_right : pre_bounds[end]

	spectrum_part = spectrum_XIC[pre_left_bound:pre_right_bound]
	part_max_scan = max_scan - pre_left_bound + 1

	# Subtract baseline from spectrum, negative values set to zero
	spectrum_part = map((x_spectr, x_base) -> (x_spectr < x_base ? 0 : x_spectr - x_base), spectrum_part, baseline)

	# Set bounds to zero closest to maximum or minimum if no zero value
	left_bound  = max_scan - findmin(spectrum_part[part_max_scan:-1:begin])[2] + 1
	right_bound = max_scan + findmin(spectrum_part[part_max_scan:end])[2] - 1

	return left_bound:right_bound
end

"""Integrates peak within bounds with baseline subtracted"""
function integrate_peak(spectrum_XIC, bounds, baseline, pre_bounds)
	spectrum_part = spectrum_XIC[pre_bounds]
	# spectrum_part .-= baseline
	bounds = bounds .- (pre_bounds[begin] - 1)

	# Subtract baseline from spectrum, negative values set to zero
	spectrum_part = map((x_spectr, x_base) -> (x_spectr < x_base ? 0 : x_spectr - x_base), spectrum_part, baseline)

	# Calculate integral within bounds
	peak_integral = sum(spectrum_part[bounds])

	return peak_integral
end


"""
	create_impurity_profile(spectrum, compound_mz_integrals, impurity_profiles)
Creates impurity profile and adds it to DataFrame "impurity_profiles"
"""
function create_impurity_profile(spectrum, compound_mz_integrals, compounds_csv, compounds_in_profile)

	# Retrieve metadata
	metadata = Array{Any,1}(undef, 3)
	metadata[1] = spectrum["Sample Name"][1]
	metadata[2] = spectrum["Filename"][1]
	metadata[3] = spectrum["Folder Name"][1]

	# No cocaine present
	if compound_mz_integrals == 0
		return metadata, zeros(length(compounds_in_profile))
	end

	# Retrieve names of major compound and internal standard, and retrieve minimum ratio
	major_IS_ratio = parse(Float16, filter(row -> row.type_int == -10, compounds_csv).ratio[1])
	major_compound = filter(row -> row.type_int == -1, compounds_csv)[1, :]
	IS_name = filter(row -> row.type_int == -10, compounds_csv).compound[1]

	# Determine integral of highest mz value of major compound
	major_highest_mz_value = maximum(parse_data(major_compound)[1])
	major_highest_mz_integral = compound_mz_integrals[Symbol(major_compound.compound)][Symbol(major_highest_mz_value)]

	# Sum integrals of all mz values of internal standard
	IS_integral = sum(values(compound_mz_integrals[Symbol(IS_name)]))


	impurity_profile_values = zeros(length(compounds_in_profile))
	# Ratio too low, no significant amount of major compound present 
	if major_highest_mz_integral / IS_integral < major_IS_ratio || major_highest_mz_integral == 0
		return metadata, zeros(length(compounds_in_profile))
	elseif IS_integral == 0 && major_highest_mz_integral != 0
		@warn "Major compound detected, but internal standard not detected"
		return metadata, zeros(length(compounds_in_profile))
	end

	# Sum all integrals of cocaine
	major_integral_sum = sum(values(compound_mz_integrals[Symbol(major_compound.compound)]))
	
	# Calculate percent ratio for all compounds
	for (i, compound_name) in enumerate(compounds_in_profile)
		mz_dict = compound_mz_integrals[Symbol(compound_name)]

		# Check if first mass has no intensity
		first_mz = parse_data(filter(row -> row.compound == compound_name, compounds_csv)[1, :])[1][1]
		if mz_dict[Symbol(first_mz)] <= 0
			impurity_profile_values[i] = 0
			continue
		end

		# Calculate percent ratio, set to zero if < 0.05
		integral_sum = sum(values(mz_dict))
		percent_ratio = round(integral_sum / major_integral_sum * 100, digits=2)
		# percent_ratio = percent_ratio > 0.05 ? percent_ratio : 0 #TEMP
		impurity_profile_values[i] = percent_ratio
	end

	return metadata, impurity_profile_values
end

# main()
