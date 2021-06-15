
"""
    compound_percentage_occurrence(compound_name, impurity_profile_csv, date_range, minimum_value)
Searches through impurity profile and determines percentage occurrence of compound
Can also search for a combination of compounds if compound_name and minimum_value are given as vector
date_range format: [yyyy-mm-dd, yyyy-mm-dd]
"""
function compound_percentage_occurrence(compound_name::String, impurity_profile_csv, date_range, minimum_value::Number)
    # Read impurity profile and filter on dates
    impurity_profile = CSV.read(impurity_profile_csv, DataFrame)
    date_range = sort(Date.(date_range))
    filter!(row -> Date(row."DateTime (y-m-d)"[1:10]) >= date_range[1] && Date(row."DateTime (y-m-d)"[1:10]) <= date_range[2] , impurity_profile)

    # Count occurrences
    occurrence_count = count(x -> x >= minimum_value, impurity_profile[!, compound_name])

    # Determine percent occurrence
    total_samples = nrow(impurity_profile)
    percent_occurrence = round(occurrence_count/nrow(impurity_profile) * 100, digits=2)
    println("Percent occurrence: $percent_occurrence, Total occurrence: $occurrence_count, Total number of samples: $total_samples")
    return percent_occurrence, occurrence_count, total_samples
end

# Searches for combination of compounds
function compound_percentage_occurrence(compound_name::Vector{String}, impurity_profile_csv, date_range, minimum_value::Vector{Number})
    if length(compound_name) != length(minimum_value)
        @error "compound_name and minimum_value not same length"
    end

    # Read impurity profile and filter on dates
    impurity_profile = CSV.read(impurity_profile_csv, DataFrame)
    date_range = sort(Date.(date_range))
    filter!(row -> Date(row."DateTime (y-m-d)"[1:10]) >= date_range[1] && Date(row."DateTime (y-m-d)"[1:10]) <= date_range[2] , impurity_profile)
    
    # Count occurrences
    occurrence_count = 0
    for row in eachrow(impurity_profile)
        meets_req = true
        for (compound, min_val) in zip(compound_name, minimum_value)
            if row[compound] < min_val
                meets_req = false
                break
            end
        end
        if meets_req
            occurrence_count += 1
        end
    end

    # Determine percent occurrence
    total_samples = nrow(impurity_profile)
    percent_occurrence = round(occurrence_count/nrow(impurity_profile) * 100, digits=2)
    println("Percent occurrence: $percent_occurrence, Total occurrence: $occurrence_count, Total number of samples: $total_samples")
    return percent_occurrence, occurrence_count, total_samples
end

# --------------------------------------------------------------------

"""
    find_similar_samples(sample, impurity_profile_csv, compound_group, max_deviation_fraction, date_range=0)
Looks for similar samples in impurity profile, samples are determined as similar if for each compound: 
{difference values / mean values <= max_deviation_fraction}

sample can be input as DataFrameRow or as index (Int) within the impurity profile

"""
function find_similar_samples(sample::DataFrameRow, impurity_profile_csv, compound_group, max_deviation_fraction, date_range=0)

    # Read impurity profile and filter dates if given
    impurity_profile = CSV.read(impurity_profile_csv, DataFrame)
    if date_range != 0
        date_range = sort(Date.(date_range))
        filter!(row -> Date(row."DateTime (y-m-d)"[1:10]) >= date_range[1] && Date(row."DateTime (y-m-d)"[1:10]) <= date_range[2] , impurity_profile)
    end

    # Read settings.json and retrieve compound names from group
    json_string = read(joinpath(@__DIR__, "settings.json"), String)
	settings_json = JSON3.read(json_string)
    compound_names = settings_json[:compound_groups][compound_group]

    # filter columns of impurity profile and sample based on compound group
    impurity_profile_filtered = select(impurity_profile, compound_names)
    select!(DataFrame(sample), compound_names)[1, :]

    # Determine similar samples
    similar_samples = DataFrame()
    for (i, row) in enumerate(eachrow(impurity_profile_filtered))
        meet_req = true
        # Compare values of both samples for all compounds
        for compound in compound_names
            samples_mean = mean((row[compound], sample[compound]))
            sample_deviation = abs(samples_mean - sample[compound]) / samples_mean
            if sample_deviation > max_deviation_fraction
                meet_req = false
                break
            end
        end
        if meet_req
            push!(similar_samples, impurity_profile[i, :])
        end
    end

    return similar_samples
end

# Sample input as index within impurity profile
function find_similar_samples(sample::Int, impurity_profile_csv, compound_group, max_deviation_fraction, date_range=0)

    # Read impurity profile and filter dates if given
    impurity_profile = CSV.read(impurity_profile_csv, DataFrame)
    if date_range != 0
        date_range = sort(Date.(date_range))
        filter!(row -> Date(row."DateTime (y-m-d)"[1:10]) >= date_range[1] && Date(row."DateTime (y-m-d)"[1:10]) <= date_range[2] , impurity_profile)
    end

    # Retrieve sample as DataFrameRow
    sample = impurity_profile[sample, :]

    # Read settings.json and retrieve compound names from group
    json_string = read(joinpath(@__DIR__, "settings.json"), String)
	settings_json = JSON3.read(json_string)
    compound_names = settings_json[:compound_groups][compound_group]

    # filter columns of impurity profile and sample based on compound group
    impurity_profile_filtered = select(impurity_profile, compound_names)
    select!(DataFrame(sample), compound_names)[1, :]

    # Determine similar samples
    similar_samples = DataFrame()
    for (i, row) in enumerate(eachrow(impurity_profile_filtered))
        meet_req = true
        # Compare values of both samples for all compounds
        for compound in compound_names
            samples_mean = mean((row[compound], sample[compound]))
            sample_deviation = abs(samples_mean - sample[compound]) / samples_mean
            if sample_deviation > max_deviation_fraction
                meet_req = false
                break
            end
        end
        if meet_req
            push!(similar_samples, impurity_profile[i, :])
        end
    end

    return similar_samples
end