"""
Functions for manual inspection of spectra.

plt: Plots spectrum for chosen mz values and retention time
plt!: Same as plt, but plots on top of previous plot
mass_ratios: Prints mz values of chosen compounds and their (normalised) ratios
display_distribution: Integrates all mz values of specific compound in spectrum and outputs to .csv file
manual_integration: Manually integrates spectrum between chosen RT and mz values

"""

# include("main.jl")

"""Plot spectrum of given mass
plt((mass), (RT), (RT_deviation))

mass: XIC at specific mass
mass = 0: Spectrum not filtered

RT given: Zooms in around given RT
RT = 0: Shows complete RT spectrum
RT_range_size: Deviation around RT (default 0.5)
"""
function plt(mass=0, RT=0, RT_range_size=0.5)


	# TODO gives error at start and end

	RT_range = [RT - RT_range_size, RT + RT_range_size]
	RT_range_index = RT_to_scans(spectrum, RT_range)

	if mass != 0
		spectrum_XIC = filter_XIC(spectrum, mass)
	else
		spectrum_XIC = filter_XIC(spectrum, 0)
	end

	if RT > 0
		plot(spectrum["Rt"][RT_range_index[1]:RT_range_index[2]], spectrum_XIC[(RT_range_index[1]):RT_range_index[2]])
	else
		plot(spectrum["Rt"], spectrum_XIC)
	end
end


"""Plot spectrum of given mass over previous spectrum
plt((mass), (RT), (RT_deviation))

mass: XIC at specific mass
mass = 0: Spectrum not filtered

RT given: Zooms in around given RT
RT = 0: Shows complete RT spectrum
RT_range_size: Deviation around RT (default 0.5)
"""
function plt!(mass=0, RT=0, RT_range_size=0.5)
	# TODO gives error at start and end

	RT_range = [RT - RT_range_size, RT + RT_range_size]
	RT_range_index = RT_to_scans(spectrum, RT_range)

	if mass != 0
		spectrum_XIC = filter_XIC(spectrum, mass)
	else
		spectrum_XIC = filter_XIC(spectrum, 0)
	end

	if RT > 0
		plot!(spectrum["Rt"][RT_range_index[1]:RT_range_index[2]], spectrum_XIC[(RT_range_index[1]):RT_range_index[2]])
	else
		plot!(spectrum["Rt"], spectrum_XIC)
	end
end

"""
    mz_integrals_to_csv(spectra, compound_name, compounds)  
Determine mz integrals for specific compound and output to csv
"""
function mz_integrals_to_csv(pathin, compound_name)
    json_string = read(joinpath(@__DIR__, "settings.json"), String)
	settings_json = JSON3.read(json_string)
	main_compound_name = settings_json[:main_settings]["main_compound"]

    spectra, metadata_headers = batch_import(pathin, settings_json)

    # Import RT and mz info of valid compounds into DataFrame
	compounds_csv = CSV.read("compounds.csv", DataFrame)
    compound_to_analyse = filter(row -> row.compound == compound_name, compounds_csv) # DataFrame
    main_compound = filter(row -> row.compound == main_compound_name, compounds_csv)[1, :] # DataFrameRow


    # Create DataFrame for storing mz integrals
    mz_values = parse_data(compound_to_analyse[1, :])[1]
    mz_integrals_dataframe = DataFrame()
	for header in metadata_headers
		insertcols!(mz_integrals_dataframe, Symbol(header) => String[])
	end
	for mz in mz_values
		insertcols!(mz_integrals_dataframe, Symbol(mz) => Float32[])
	end
	
    # Analyse all spectra
	for i=1:length(spectra)
		println("Analysing spectrum $i...")
        spectrum = spectra[i]["MS1"]

        # Determine mz integrals, retrieve metadata, and store in DataFrame
        mz_integrals = analyse_spectrum(spectrum, compounds_csv, main_compound, compound_to_analyse)
        sample_metadata = Vector{Any}([spectrum[header] for header in metadata_headers])
        if mz_integrals == 0
            push!(mz_integrals_dataframe, append!(sample_metadata, zeros(length(mz_values))))
        else
            mz_integrals = mz_integrals[Symbol(compound_name)]
            push!(mz_integrals_dataframe, append!(sample_metadata, [mz_integrals[Symbol(mz)] for mz in mz_values]))
        end
        
        
    end

    CSV.write(joinpath(@__DIR__, "$(compound_name)_mz integrals.csv"), mz_integrals_dataframe)
end



"""Manually integrate spectrum between RT range at chosen mz values"""
function manual_integration(spectrum, RT_range, mz_values)
    index_left, index_right = RT_to_scans(spectrum, RT_range)
    spectrum_XIC = filter_XIC(spectrum, mz_values)

    return sum(spectrum_XIC[index_left : index_right])
end

"""Visualize determined peak range"""
function visualize_peak_range(spectrum, compound_name, mz, secondary_ylims=0, predicted_RT=0)
    PLOT_EXTRA_SCANS = 60

    json_string = read(joinpath(@__DIR__, "settings.json"), String)
	settings_json = JSON3.read(json_string)
	main_compound_name = settings_json[:main_settings]["main_compound"]

    # Read compounds.csv
	compounds_csv = CSV.read("compounds.csv", DataFrame)
	filter!(row -> !(any(ismissing, (row.RT, row.mz)) || any((row.RT, row.mz) .== 0)), compounds_csv)

    # Retrieve row of compound and main compound from compounds.csv
    compound = filter(row -> row.compound == compound_name, compounds_csv)[1, :]
    main_compound = filter(row -> row.compound == main_compound_name, compounds_csv)[1, :]

    RT_modifier = determine_RT_modifier(spectrum, main_compound)

    # Cocaine not found
    if RT_modifier == -1
        return @error "Cocaine peak not found"
    end
    # Determine predicted RT
    if predicted_RT == 0
        predicted_RT = compound.RT * RT_modifier
    end    

    spectrum_XIC = filter_XIC(spectrum, mz)
    peak_exists, max_scan = search_peak(spectrum_XIC, predicted_RT, spectrum)


    if peak_exists == false
        _x_lims = [i for i in RT_to_scans(spectrum, predicted_RT)]
        _x_lims[1] -= PLOT_EXTRA_SCANS; _x_lims[2] += PLOT_EXTRA_SCANS
        return plot(spectrum_XIC, xlims=_x_lims, ylims=(0,
                maximum(spectrum_XIC[_x_lims[1]:_x_lims[2]]) * 1.1), 
                        xlabel="scan number", ylabel="intensity", label=compound_name)
    end
    
    # Determine baseline and bounds
    pre_bounds = max_scan - MAX_SCANS_PEAK_LEFT : max_scan + MAX_SCANS_PEAK_RIGHT
    overlap_RT_vals = parse_data(compound)[2] .* RT_modifier
    overlap_max_left, overlap_max_right = determine_overlap(spectrum_XIC, max_scan, overlap_RT_vals, spectrum)
    baseline = determine_baseline(spectrum_XIC, max_scan, overlap_max_left, overlap_max_right, pre_bounds)
    bounds = determine_bounds(spectrum_XIC, max_scan, overlap_max_left, overlap_max_right, baseline, pre_bounds)
    left_bound, right_bound = bounds[begin], bounds[end]

    
    if secondary_ylims == 0
        secondary_ylims = max(spectrum_XIC[left_bound] * 5, spectrum_XIC[right_bound] * 5, 4000)
    end
    _xlims = (left_bound - PLOT_EXTRA_SCANS, right_bound + PLOT_EXTRA_SCANS)

    plot(spectrum_XIC, xlims=_xlims, ylims=(0, maximum(spectrum_XIC[left_bound : right_bound]) * 1.1), 
               xlabel="scan number", ylabel="intensity", label=compound_name, left_margin = 5Plots.mm, 
                                                legend=:topleft,right_margin = 18Plots.mm, grid=:off)
    plot!(title="Filename: $(spectrum["Filename"][1]) \nItem Number: $(spectrum["item_number"][1]) \nmz: $mz", titlefontsize=7)
    vline!([left_bound, right_bound], linestyle=:dash, label="peak cutoffs")
    plot!(pre_bounds, baseline, linestyle=:dash, label="baseline ($(mean(baseline)))", seriescolor=:green)

    if maximum(spectrum_XIC[left_bound : right_bound]) / secondary_ylims > 2.5

        plot!(twinx(), pre_bounds, baseline, xlims=_xlims, linestyle=:dot, ylims=(0, secondary_ylims), label="", seriescolor=:green, xaxis=:false)

        plot!(twinx(), spectrum_XIC, xlims=_xlims,
        ylims=(0, secondary_ylims), linestyle=:dot, ylabel="intensity (zoomed)", label="", grid=:off)
        # plot!(linestyle=:dot)
    end
    #, seriescolor=RGB(160/255, 193/255, 249/255)
    return plot!()
end


"""
    function sort_spectra(spectra, sort_by)
Sorts spectra by chosen header
"""
function sort_spectra(spectra, sort_by)
    sort!(spectra, by = x -> x["MS1"][sort_by])
end

"""
    retrieve_spectrum(spectra, header, value)
Returns first spectrum found that has chosen value under chosen header
"""
function retrieve_spectrum(spectra, header, value)
    for spectrum in spectra
        if spectrum["MS1"][header] == value
            return spectrum["MS1"]
        end
    end
    
    println("Spectrum Not Found")
    return 0
end