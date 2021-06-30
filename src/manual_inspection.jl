"""
    mz_integrals_to_csv(spectra, compound_name, compounds)  
Determine mz integrals for specific compound and output to csv
"""
function mz_integrals_to_csv(pathin, compound_name)
    json_string = read(joinpath(SETTINGS_JSON_LOCATION, "settings.json"), String)
	settings_json = JSON3.read(json_string)
	main_compound_name = settings_json[:main_settings]["main_compound"]
    IS_name = settings_json[:main_settings]["internal_standard"]

    spectra, metadata_headers = batch_import(pathin, settings_json)

    # Import RT and mz info of valid compounds into DataFrame
	compounds_csv = CSV.read(joinpath(COMPOUNDS_CSV_LOCATION, "compounds.csv"), DataFrame)
    compound_to_analyse = filter(row -> row.compound == compound_name, compounds_csv) # DataFrame
    main_compound = filter(row -> row.compound == main_compound_name, compounds_csv)[1, :] # DataFrameRow
    internal_standard = filter(row -> row.compound == IS_name, compounds_csv)[1, :] # DataFrameRow

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
        mz_integrals = analyse_spectrum(spectrum, compounds_csv, main_compound, internal_standard)
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


"""
    manual_integration(spectrum, RT_range, mz_values)
Manually integrate spectrum between RT range at chosen mz values
"""
function manual_integration(spectrum, RT_range, mz_values)
    index_left, index_right = RT_to_scans(spectrum, RT_range)
    spectrum_XIC = filter_XIC(spectrum, mz_values)

    return sum(spectrum_XIC[index_left : index_right])
end


"""
    visualize_peak_range(spectrum, compound_name, mz; secondary_ylims=0, predicted_RT=0, hidetitle=false, transparent=false)
Shows plot of peak with determined bounds and baseline. Will plot secondary zoomed in plot if peak has high intensity.
Set secondary_ylims to enforce a specific y max for secondary plot.
Set predicted_RT to enforce specific retention time, instead of using the retention time from compounds.csv
Setting hidetitle and transparent to true will not show a title and set the background to be transparent, respectively.
"""
function visualize_peak_range(spectrum, compound_name, mz; secondary_ylims=0, predicted_RT=0, hidetitle=false, transparent=false)
    PLOT_EXTRA_SCANS = 60
    H_V_LINEWIDTH = 3

    if transparent
        default(background_color=:transparent, foreground_color=:black)
    end

    json_string = read(joinpath(SETTINGS_JSON_LOCATION, "settings.json"), String)
	settings_json = JSON3.read(json_string)
	main_compound_name = settings_json[:main_settings]["main_compound"]
    IS_name = settings_json[:main_settings]["internal_standard"]

    # Read compounds.csv
	compounds_csv = CSV.read(joinpath(COMPOUNDS_CSV_LOCATION, "compounds.csv"), DataFrame)
	filter!(row -> !(any(ismissing, (row.RT, row.mz)) || any((row.RT, row.mz) .== 0)), compounds_csv)

    # Retrieve row of compound and main compound from compounds.csv
    compound = filter(row -> row.compound == compound_name, compounds_csv)[1, :]
    main_compound = filter(row -> row.compound == main_compound_name, compounds_csv)[1, :]
    internal_standard = filter(row -> row.compound == IS_name, compounds_csv)[1, :]

    RT_modifier = determine_RT_modifier(spectrum, main_compound, internal_standard)

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
        _plot = plot(spectrum_XIC, xlims=_x_lims, ylims=(0,
                maximum(spectrum_XIC[_x_lims[1]:_x_lims[2]]) * 1.1), 
                        xlabel="scan number", ylabel="intensity", label=spectrum)
        default(background_color=:white, background_color_outside=:white)
        return _plot
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

    _plot = plot(_xlims[1]:_xlims[2], spectrum_XIC[_xlims[1]:_xlims[2]], xlims=_xlims, ylims=(0, maximum(spectrum_XIC[left_bound : right_bound]) * 1.1), 
               xlabel="scan number", ylabel="intensity", label="spectrum", legend=:topleft)
    if !hidetitle
        try
            plot!(title="Filename: $(spectrum["Filename"]) \nItem Number: $(spectrum["item_number"]) \nmz: $mz", titlefontsize=7)
        catch KeyError
            @warn "No metadata for title"
        end
    end

    vline!([left_bound, right_bound], linestyle=:dash, label="peak cutoffs", linewidth=H_V_LINEWIDTH)
    plot!(pre_bounds, baseline, linestyle=:dash, label="baseline ($(round(Int, mean(baseline))))", seriescolor=:green, linewidth=H_V_LINEWIDTH)

    if maximum(spectrum_XIC[left_bound : right_bound]) / secondary_ylims > 2.5
        
        plot!(twinx(), pre_bounds, baseline, xlims=_xlims, linestyle=:dot, ylims=(0, secondary_ylims), label="", 
              seriescolor=:green, xaxis=:false, linewidth=1, right_margin=35mm)

        plot!(twinx(), _xlims[1]:_xlims[2], spectrum_XIC[_xlims[1]:_xlims[2]], xlims=_xlims,
        ylims=(0, secondary_ylims), linestyle=:dot, ylabel="intensity (zoomed)", label="", grid=:off,
        linewidth=1)
        # plot!(linestyle=:dot)
    end
    default(background_color=:white, background_color_outside=:white)
    return _plot
end


"""
    visualize_peak_determination_steps(spectrum, compound_name, mz; pathout="", hidetitle=false, transparent=false)
Visualizes different steps in determining bounds and baseline with plots.
Setting a pathout will save the plots to this path.
Setting hidetitle and transparent to true will not show a title and set the background to be transparent, respectively.
"""
function visualize_peak_determination_steps(spectrum, compound_name, mz; transparent=false, pathout="", hidetitle=false)
    if transparent
        default(background_color=:transparent, foreground_color=:black)
    end
    # default(linewidth=4, labelfontsize=20, tickfontsize=18)
    H_V_LINEWIDTH = 3
    plots = []

	# Read settings.json
	json_string = read(joinpath(SETTINGS_JSON_LOCATION, "settings.json"), String)
	settings_json = JSON3.read(json_string)
	main_compound_name = settings_json[:main_settings]["main_compound"]
	IS_name = settings_json[:main_settings]["internal_standard"]

	# Import RT and mz info of valid compounds into DataFrame
	compounds_csv = CSV.read(joinpath(COMPOUNDS_CSV_LOCATION, "compounds.csv"), DataFrame)
	filter!(row -> !(any(ismissing, (row.RT, row.mz)) || any((row.RT, row.mz) .== 0)), compounds_csv)
	main_compound = filter(row -> row.compound == main_compound_name, compounds_csv)[1, :]
	internal_standard = filter(row -> row.compound == IS_name, compounds_csv)[1, :]
    compound_row = filter(row -> row.compound == compound_name, compounds_csv)[1, :]
    
    RT_modifier = determine_RT_modifier(spectrum, main_compound, internal_standard)
    @info "" RT_modifier
    if RT_modifier == -1
        @error "Main compound and/or internal standard not detected"
    end

    spectrum_XIC = filter_XIC(spectrum, mz)
    RT = compound_row.RT * RT_modifier
    overlap_RT_vals = parse_data(compound_row)[2]
    overlap_RT_vals .*= RT_modifier

    # --- Function search_peak --- #
    RT_range = (RT - MAX_RT_SHIFT, RT + MAX_RT_SHIFT)
	scan_range = collect(RT_to_scans(spectrum, RT_range))

    push!(plots, plot(1:length(spectrum_XIC), spectrum_XIC, legend=false, xlabel="Scan Number", ylabel="Intensity", title="First"))
    push!(plots, plot(deepcopy(plots[1]), title="Range to search max"))
    vline!(scan_range, linestyle=:dot, linewidth=H_V_LINEWIDTH)

    max_intensity, max_scan_number = findmax(spectrum_XIC[scan_range[1] : scan_range[2]])
	max_scan_number += scan_range[1] - 1
    if max_intensity <= 0
        @warn "$compound_name at m/z: $mz not detected"
        return plots
    end

    push!(plots, plot(deepcopy(plots[2]), title="Found max"))
    vline!([max_scan_number], linestyle=:dot, linewidth=H_V_LINEWIDTH, seriescolor=:black)

	# Bounds used for determining whether peak is noise or not
	left_scan = max_scan_number - round(Int, MAX_SCANS_PEAK_SEARCH / 2)
	right_scan = max_scan_number + round(Int, MAX_SCANS_PEAK_SEARCH / 2)
    push!(plots, plot(deepcopy(plots[1]), title="Range to determine noise"))
    vline!([left_scan, right_scan], linestyle=:dot, linewidth=H_V_LINEWIDTH)
    push!(plots, plot(deepcopy(plots[4]), title="Median"))
    hline!([median(spectrum_XIC[left_scan:right_scan])], linestyle=:dot, linewidth=H_V_LINEWIDTH)

    
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
        @warn "Noise detected instead of peak"
		return plots
    end
    # --- End of function --- #

    max_scan = max_scan_number
    pre_bounds = max_scan - MAX_SCANS_PEAK_LEFT : max_scan + MAX_SCANS_PEAK_RIGHT
    overlap_max_left, overlap_max_right = determine_overlap(spectrum_XIC, max_scan, overlap_RT_vals, spectrum)
	push!(plots, plot(deepcopy(plots[1]), title="Show overlap"))
    vline!([overlap_max_left, overlap_max_right], linestyle=:dot, linewidth=H_V_LINEWIDTH, seriescolor=:black)

    # --- Function determine_baseline --- #
    # Noise end as minimum closest to peak
	noise_end = max_scan - findmin(spectrum_XIC[max_scan:-1:max_scan - 30])[2] + 1
	spectrum_part = spectrum_XIC[noise_end - BASELINE_SCANS:noise_end]

	# Calculate median and set median as baseline
	spectrum_part_median = round(Int, median(spectrum_part))
	baseline = fill(spectrum_part_median, length(pre_bounds))

    push!(plots, plot(deepcopy(plots[1]), title="Show baseline"))
    vline!([noise_end - BASELINE_SCANS, noise_end], linestyle=:dot, linewidth=H_V_LINEWIDTH)
    plot!(pre_bounds, baseline, linestyle=:dot, linewidth=H_V_LINEWIDTH)
    # --- End of function --- #

    push!(plots, plot(deepcopy(plots[6]), title="Baseline with overlap values"))
    plot!(pre_bounds, baseline, linestyle=:dot, linewidth=H_V_LINEWIDTH)

    # --- Function determine_bounds --- #
    peak_bounds = determine_bounds(spectrum_XIC, max_scan, overlap_max_left, overlap_max_right, baseline, pre_bounds)
	pre_left_bound = overlap_max_left > 0 && overlap_max_left > pre_bounds[begin] ? overlap_max_left : pre_bounds[begin]
	pre_right_bound = overlap_max_right > 0 && overlap_max_right < pre_bounds[end] ? overlap_max_right : pre_bounds[end]
    # --- End of function --- #

    push!(plots, plot(deepcopy(plots[1]), title="Determining bounds"))
    vline!([peak_bounds[begin], peak_bounds[end]], linestyle=:dashdot, linewidth=H_V_LINEWIDTH)
    plot!([pre_bounds[begin]:pre_bounds[end]], baseline, linestyle=:dot, linewidth=H_V_LINEWIDTH)
    
    # Fill integral
    fill_bounds = peak_bounds[begin]:peak_bounds[end]
    baseline_bounds = (peak_bounds[1] - pre_bounds[begin] + 1):((peak_bounds[1] - pre_bounds[begin] + 1) + length(fill_bounds) - 1)
    plot_ = plot(deepcopy(plots[9]))
    plot!(plot_, fill_bounds, spectrum_XIC[fill_bounds], fillrange=baseline[baseline_bounds], fillalpha = 0.5, fillcolor = :grey, linewidth=0)
    vline!([pre_left_bound, pre_right_bound], seriescolor=:grey, linewidth=H_V_LINEWIDTH, linestyle=:dot)
    push!(plots, plot_)
    
    for (i, plot) in enumerate(plots)
        xlims!(plot, max_scan_number - (MAX_SCANS_PEAK_LEFT + 10), max_scan_number + (MAX_SCANS_PEAK_RIGHT + 10))
        ylims!(plot, 0, max_intensity*1.1)
        plot!(right_margin=15mm)
        if hidetitle
            plot!(title="")
        end
        if pathout != ""
            savefig(plot, joinpath(pathout, "plot_$i"))
        end
    end
    default(background_color=:white, background_color_outside=:white)
    return plots
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