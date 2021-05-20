"""
Functions for manual inspection of spectra.

plt: Plots spectrum for chosen mz values and retention time
plt!: Same as plt, but plots on top of previous plot
mass_ratios: Prints mz values of chosen compounds and their (normalised) ratios
display_distribution: Integrates all mz values of specific compound in spectrum and outputs to .csv file
manual_integration: Manually integrates spectrum between chosen RT and mz values

"""

# include("main.jl")

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
		plot(spectrum["Rt"][RT_range_index[1]:RT_range_index[2]], spectrum_XIC[(RT_range_index[1]):RT_range_index[2]])
	else
		plot(spectrum["Rt"], spectrum_XIC)
	end
end

function plt!(mass=0, RT=0, RT_range_size=0.5)
	"""Plot spectrum of given mass over previous spectrum
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
		plot!(spectrum["Rt"][RT_range_index[1]:RT_range_index[2]], spectrum_XIC[(RT_range_index[1]):RT_range_index[2]])
	else
		plot!(spectrum["Rt"], spectrum_XIC)
	end
end




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
            mass_intensity = integrate_mz_values(spectrum, RT * RT_modifier, mass_vals)
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
    
        # println("Normalized list second mass:")
        # for i in normalized_list
        # 	@printf("%4i\n", i)
        # end
    
    end


	
function display_distribution(compound)
    data_folder = joinpath(@__DIR__, "data")
    subfolder = "coca_caf_cal"
    pathin = joinpath(data_folder, subfolder)
    csvout = joinpath(pathin, "distribution.csv")

    spectra = batch_import(pathin)

    compounds = CSV.read("compounds.csv", DataFrame)
    filter!(row -> !(any(ismissing, (row.RT, row.mz)) || any((row.RT, row.mz) .== 0)), compounds)

    compound_row = compounds[findfirst(comp -> comp == compound, compounds.compound), :]
    

    # Read mz values
    mass_str = split(compound_row.mz, ";")
    mass_values = Array{Any,1}(undef, length(mass_str))
    for (i, mass) in enumerate(mass_str)
        if startswith(mass, "(")
            mass = split(mass[2:end - 1], ",")
            mass_values[i] = Tuple(([parse(Float32, sub_mass) for sub_mass in mass]))
        else
            mass_values[i] = parse(Float32, mass)
        end	
    end

    distribution = DataFrame()
    insertcols!(distribution, :item_number => String[])
    insertcols!(distribution, :filename => String[])
    insertcols!(distribution, :folder => String[])
    for mz in mass_values
        insertcols!(distribution, Symbol(mz) => Int32[])
    end

    for i=1:length(spectra)
        println("Analysing spectrum $i...")
        spectrum = spectra[i]["MS1"]

        RT_modifier = process_major(compounds, spectrum)[2]
        RT = compound_row.RT * RT_modifier

        # Save metadata
        sample_metadata = Array{Any,1}(undef, 3)
        sample_metadata[1] = spectrum["Sample Name"][1]
        sample_metadata[2] = spectrum["Filename"][1]
        sample_metadata[3] = subfolder

        # Integrate peaks
        mass_integral = integrate_mz_values(spectrum, RT, mass_values)
        push!(distribution, append!(sample_metadata, zeros(length(mass_values))))
        for j=1:length(mass_values)
            distribution[i, Symbol(string(mass_integral[j, 1]))] = mass_integral[j, 2]
        end

    end

    CSV.write(csvout, distribution)
end 

function manual_integration(spectrum, RT_range, mz_values)
    """Manually integrate spectrum between RT range at chosen mz values"""
    index_left, index_right = RT_indices(spectrum, RT_range)
    spectrum_XIC = filter_XIC(spectrum, mz_values)

    return sum(spectrum_XIC[index_left : index_right])
end


function visualize_peak_range(spectrum, compound_name, mz, secondary_ylims=0, predicted_RT=0)
    """Visualize determined peak range"""
    PLOT_EXTRA_SCANS = 60

    # Read compounds.csv
	compounds = CSV.read("compounds.csv", DataFrame)
	filter!(row -> !(any(ismissing, (row.RT, row.mz)) || any((row.RT, row.mz) .== 0)), compounds)

    # Retrieve row of compound in compounds.csv
    compound = filter(row -> row.compound == compound_name, compounds)

    RT_modifier = process_major(compounds, spectrum)[2]

    # Determine predicted RT
    if predicted_RT == 0
        predicted_RT = compound.RT[1] * RT_modifier
    end

    # Read overlap_RT and determine RT (index) range
    overlap_RT = ismissing(compound.overlap[1]) ? 0 : compound.overlap[1] * RT_modifier # TODO support 2 overlap RT's or return error
    RT_range = (predicted_RT - MAX_RT_DEVIATION, predicted_RT + MAX_RT_DEVIATION)
	RT_range_index = RT_indices(spectrum, RT_range)

    # Determine peak range and baseline (noise median)
    spectrum_XIC = filter_XIC(spectrum, mz)
    left_index, right_index, noise_median = determine_peak_range(spectrum_XIC, RT_range_index, compound.RT[1], overlap_RT)

    if left_index == -1 || maximum(spectrum_XIC[left_index : right_index]) <= 0
        return plot(spectrum_XIC, xlims=(RT_range_index[1], RT_range_index[2]), ylims=(0,
                     maximum(spectrum_XIC[RT_range_index[1] : RT_range_index[2]]) * 1.1), 
                          xlabel="scan number", ylabel="intensity", label=compound_name)
    end
    if secondary_ylims == 0
        secondary_ylims = max(spectrum_XIC[left_index] * 5, spectrum_XIC[right_index] * 5, 4000)
    end
    _xlims = (left_index - PLOT_EXTRA_SCANS, right_index + PLOT_EXTRA_SCANS)

    plot(spectrum_XIC, xlims=_xlims, ylims=(0, maximum(spectrum_XIC[left_index : right_index]) * 1.1), 
               xlabel="scan number", ylabel="intensity", label=compound_name, left_margin = 5Plots.mm, 
                                                legend=:topleft,right_margin = 18Plots.mm, grid=:off)
    plot!(title="Filename: $(spectrum["Filename"][1]) \nSample Name: $(spectrum["Sample Name"][1]) \nmz: $mz", titlefontsize=7)
    vline!([left_index, right_index], linestyle=:dash, label="peak cutoffs")
    hline!([noise_median], linestyle=:dash, label="baseline ($noise_median)", seriescolor=:green)

    if maximum(spectrum_XIC[left_index : right_index]) / secondary_ylims > 2.5

        hline!(twinx(), [noise_median], xlims=_xlims, linestyle=:dot, ylims=(0, secondary_ylims), label="", seriescolor=:green, xaxis=:false)

        plot!(twinx(), spectrum_XIC, xlims=_xlims,
        ylims=(0, secondary_ylims), linestyle=:dot, ylabel="intensity (zoomed)", label="", grid=:off)
        # plot!(linestyle=:dot)
    end
    #, seriescolor=RGB(160/255, 193/255, 249/255)
    return plot!()
end

