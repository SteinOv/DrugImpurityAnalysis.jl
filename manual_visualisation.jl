"""
Functions for visualising spectra manually.

plt: Plots spectrum for chosen mz values and retention time
plt!: Same as plt, but plots on top of previous plot
mass_ratios: Prints mz values of chosen compounds and their (normalised) ratios

"""

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
		add_left = 100
		add_right = 100
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
		add_left = 100
		add_right = 100
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
            mass_intensity = integrate_peaks(spectrum, RT * RT_modifier, mass_vals)
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