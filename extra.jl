include("main.jl")
include("manual_visualisation.jl")
include("helpers.jl")



function display_distribution(compound)
    data_folder = joinpath(@__DIR__, "data")
    subfolder = "tmp"
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
        mass_integral = integrate_peaks(spectrum, RT, mass_values)
        push!(distribution, append!(sample_metadata, zeros(length(mass_values))))
        for j=1:length(mass_values)
            distribution[i, Symbol(string(mass_integral[j, 1]))] = mass_integral[j, 2]
        end

    end

    CSV.write(csvout, distribution)
end 