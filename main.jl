using MS_Import
using Plots
using Printf

function main()

	# Specify path
	folder = "P:/Git/bachelor_project"
	mzxml_path = "data/mzxml"
	pathin = joinpath(folder, mzxml_path)

	# Import spectra
	spectra = batch_import(pathin)



	# spectrum = spectra[6]["MS1"]


	normalized_list = zeros(Float32, length(spectra))

	# Print intensities of mz peaks of cocaine for each spectrum
	for i=1:length(spectra)
		spectrum = spectra[i]["MS1"]
		println("Spectrum $i: \t Mass \t   Intensity     Normalised Intensity")

		mass_intensity = cocaine_intensity(spectrum)
		norm = mass_intensity[1, 2]
		
		for (mass, intensity) in zip(mass_intensity[:, 1], mass_intensity[:, 2])
			# println("$mass: $intensity")
			@printf("\t\t%3.2f\t| %10.3E  |  %1.3f\n", mass, intensity, intensity/norm)

			if mass == 182.12
				normalized_list[i] = intensity/norm
			end

		end

		println("----------------------------")
	end

	println("Normalized list:")
	for i in normalized_list
		@printf("%1.3f", i)
	end

end


function cocaine_intensity(spectrum)
	"""Returns intensity of cocaine peak if present, else 0"""

	# Cocaine approximate peak locations
	RT = 6.66 # minutes
	mass_vals = [82.07, 182.12, 83.07, 94.07, 77.04, 105.03, 96.08, 42.03, 303.15]
	
	return integrate_peak(spectrum, RT, mass_vals)

end


"""
TODO
- Improve overlapping peak integration by predicting actual integral after peak overlap
- Noise cut off determined dynamically

"""
function integrate_peak(spectrum, RT, mass_vals)
	"""Integrates peak of specific compound"""
	noise_cutoff = 5000 # TODO Maybe determined dynamically

	max_RT_deviation = 0.1 # minutes
	max_mass_deviation = 0.5 
	
	RT_range = [RT - max_RT_deviation, RT + max_RT_deviation]
	RT_range_index = RT_indices(spectrum, RT_range)

	mass_integral = zeros(Float64, length(mass_vals), 2)

	# Loop through every mass
	for (i, mass) in enumerate(mass_vals)
		mass_range = [mass - max_mass_deviation, mass + max_mass_deviation]
		spectrum_XIC = filter_XIC(spectrum, mass_range)

		# Determine maximum within RT and mass range
		(max_intensity, max_index) = findmax(spectrum_XIC[RT_range_index[1]:RT_range_index[2]])
		max_index += RT_range_index[1] - 1


		# Below noise cutoff, set intensity to zero
		if max_intensity < noise_cutoff
			mass_integral[i,:] = [mass, 0]
			continue
		end


		# Used for integrating peak
		integral = max_intensity
		index_range = zeros(Int16, 2)

		# Look for left and right side of peak
		for (j, direction) in enumerate([-1, 1])
			current_intensity = min_intensity = max_intensity
			current_index = min_index = max_index
			count_intensity_incr = 0

			# Continue until below noise_cutoff or begin/end of spectrum reached
			while current_intensity > noise_cutoff & !(current_index in [1, length(spectrum_XIC)])
				current_index += direction
				current_intensity = spectrum_XIC[j]

				if current_intensity < min_intensity
					min_intensity = current_intensity 
					min_index = current_index
					count_intensity_incr = 0
				elseif count_intensity_incr < 5
					count_intensity_incr += 1
				# Peaks overlap, found increase in intensity 5 consecutive times
				else
					current_index = min_index
					println("WARNING: Peak overlap at RT: $(spectrum["RT"][current_index]) and index: $current_index")
					break
				end
			end

			index_range[j] = current_index
		end	

		integral = sum(spectrum_XIC[index_range[1]:index_range[2]])

		mass_integral[i,:] = [mass, integral]
		
	end

	return mass_integral

end


function plt(mass, mass2=0, RT=6.66, RT_deviation=0.1)
	"""Plot spectrum of given mass
	plt(mass, (mass2), (RT), (RT_deviation) )

	only mass given: XIC at specific mass
	mass and mass2: XIC between range
	mass = 0: Spectrum not filtered

	RT given: Zooms in around given RT (6.66 default)
	RT = 0: Shows complete RT spectrum
	RT_deviation: Deviation around RT (default 0.1)
	"""

	# TODO gives error at start and end

	max_mass_deviation = 0.05 
	RT_range = [RT - RT_deviation, RT + RT_deviation]
	RT_range_index = RT_indices(spectrum, RT_range)

	if mass > 0 & mass2 <= 0
		mass_range = [mass - max_mass_deviation, mass + max_mass_deviation]
	elseif mass > 0
		mass_range = [mass, mass2]
	else
		mass_range = [0, Inf]
	end

	spectrum_XIC = filter_XIC(spectrum, mass_range)

	if RT > 0
		add_left = 100
		add_right = 100
		plot(spectrum["Rt"][RT_range_index[1] - add_left:RT_range_index[2] + add_right], spectrum_XIC[(RT_range_index[1] - add_left):(RT_range_index[2] + add_right)])
	else
		plot(spectrum["Rt"], spectrum_XIC)
	end
end

function RT_indices(spectrum, RT)
	"""Returns index range of given retention time range"""

	# Only possible if Rt is sorted
	if !issorted(spectrum["Rt"])
		error("Error: Retention time array is not in order")
	end
	
	return [findfirst(i -> i >= RT[1], spectrum["Rt"]), findlast(i -> i <= RT[2], spectrum["Rt"])]
end

function filter_XIC(spectrum, mass)
	"""Returns filtered spectrum based on mass range"""
	
	# Store which indices are between mass range
	filter = (spectrum["Mz_values"] .>= mass[1]) .& (spectrum["Mz_values"] .<= mass[2])

	# For every scan (Rt), sum mz intensities at chosen mz values
	spectrum_filtered = zeros(Int, size(spectrum["Mz_values"][:, 1]))
	for i=1:length(spectrum_filtered)
		for j=1:length(filter[i, :])
			if filter[i, j] == true
				spectrum_filtered[i] += spectrum["Mz_intensity"][i, j]
			end
		end
	end

	return spectrum_filtered

end


function batch_import(pathin)
	"""Returns list containing all imported spectra from folder"""
	mz_thresh = [0, 0]
	Int_thresh = 0
	files = readdir(pathin)
	files_supported = zeros(Bool, length(files),1)

	# Determine supported files
	for i=1:length(files)
		if lowercase(files[i][end-4:end]) == "mzxml" || lowercase(files[i][end-3:end]) == "cdf"
			files_supported[i] = true
		end
	end

	# Save indices of supported files
	supported_indices=findall(x -> x == true, files_supported)
	num_of_spectra = length(supported_indices)

	# Create array for storing spectra
	spectra = Array{Any}(undef,size(supported_indices,1))

	# Load files into spectra array
	for i in supported_indices
		filenames = [files[i]]
		spectra[i] = import_files(pathin,filenames,mz_thresh,Int_thresh)
		println("Read spectra $(i[1]) of $num_of_spectra")
	end

	return spectra
end

main()

