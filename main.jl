using MS_Import
function main()
	# filenames=["B039695.mzXML", "B039699.mzXML"]

	folder = "P:/Git/bachelor_project"
	mzxml_path = "data/mzxml"
	pathin = joinpath(folder, mzxml_path)

	spectra = batch_import(pathin)

	spectra[2]


end


function batch_import(pathin)
	"""Returns list containing all spectra in folder"""
	mz_thresh = [0, 0]
	Int_thresh = 0
	files = readdir(pathin)
	files_supported = zeros(length(files),1)

	# Determine supported files
	for i=1:length(files)
		if lowercase(files[i][end-4:end]) == "mzxml" || lowercase(files[i][end-3:end]) == "cdf"
			files_supported[i] = 1
		end
	end

	# Save indices of supported files
	supported_indices=findall(x -> x > 0, files_supported)
	
	# Create array for storing spectra
	spectra = Array{Any}(undef,size(supported_indices,1))

	# Load files into spectra array
	for i in supported_indices
		filenames = [files[i]]
		spectra[i]=import_files(pathin,filenames,mz_thresh,Int_thresh)
	end

	return spectra
end

main()

