#=
This is a modified version of MS_Import by Saer Samanipour: https://bitbucket.org/SSamanipour/ms_import.jl/src/master/

Copyright (c) 2020 Saer Samanipour, PhD, Computational Mass Spec Lab (CMSL)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
=#

using LightXML
using Codecs

#############################################################################
# import_files_light

function import_files_light(pathin,filenames,mz_thresh=[0,0],Int_thresh=0)

    m=split(filenames[1],".")

    if m[end] == "mzxml" || m[end] == "mzXML" || m[end] == "MZXML"

        if isa(filenames,Array)==1
            path_in=joinpath(pathin,filenames[])
        else
            path_in=joinpath(pathin,filenames)
        end

        chrom = mzxml_read_light(path_in,mz_thresh,Int_thresh)
    end

    return (chrom)

end

#############################################################################
# mzxml_read_light

function mzxml_read_light(path2mzxml,mz_thresh,Int_thresh)
        xdoc = parse_file(path2mzxml);
        xroot = root(xdoc);

        if LightXML.name(xroot) != "mzXML"
            error("Not an mzXML file")
        end
        # Find the msRun node
        msRun = find_element(xroot, "msRun");
        mslevels = ["1"]
        row_length, col_length = determine_matrix_size(msRun,mslevels)
        retentionTime,mz,mz_int=read_scan_light(msRun,mz_thresh,Int_thresh, row_length, col_length, mslevels)


    chrom=Dict("MS1" => Dict{String, Any}(
        "Rt" => retentionTime,
        "Mz_values" => mz,
        "Mz_intensity" => mz_int

    ))
    free(xdoc)
    return chrom
    
end # function



    
###############################################################################
# Reading the scans

function read_scan_light(msRun,mz_thresh,Int_thresh, row_length, col_length, mslevels)
    retentionTime = zeros(Float32, row_length)
    mz, mz_int = zeros(Float32, row_length, col_length), zeros(Int64, row_length, col_length)
    i = 0
    for line in child_elements(msRun)
        if LightXML.name(line) == "scan" && attribute(line, "msLevel") in mslevels
            i += 1
            retentionTime1,mz1,I=read_scan_info_light(line,mz_thresh,Int_thresh)
            # if retentionTime1 == false #TEMP
            #     continue
            # end
            retentionTime[i] = round(retentionTime1 / 60, digits=2)
            mz[i, 1:length(mz1)] = mz1
            mz_int[i, 1:length(I)] = I
            # retentionTime=vcat(retentionTime,retentionTime1) #TEMP
            # mz=vcat(mz,[mz1])
            # mz_int=vcat(mz_int,[I])
        end

    end
    return(retentionTime,mz,mz_int)

end # function




###############################################################################
# Extracting info from the scans

function read_scan_info_light(line,mz_thresh,Int_thresh)

    precisiondict = Dict("32" => (Float32, Float32, Float32), "64" => (Float64, Float64, Float64))
    # line = find_element(msRun,"scan")
    # println(line)

    tstr=attribute(line, "retentionTime")
    retentionTime = parse(Float64, tstr[3:end-1])

    # # Only read ms level 1 #TEMP
    # if attribute(line, "msLevel") != "1"
    #     return false, 0, 0
    # end

    peak = find_element(line, "peaks")
    # 1840-element Array{UInt8,1}:
    if attribute(peak, "compressionType") ==  "zlib"
        data = decode(Zlib, decode(Base64, content(peak)))
    else
        data = decode(Base64, content(peak))
    end

    # data = decode(Zlib, decode(Base64, content(peak)))
    # println(bytestring(data))
    TI, T, nochildren = precisiondict[attribute(peak, "precision")]
    A = reinterpret(TI, data)

    bo = attribute(peak, "byteOrder")
    if bo == "network"
        ntoh!(A)
    else
        error("Don't know what to do with byteOrder $bo")
    end
    I = A[2:2:end]

    mz = reinterpret(T, A[1:2:end])

    
    if Int_thresh > 0
        mz1=mz[findall(x -> Int_thresh<= x, I)]
        I1=I[findall(x -> Int_thresh<= x, I)]
    else
        mz1=mz
        I1=I
    end

    if mz_thresh[2]>0
        mz2=mz1[findall(x -> mz_thresh[1]<= x <= mz_thresh[2], mz1)]
        I2=I1[findall(x -> mz_thresh[1]<= x <= mz_thresh[2], mz1)]
    elseif mz_thresh[1] == 0
        mz2=mz1[mz1 .>= mz_thresh[1]]
        I2=I1[mz1 .>= mz_thresh[1]]

    end


    return(retentionTime,mz2,I2)
end # function



###############################################################################
# Modified for cocaine batch processing
# Determines scan with highest peak count and number of scans
# This determines spectrum mz and intensity matrix size

function determine_matrix_size(msRun,mslevels)

    row_length, col_length = 0, 0
    for line in child_elements(msRun)
        if LightXML.name(line) == "scan" && attribute(line, "msLevel") in mslevels
            row_length += 1
            peakscount = parse(Int, attribute(line, "peaksCount"))

            # If peakscount is longer than current length, change to peakscount
            col_length = peakscount > col_length ? peakscount : col_length
        end
    end

    return row_length, col_length

end # function

##############################################################################
# XML encoding


function ntoh!(A)
    for i = 1:length(A)
        A[i] = ntoh(A[i])
    end
end


##############################################################################

