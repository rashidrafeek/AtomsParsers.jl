"""
    XSF

Type containing the data in an XSF file as defined by the [XSF specification](http://www.xcrysden.org/doc/XSF.html)
"""
struct XSF end

function read_xsf(filename)
    filedata = readlines(filename)
    alldata = Dict{String, Any}()
    nlines = length(filedata)

    # Find type of structure
    sectype, curline = _xsf_next_section(filedata, 0)

    if sectype == "CRYSTAL" # Periodic crystal structure
        while(curline < nlines)
            curline = _xsf_read_crystalsection!(alldata, filedata, curline)
        end
    else
        error("Currently only periodic crystal structures are implemented.")
    end

    alldata, filedata, curline
end

# Read the next section title excluding comment lines and empty lines
function _xsf_next_section(lines, curline)
    curline += findfirst(!startswith(r"^ *#|^$"), lines[curline+1:end])
    sectiontype = strip(lines[curline])

    return sectiontype, curline
end

# Read various crystal structure data
function _xsf_read_crystalsection!(alldata, lines, curline)
    sectiontype, curline = _xsf_next_section(lines, curline)

    if sectiontype == "PRIMVEC" || sectiontype == "CONVVEC"
        alldata[sectiontype] = _xsf_read_lattice(lines[curline+1:curline+3])

        curline += 3
    elseif sectiontype == "PRIMCOORD"
        natoms, _ = parse.(Int, split(lines[curline+1]))

        atomtypes, prim_coord = _xsf_read_coords(
            lines[curline+2:curline+(natoms+1)], natoms
        )
        alldata[sectiontype] = (atomtypes, prim_coord)

        curline += natoms+1
    elseif sectiontype == "BEGIN_BLOCK_DATAGRID_3D"
        if !haskey(alldata, sectiontype)
            alldata[sectiontype] = Dict{String, Any}()
        end
        blockname = strip(lines[curline+1])
        curline += 1

        datanameraw = lines[curline+1]
        while !contains(datanameraw, "END_BLOCK")
            dataname = split(datanameraw, "DATAGRID_3D_")[end]
            nx, ny, nz = parse.(Int, split(lines[curline+2]))
            dataendline = findfirst(
                contains("END_DATAGRID_3D"), 
                lines[curline+1:end]
            )

            origin, spanvec, datagrid = _xsf_read_datagrid3d(
                lines[curline+3:curline+dataendline], nx, ny, nz
            )

            alldata[sectiontype][blockname] = (;
                size = (nx,ny,nz),
                origin, spanvec, datagrid
            )

            curline += dataendline
            datanameraw = lines[curline+1]
        end
        
        curline += 1
    else
        error("Unknown section \"$(sectiontype)\"")
    end

    return curline
end

# Read lattice data
function _xsf_read_lattice(lines)
    lattice_vecs = Vector{Vector{Float64}}(undef, 3)
    
    for (i,line) in pairs(lines)
        lattice_vecs[i] = parse.(Float64, split(line))
    end

    return lattice_vecs
end

# Read coordinates data
function _xsf_read_coords(lines, natoms)
    atomtypes = Vector{Union{Symbol, Int}}(undef, natoms)
    coordinates = Array{Float64}(undef, natoms, 3)
    
    for (i,line) in pairs(lines)
        linedata = split(line)
        
        atnum = tryparse(Int, linedata[1])
        if isnothing(atnum)
            atomtypes[i] = Symbol(linedata[1])
        else
            atomtypes[i] = atnum
        end

        coordinates[i,:] .= parse.(Float64, linedata[2:4])
    end

    return atomtypes, coordinates
end

# Read 3D data grid
function _xsf_read_datagrid3d(lines, nx, ny, nz)
    orig = parse.(Float64, split(lines[1]))
    spanvec = Vector{Vector{Float64}}(undef, 3)
    datagrid = Array{Float64}(undef, nx, ny, nz)

    for i in 1:3
        spanvec[i] = parse.(Float64, split(lines[i+1]))
    end

    counter = 0
    for line in lines[5:end-1]
        vals = parse.(Float64, split(line))
        nvals = length(vals)

        datagrid[counter+1:(counter+nvals)] .= vals

        counter += nvals
    end

    return orig, spanvec, datagrid
end
