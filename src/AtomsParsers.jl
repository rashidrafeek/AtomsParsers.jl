module AtomsParsers

"""
    AtomsFile

Abstract type to be extended by any file format representing atomic data.
"""
abstract type AtomsFile end

export XSF, read_xsf
include("xsf.jl")

end
