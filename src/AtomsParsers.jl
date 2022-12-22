module AtomsParsers

export AtomsFile, parse_file

using Unitful, AtomsBase

"""
    AtomsFile

Abstract type to be extended by any file format representing atomic data.
"""
abstract type AtomsFile end

parse_file(::T) where T <: AtomsFile = error(
    "`parse_file` not implemented for type `$T`. Add a method `parse_file(::$T)` for this to work."
)


export XSF, read_xsf
include("xsf.jl")

end
