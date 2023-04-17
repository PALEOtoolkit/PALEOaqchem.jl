

using Test
using Documenter


@testset "PALEOaqchem all" begin

include("PALEOcarbchem/runtests.jl")

include("runcarbchemtests.jl")

doctest(PALEOaqchem; manual=false)  

end