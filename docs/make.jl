using Documenter

import PALEOaqchem

using DocumenterCitations

bib = CitationBibliography(
    joinpath(@__DIR__, "src/paleoaqchem_references.bib");
    style=:authoryear,
)

makedocs(;
    sitename = "PALEOaqchem Documentation", 
    pages = [
        "index.md",
        # "Examples and Tutorials" => examples_pages,
        # no Design docs yet
        # "Design" => [
        #     "COPSE_Domains.md",
        # ],
        # no HOWTO docs yes
        "Reference" => [
            "Organic matter and remineralization.md",
            "Generic Chemistry.md",
            "Secondary Redox.md",
            "Carbonate chemistry.md",
            "Molecular Diffusivity.md",
            "Isotopes.md",
            "functions.md",
            
        ],
        "References.md",
        "indexpage.md",
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    plugins = [bib],
)

@info "Local html documentation is available at $(joinpath(@__DIR__, "build/index.html"))"

deploydocs(
    repo = "github.com/PALEOtoolkit/PALEOaqchem.jl.git",
)
