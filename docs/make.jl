using Documenter

import PALEOaqchem

using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "src/paleoaqchem_references.bib"))

makedocs(bib, sitename="PALEOaqchem Documentation", 
    pages = [
        "index.md",
        # "Examples and Tutorials" => examples_pages,
        # no Design docs yet
        # "Design" => [
        #     "COPSE_Domains.md",
        # ],
        # no HOWTO docs yes
        "Reference" => [
            "PALEOaqchem Reactions.md",
            "PALEOaqchem functions.md",
            "PALEOcarbchem.md",
        ],
        "References.md",
        "indexpage.md",
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
)

@info "Local html documentation is available at $(joinpath(@__DIR__, "build/index.html"))"

deploydocs(
    repo = "github.com/PALEOtoolkit/PALEOaqchem.jl.git",
)
