using Documenter, DiskJockey

makedocs(
    modules = [DiskJockey],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == true, analytics = "UA-5472810-7"),
    sitename = "DiskJockey.jl Documentation",
    authors = "Ian Czekala and contributors.",
    pages = ["Home" => "index.md",
    "Introduction" => "dynamical_mass_intro.md",
    "Installation" => "installation.md",
    "Conventions" => "conventions.md",
    "Models" => "models.md",
    "Formats" => "formats.md",
    "RADMC Setup" => "RADMC3D_setup.md",
    "Cookbook" => "cookbook.md",
    "Priors" => "priors.md",
    "API" => "api.md",
    "Changelog" => "changelog.md"],
)

deploydocs(
    repo = "github.com/iancze/DiskJockey.git",
)
