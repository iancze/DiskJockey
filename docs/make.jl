using Documenter, DiskJockey

makedocs(
    modules=[DiskJockey],
    format = :html,
    assets = ["assets/dust.mp4", "assets/vis.mp4", "assets/gas.mp4", "assets/uv_spacings.png", "assets/disk.svg"],
    sitename = "DiskJockey.jl Documentation",
    authors = "Ian Czekala and contributors.",
    analytics = "UA-5472810-7",
    pages = ["Home" => "index.md",
    "Introduction" => "dynamical_mass_intro.md",
    "Installation" => "installation.md",
    "Models" => "models.md",
    "Formats" => "formats.md",
    "RADMC Setup" => "RADMC3D_setup.md",
    "Cookbook" => "cookbook.md",
    "Priors" => "priors.md",
    "API" => "api.md",
    "Changelog" => "changelog.md"],
    # html_prettyurls = !("local" in ARGS),
)

deploydocs(
    repo = "github.com/iancze/DiskJockey.git",
    julia = "0.6",
    branch = "gh-pages",
    latest = "master",
    osname = "linux",
    target = "build",
    deps = nothing,
    make = nothing,
)
