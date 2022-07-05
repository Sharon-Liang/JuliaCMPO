using Documenter
using JuliaCMPO

makedocs(
    sitename = "JuliaCMPO",
    format = Documenter.HTML(),
    modules = [JuliaCMPO]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
