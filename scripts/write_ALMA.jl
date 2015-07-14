# Generate a set of model visibilities and write them to HDF5
# Rather than run everything in parallel, like `burma_shave.jl`, do this
# all in serial, here.

# The follow up step is to use a the python script `write_ALMA.py` to convert
# these from HDF5 and pack into numpy arrays.

using constants
using image
using model
using HDF5
using visibilities
using gridding

using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "--run_index", "-r"
    help = "Output run index"
    arg_type = Int
    "config"
    help = "a YAML configuration file"
    required = true
end

parsed_args = parse_args(ARGS, s)

import YAML
config = YAML.load(open(parsed_args["config"]))

# read the wavelengths for all channels
fid = h5open(config["data_file"], "r")
lams = read(fid["lams"]) # [μm]
close(fid)
nchan = length(lams)

# Read the parameters from the config file
pp = config["parameters"]

a = (st)->pp[st][1]

# The parameters we'll be using
pars = Parameters(a("M_star"), a("r_c"), a("T_10"), a("q"), a("gamma"), 10.^a("logM_CO"), a("ksi"), a("dpc"), a("incl"), a("PA"), a("vel"), a("mu_RA"), a("mu_DEC"))

# Create the model grid
grd = config["grid"]
global const grid = Grid(grd["nr"], grd["ntheta"], grd["r_in"], grd["r_out"], true)


global const basedir = ""
# Regenerate all of the static files (e.g., amr_grid.inp)
# so that they may be later copied
write_grid(basedir, grid)

vel = pars.vel # [km/s]
# RADMC conventions
incl = pars.incl # [deg]
PA = pars.PA # [deg] Position angle runs counter clockwise, due to looking at sky.
npix = config["npix"] # number of pixels

# Doppler shift the dataset wavelength to rest-frame wavelength
beta = vel/c_kms # relativistic Doppler formula
shift_lams =  lams .* sqrt((1. - beta) / (1. + beta)) # [microns]

write_model(pars, basedir, grid)
write_lambda(shift_lams, basedir)

# files = ["lines.inp", "molecule_co.inp", "wavelength_micron.inp"]
# for file in files
#     cp("../" * file, file)
# end

cp("radmc3d.inp.gas", "radmc3d.inp")


run(`radmc3d image incl $incl posang $PA npix $npix loadlambda`)

im = imread()
skim = imToSky(im, pars.dpc)
corrfun!(skim, pars.mu_RA, pars.mu_DEC) # alpha = 1.0

dvarr = DataVis(config["data_file"])
mvarr = Array(DataVis, nchan)
chi2s = Array(Float64, nchan)

for i=1:nchan

    dv = dvarr[i]

    # FFT the appropriate image channel
    vis_fft = transform(skim, i)

    # Apply the phase correction here, since there are fewer data points
    phase_shift!(vis_fft, pars.mu_RA, pars.mu_DEC)
    # Interpolate the `vis_fft` to the same locations as the DataSet
    mvis = ModelVis(dv, vis_fft)

    dvis = visibilities.ModelVis2DataVis(mvis)

    # Complex conjugate for SMA convention
    visibilities.conj!(dvis)
    mvarr[i] = dvis

    chi2s[i] = visibilities.chi2(dv, mvis)
end

println("Chi2s", chi2s)
N = nchan * 2 * length(dvarr[1].VV)
println("Chi_R", sum(chi2s)/N)

visibilities.write(mvarr, "data/AKSco/AKSco_model.hdf5")
