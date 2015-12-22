#File Formats

Interferometric data from reduction software like CASA or MIRIAD is generally stored in measurement sets (`*.ms`) or FITS files (`*.fits`), and usually contains a lot of ancillary data unnecessary for a dynamical mass measurement (or fitting in the UV plane in general). To reduce dependency on these outside programs and overall memory footprint, we have made it so that all of the Julia code only interfaces with the minimal necessary visibility data stored in an HDF5 file.

The format of this datafile is as follows:

The key Julia code used to read this file is provided in `/src/visibilities.jl`

Which reads the HDF5 file and stores each channel in an array of `DataVis` instances. For example, this `data.hdf5` file contains the following datasets in arrays of (nrows, ncols) form. `nchan` is the number of channels in the dataset, `nvis` is the number of complex visibilities in the measurement set.

    data.hdf5
      lams # [μm] (nchan) Wavelength (in microns) corresponding to channel
      uu # [kλ] (nchan, nvis) Vectors of the u locations in kilolambda
      vv # [kλ] (nchan, nvis) Vectors of the v locations in kilolambda
      real # [Jy] (nchan, nvis) real component of the complex visibilities
      imag # [Jy] (nchan, nvis) imaginary component of the complex visibilities
      invsig # [1/Jy] (nchan, nvis) the inverse of sigma for each visibility (1/sigma)


and an example datafile is provided for AK Sco, available for download here (TBD). The invsig field may seem like an non-traditional way to store visibility weights (typically measured in `[1/Jy^2]`), but enables quick computation of the `chi^2` statistic via Julia's `sumabs2` routine.
