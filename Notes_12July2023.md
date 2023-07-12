### Building

### Steps in running
1. Get input files from analysis team.

   Input files should contain oscillated spectra for the grid point we want to calculate. Format should be 2D matrix: row is the grid point, column the bins of the spectra for each sub-component for each channel (so decide on bin numbering unique across all channels).

   Could also contain data spectrum if that's already there, and the null though that should be a grid point?
   
   Input files should also contain the full fractional covariance matrix in the same dimensions as the oscillated spectra, not including any (data) statistical uncertainty.

   These should be one HDF5 file that has all of this.   
2. Make sure the function doing the collapsing works (maybe there are some optimizations to be done there)
3. 


### Output file
`grid_point_1d = n_pts_z*n_pts_y*idx_x + n_pts_z*idx_y + idx_z`

  * `last_chi_min`: Min chi2 per grid point per universe (will be all same for wilks scan)
  * `delta_chi`:Delta chi2 per grid point per univers
  * `grid_x`, `grid_y`, `grid_z`: Grid points (x, y, z)
  * `best_grid_pointx`, `best_grid_pointy`, `best_grid_pointz`:Best fit point for minimizer
  * `best_grid_point`: Best fit point in 1d grid point for grid scan
  * `best_start_pointx`: For minimizer, save best fit's starting point `(x,y,z)`
  * `this_chi_mfa`: chi2 with respect to null using MFA spectrum
  * `this_chi_gs`: chi2 with respect to null spectrum from input file (at grid point)
  * `this_chi_raw`: chi2 of MFA spectrum with respect to truth spectrum from grid
  * `speccoll_gs`: collapsed spectrum from input file (at grid point under calc) -- wilks only
  * `speccold_mfa`: collapsed spectrum from mfa (at grid point under calc) -- wilks only
  * `n_iter`: number of iterations
  * `warn_flag`: warning flag on if the best fit hasn't converged
  * `num_iters`: same as `n_iter`?
  * `i_grid`: index of grid under calc (see formula above)
  * `i_univ`: index of universe

### Notes on functions 
  * `writeGrid`: determines the values of the osc parameters (not needed if done already, like in WC analysis)

### Potential improvments
  * Have a mode where we can save the MFA-produced spectra to allow for later comparison
    
    (Previously was done with just a hard-coded dedicated run.)
  * Remove ROOT entirely from container? We can translate from ROOT to HDF5 outside of the container. WC pieces definitely not using ROOT in calcs. Need to check if SBNFit under the hood is, but likely could be replaced if so.
  * Standardize how the inputs are provided and optimizie the collapsing based on that (e.g. do loops? do matrix multiplication? etc.)
