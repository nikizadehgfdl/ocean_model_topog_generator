# A tool set for creating model topography
This tool creates a model topography for a given horizontal model grid and a given observation dataset for height. It uses a recursive refine-sample-coarsen scheme to interpolate the observation dataset onto the model grid. 
This tool is based on the thin wall interpolation algorithm of Adcroft et.al. for interpolating high-resolution topography for use in numerical circulation models of the ocean and atmosphere which is hosted at https://github.com/adcroft/thin-wall-topography. 

# Citation
Ocean Modelling (Elsevier), http://www.journals.elsevier.com/ocean-modelling/. If you are publishing data or research that uses this software or algorithm, please cite the original algorithm as follows:
 Adcroft, A., 2013: Representation of topography by porous barriers and objective interpolation of topographic data. Ocean Modelling, doi:10.1016/j.ocemod.2013.03.002.(http://dx.doi.org/10.1016/j.ocemod.2013.03.002)

# License
Distributed under the Gnu General Public License, version 3. See LICENSE.txt for details or http://www.gnu.org/licenses/gpl.html.

# Sample usage
create_topog_refinedSampling.py --hgridfilename PATH_TO/ocean_hgrid.ncSO.nc --outputfilename topog_SO_p25.nc --source_file PATH_TO/GEBCO_2020.nc --source_lon lon --source_lat lat --source_elv elevation

