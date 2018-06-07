% Topographic Analysis Kit
% Version 0.1 pre-release  1-Jun-2018
%
% These are a series of matlab functions that build upon the functionality 
% of TopoToolbox [https://github.com/wschwanghart/topotoolbox]. Each function
% contains a header with basic functionality info along with expected inputs 
% and possible outputs. This readme compiles some of that information and lays
% out possible workflows, see also 'Suggested_Workflows.pdf' within the repository 
% for complete flowchart. 
%
% If you encounter errors or have suggestions please contact:
%
% Adam M. Forte
% aforte8 'at' lsu.edu
%
% Publication of tools is pending, please contact Adam if you are using the codes
%
%
% COMPLETE FUNCTION LIST
%
% ----------------------------------------------------------------------------------------------------------------------------
% Basin2Raster
%
% Function takes outputs from 'ProcessRiverBasins' function and produces a single GRIDobj with individual drainage
% basins (as selected by 'ProcessRiverBasins' and 'SubDivideBigBasins') assinged various values
%
% Required Inputs:
%	DEM - GRIDobj of full extent of datasets
% 	valueOI - value to assign to basins, acceptable inputs are:
%		'ksn' - mean ksn value of basin
%		'gradient' - mean gradient of basin
%		'elevation' - mean elevation of basin
%		'chir2' - R^2 value of chi-z fit (proxy for disequilibrium)
%		'drainage_area' - drainage area in km2 of basin
%		'id' - basin ID number (i.e third column RiverMouth output)
%   	'theta' - best fit concavity resultant from the topo toolbox chiplot function (will not work if method for 
%			ProcessRiverBasin was 'segment')
%		'NAME' - where name is the name provided for an extra grid (i.e. entry to second column of 'add_grid' or entry to 
%			third column of 'add_cat_grid')
%	location_of_data_files - full path of folder which contains the mat files from 'ProcessRiverBasins' as a string
%
% Optional Inputs:
%	file_name_prefix ['basins'] - prefix for outputs, will automatically append the type of output, i.e. 'ksn', 'elevation', etc
%	location_of_subbasins ['SubBasins'] - name of folder that contains subbasins of interest (if you created subbasins using
%		"SubDivideBigBasins"), expected to be within the main Basin folder provided with "location_of_data_files".
%	method ['subdivided'] - method used for subdividing watersheds. If you used 'ProcessRiversBasins' and then
%		'SubDivideBigBasins' or if you only used 'ProcessRiverBasins' but did not pick any nested catchments, i.e.
%		none of the river mouths supplied to 'ProcessRiverBasins' were within the catchment boundaries of other 
%		watersheds for which you provided river mouths, then you should use use 'subdivided' which is the default
%		so you do not need to specify a value for this property. If you picked nested catchments manually and then
%		ran 'ProcessRiverBasins' you should use 'nested'.
%
%
% ----------------------------------------------------------------------------------------------------------------------------
% Basin2Shape
%
% Function to take the outputs from 'ProcessRiverBasins' and 'SubDivideBigBasins' and produce a single shapefile showing the outlines of polygons
% 	and with commonly desired attributes from the results of 'ProcessRiverBasins' etc. See below for a full list of fields that the output shapefile
% 	will include. If additional grids were provided to 'ProcessRiverBasins', mean and standard error values for those grids will be auto-populated in
% 	the shapefile and the name of the fields will be the character array provided in the second column of additional grids input. This function also
% 	allows you to input a list of additional fields you wish to include (see Optional Inputs below). If you would rather create a GRIDobj with specified
% 	values, use 'Basin2Raster'.
%
% Required Inputs:
%		DEM - GRIDobj of the DEM originally used as input for 'ProcessRiverBasins' for the basins of interest.
% 		location_of_data_files - full path of folder which contains the mat files from 'ProcessRiverBasins'
%
% Optional Inputs:
%		location_of_subbasins ['SubBasins'] - name of folder that contains subbasins of interest (if you created subbasins using
%			"SubDivideBigBasins"), expected to be within the main Basin folder provided with "location_of_data_files". Note that if you do not provide
%			the correct directory name for the location of the subbasins, subbasin values will not be included in the output regardless of your choice
%			for the "include" parameter.
%		shape_name ['basins'] - name for the shapefile to be export, must have no spaces to be a valid name for ArcGIS and should NOT include the '.shp';
%		include ['all'] - parameter to specify which basins to include in building the shapfile. The default 'all' will include all basin mat files in the 
%			folder you specify. Providing 'subdivided' will check to see if a given main basin was subdivided using 'SubdivideBigBasins' and then only include 
%			the subdivided versions of that basin (i.e. the original main basin for those subbasins will not be included in the shapefile). Providing 'bigonly'
%			will only include the original basins produced by 'ProcessRiverBasins' even if 'SubDivideBigBasins' was run. If 'SubDivideBigBasins' was never run,
%			result of 'all' and 'bigonly' will be the same.
%		extra_field_values [] - cell array of extra field values you wish to include. The first column in this cell array must be the river basin number
%			(i.e. the identifying number in the third column of the RiverMouth input to ProcessRiverBasins or the number generated for the basin in
%			SubDivideBigBasins). Only one row per river basin number is allowed and ALL river basin numbers in the basins being processed must have a value
%			associated with them. Additional columns are interpreted as the values with which you wish to populate the extra fields. These can either be character
%			arrays or numbers, other values will results in an error. 
%		extra_field_names [] - a 1 x m cell array of field names, as characters (no spaces as this won't work with shapefile attributes), associated with the field values. 
%			These must be in the same order as values are given in extra_field_values. If for example your extra_field_values cell array is 3 columns with the river number, 
%			sample name, and erosion rate then your extra_field_names cell array should include entries for 'sample_name' and 'erosion_rate' in that order. 
%		uncertainty ['se'] - parameter to control which measure of uncertainty is included, expects 'se' for standard error (default), 'std' for standard deviation, or 'both'
%			to include both standard error and deviation.
%		populate_categories [false] - logical flag to add entries that indicate the percentage of a watershed occupied by each category from a categorical grid, e.g. if you
%			provided an entry for 'add_cat_grids' to ProcessRiverBasins that was a geologic map that had three units, 'Q', 'Mz', and 'Pz' and you set 'populate_categories' 
%			to true there will be field names in the resulting shapefile named 'Q', 'Mz', and 'Pz' and the values stored in those columns will correspond to the percentage 
%			of each basin covered by each unit for each basin. Setting populate_categories to true will not have any effect if no entry was provided to 'add_cat_grids' when
%			running ProcessRiverBasins.
%
% Output:
%		Outputs a mapstructure (MS) and saves a shapefile with the following default fields:
%			river_mouth - river mouth number provided to ProcessRiverBasins
%			drainage_area - drainage area of basin in km^2
%			center_x - x coordinate of basin in projected coordinates
%			center_y - y coordinate of basin in projected coordinates
%			outlet_elevation - elevation of pour point in m
%			mean_el - mean elevation of basin in meters
%			max_el - maximum elevation of basin in meters
%			mean_ksn - mean channel steepenss
%			mean_gradient - mean gradient
%		Either standard errors, standard deviations or both will be populated for elevation, ksn, and gradient depending on value of 'uncertainty'
%		Mean and standard error / standard deviation / both values will be populated for any additional grids
%
% ----------------------------------------------------------------------------------------------------------------------------
% BasinPicker
%
% Function takes results of makes streams and allows for interactive picking of basins (watersheds). Function was
%   designed intially for choosing basins suitable for detrital analyses (e.g. Be-10 cosmo). Displays two panel figure 
%   with topography colored by elevation and local relief on which to pick individual basins. After the figure displays,
%   it will wait until you press enter to begin the watershed picking process. This is to allow you to zoom, pan, etc to 
%   find a stream you are interested in. When you click enter, cross hairs will appear in the elevation map so you can 
%   select a pour point. Once you select a pour point, a new figure will display this basin and stream to confirm that's 
%   the watershed you wanted (it will also display the drainage area). You can either accept this basin or reject it if 
%   it was misclick. If you accept it will then display a new figure with the chi-z and longitudinal profiles for that basin. 
%   It will then give you a choice to either save the choice or discard it. Finally it will ask if you want to keep picking 
%   streams, if you choose yes (the default) it will start the process over. Note that any selected (and saved) pour point 
%   will be displayed as a white circle on the main figure. As you pick basins the funciton saves a file called 'Outlets.mat'
%   that contains the outlets you've picked so far. If you exit out of the function and restart it later, it looks for this 
%   Outlets file in the current working directory so you can pick up where you left off.
%
%
% Required Inputs:
%       DEM - GRIDobj of the DEM
%       FD - FLOWobj from the supplied DEM
%       A - Flow accumulation grid (GRIDobj)
%       S - STREAMobj derived from the DEM
%
% Optional Inputs:
%       ref_concavity [0.50]- reference concavity for chi-Z plots
%       rlf_radius [2500] - radius in map units for calculating local relief OR
%       rlf_grid [] - if you already have a local relief grid that you've calculated, you can provide it with 'rlf_grid', it must be a GRIDobj and must be the same
%           coordinates and dimensions as the provided DEM.
%       extra_grid [] - sometimes it can be useful to also view an additional grid (e.g. georeferenced road map, precipitation grid, etc) along with the DEM and relief.
%           This grid can be a different size or have a different cellsize than the underlying dem (but still must be the same projection and coordinates system!), it will be
%           resampled to match the provided DEM. 
%       cmap ['jet'] - colormap to use for the displayed maps. Input can be the name of a standard colormap or a nx3 array of rgb values
%           to use as a colormap.  
%       conditioned_DEM [] - option to provide a hydrologically conditioned DEM for use in this function (do not provide a conditoned DEM
%           for the main required DEM input!) which will be used for extracting elevations. See 'ConditionDEM' function for options for making a 
%           hydrological conditioned DEM. If no input is provided the code defaults to using the mincosthydrocon function.
%       interp_value [0.1] - value (between 0 and 1) used for interpolation parameter in mincosthydrocon (not used if user provides a conditioned DEM)
%       plot_type ['vector'] - expects either 'vector' or 'grid', default is 'vector'. Controls whether all streams are drawn as individual lines ('vector') or if
%           the stream network is plotted as a grid and downsampled ('grid'). The 'grid' option is much faster for large datasets, 
%           but can result in inaccurate site selections. The 'vector' option is easier to see, but can be very slow to load and interact with.
%       threshold_area [1e6] - used to redraw downsampled stream network if 'plot_type' is set to 'grid'
%
% 
% Outputs:
%       Outlets - n x 3 matrix of sample locations with columns basin number, x coordinate, and y coordinate (valid input to 'ProcessRiverBasins'
%           as 'river_mouths' parameter)
%
% Examples:
%       [Outs]=DetritalSamplePicker(DEM,FD,A,S);
%       [Outs]=DetritalSamplePicker(DEM,FD,A,S,'rlf_radius',5000);
%       [Outs]=DetritalSamplePicker(DEM,FD,A,S,,'rlf_grid',RLF);
%
% ----------------------------------------------------------------------------------------------------------------------------
% BasinStatsPlots
%
% Function to take the complied outputs from 'ProcessRiverBasins' and 'SubDivideBigBasins' and produce various plots
%	of aggregated basin values. 
%
% Required inputs:
%	basin_table - Table output from 'CompileBasinStats'
%	plots - Type of plot you want to produce, valid inputs are:
%		'grd_ksn' - plot of mean basin gradient vs mean basin channel steepness (e.g. see Forte et al, 2016, Earth and Planetary 
%					Science Letters for discussion of use of these plots)
%		'grd_rlf' - similar to 'grd_ksn' but uses local relief instead of ksn, requires that relief was calculated when running
%					ProcessRiverBasins. Assumes relief radius is 2500 (can set alternative radii with 'rlf_radius' optional parameter)
%		'rlf_ksn' - plot of mean basin relief vs mean basin channel steepness
%		'compare_filtered' - plot comparing mean values vs filtered mean values if you ran 'CompileBasinStats' and filtered by a category
%		'category_mean_hist' - if you calculated 'means_by_category' when running 'CompileBasinStats', you can plot distributions of the
%					means by category as histograms using this option. Requires an input to 'cat_mean1'
%		'category_mean_compare' -if you calculated 'means_by_category' for more than one value (e.g. both gradient and ksn), you can compare
%					the mean values by category using this plot. Requires inputs to both 'cat_mean1' (value that will be plotted on x axis) 
%					and 'cat_mean2' (value that will be plotted on y axis)
%		'stacked_hypsometry' - plot hypsometries for the basins
%		'xy' - generic plot, requires entries to optional 'xval' and 'yval' inputs
%
% Optional Inputs:
%	uncertianty ['se'] - uncertainty to value use for plots, valid options are 'se' (standard error), 'std' (standard deviation), or 'none'. 
%		Providing 'none' indicates you do not want to plot errorbars. Behavior of this option will depend on how you ran ProcessRiverBasins, 
%		e.g. if you only calculated standard deviations when running ProcessRiverBasins but supply 'se'	here, the code will ignore your choice
%		and use the standard deviation values.
%	use_filtered [false] - logical flag to use filtered values for 'grd_ksn', 'grd_rlf', or 'rlf_ksn'. Will only work if you calculated filtered
%		values when running 'CompileBasinStats'.
%	color_by [] - value to color points by, valid for 'grd_ksn','grd_rlf','rlf_ksn', and 'xy'
%	cmap [] - colormap to use if an entry is provided to 'color_by', can be the name of a standard colormap or a nx3 array of rgb values
%		to use as a colormap. 
%	xval [] - value to plot on x axis (name of column as it appears in the provided table) for plot type 'xy'
%	yval [] - value to plot on y axis (name of column as it appears in the provided table) for plot type 'xy'
%	rlf_radius [2500] - radius of relief used when plotting relief related values
%	cat_mean1 [] - category to use for plotting, see 'category_mean_hist' or 'category_mean_compare', valid inputs are 'ksn', 'rlf', 'gradient', or 
%		the name of an additional grid provided to ProcessRiverMeans.
%	cat_mean2 [] - category to use for plotting, 'category_mean_compare' , valid inputs are 'ksn', 'rlf', 'gradient', or the name of an additional grid 
%		provided to ProcessRiverMeans.
%	only_positive [false] - filter out negative values when using either 'category_mean_hist' or 'category_mean_compare'	
%	save_figure [false] - logical flag to save pdfs of all figures produced
%
% ----------------------------------------------------------------------------------------------------------------------------
% CatPoly2GRIDobj
%
% Function to convert a categorical polygon shape file (e.g. a digitzed geologic map) to a GRIDobj. Can
%	be useful for use in 'ProcessRiverBasins'
%
% Required Inputs:
% 	DEM - DEM that you want the output to match 
%	poly_shape - name or path to shapefile containing the categorical data
%	field - field name of categorical data within the shapefile 
%
% Outputs:
%	OUT - GRIDobj of the same size as DEM where values correspond to categorical data
%		as defined in the look_table
%	look_table - nx2 table with columns Numbers and Categories that serves as a lookup table 
%		to convert between the numbers and the original categories.
%
% ----------------------------------------------------------------------------------------------------------------------------
% CheckTAKDepedencies
%
% Function to check required toolboxes for running Topographic Analysis Kit (TAK)
%
% ----------------------------------------------------------------------------------------------------------------------------
% ClassifyKnicks
%
% Function to iterate through a set of bounds (i.e. knickpoints) selected while running 'KsnProfiler'. The function 
%	will display a long profile and chi - elevation plot for individual stream segments and will iterate through each
%	bound point you selected in KsnProfiler. The code expects you to input a number or character (at the command prompt)
%	to categorize the knickpoint higlighted in red. You must be consistent in your choice (i.e. you must either use 
%	numbers for all of the classifications or characters for all the classifications within a given run), mixing numbers 
%	and characters will result in an error at the end of the run. For entering characters, it's recommended you keep these 
%	short strings without spaces (i.e. entries supported into a shapefile a attribute table), e.g. knick or bound  
%
% Required Inputs:
%	DEM - Digital Elevation used as input to KsnProfiler
%	FD - Flow direction used as input to KsnProfiler
%	A - Flow accumulation used as input to KsnProfiler
%	Sc - Stream network output from KsnProfiler
%	ksn_master - cell array of selected channels output from KsnProfiler
%	bnd_list - matrix of bounds (i.e. knickpoints) output from KsnProfiler
%
% Optional Inputs:
%	shape_name ['ksn'] - name for the shapefile to be export, must have no spaces to be a valid name for ArcGIS and should NOT include the '.shp'
%
% ----------------------------------------------------------------------------------------------------------------------------
% CompileBasinStats
%
% Function to take the outputs from 'ProcessRiverBasins' and 'SubDivideBigBasins' and produce a single shapefile showing the outlines of polygons
% 	and with commonly desired attributes from the results of 'ProcessRiverBasins' etc. See below for a full list of fields that the output shapefile
% 	will include. If additional grids were provided to 'ProcessRiverBasins', mean and standard error values for those grids will be auto-populated in
% 	the shapefile and the name of the fields will be the character array provided in the second column of additional grids input. This function also
% 	allows you to input a list of additional fields you wish to include (see Optional Inputs below). If you would rather create a GRIDobj with specified
% 	values, use 'Basin2Raster'.
%
% Required Inputs:
% 		location_of_data_files - full path of folder which contains the mat files from 'ProcessRiverBasins'
%
% Optional Inputs:
%		location_of_subbasins ['SubBasins'] - name of folder that contains subbasins of interest (if you created subbasins using
%			"SubDivideBigBasins"), expected to be within the main Basin folder provided with "location_of_data_files". Note that if you do not provide
%			the correct directory name for the location of the subbasins, subbasin values will not be included in the output regardless of your choice
%			for the "include" parameter.
%		include ['all'] - parameter to specify which basins to include in building the shapfile. The default 'all' will include all basin mat files in the 
%			folder you specify. Providing 'subdivided' will check to see if a given main basin was subdivided using 'SubdivideBigBasins' and then only include 
%			the subdivided versions of that basin (i.e. the original main basin for those subbasins will not be included in the table). Providing 'bigonly'
%			will only include the original basins produced by 'ProcessRiverBasins' even if 'SubDivideBigBasins' was run. If 'SubDivideBigBasins' was never run,
%			result of 'all' and 'bigonly' will be the same.
%		extra_field_values [] - cell array of extra field values you wish to include. The first column in this cell array must be the river basin number
%			(i.e. the identifying number in the third column of the RiverMouth input to ProcessRiverBasins or the number generated for the basin in
%			SubDivideBigBasins). Only one row per river basin number is allowed and ALL river basin numbers in the basins being processed must have a value
%			associated with them. Additional columns are interpreted as the values with which you wish to populate the extra fields. These can either be character
%			arrays or numbers, other values will results in an error. 
%		extra_field_names [] - a 1 x m cell array of field names, as characters (no spaces as this won't work with shapefile attributes), associated with the field values. 
%			These must be in the same order as values are given in extra_field_values. If for example your extra_field_values cell array is 3 columns with the river number, 
%			sample name, and erosion rate then your extra_field_names cell array should include entries for 'sample_name' and 'erosion_rate' in that order. 
%		uncertainty ['se'] - parameter to control which measure of uncertainty is included, expects 'se' for standard error (default), 'std' for standard deviation, or 'both'
%			to include both standard error and deviation.
%		filter_by_category [false] - logical flag to recalculate selected mean values based on filtering by particular categories within a categorical grid (provided to
%			ProcessRiverBasins as 'add_cat_grids'). Requires entries to 'filter_type', 'cat_grid', and 'cat_values'. Will produce filtered values for channel steepness, gradient,
%			and mean  elevation by default along with any additonal grids present (i.e. grids provided with 'add_grids' to ProcessRiverBasins).
%		filter_type ['exclude'] - behavior of filter, if 'filter_by_categories' is set to true. Valid inputs are 'exclude', 'include', or 'mode'. If set to 'exclude', the filtered 
%			means will be calculated excluding any portions of grids have the values of 'cat_values' in the 'cat_grid'. If set to 'include', filtered means will only be calculated 
%			for portions of grids that are within specified categories. If set to 'mode', filtered means will be calculated based on the modal value of the categorical grid by basin,
%			e.g. if the mode of basin 1 is 'grMz' and the mode of basin 2 is 'T', then the filtered mean will be calculated based on nodes that are 'grMz' in basin 1 and are 'T' in 
%			basin 2. The idea behind this filter is if you wish to find characteristic stats for particular categories. If filter type is 'mode' then an entry for 'cat_values' is not
%			required.
%		cat_grid [] - name of categorical grid to use as filter, must be the same as the name provided to ProcessRiverBasins (i.e. third column in the cell array provided to
%			'add_cat_grids').
%		cat_values [] - 1xm cell array of categorical values of interest to use in filter. These must match valid categories in the lookup table as output from CatPoly2GRIDobj
%			(i.e. second colmun in cell array provided to 'add_cat_grids')
%		populate_categories [false] - logical flag to add entries that indicate the percentage of a watershed occupied by each category from a categorical grid, e.g. if you
%			provided an entry for 'add_cat_grids' to ProcessRiverBasins that was a geologic map that had three units, 'Q', 'Mz', and 'Pz' and you set 'populate_categories' 
%			to true there will be field names in the resulting shapefile named 'Q', 'Mz', and 'Pz' and the values stored in those columns will correspond to the percentage 
%			of each basin covered by each unit for each basin. Setting populate_categories to true will not have any effect if no entry was provided to 'add_cat_grids' when
%			running ProcessRiverBasins.
%		means_by_category [] - method to calculate means of various continuous values within by categories. Requires that a categorical grid(s) was input to ProcessRiverBasins.
%			Expects a cell 1 x m cell array where the first entry is the name of the category to use (i.e. name for categorical grid you provided to ProcessRiverBasins) and
%			following entries are names of grids you wish to use to find means by categories, e.g. an example array might be {'geology','ksn','rlf2500','gradient'} if you 
%			were interested in looking for patterns in channel steepness, 2.5 km^2 relief, and gradient as a function of rock type/age. Valid inputs for the grid names are:
%				'ksn' - uses channel steepness map structure with user provided reference concavity
%				'gradient' - uses gradient grid
%				'rlf####' - where #### is the radius you provided to ProcessRiverBasins (requires that 'calc_relief' was set to true when running ProcessRiverBasins
%				'NAME' - where NAME is the name of an additional grid provided with the 'add_grids' option to ProcessRiverBasins
%
% Output:
%		Outputs a table (T) with the following default fields:
%			river_mouth - river mouth number provided to ProcessRiverBasins
%			drainage_area - drainage area of basin in km^2
%			out_x - x coordinate of basin mouth
%			out_y - y coordinate of basin mouth
%			center_x - x coordinate of basin in projected coordinates
%			center_y - y coordinate of basin in projected coordinates
%			outlet_elevation - elevation of pour point in m
%			mean_el - mean elevation of basin in meters
%			max_el - maximum elevation of basin in meters
%			mean_ksn - mean channel steepenss
%			mean_gradient - mean gradient
%		Either standard errors, standard deviations or both will be populated for elevation, ksn, and gradient depending on value of 'uncertainty'
%		Mean and standard error / standard deviation / both values will be populated for any additional grids
%
% ----------------------------------------------------------------------------------------------------------------------------
% ConditionDEM
%
% Wrapper around the variety of methods provided by TopoToolbox for smoothing a stream profile. With the exception of
% 	'quantc_grid' and 'mingrad' these methods will only modify elevations along the stream network provided to the code.
% 	See the relevant parent functions for a more in depth description of the behavior of these individual methods. Produces
%	one figure that compares the long profile of the longest stream within the dataset using the uncondtioned and 
%	conditioned DEM to provide a quick method of evaluating the result. These methods vary in their complexity and processing
%	times so it is recommended you understand your choice. Using the 'mincost' method is a good starting place before
%	exploring some of the more complicated methods. Note that the majority of methods (other than 'mincost' and 'mingrad')
%	require the Optimization Toolbox to run and will also process more quickly if you have the Parallel Processing Toolbox.
%	If you do not have access to the Optimization Toolbox, consider using the compiled version of this function and converting
%	the ascii output to a GRIDobj.
%
% Required Inputs:
%	DEM - Digital Elevation as a GRIDobj, assumes unconditioned DEM (e.g. DEMoc from ProcessRiverBasins)
%	FD - Flow direction as FLOWobj
%	S - Stream network as STREAMobj
%	method - method of conditioning, valid inputs are as follows:
%		'mincost' - uses the 'mincosthydrocon' function, valid optional inputs are 'mc_method' and 'fillp'.
%		'mingrad' - uses the 'imposemin' function, valid optional inputs are 'ming'. Note that providing a large minimum gradient
%			to this code can carve the stream well below the topography.
%		'quantc' - uses the 'quantcarve' function (for STREAMobjs), valid optional inputs are 'tau','ming', and 'split'. Requires the
%			Optimization Toolbox and if 'split' is set to true, requires Parallel Processing Toolbox.
%		'quantc_grid' - uses the 'quantcarve' function for (GRIDobjs), valid optional inputs are 'tau'. Requires the Optimization Toolbox.
%			This is a computationally expensive calculation and because it operates it on the whole grid, it can take a long time and/or
%			fail on large grids. The 'quantc' method which only operates on the stream network is significantly fasters and less prone
%			to failure.
%		'smooth' - uses the 'smooth' function, valid optional inputs are 'sm_method','split','stiffness','stiff_tribs', and 'positive' depending 
%			on inputs to optional parameters may require Optimization Toolbox ('sm_method'='regularization' and 'positive'=true) and Parallel
%			Processing Toolbox ('split'=true).
%		'crs' - uses the 'crs' function, valid optional inputs are 'stiffness', 'tau', 'ming', 'stiff_tribs', 'knicks', and 'split'. Requires
%			Optimization Toolbox.
%		'crslin' - uses the 'crslin' function, valid optional inputs are 'stiffness', 'stiff_tribs', 'ming', 'imposemin', 'attachtomin', 
%			'attachheads', 'discardflats','precisecoords'
%			
% Optional Inputs:
%	mc_method [interp] - method for 'mincost', valid inputs are 'minmax' or 'interp'
%	fillp [0.1] - scalar value between 0 and 1 controlling the ratio of carving to filling for 'mincost'
%	ming [0] - minimum gradient [m/m] in downslope direction, used in 'mingrad','quantc','crs','crslin'
% 	tau [0.5] - quantile for carving, used in 'quantc', 'quantc_grid', 'crs'.
%	split [true] - logical flag to utilized parallel processing to independently process tributaries, used in 'quantc_grid', 'smooth', and 'crs'
%	sm_method ['regularization'] - method for 'smooth', valid inputs are 'regularization' and 'movmean'. 
%	stiffness [10] - scalar positive value for stiffness penalty, used in 'smooth', 'crs', and 'crslin'
%	stiff_tribs [true] - logical flag to relax the stiffness penalty at tributary junctions, used in 'smooth', 'crs', and 'crslin'
% 	knicks [] - nx2 matrix of x and y locations of knickpoints where stiffness penalty should be relaxed, used in 'crs' and 'crslin'
%	imposemin [false] -logical flag to preprocess DEM with imposemin during crslin
%	attachtomin [false] - logical flag to prevent elevations from going below profile minima, used in crslin
%	attachheads [false] - logical flag to fix the channel head elevations, used in crslin
%	discardflats [false] - logical flag to discard flat portions of profiles, used in crslin
%	maxcurvature [] - maximum convex curvature at any vertices along profile, used in crslin
%	precisecoords [] - nx3 matrix with x, y, and z coordinates of points that the smoothed profile must past through, used in crslin
%
% ----------------------------------------------------------------------------------------------------------------------------
% DippingBedFinder
%
% Function to determine the expected location of a planar dipping bed within a landscape based on an input coordinate
%
% Required Inputs:
% 	DEM - DEM GRIDobj
%   xy - 1 x 2 vector with the x and y coordinate (i.e. easting and northing) of the location of interest, if you provide an empty vector
%        you will be given the opportunity to pick a location on the DEM
% 	hght_abv_base - height of the outcrop of interest above the base of the bed of interest (i.e. positin of the outcrop in the section)
%   thickness - thickness of the bed (hght_abv_base must be smaller than total thickness)
% 	strike - strike of bed, report with right hand rule
% 	dip - dip of bed
%
% Output:
% 	Code will produce a figure showing expected location of bed and will produce a binary GRIDobj with expected location of the bed 
%   (1 where the bed should appear, 0 where it should not)
%
% ----------------------------------------------------------------------------------------------------------------------------
% FindBasinKnicks
%
% Function for manually selecting knickpoints within a Basin_Data_File (i.e. result of ProcessRiverBasins). 
% 	Choose knickpoints on Chi-Elevation plot with mouse clicks and press return when you have selected
% 	all the knickpoints for a given stream segment. As you progress through, knickpoints you have already picked 
% 	(i.e. on shared portions of river profiles) will be displayed as red dots. If you're interested in trying out
%	an automated method of finding knickpoints, try 'knickpointfinder' included with TopoToolbox.
%
% Required Inputs:
% 	Basin_Data_File - full file path to a saved result from the ProcessRiverBasins script
% 	plot_result - logical flag to either plot the results (true) or not (false) 
%
% Optional Inputs
% 	theta_ref [0.5] - reference concavity for chi calculation
%	save_mat [true] - logical flag to save output mat file containing the KnickPoints array. The name of the file will 
%		be 'Knicks_NUM.mat' where NUM is the river number. Do not change the file name if you want to plot knickpoints
%		using 'MakeCombinedSwath'.
%	shape_name [] - character string to name output shapefile (without .shp), if no input is provided then
%		no shapefile is output
%
% 
% ----------------------------------------------------------------------------------------------------------------------------
% FindCentroid
%
% Function to find centroid of drainage basin
% Input is DEM GRID object 
% Output is two 1x1 matrices, Cx and Cy with x and y positions of centroid in same coordinates as DEM
%
% ----------------------------------------------------------------------------------------------------------------------------
% FitThreshold
%
% Function to interactively select an appopriate threshold area for a given stream
%	network. Function will have you iterate through a number of single streams 
%	(controlled by the number passed to 'num_streams') extracted from the drainage divide.
%	You can use either chi-elevation or slope-area plots (both will be displayed regardless
%	of choice) to visually select where channels begin. Code will find the average of 
%	the selections and generate a new STREAMobj using this new average threshold area. Code
%	also outputs the lists of selected threshold areas and distance from channel head to divide.
%
% Required Inputs:
%	DEM - Digital Elevation as a GRIDobj, assumes unconditioned DEM (e.g. DEMoc from ProcessRiverBasins 
%		or output from MakeStreams)
%	FD - Flow direction as FLOWobj
%	A - Flow accumulation GRIDobj
%	S - Stream network as STREAMobj (doesn't really matter what the threshold area used to calculate this was originally)
%	num_streams - Number of stream profiles to view and select threshold areas
%	pick_method - Type of plot you wish to choose the threshold area on, valid options are:
%		'chi' - Choose threshold areas on a chi elevation plot
%		'slopearea' - Choose threshold areas on slope-area plot
%
% Optional Inputs;
%	ref_concavity [0.50] - refrence concavity used to generate the chi-elevation plot
%
% Outputs:
%	Sn - New version of the STREAMobj using the mean threshold area to define streams
%	thresh_list - list of chosen threshold areas
%	xd_list - list of chosen distances from channel head to divide
%
% ----------------------------------------------------------------------------------------------------------------------------
% KsnChiBatch
%
% Function to produce channel steepness, chi maps or chi grids for all channels within a DEM
% 
% Reqiured Inputs:
% 	DEM - DEM Grid Object (assumes unconditioned DEM)
% 	FD - FLOW object
% 	A - GRID object of flow accumulations
%	S - STREAM object
% 	product - switch to determine which products to produce
%		'ksn' - ksn map as a shapefile
%		'ksngrid' - ascii file with ksn interpolated at all points in a grid
%		'chimap' - ascii file with chi calculated in channel networks
%		'chigrid' - ascii file with chi calculate at all points in a grid
%		'chi' - results for both chimap and chigrid
%		'all' - ksn, ksngrid, chimap, and chigrids
%
% Optional Inputs:
%	conditioned_DEM [] - option to provide a hydrologically conditioned DEM for use in this function (do not provide a conditoned DEM
%		for the main required DEM input!) which will be used for extracting elevations. See 'ConditionDEM' function for options for making a 
%		hydrological conditioned DEM. If no input is provided the code defaults to using the mincosthydrocon function.
%	file_name_prefix ['batch'] - prefix for outputs, will append the type of output, i.e. 'ksn', 'chimap', etc
% 	segment_length [1000] - length of segments in map units for smoothing ksn values, equivalent to smoothing in Profiler
% 	ref_concavity [0.50] - reference concavity (as a positive value) for calculating ksn
% 	output [false]- switch to either output matlab files to the workspace (true) or to not only save the specified files
%		without any workspace output (false)
%	ksn_method [quick] - switch between method to calculate ksn values, options are 'quick' and 'trib', the 'trib' method takes 3-4 times longer 
%		than the 'quick' method. In most cases, the 'quick' method works well, but if values near tributary junctions are important, then 'trib'
%		may be better as this calculates ksn values for individual channel segments individually
%	outlet_level_method [] - parameter to control how stream network outlet level is adjusted. Options for control of outlet elevation are:
%			'elevation' - extract streams only above a given elevation (provided by the user using the 'min_elevation' parameter) to ensure that base level
%				elevation for all streams is uniform. If the provided elevation is too low (i.e. some outlets of the unaltered stream network are above this
%				elevation) then a warning will be displayed, but the code will still run.
%			'max_out_elevation' - uses the maximum elevation of all stream outlets to extract streams only above this elevation, only valid for options that operate
%				on streamlines only (i.e. will not work with 'ksngrid' or 'chigrid').
%	min_elevation [] - parameter to set minimum elevation for outlet level, required if 'outlet_level_method' is set to 'elevation'
%	complete_networks_only [true] - if true (default) the code will only populate portions of the stream network that are complete. Generally, this
%			option should probably be left as true (i.e. chi will not be accurate if drainage area is not accurate), but this can be overly agressive
%			on certain DEMs and when used in tandem with 'min_elevation', it can be slow to calculate as it requires recalculation of the FLOWobj.
%	interp_value [0.1] - value (between 0 and 1) used for interpolation parameter in mincosthydrocon (not used if user provides a conditioned DEM)
%
% ----------------------------------------------------------------------------------------------------------------------------
% KsnColor
%
% Function to generate colormap for channel steepness data
% 
% Code based off colormap functions provided with TopoToolbox
%
% Colormap originally generated from function 'cbrewer2' using
% 'RdYlGn'
%
%
% ----------------------------------------------------------------------------------------------------------------------------
% KsnProfiler
%
% Function to interactively select channel heads and define segements over which to calculate channel steepness values.
% 	This function is designed to be similar to the operation of Profiler_51, with some improvements. Function will display map
%	with the stream network and expects the user to select a location near a channel head of interest. The user will be then 
%	prompted to confirm that the defined stream is the desired choice. Finally, displays of the chi-z and longitudinal profile 
%	of the selected river will appear and the user is expected to define (with mouse clicks) any obvious segments with different 
%	channel steepness (or concavity) on either the chi-z plot or the stream profile (see 'pick_method' option). When done selecting 
%	press enter/return. The user will be prompted whether they wish to continue picking streams or if they are done. When done 
%	picking streams, the function will output three different products (see below) and produce a shapefile of the selected streams 
%	with ksn, concavity, area, and gradient.
%	
% Required Inputs:
%	DEM - Digital Elevation as a GRIDobj, assumes unconditioned DEM (e.g. DEMoc from ProcessRiverBasins or output from MakeStreams)
%	FD - Flow direction as FLOWobj
%	S - Stream network as STREAMobj
%	A - Flow accumulation GRIDobj
%
%%%%%%%%%%%%%%%%%%
% Optional Inputs:
%
%%% Main Options
% 	input_method ['interactive'] - parameter which controls how streams of interest are supplied:
%		'interactive' - user picks streams of interest by selecting channelheads on a map, this option will also iteratively build a 
%			channel steepness map as the user picks more streams.
%		'all_streams' - will use the supplied STREAMobj and iterate through all channel heads. There is an internal parameter to avoid 
%			selecting streams that are too short to properly fit (mostly relevant if 'junction method' is set to 'check'). The default 
%			value is ~4 * the DEM cellisze, the user can change this value by providing an input for the optional parameter 'min_channel_length',
%			input should be in map units and greater than the default. You can use a code like 'SegmentPicker' to select portions of a STREAMobj
%		'stream_length' - will use supplied STREAMobj and entry to 'min_length_to_extract' to iterate through all streams that are longer than
%			the length provided to 'min_length_to_extract'. There is an internal parameter to avoid selecting streams that are too short to fit
%			(mostly relevant if 'junction method' is set to 'check'). The default value is ~4 * the DEM cellsize, the user can change this value
%			by providing an input for the optional parameter 'min_channel_length', input should be in map units and greater than the default.
%		'channel_heads' - will use a supplied list of coordinates of channel heads to select and iterate through streams of interest. If this
%			option is used, the user must provide an input for the optional 'channel_head_list' parameter.
%	pick_method ['chi'] - choice of how you want to pick stream segments. The diagram within which to pick based on your selection will be 
%			outline in red. Valid inputs are:
%		'chi' - select segments on a chi - z plot (recommended and default)
%		'stream' - select segments on a longitudinal profile
%		'slope_area' - select segments on a slope area plot
%	junction_method ['check'] - choice of how to deal with stream junctions:
%		'check' - after each choice, will check whether downstream portions of the selected stream have already been fit, and if it has,
%			the already fit portion of the stream will not be displayed or refit
%		'ignore' - each stream will be displayed from its head to mouth independent of whether portions of the same stream network have 
%			been fit
%	concavity_method ['ref']- options for concavity:
%		'ref' - uses a reference concavity, the user can specify this value with the reference concavity option (see below)
%		'auto' - function finds a best fit concavity for each selected stream, if used in conjunction with 'junction_method','check'
%			this means that short sections of streams picked will auto fit concavity that may differ from downstream portions of the same
%			streams
%
%%% Input Method Options 
%	min_channel_length [] - minimum channel length for consideration when using the 'all_streams' method of input, provide in map units.
%	channel_head_list [] - m x 2 array of x and y coordinates of channel heads OR the name / location of a point shapefile of channel heads, 
%			one of these is required when using 'channel_heads' method of input, must be in the same coordinate system as the input DEM etc. 
%			The code will attempt to find the nearest channel head to the coordinates you provided, so the closer the provided user coordinates
%			are to channel heads, the more accurate this selection method will be.
%	min_length_to_extract [] - minimum stream length (in map units) to extract streams if 'input_method' is set to 'stream_length'.
%
%%% Redefine Threshold Area Options
%	redefine_threshold [false] - logical flag to initiate an extra step for each stream where you manually define the hillslope-fluvial 
%			transition (this will result in overriding the threshold area you used to generate the supplied STREAMobj, and it will also produce
%			a STREAMobj with a variable threshold area for channel definition). See additional optional input 'rd_pick_method'.
%	rd_pick_method ['slope_area'] - plot to use to choose new threshold area if 'redefine_threshold' is set to true. Valid inputs are 
%			'slopearea' and 'chi'.
%
%%% Stream Network Modification Options
%	complete_networks_only [false] - if true, the code will filter out portions of the stream network that are incomplete prior to choosing
%			streams
%	min_elev [] - minimum elevation below which the code stops extracting channel information (no action if left empty)
%	max_area [] - maximum drainage area above which the code stops extracting channel information (in square map units, no action if left empty)
%
%%% Hydrological Conditioning Options
%	conditioned_DEM [] - option to provide a hydrologically conditioned DEM for use in this function (do not provide a conditoned DEM
%			for the main required DEM input!) which will be used for extracting elevations. See 'ConditionDEM' function for options 
%			for making a hydrological conditioned DEM. If no input is provided the code defaults to using the mincosthydrocon function.
%	interp_value [0.1] - value (between 0 and 1) used for interpolation parameter in mincosthydrocon (not used if user provides a 
%			conditioned DEM). Values closer to 0 tend to 'carve' more, whereas values closer to 1 tend to fill. See info for 
%			'mincosthydrocon'
%
%%% Display Options
%	display_slope_area [false] - logical flag to display slope area plots. Some people love slope area plots (like one of the authors of
%			the supporting paper), some people hate slope area plots (like the other author of the supporting paper), so you can either 
%			not draw them at all (false - default) or include them (true). This will automatically be set to true if you select 'slope_area'
%			as the 'pick_method'.
%	plot_type ['vector'] - expects either 'vector' or 'grid', default is 'vector'. Controls whether all streams are drawn as individual 
%			lines ('vector') or if the stream network is plotted as a grid and downsampled ('grid'). The 'grid' option is much faster on 
%			large datasets, but can result in inaccurate channel head selection. The 'vector' option is easier to see, but can be very 
%			slow to load and interact with on large datasets.	
%
%%% Constants
%	ref_concavity [0.50] - refrence concavity used if 'theta_method' is set to 'ref'
%	smooth_distance [1000] - distance in map units over which to smooth ksn measures when converting to shapefile
%	max_ksn [250] - maximum  ksn used for the color scale, will not effect actual results, for display purposes only
%	threshold_area [1e6] - used to redraw downsampled stream network if 'plot_type' is set to 'grid' 
%
%%% Output Options
%	shape_name ['ksn'] - name for the shapefile to be export, must have no spaces to be a valid name for ArcGIS and should NOT include the '.shp'
%	save_figures [false] - logical flag to either save figures showing ksn fits (true) or to not (false - default)	
%
%%%%%%%%%%	
% Outputs:
%	knl - n x 11 matrix of node list for selected stream segments, columns are x coordinate, y coordinate, drainage area, ksn, negative ksn error,
%		positive ksn error, reference concavity, best fit concavity, mininum threshold area, gradient, and an identifying number. Note that if using the 
%		code in 'concavity_method','auto' mode then the reference concavity and best fit concavity columns will be the same.
%	ksn_master - identical to knl but as a cell array where individual cells are individual selected channels
%	bnd_list - n x 4 matrix of selected bounds for fitting ksn, columns are x coordinate, y coordinate, elevation, and the stream identifying number 
%		(this could be thought of as a list of knickpoints), also output as a seperate shapefile. If x y and z values appear as NaN, this indicates
%		that bounds for this stream were not selected. 
%	Sc - STREAMobj of selected streams
%
% ----------------------------------------------------------------------------------------------------------------------------
% MakeCombinedSwath
%
% Function to plot various additional data onto a swath profile.
%
% Required Inputs:
% 	DEM - DEM Grid Object with which to make topo swath
% 	points - n x 2 matrix containing x,y points for swath, minimum are two points (start and end points).
%		First row contains starting point and proceeds down rows, additional points besides a start and end are
%		treated as bends in the swath. Coordinates for points must be in the same coordinate system as DEM and must
%		lie within the DEM (cannot be coordinates on the very edge of the DEM)
% 	width - width of swath in map units
% 	data_type - the type of additional data you are providing to plot along with the swath, supported inputs are:
%		'points3' - generic point dataset, expects a n x 3 matrix with values of x, y, and z
%		'points4' - generic point dataset, expects a n x 4 matrix with values of x, y, z, and extra value. Dots will 
%					be sclaed by this extra value
%		'eqs' - earthquakes, expects a n x 4 matrix with x, y, depth, and magnitude. Points will be scaled by magnitude 
%					and colored by distance from swath line. Expects depth to be positive.
%		'STREAMobj' - will project portions of selected stream profiles (as points) onto a swath. Expects a STREAMobj 
%					that was generated from the provided DEM.
%		'ksn_chandata' - will plot swath through ksn values, expects a chandata file as output from old Profiler51 code 
%					(just in case you have some sitting around)
%		'ksn_batch' - will plot swath through ksn values, expects the map structure output from 'KsnChiBatch' function 
%					(i.e. run 'KsnChiBatch' with product set to 'ksn' and 'output' set to true, and provide the second 
%					output here as 'data') 
%		'ksn_profiler' - will plot swath through ksn values, expects the 'knl' output from 'KsnProfiler' function
%		'basin_stats' - will plot swath through selected mean basin values as calculated from 'ProcessRiverBasins', 
%					expects output from 'CompileBasinStats' and requires an entry to optional input 'basin_value' and  
%					accepts optional input to 'basin_scale'. Will place point for basin at mean elevation and projected  
%					location of the basin centroid, will color by value provided to 'basin_value' and will optionall scale  
%					the point by the value provided to 'basin_scale'
%		'basin_knicks' - will plot swath through knickpoints as chosen by 'FindBasinKnicks'. For 'data' provide name of folder 
%					(or file path) where to find knickpoint files saved as a result of running 'FindBasinKnicks' on a series of  
%					basins selected from 'ProcessRiverBasins'
%	data - input data, form varies depending on choice of data_type
% 	data_width - width in map units of swath through provided data. Values greater than data_width/2 from the center line 
%					of the toposwath will not be plotted
%
% Optional Inputs:
% 	sample [] - resampling distance along topographic swath in map units, if no input is provided, code will use the cellsize 
%				of the DEM which results in no resampling.
% 	smooth [0] - smoothing distance, width of filter in map units over which to smooth values, default (0) results in no smoothing
%	vex [10] - vertical exaggeration for the topographic swath. Note that because matlabs controls on physical axis dimensions are
%				problematic, the vertical exaggeration controls don't work on plots that have two panels (e.g. 'ksn_batch', 'ksn_profiler',
%				'ksn_chandata', and 'eqs')
%	basin_value [] - required for option 'basin_stats', name (as it appears in the provided table provided to 'data') of the value 
%				you wish to color points by
%	basin_scale [] - optional input for option 'basin_stats', name (as it appears in the provided table provided to 'data') of the 
%				value you wish to scale points by
%	plot_map [true] - logical flag to plot a map displaying the location of the topographic swath and the additional data included 
%				in the swath (red dots) and those not (white dots) based on the provided data_width parameter.
%	cmap ['parula'] - valid name of colormap (e.g. 'jet') or a nx3 colormap array to use to color points.
%	save_figure [false] - logical flag to save the swath figure as a pdf
%
% Outputs:
%	SW - TopoToolbox Swath object, contains various information as a structure. Can plot path and box of swath with plot(SW) and
%		plot version of swath profile with plotdz(SW);
% 	SwathMat - n x 4 matrix containing distance along the swath, min elevation, mean elevation, max elevation
% 	xypoints - n x 2 matrix containing x,y points of each swath sample point, along swath center line
% 	outData - data for plotting the swath through the provided data, distances that area 'NaN' indicate those data do not
%			fall on the swath line provided. Form of output depends on data_type:
%		'points3' - distances, elevation, distance from base line, x coordinate, y coordinate
%		'points4' - distances, elevation, value, distance from base line, x coordinate, y coordinate
%		'eqs' - distances, depth, magnitude, distance from base line, x coordinate, y coordinate
%		'STREAMobj' - distances, elevation, distance from base line, x coordinate, y coordinate
%		'ksn_chandata' - distances, elevation, ksn, distance from base line, x coordinate, y coordinate
%		'ksn_batch' - distances, ksn, distance from base line, x coordinate, y coordinate
%		'ksn_profiler' - distances, ksn, distance from base line, x coordinate, y coordinate
%		'basin_stats' - distances, mean basin elevation, 'basin_value', 'basin_scale' (if provided), distance from base line, 
%						x coordinate, y coordinate
%
% ----------------------------------------------------------------------------------------------------------------------------
% MakeStreams
%
% Function takes a dem and outputs the necessary base datasets for use in other TopoToolbox functions.
% 	Input DEMs with grid resolutions (i.e. cellsizes) that are not whole numbers sometimes cause issues
% 	in companion functions. If the provided DEM has a non-whole number for a cellsize, the code will
% 	warn the user (but not do anything). If you want to fix the cellsize issue, you can either reproject
% 	in a GIS program or you can use this code (with 'resample_grid' set to true) to do it for you.
%
% Required Inputs:
% 	dem - either full path of dem file as either an ascii text file (recommended) or geotiff OR 
%		a GRIDobj of a DEM
% 	threshold_area - minimum accumulation area to define streams in meters squared
%
% Optional Inputs:
%	file_name [] - name for matfile containing the DEM, FD, A, and S and the shapfile of the stream network.
%		If file_name is not provided, the function assumes the user does not wish to save the results to a
%		mat file (results will still appear in the workspace) or shapefile.
%	precip_grid [] - optional input of a GRIDobj of precipitation. If you provide an argument for this, the code will use
%		this to produce a weighted flow accumulation grid.
%	rr_grid [] - optional input of a GRIDobj of runoff ratios. If you provide an argument for this, the code will use
%		this, along with the input to 'precip_grid' to produce a weighted flow accumulation grid.
%	no_data_exp [] - input to define no data conditions. Expects a string that defines a valid equality using
%		the variable DEM OR 'auto'. E.g. if you wish to define that any elevation less that or equal to 0 should 
%		be set to no data, you would provide 'DEM<=0' or if you wanted to set elevations less than 500 and greater  
%		than 1000 ot no data, you would provide 'DEM<500 | DEM>1000'. If the expression is not valid the user will be
%		warned, but the code will continue and ignore this continue. If you provide 'auto' the code will use the log 
%		of the gradient to identify true connected flats and set these to nan. If you want more control on removing flat 
%		ares that are at multiple elevations (e.g. internally drained basins), consider using 'RemoveFlats'. 
%	min_flat_area [1e5] - minimum area (in m^2) for a portion of the DEM to be identified as flat (and set to nan) if 'no_data_exp'
%		is set to 'auto'. If 'no_data_exp' is not called or a valid logical expression is provided, the input to 'min_flat_area'
%		is ignored.
%	resample_grid [false] - flag to resample the grid. If no input is provided for new_cellsize, then the
%		grid will be resampled to the nearest whole number of the native cellsize.
%	new_cellsize [] - value (in map units) for new cellsize.
%
% Outputs:
% 	DEM - GRIDobj of the DEM
% 	FD - FLOWobj from the supplied DEM
% 	A - Flow accumulation grid (GRIDobj)
% 	S - STREAMobj derived from the DEM
%
% ----------------------------------------------------------------------------------------------------------------------------
% MakeTopoSwath
%
% Wrapper around TopoToolbox SWATHobj functionality
%
% Required Inputs:
% 	DEM - DEM Grid Object with which to make topo swath
% 	points - n x 2 matrix containing x,y points for swath, minimum are two points (start and end points).
%		First row contains starting point and proceeds down rows, additional points besides a start and end are
%		treated as bends in the swath. Coordinates for points must be in the same coordinate system as DEM and must
%		lie within the DEM (cannot be coordinates on the very edge of the DEM)
% 	width - width of swath in map units
%
% Optional Inputs:
% 	sample [] - resampling distance along swath in map units, if no input is provided, code will use the cellsize of the DEM 
%		which results in no resampling.
% 	smooth [0] - smoothing distance, width of filter in map units over which to smooth values, default (0) results in no smoothing
%	vex [10] - vertical exaggeration for displaying plot.
% 	plot_figure [false] - logical flag to plot result. 
%	plot_as_points [false] - logical flag to switch plot type to distributions of points
%	plot_as_heatmap [false] - logical flag to switch plot type to a heat map
%
% Outputs:
% 	SW - TopoToolbox Swath object, contains various information as a structure. Can plot path and box of swath with plot(SW) and
%		plot version of swath profile with plotdz(SW);
% 	SwathMat - n x 4 matrix containing distance along the swath, min elevation, mean elevation, max elevation
% 	xypoints - n x 2 matrix containing x,y points of each swath sample point, along swath center line
% 	bends - distances along swath of any bends, 0 if no bends
%
% ----------------------------------------------------------------------------------------------------------------------------
% Mat2Arc
%
% Function converts all valid topotoolbox files contained within a mat file
%	to Arc compatible outputs. Specifically converts any GRIDobjs to
%	ascii files, any STREAMobjs to shapefiles, any FLOWobjs to ArcGIS 
%	flow direction grids saved as an ascii file, and any valid mapstructures
%	to shapefiles.
%
% Input:
%	mat_file - name or path to matfile of interest
%	file_prefix - characters to add to the front of all output files
%
% ----------------------------------------------------------------------------------------------------------------------------
% PlotIndividualBasins
%
% Function takes outputs from 'ProcessRiverBasins' function and makes and saves plots for each basin with stream profiles, chi-z, and slope area
%
% Required Inputs:
% 	location_of_data_files - full path of folder which contains the mat files from 'ProcessRiverBasins'
%
% Optional Inputs:
%	location_of_subbasins ['SubBasins'] - name of folder that contains subbasins of interest (if you created subbasins using
%		"SubDivideBigBasins"), expected to be within the main Basin folder provided with "location_of_data_files"
%
% ----------------------------------------------------------------------------------------------------------------------------
% ProcessRiverBasins
%
% Function takes grid object outputs from MakeStreams script (DEM,FD,A,S), a series of x,y coordinates of river mouths,
% and outputs clipped dem, stream network, variout topographic metrics, and river values (ks, ksn, chi)
%
% Required Inputs:
% 		DEM - GRIDobj of the digital elevation model of your area loaded into the workspace
% 		FD - FLOWobj of the flow direction of your area loaded into the workspace
% 		A - GRID object of flow accumulation of your ara loaded into the workspace
% 		S - STREAMobj of the stream network of your area loaded into the workspace	
% 		river_mouths - locations of river mouths (i.e. pour points) above which you wish to extract basins, can take one of three forms:
%			1) nx3 matrix of river mouths with x, y, and a number identifying the stream/basin of interest (must be same projection as DEM).
%			2) a single value that will be interpreted as an elevation that  the code will use this to autogenerate river mouths at this elevation.
%			3) point shapefile with one numeric user input field (e.g. the default 'ID' field generated by ArcGIS) that will be used as the
%				river mouth ID (must be same projection as DEM).
%		basin_dir - location of folder to store basin files (if specified folder does not exist, code will create it)
%
% Optional Inputs:
%		conditioned_DEM [] - option to provide a hydrologically conditioned DEM for use in this function (do not provide a conditoned DEM
%			for the main required DEM input!) which will be used for extracting elevations. See 'ConditionDEM' function for options for making a 
%			hydrological conditioned DEM. If no input is provided the code defaults to using the mincosthydrocon function.
%		interp_value [0.1] - value (between 0 and 1) used for interpolation parameter in mincosthydrocon (not used if user provides a conditioned DEM)
% 		threshold_area [1e6] - minimum accumulation area to define streams in meters squared
% 		segment_length [1000] - smoothing distance in meters for averaging along ksn, suggested value is 1000 meters
% 		ref_concavity [0.5] - reference concavity for calculating ksn, suggested value is 0.45
%		ksn_method [quick] - switch between method to calculate ksn values, options are 'quick' and 'trib', the 'trib' method takes 3-4 times longer 
%			than the 'quick' method. In most cases, the 'quick' method works well, but if values near tributary junctions are important, then 'trib'
%			may be better as this calculates ksn values for individual channel segments individually
% 		write_arc_files [false] - set value to true to output a ascii's of various grids and a shapefile of the ksn, false to not output arc files
%		add_grids [] - option to provide a cell array of additional grids to clip by selected river basins. The expected input is a nx2 cell array,
%			where the first column is a GRIDobj and the second column is a string identifying what this grid is (so you can remember what these grids
%			are when looking at outputs later, but also used as the name of field values if you use 'Basin2Shape' on the output basins so these should be short 
%			strings with no spaces). The code will perform a check on any input grid to determine if it is the same dimensions and cellsize as the input DEM, if
%			it is not it will use the function 'resample' to transform the input grid. You can control the resampling method used with the 'resample_method' optional
%			parameter (see below), but this method will be applied to all grids you provide, so if you want to use different resampling methods for different grids
%			it is recommnended that you use the 'resample' function on the additional grids before you supply them to this function.
%		add_cat_grids [] - option to provide a cell array of additional grids that are categoricals (e.g. geologic maps) as produced by the 'CatPoly2GRIDobj' function.
%			The expected input is a nx3 cell array where the first column is the GRIDobj, the second column is the look_table, and the third column is a string identifying
%			what this grid is. It is assumed that when preprocessing these grids using 'CatPoly2GRIDobj' you use the same DEM GRIDobj you are inputing to the main function
%			here. These grids are treated differently that those provided to 'add_grids' as it is assumed because they are categorical data that finding mean values is 
%			not useful. Instead these use the 'majority' as the single value but also calculate statistics on the percentages of each clipped watershed occupied by each
%			category.
%		resample_method ['nearest'] - method to use in the resample function on additional grids (if required). Acceptable inputs are 'nearest', 'bilinear', 
%			or 'bicubic'. Method 'nearest' is appropriate if you do not want the resampling to interpolate between values (e.g. if an additinal grid has specific values
%			that correlate to a property like rock type) and either 'bilinear' or 'bicubic' is appropriate if you want smooth variations between nodes. 
%		gradient_method ['arcslope'] - function used to calculate gradient, either 'arcslope' (default) or 'gradient8'. The 'arcslope' function calculates
%			gradient the same way as ArcGIS by fitting a plane to the 8-connected neighborhood and 'gradient8' returns the steepest descent for the same
%			8-connected neighborhood. 'gradient8' will generally return higher values than 'arcslope'.
%		calc_relief [false] - option to calculate local relief. Can provide an array of radii to use with 'relief_radii' option.
%		relief_radii [2500] - a 1d vector (column or row) of radii to use for calculating local relief, values must be in map units. If more than one value is provided
%			the function assumes you wish to calculate relief at all of these radii. Note, the local relief function is slow so providing multiple radii will
%			slow code performance. Saved outputs will be in a m x 2 cell array, with the columns of the cell array corresponding to the GRIDobj and the input radii.
%
% ----------------------------------------------------------------------------------------------------------------------------
% ProjectOntoSwath
%
% Function projects points on a SWATHobj (SW) and finds distance along (ds) and 
%	from center line (db) of the SWATHobj, used in 'MakeCombinedSwath'
% 
% Values for points that do not project onto swath are set to NaN
%
% ----------------------------------------------------------------------------------------------------------------------------
% RemoveFlats
%
% Function takes DEM and attempts a semi-automated routine to remove flat areas with some input from 
% 	the user to select areas considred to be flat. This function sometimes works reliably, but will 
% 	never produce as clean a result as manually clipping out flat areas in gis software (but it's
%	a lot faster!)
%
% Required Inputs:
% 	dem - either full path of dem file as either an ascii text file or geotiff OR 
%			the name of a GRIDobj of a DEM stored in the workspace
%		strength - integer value between 1 and 4 that controls how aggressively the function defines
%			flat areas, specifically related to the size of the neighborhood the function uses to
%			connect ares of similar elevation. A strength of 1 = a 3x3 neighborhood, 2=5x5, 3=7x7, and 
%			4=9x9. If the results of the function do not capture enough of the flat areas in the MASK,
%			increase the strength and rerun. Similarly, if the function erroneously includes areas that
%			are not part of what you consider the flats, try decreasing the strength.
%
% Outputs:
%		DEMn - Version of the DEM with idenitifed flat areas masked out (values set to nan)
%		MASK - Logical GRIDobj, true where area was identified as a flat.
%
% ----------------------------------------------------------------------------------------------------------------------------
% SegmentPicker
%
% Function to select a segment of a stream network from the top of the stream, and plot the long profile 
% 	and chi-Z relationship of that segment,also outputs the extraced portion of the stream network and chi structure 
% 	(out of 'chiplot'). Allows user to iteratively select different parts of the stream network and display. 
% 	Keeps running dataset of all the streams you pick and accept.
%
% Required Inputs:
%	DEM - Digital Elevation as a GRIDobj, assumes unconditioned DEM (e.g. DEMoc from ProcessRiverBasins)
%	FD - Flow direction as FLOWobj
%	A - Flow accumulation GRIDobj
%	S - Stream network as STREAMobj
%	basin_num - basin number from process river basins for output name or other identifying number for the set of streams you will pick
%
% Optional Inputs
%	conditioned_DEM [] - option to provide a hydrologically conditioned DEM for use in this function (do not provide a conditoned DEM
%		for the main required DEM input!) which will be used for extracting elevations. See 'ConditionDEM' function for options for making a 
%		hydrological conditioned DEM. If no input is provided the code defaults to using the mincosthydrocon function.
%	direction ['down'] - expects either 'up' or 'down', default is 'down', if 'up' assumes individual selections are points above
%		which you wish to extract and view stream profiles (i.e. a pour point), if 'down' assumes individual
%		selections are channel heads if specific streams you wish to extract and view stream profiles. 
%	method ['new_picks'] - expects either 'new_picks' or 'prev_picks', default is 'new_picks' if no input is provided. If 'prev_picks' is
%			 given, the user must also supply an input for the 'picks' input (see below)
%	plot_style ['refresh'] - expects either 'refresh' or 'keep', default is 'refresh' if no input is provided. If 'refresh' is given, the plots reset
%			after each new stream pick, but if 'keep' is given, all selected streams remain on both the map (as thick red lines) and the
%			chi-z/longitudinal profile/slope-area plots.
%	plot_type ['vector'] - expects either 'vector' or 'grid', default is 'vector'. Controls whether all streams are drawn as individual lines ('vector') or if
%		   the stream network is plotted as a grid and downsampled ('grid'). The 'grid' option is much faster with large datasets, 
%			but can result in inaccurate choices. The 'vector' option is easier to see, but can be very slow to load and interact with.
%	complete_networks_only [false] - if true, the code will filter out portions of the stream network that are incomplete prior to choosing
%			streams
%	picks - expects a m x 3 matrix with columns as an identifying number, x coordinates, and y coordinates OR the name of a point shapefile with a single value column of 
%			identifying numbers. Will interpret this input as a list of channel heads if 'direction' is 'down' and a list of channel outlets if 'direction' is 'up'.
%	ref_concavity [0.50] - reference concavity for calculating Chi-Z, default is 0.50
%	min_elev [] - minimum elevation below which the code stops extracting channel information, only used if 'direction'
%			   is 'down'
%	max_area [] - maximum drainage area above which the code stops extracting channel information, only used if 'direction'
%			   is 'down'
%	recalc [false] - only valid if either min_elev or max_area are specified. If recalc is false (default) then extraction of 
%			 streams stops downstream of the condition specified in either min_elev or max_area, but chi is not recalculated 
%			 and distances will remain tied to the original stream (i.e. distances from the outlet will be relative to the outlet
%			 of the stream if it continued to the edge of the DEM, not where it stops extracting the stream profile). If recalc is true, 
%	     	 then chi and distance are recalculated (i.e. the outlet as determined by the min_elev or max_area condition
%			 will have a chi value of zero and a distance from mouth value of zero).
%	threshold_area [1e6] - used to redraw downsampled stream network if 'plot_type' is set to 'grid'
%	interp_value [0.1] - value (between 0 and 1) used for interpolation parameter in mincosthydrocon (not used if user provides a conditioned DEM)
%
% Outputs:
% 	Sc - STREAMobj containing all the stream segments chosen.
%
%	Saves an output called 'PickedSegements_*.mat' with the provided basin number containing these results:
%		StreamSgmnts - Cell array of selected stream segments as STREAMobj
%		ChiSgmnts - Cell array of selected chi structures 
%		Sc - Single STREAMobj containing all the streams chosen.
%		and if 'down' is selected:
%			Heads - nx3 matrix of channel heads you picked with x cooord, y coordinate, and pick number as the columns
%		and if 'up' is selected:
%			Outlets - nx3 matrix of outlets you picked with x cooord, y coordinate, and pick number as the columns (valid input to 'ProcessRiverBasins'
%			as 'river_mouths' parameter)
%
% ----------------------------------------------------------------------------------------------------------------------------
% SegmentPlotter
%
% Function to plot all of the chi-Z relationships from a series of picked segments of river networks that result
% from the 'SegmentPicker' function. 
%
% Required Input:
% 	basin_nums - row or column vector of basin numbers used for the SegmentPicker you wish to plot together.
% 
% Optionl Input:
%	separate [false] - logical flag to plot all segments as separate figures 
%	subset [] - list of specific river numbers (i.e. the first column of either the 'Heads' or the 'Outlets' variable) that you wish to include
%				in the plot. Only valid if you have only provided a single basin number for 'basin_nums'.
%	label [false]  - logical flag to either label individual streams with the river number (true) or not label them (false, default). If 'separate' flag
%					is true then the input for label is ignored as the stream number will be in the title of the plots
%	names [] - option to add an identifying name for streams in a basin when 'label' is set to true. If you have input one basins data, the argument provided
%				to 'names' should be a string, if you've provided several basins, it should be a cell array of strings in the order of the basins
%
% ----------------------------------------------------------------------------------------------------------------------------
% SegmentProjector
%
% Function to interactively select segments of a channel profile you wish to project (e.g. projecting a portion of the profile with a different ksn).
%	You can use the 'SegmentPicker' function to interactively choose channels to provide to the StreamProjector function. If the STREAMobj has more than 
%	one channel head, this code will iterate through all channel heads (i.e. make sure you're only providing it stream you want to project, not an entire
%	network!). It calculates and will display 95% confidence bounds on this fit.
%	
% Required Inputs:
%	DEM - Digital Elevation as a GRIDobj, assumes unconditioned DEM (e.g. DEMoc from ProcessRiverBasins)
%	FD - Flow direction as FLOWobj
%	A - Flow accumulation GRIDobj
%	S - Streams you wish to project saved as a STREAMobj
%
% Optional Inputs:
%	conditioned_DEM [] - option to provide a hydrologically conditioned DEM for use in this function (do not provide a conditoned DEM
%		for the main required DEM input!) which will be used for extracting elevations. See 'ConditionDEM' function for options for making a 
%		hydrological conditioned DEM. If no input is provided the code defaults to using the mincosthydrocon function.
%	theta_method ['ref']- options for concavity
%		'ref' - uses a reference concavity, the user can specify this value with the reference concavity option (see below)
%		'auto' - function finds a best fit concavity for the provided stream
%	pick_method ['chi'] - choice of how you want to pick the stream segment to be projected:
%		'chi' - select segments on a chi - z plot
%		'stream' - select segments on a longitudinal profile
%	ref_concavity [0.50] - refrence concavity used if 'theta_method' is set to 'auto'
%	refit_streams [false] - option to recalculate chi based on the concavity of the picked segment (true), useful if you want to try to precisely 
%		match the shape of the picked segment of the profile. Only used if 'theta_method' is set to 'auto'
%	save_figures [false] - option to save (if set to true) figures at the end of the projection process
%	interp_value [0.1] - value (between 0 and 1) used for interpolation parameter in mincosthydrocon (not used if user provides a conditioned DEM)
%
% Output:
%	Produces a 2 x n cell array with a column for each stream segment provided (or channel head if a network is provided). The first row is the x-y
%	coordinate of the channel head for that stream. The second row is an array containing the following information about the segment of interest:
%	x coordinate, y coordinate, distance from mouth, drainage area, chi, concavity, true elevation, projected elevation, upper bound of projected elevation,
%	and lower bound of projected elevation. This output is also saved in a mat file called 'ProjectedSegments.mat'.
%
% ----------------------------------------------------------------------------------------------------------------------------
% SubDivideBigBasins
%
% Function takes outputs from 'ProcessRiverBasins' function and subdvides any basin with a drainage area above a specified size and
% outputs clipped dem, stream network, variout topographic metrics, and river values (ks, ksn, chi)
%
% Required Inputs:
% 		basin_dir - full path of folder which contains the mat files from 'ProcessRiverBasins'
% 		max_basin_size - size above which drainage basins will be subdivided in square kilometers
%		divide_method - method for subdividing basins, options are ('confluences' and 'up_confluences' is NOT recommended large datasets):
%			'order' - use the outlets of streams of a given order that the user can specify with the optional 's_order' parameter 
%			'confluences' - use the locations of confluences (WILL PRODUCE A LOT OF SUB BASINS!). There is an internal parameter to remove
%				extremely short streams that would otherwise result in the code erroring out.
%			'up_confluences' - use locations just upstream of confluences (WILL PRODUCE A LOT OF SUB BASINS!). There is an internal parameter
%				to remove extremely short streams that otherwise result in the code erroring out.
%			'filtered_confluences' - use locations of confluences if drainage basin above confluence is of a specified size that the user
%				can specify with the optional 'min_basin_size'  
%			'p_filtered_confluences' - similar to filtered confluences, but the user defines a percentage of the main basin area
%				with the optional 'min_basin_size'
%			'trunk' - uses the tributary junctions with the trunk stream within the main basin as pour points for subdivided basins. There is
%				an internal parameter to remove extremely short streams that would otherwise result in the code erroring out.
%			'filtered_trunk' - same as 'trunk' but will only include basins that are greater than the min_basin_size
%			'p_filtered_trunk' - same as 'filtered_trunk' but 'min_basin_size' is interpreted as a percentage of the main basin area
%
% Optional Inputs:
%		SBFiles_Dir ['SubBasins'] - name of folder (within the main Basins folder) to store the subbasin files. Subbasin files are now stored in
%			a separate folder to aid in the creation of different sets of subbasins based on different requirements. 
%		recursive [true] - logical flag to ensure no that no subbasins in the outputs exceed the 'max_basin_size' provided. If 'divide_method' is 
%			one of the trunk varieties the code will continue redefining trunks and further split subbasins until no extracted basins are greater
%			than the 'max_basin_size'. If the 'divide_method' is one of the confluence varities, subbasins greater than 'max_basin_size' will simply
%			no be included in the output. The 'recursive' check is not implemented for the 'order' method.
% 		threshold_area [1e6] - minimum accumulation area to define streams in meters squared
% 		segment_length [1000] - smoothing distance in meters for averaging along ksn, suggested value is 1000 meters
% 		ref_concavity [0.5] - reference concavity for calculating ksn
% 		write_arc_files [false] - set value to true to output a ascii's of various grids and a shapefile of the ksn, false to not output arc files
%		s_order [3] - stream order for defining stream outlets for subdividing if 'divide_method' is 'order' (lower number will result in more sub-basins)
%		min_basin_size [10] - minimum basin size for auto-selecting sub basins. If 'divide_method' is 'filtered_confluences' this value is
%			interpreted as a minimum drainage area in km^2. If 'divide_method' is 'p_filtered_confluences', this value is interpreted as
%			the percentage of the input basin drainage area to use as a minimum drainage area, enter a value between 0 and 100 in this case.
%		no_nested [false] - logical flag that when used in conjunction with either 'filtered_confluences' or 'p_filtered_confluences' will only extract
%			subbasins if they are the lowest order basin that meets the drainage area requirements (this is to avoid producing nested basins)
























