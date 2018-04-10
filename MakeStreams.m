function [DEM,FD,A,S]=MakeStreams(dem,threshold_area,varargin)
	% Function takes a dem and outputs the necessary base datasets for use in other TopoToolbox functions.
	% Input DEMs with grid resolutions (i.e. cellsizes) that are not whole numbers sometimes cause issues
	% in companion functions. If the provided DEM has a non-whole number for a cellsize, the code will not
	% warn the user (but not do anything). If you want to fix the cellsize issue, you can either reproject
	% in a GIS program or you can use this code (with 'resample_grid' set to true) to do it for you.
	%
	% Required Inputs:
	% 	dem - either full path of dem file as either an ascii text file (recommended) or geotiff OR 
	%		a GRIDobj of a DEM
	% 	threshold_area - minimum accumulation area to define streams in meters squared
	%
	% Optional Inputs:
	%	file_name [] - name for matfile containing the DEM, FD, A, and S and the shapfile of the stream network.
	%		If file_name is not provided, the function assumes the user does not wish to save the results to a
	%		mat file (results will still appear in the workspace).
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
	%
	% examples [DEM,FD,A,S]=MakeStreams('/Users/forte/topo/dem.tif',1e6);
	%		   [DEM,FD,A,S]=MakeStreams('/Users/forte/topo/dem.tif',1e6,'file_name','AreaFiles');
	%		   [DEM,FD,A,S]=MakeStreams(DEMgrid,1e6,'resample_grid',true); %Where DEMgrid is a GRIDobj	
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Last Revised Spring 2018 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Parse Inputs
	p = inputParser;         
	p.FunctionName = 'MakeStreams';
	addRequired(p,'dem',@(x) isa(x,'GRIDobj') | ischar(x));
	addRequired(p,'threshold_area', @(x) isscalar(x));

	addParamValue(p,'file_name',[],@(x) ischar(x));
	addParamValue(p,'resample_grid',false,@(x) isscalar(x) && islogical(x));
	addParamValue(p,'new_cellsize',[],@(x) isscalar(x) && isnumeric(x));

	parse(p,dem,threshold_area,varargin{:});
	dem=p.Results.dem;
	threshold_area=p.Results.threshold_area;

	file_name=p.Results.file_name;
	resample_grid=p.Results.resample_grid;
	new_cellsize=p.Results.new_cellsize;


	% Check for filename
	if isempty(file_name)
		save_output=false;
	else
		save_output=true;
	end

	% Check type of input
	if isa(dem,'GRIDobj');
		DEM=dem;
	elseif ischar(dem);
		disp('Loading and processing DEM')
		DEM=GRIDobj(dem);
	else
		error('Input for dem not recognized as either a GRIDobj or character')
	end

	% Resample grid if flag is thrown
	if resample_grid & isempty(new_cellsize)
		disp('Resampling DEM - May take some time, please be patient')
		DEM=resample(DEM,ceil(DEM.cellsize),'bicubic');
	elseif resample_grid & ~isempty(new_cellsize)
		disp('Resampling DEM - May take some time, please be patient')
		DEM=resample(DEM,new_cellsize,'bicubic');		
	end

	% Check resolution of DEM
	if mod(DEM.cellsize,1)~=0 & ~resample_grid
		warning('Grid Cellsize is not a whole number, this may cause problems in some TopoToolbox functions, consider using resample_grid option')
	end

	% Clippin gout NaNs, max and min elev set arbitrarily large but below (or above) the (-)32,768 internal NaN value.
	disp('Cleaning Up DEM')
	max_elev=10000;
	min_elev=-1000;
	IDX=DEM<max_elev & DEM>min_elev;
	DEM=crop(DEM,IDX,nan);

	if save_output
        fileNameBase=file_name;
        MatFileName=[fileNameBase '.mat'];
        ShpFileName=[fileNameBase '.shp'];
		save(MatFileName,'DEM','-v7.3');
	end

	disp('Calculating Flow Direction')
	FD=FLOWobj(DEM,'preprocess','carve','verbose',true);

	if save_output
		save(MatFileName,'FD','-append');
	end

	disp('Calculating Flow Accumulation')
	A=flowacc(FD);

	disp('Extracting total stream network')
	DEM_res=DEM.cellsize;
	min_area=floor(threshold_area/(DEM_res*DEM_res));
	isstream=A>min_area;
	S=STREAMobj(FD,isstream);

	if save_output
		save(MatFileName,'A','S','-append');
		MS=STREAMobj2mapstruct(S);
		shapewrite(MS,ShpFileName);
	end

end
