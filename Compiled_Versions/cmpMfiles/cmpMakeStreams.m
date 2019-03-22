function cmpMakeStreams(wdir,dem,threshold_area,file_name,varargin)
	% Function takes a dem and outputs the necessary base datasets for use in other TopoToolbox functions.
	% 	Input DEMs with grid resolutions (i.e. cellsizes) that are not whole numbers sometimes cause issues
	% 	in companion functions. If the provided DEM has a non-whole number for a cellsize, the code will
	% 	warn the user (but not do anything). If you want to fix the cellsize issue, you can either reproject
	% 	in a GIS program or you can use this code (with 'resample_grid' set to true) to do it for you.
	%
	% Required Inputs:
	%	wdir - full path of working directory
	% 	dem - name of dem file as either an ascii text file (recommended) or geotiff 
	% 	threshold_area - minimum accumulation area to define streams in meters squared
	%	file_name [] - name for matfile containing the DEM, FD, A, and S (for use in all the other compiled versions
	%		of the TAK codes) and the shapfile of the stream network, do not include a file type suffix, this will be
	%		added by the function.
	%
	% Optional Inputs:
	%	precip_grid [] - optional input of name of a precipitation raster (ascii or geotiff). If you provide an argument for this, 
	%		the code will use this to produce a weighted flow accumulation grid.
	%	rr_grid [] - optional input of nameof a runoff ratio raster (ascii or geotiff). If you provide an argument for this, the code 
	%		will use this, along with the input to 'precip_grid' to produce a weighted flow accumulation grid.
	%	no_data_exp [] - input to define no data conditions. Expects a string that defines a valid equality using
	%		the variable DEM OR 'auto'. E.g. if you wish to define that any elevation less that or equal to 0 should 
	%		be set to no data, you would provide 'DEM<=0' or if you wanted to set elevations less than 500 and greater  
	%		than 1000 ot no data, you would provide 'DEM<500 | DEM>1000'. If the expression is not valid the user will be
	%		warned, but the code will continue and ignore this continue. If you provide 'auto' the code will use the log 
	%		of the gradient to identify true connected flats and set these to nan. If you want more control on removing flat 
	%		ares that are at multiple elevations (e.g. internally drained basins), consider using 'RemoveFlats'. 
	%	min_flat_area [1e8] - minimum area (in m^2) for a portion of the DEM to be identified as flat (and set to nan) if 'no_data_exp'
	%		is set to 'auto'. If 'no_data_exp' is not called or a valid logical expression is provided, the input to 'min_flat_area'
	%		is ignored.
	%	resample_grid [false] - flag to resample the grid. If no input is provided for new_cellsize, then the
	%		grid will be resampled to the nearest whole number of the native cellsize.
	%	new_cellsize [] - value (in map units) for new cellsize.
	%
	% Outputs:
	% 	Saved matfile containing (for use with other 'cmp*' codes):
	%		DEM - GRIDobj of the DEM
	% 		FD - FLOWobj from the supplied DEM
	% 		A - Flow accumulation grid (GRIDobj)
	% 		S - STREAMobj derived from the DEM
	%	Shapefile of streams
	%
    % Examples if running for the command line, minus OS specific way of calling main TAK function:
    %   MakeStreams /path/to/wdir dem.txt 1e6 Topo
    %   MakeStreams /path/to/wdir dem.txt 1e6 Topo no_data_exp auto min_flat_area 1e6
   	%   MakeStreams /path/to/wdir dem.txt 1e6 Topo no_data_exp DEM<=0
	% 
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 07/02/18 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if isdeployed
		if ~isempty(varargin)
			varargin=varargin{1};
		end
	end

	% Parse Inputs
	p = inputParser;         
	p.FunctionName = 'cmpMakeStreams';
	addRequired(p,'wdir',@(x) ischar(x));
	addRequired(p,'dem',@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.txt'))) || ~isempty(regexp(x,regexptranslate('wildcard','*.tif'))));
	addRequired(p,'threshold_area', @(x) isscalar(x));
	addRequired(p,'file_name',@(x) ischar(x));
	
	addParameter(p,'no_data_exp',[],@(x) ischar(x));
	addParameter(p,'min_flat_area',1e8,@(x) isnumeric(x) && isscalar(x));
	addParameter(p,'resample_grid',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'new_cellsize',[],@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'precip_grid',[],@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.txt'))) || ~isempty(regexp(x,regexptranslate('wildcard','*.tif'))));
	addParameter(p,'rr_grid',[],@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.txt'))) || ~isempty(regexp(x,regexptranslate('wildcard','*.tif'))));

	parse(p,wdir,dem,threshold_area,file_name,varargin{:});
	wdir=p.Results.wdir;
	dem=p.Results.dem;
	threshold_area=p.Results.threshold_area;
	file_name=p.Results.file_name;

	no_data_exp=p.Results.no_data_exp;
	min_flat_area=p.Results.min_flat_area;
	resample_grid=p.Results.resample_grid;
	new_cellsize=p.Results.new_cellsize;
	precip_grid=p.Results.precip_grid;
	rr_grid=p.Results.rr_grid;

	% Process DEM
	disp('Loading and processing DEM')
	DEM=GRIDobj(fullfile(wdir,dem));

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

	% Optional cleaning step depending on user input
	disp('Cleaning Up DEM')
	if ~isempty(no_data_exp) & ~strcmp(no_data_exp,'auto')
		try 
			IDX=eval(no_data_exp);
			DEM.Z(IDX.Z)=nan;
			% Remove any borders of nans
			DEM=crop(DEM);
		catch
			warning('Provided "no_data_exp" was not a valid expression, proceeding without this no data condition');
		end
	elseif strcmp(no_data_exp,'auto')
		[DEM]=AutoFlat(DEM,min_flat_area);
	else
		% Remove any borders of nans
		DEM=crop(DEM);
	end



    fileNameBase=fullfile(wdir,file_name);
    MatFileName=[fileNameBase '.mat'];
    ShpFileName=[fileNameBase '.shp'];
	save(MatFileName,'DEM','-v7.3');


	disp('Calculating Flow Direction')
	FD=FLOWobj(DEM,'preprocess','carve','verbose',true);

	save(MatFileName,'FD','-append','-v7.3');

	disp('Calculating Flow Accumulation')
	if isempty(precip_grid)
		A=flowacc(FD);
	elseif ~isempty(precip_grid) & isempty(rr_grid)
		precipLoc=fullfile(wdir,precip_grid);
		PG=GRIDobj(precipLoc);
		if ~validatealignment(PG,DEM);
			PG=resample(PG,DEM,'nearest');
		end

		A=flowacc(FD,PG);
	elseif isempty(precip_grid) & ~isempty(rr_grid)
		precip_grid=GRIDobj(DEM);
		precip_grid.Z=ones(DEM.size);

		rrLoc=fullfile(wdir,rr_grid);
		RRG=GRIDobj(rr_grid);

		if ~validatealignment(RRG,DEM);
			RRG=resample(RRG,DEM,'nearest');
		end

		A=flowacc(FD,precip_grid,RRG);			
	elseif ~isempty(precip_grid) & ~isempty(rr_grid)
		precipLoc=fullfile(wdir,precip_grid);
		PG=GRIDobj(precipLoc);
		rrLoc=fullfile(wdir,rr_grid);
		RRG=GRIDobj(rrLoc);
		if ~validatealignment(PG,DEM);
			PG=resample(PG,DEM,'nearest');
		end

		if ~validatealignment(RRG,DEM);
			RRG=resample(RRG,DEM,'nearest');
		end

		A=flowacc(FD,PG,RRG);
	end

	disp('Extracting total stream network')
	S=STREAMobj(FD,'unit','mapunits','minarea',min_area);

	save(MatFileName,'A','S','-append','-v7.3');
	MS=STREAMobj2mapstruct(S);
	shapewrite(MS,ShpFileName);

end

function [DEMn] = AutoFlat(DEM,min_area)

    num_pix=round(min_area/(DEM.cellsize^2));

    LG=log10(gradient8(DEM));
    BW=isnan(LG.Z) | isinf(LG.Z);
    CC=bwconncomp(BW);
    FLATS=GRIDobj(DEM,'logical');

    for ii=1:numel(CC.PixelIdxList)
        if numel(CC.PixelIdxList{ii})>=num_pix
            idx=CC.PixelIdxList{ii};
            FLATS.Z(idx)=true;
        end
    end

    DEMn=DEM;
    DEMn.Z(FLATS.Z)=nan;

end
