function Mat2Arc(mat_file,file_prefix,varargin)
	%
	% Usage:
	%	Mat2Arc(mat_file,file_prefix);
	%
	% Description:
	% 	Function converts all valid topotoolbox files contained within a mat file
	%	to Arc compatible outputs. Specifically converts any GRIDobjs to
	%	ascii files, any STREAMobjs to shapefiles, any FLOWobjs to ArcGIS 
	%	flow direction grids saved as an ascii file, and any valid mapstructures
	%	to shapefiles.
	%
	% Required Inputs:
	%	mat_file - name or path to matfile of interest
	%	file_prefix - characters to add to the front of all output files
	%
	% Optional Inputs:
	%	raster_type ['ascii'] - option to specify the format of the raster export, valid
	%		inputs are 'tif' or 'ascii'
	%
	% Example:
	%	Mat file 'Data.mat' containing:
	%		DEM - GRIDobj
	%		FD - FLOWobj
	%		S - STREAMobj
	%		MS_KSN - Mapstructure with ksn data
	%		
	%	Mat2Arc('Data.mat','Example'); would produce:
	%		Example_DEM.txt - ascii of the DEM
	%		Example_FD.txt - ascii of the Arc flow direction grid
	%		Example_S.shp - shapefile of the streams in STREAMobj
	%		Example_MS_KSN.shp - shapefile of the stream data with ksn
	%
	%	Mat2Arc('Data.mat','Example','raster_type','tif'); would produce:
	%		Example_DEM.tif - GeoTiff of the DEM
	%		Example_FD.tif - GeoTiff of the Arc flow direction grid
	%		Example_S.shp - shapefile of the streams in STREAMobj
	%		Example_MS_KSN.shp - shapefile of the stream data with ksn
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 06/18/18 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'Mat2Arc';
	addRequired(p,'mat_file',@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.mat'))));
	addRequired(p,'file_prefix',@(x) ischar(x));

	addParameter(p,'raster_type','ascii',@(x) ischar(validatestring(x,{'ascii','tif'})));
	addParameter(p,'out_dir',[],@(x) isdir(x));

	parse(p,mat_file,file_prefix,varargin{:});
	mat_file=p.Results.mat_file;
	file_prefix=p.Results.file_prefix;

	raster_type=p.Results.raster_type;
	out_dir=p.Results.out_dir;

	if isempty(out_dir)
		out_dir=pwd;
	end

	file_prefix=[out_dir filesep file_prefix];

	d=whos('-file',mat_file);

	num_vars=numel(d);

	for ii=1:num_vars

		classOI=d(ii,1).class;
		varNM=d(ii,1).name;

		if strcmp(classOI,'GRIDobj')
			varOI=load(mat_file,varNM);
			switch raster_type
			case 'ascii'
				out_name=[file_prefix '_' varNM '.txt'];
				GRIDobj2ascii(varOI.(varNM),out_name);
			case 'tif'
				out_name=[file_prefix '_' varNM '.tif'];
				GRIDobj2geotiff(varOI.(varNM),out_name);
			end
		elseif strcmp(classOI,'STREAMobj')
			varOI=load(mat_file,varNM);
			out_name=[file_prefix '_' varNM '.shp'];
			MS=STREAMobj2mapstruct(varOI.(varNM));
			shapewrite(MS,out_name);
		elseif strcmp(classOI,'FLOWobj')
			varOI=load(mat_file,varNM);
			G=FLOWobj2GRIDobj(varOI.(varNM));
			switch raster_type
			case 'ascii'
				out_name=[file_prefix '_' varNM '.txt'];
				GRIDobj2ascii(G,out_name);
			case 'tif'
				out_name=[file_prefix '_' varNM '.tif'];
				GRIDobj2geotiff(G,out_name);
			end
		elseif strcmp(classOI,'struct')
			varOI=load(mat_file,varNM);
			st=varOI.(varNM);
			fn=fieldnames(st);
			if any(strcmp(fn,'Geometry')) & any(strcmp(fn,'X')) & any(strcmp(fn,'Y'))
				out_name=[file_prefix '_' varNM '.shp'];
				shapewrite(st,out_name);
			end
		end
	end
end
		