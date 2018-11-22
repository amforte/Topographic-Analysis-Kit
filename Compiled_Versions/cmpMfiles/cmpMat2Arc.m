function cmpMat2Arc(wdir,MatFile,file_prefix,varargin)
	% Description:
	% 	Function converts all valid topotoolbox files contained within a mat file
	%	to Arc compatible outputs. Specifically converts any GRIDobjs to
	%	ascii files, any STREAMobjs to shapefiles, any FLOWobjs to ArcGIS 
	%	flow direction grids saved as an ascii file, and any valid mapstructures
	%	to shapefiles.
	%
	% Required Inputs:
	%	wdir - full path of working directory
	%	MatFile - full path to matfile of interest
	%	file_prefix - characters to add to the front of all output files
	%
	% Optional Inputs:
	%	raster_type ['ascii'] - option to specify the format of the raster export, valid
	%		inputs are 'tif' or 'ascii'
	%
   	% Examples if running for the command line, minus OS specific way of calling main TAK function:
    %   Mat2Arc /path/to/wdir Topo.mat outputs
    %	Mat2Arc /path/to/wdir Topo.mat outputs raster_type tif
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
	p.FunctionName = 'cmpMat2Arc';
	addRequired(p,'wdir',@(x) ischar(x));
	addRequired(p,'MatFile',@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.mat'))));
	addRequired(p,'file_prefix',@(x) ischar(x));

	addParameter(p,'raster_type','ascii',@(x) ischar(validatestring(x,{'ascii','tif'})));

	parse(p,wdir,MatFile,file_prefix,varargin{:});
	wdir=p.Results.wdir;
	MatFile=p.Results.MatFile;
	file_prefix=p.Results.file_prefix;

	raster_type=p.Results.raster_type;


	MatFile=fullfile(wdir,MatFile);
	d=whos('-file',MatFile);

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
			out_name=fullfile(wdir,[file_prefix '_' varNM '.shp']);
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
				out_name=fullfile(wdir,[file_prefix '_' varNM '.shp']);
				shapewrite(st,out_name);
			end
		end
	end
end
		