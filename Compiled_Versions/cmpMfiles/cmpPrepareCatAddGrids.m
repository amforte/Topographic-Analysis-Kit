function cmpPrepareAddCatGrids(wdir,out_file_name,MakeStreamsMat,varargin)
	% Function to prepare categorical grids for input to 'cmpProcessRiverBasins'
	%
	% Inputs:
	%	out_file_name - name for the mat file to be produced, do not include the '.mat'
	%	MakeStreamMat - full path of the output produced by 'cmpMakeStreams' and that you will be
	%		using as an input to 'cmpProcessRiverBasins'
	%
	%	Additional inputs must be given in groups of threes, and be in the order: 
	%		1) name of the shapefile containing the field you want to convert into a categorical grid
	%		2) the field name within the shapefile you want to convert into a categorical grid
	%		3) a reference name for the produced grid
	%
   	% Examples if running for the command line, minus OS specific way of calling main TAK function:
    %   PrepareAddCatGrids /path/to/wdir AddCatGrids Topo.mat geo_polygons.shp RTYPE rock_type
    %   PrepareAddCatGrids /path/to/wdir AddCatGrids Topo.mat geo_polygons.shp RTYPE rock_type geo_polygons.shp UNIT unit_name  
    %
    %
	if isdeployed
		if ~isempty(varargin)
			varargin=varargin{1};
		end
	end

	if ~isempty(regexp(MakeStreamsMat,regexptranslate('wildcard','*.mat')));
		load(fullfile(wdir,MakeStreamsMat),'DEM');
	else
		error('Must provide a valid .mat file to "MakeStreamsMat"');
	end

	if isempty(varargin) | mod(numel(varargin),3)~=0
		error('You must provide three entries for each additional grid you wish to prepare');
	end

	num_grids=numel(varargin)/3;

	% Generate indices
	ix1=[1:3:numel(varargin)];
	ix2=[2:3:numel(varargin)];
	ix3=[3:3:numel(varargin)];

	% Generate empty cell
	ACG=cell(num_grids,3);

	for ii=1:num_grids

		if isempty(regexp(varargin{ix1},regexptranslate('wildcard','*.shp')))
			error('First entry not recognized as a valid shapefile');
		end

		if ~ischar(varargin{ix2})
			error('Second entry not recognized as a valid character entry');
		end

		[ACG{ii,1},ACG{ii,2}]=CatPoly2GRIDobj(DEM,fullfile(wdir,varargin{ix1}),varargin{ix2});

		if ischar(varargin{ix3})
			ACG{ii,3}=varargin{ix3};
		else
			error('Third entry not recognized as a valid character entry');
		end
	end

	fn=fullfile(wdir,[out_file_name '.mat']);
	save(fn,'ACG');
end