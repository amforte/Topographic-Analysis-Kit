function [varargout]=KSN_Chi_Batch(DEM,FD,A,S,product,varargin)
	% Function to produce channel steepness, chi maps or chi grids for all channels within a DEM
	% 
	% Reqiured Inputs:
	% 	DEM - DEM Grid Object (assumes unconditioned DEM)
	% 	FD - FLOW object
	% 	A - GRID object of flow accumulations
	%	S - STREAM object
	% 	product - switch to determine which products to produce
	%		'ksn' - ksn map as a shapefile
	%		'chimap' - ascii file with chi calculated in channel networks
	%		'chigrid' - ascii file with chi calculate at all points in a grid
	%		'chi' - results for both chimap and chigrid
	%		'all' - ksn, chimap, and chigrids
	%
	% Optional Inputs:
	% 	smooth_distance [1000] - smoothing distance in map units for smoothing ksn values, equivalent to smoothing in Profiler
	% 	min_ksn [1] - minimum ksn for calculating minimum gradients at different drainage areas 
	% 	ref_concavity [0.45] - reference concavity (as a positive value) for calculating ksn
	% 	output [false]- switch to either output matlab files to the workspace (true) or to not only save the specified files
	%		without any workspace output (false)
	%	adjust_base_level [false] - flag to adjust base level of stream network. Chi values are sensitive to choice of baselevel so it is recommended 
	%		that you (1) either clip the DEM prior to running 'MakeStreams' so that extracted streams will have what you deem to be appropriate base levels
	%		 or (2) use one of the provided options with this flag to set base levels (see 'base_level_method').
	%	base_level_method [] - parameter to control how stream network base level is adjusted, only valid if 'adjust_base_level' is true. Options for control
	%		 of base level are:
	%			'elevation' - extract streams only above a given elevation (provided by the user using the 'min_elevation' parameter) to ensure that base level
	%				elevation for all streams is uniform. If the provided elevation is too low (i.e. some outlets of the unaltered stream network are above this
	%				elevation) then a warning will be displayed, but the code will still run.
	%			'drain_area' - extract streams only below a given maximum drainage area (provided by the user using the 'max_drainage_area' parameter) to ensure
	%				that the outlets of all extracted streams have the same drainage areas. If the provided maximum drainage area is too large (i.e. some outlets
	%				have drainage areas smaller than this maximum) then a warning will be displayed, but the code will still run.
	%			'max_out_elevation' - uses the maximum elevation of all stream outlets to extract streams only above this elevation
	%			'min_out_drain_area' - uses the minimum drainage area of all stream outlets to extract streams only below this drainage area
	%	min_elevation [] - parameter to set minimum elevation for base level, required if 'base_level_method' is set to 'elevation'
	%	max_drainage_area [] - parameter to set maximum drainage area for base level, required if 'base_level_method' is set to 'drain_area'
	%
	% Notes:
	%	Please be aware that the production of the chigrid can be time consuming, so be patient...
	%
	% Example:
	%	KSN_Chi_Batch(DEM,FD,A,S,'ksn');
	%	[KSN,ChiMap,ChiGrid]=KSN_Chi_Batch(DEM,FD,A,S,'output',true,'ref_concavity',0.55);
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Last Revised Winter 2017 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Parse Inputs
	p = inputParser;         
	p.FunctionName = 'SetBaseLevel';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'FD', @(x) isa(x,'FLOWobj'));
	addRequired(p,'A', @(x) isa(x,'GRIDobj'));
	addRequired(p,'S',@(x) isa(x,'STREAMobj'));
	addRequired(p,'product',@(x) ischar(validatestring(x,{'ksn','chimap','chigrid','all'})));

	addParamValue(p,'smooth_distance',1000,@(x) isscalar(x) && isnumeric(x));
	addParamValue(p,'min_ksn',0.1,@(x) isscalar(x) && isnumeric(x));
	addParamValue(p,'ref_concavity',0.45,@(x) isscalar(x) && isnumeric(x));
	addParamValue(p,'output',false,@(x) isscalar(x) && islogical(x));
	addParamValue(p,'adjust_base_level',false,@(x) isscalar(x) && islogical(x));
	addParamValue(p,'base_level_method',[],@(x) ischar(validatestring(x,{'elevation','drain_area','max_out_elevation','min_out_drain_area'})));
	addParamValue(p,'min_elevation',[],@(x) isnumeric(x));
	addParamValue(p,'max_drainage_area',[],@(x) isnumeric(x));

	parse(p,DEM,FD,A,S,product,varargin{:});
	DEM=p.Results.DEM;
	FD=p.Results.FD;
	A=p.Results.A;
	S=p.Results.S;
	product=p.Results.product;

	smooth_distance=p.Results.smooth_distance;
	min_ksn=p.Results.min_ksn;
	ref_concavity=p.Results.ref_concavity;
	output=p.Results.output;
	abl=p.Results.adjust_base_level;
	blm=p.Results.base_level_method;
	me=p.Results.min_elevation;
	ma=p.Results.max_drainage_area;

	% Check that cut off values have been provided if base level
	if abl & isempty(blm)
		error('Base level adjustment set to true, but no base level adjustment method was provided');
	elseif abl & strcmp(blm,'elevation') & isempty(me)
		error('Selected base level adjust method "elevation" requires that you provide an input for parameter "min_elevation"');
	elseif abl & strcmp(blm,'drain_area') & isempty(ma)
		error('Selected base level adjust method "drain_area" requires that you provide an input for parameter "max_drainage_area"');
	end

	switch product
	case 'ksn'

		disp('Calculating channel steepness')
		[G,Gmin] = CompGrad(DEM,FD,A,min_ksn,ref_concavity);

		if abl & strcmp(blm,'elevation')
			[S]=SetBaseLevel(DEM,FD,A,S,blm,'min_elevation',me);
		elseif abl & strcmp(blm,'drain_area');
			[S]=SetBaseLevel(DEM,FD,A,S,blm,'max_drainage_area',ma);
		elseif abl
			[S]=SetBaseLevel(DEM,FD,A,S,blm);
		end	

		ksn=G./(A.*(A.cellsize^2)).^(-ref_concavity);
		KSN=STREAMobj2mapstruct(S,'seglength',smooth_distance,'attributes',...
			{'ksn' ksn @mean 'uparea' (A.*(A.cellsize^2)) @mean 'gradient' G @mean 'gmin' Gmin @mode});

		disp('Writing ARC files')
		shapewrite(KSN,'ksn.shp');

		switch output
		case true
			KSNG=GRIDobj(DEM);
			for ii=1:numel(KSN)
				ix=coord2ind(DEM,KSN(ii).X,KSN(ii).Y);
				KSNG.Z(ix)=KSN(ii).ksn;
			end
			varargout{1}=KSNG;
		end


	case 'chimap'

		if abl & strcmp(blm,'elevation')
			[S]=SetBaseLevel(DEM,FD,A,S,blm,'min_elevation',me);
		elseif abl & strcmp(blm,'drain_area');
			[S]=SetBaseLevel(DEM,FD,A,S,blm,'max_drainage_area',ma);
		elseif abl
			[S]=SetBaseLevel(DEM,FD,A,S,blm);
		end	

	    disp('Calculating chi map');
		[ChiMap]=MakeChiMap(DEM,FD,A,S,ref_concavity);

		disp('Writing ARC files')
		GRIDobj2ascii(ChiMap,'chimap.txt');

		switch output
		case true
			varargout{1}=ChiMap;
		end

	case 'chigrid'

	    disp('Calculating chi grid');
		if abl & strcmp(blm,'elevation')
		    [ChiGrid]=MakeChiGrid(DEM,FD,A,ref_concavity,abl,blm,me);
		elseif abl & strcmp(blm,'drain_area');
		    [ChiGrid]=MakeChiGrid(DEM,FD,A,ref_concavity,abl,blm,ma);
		elseif abl
		    [ChiGrid]=MakeChiGrid(DEM,FD,A,ref_concavity,abl,blm);
		else
			[ChiGrid]=MakeChiGrid(DEM,FD,A,ref_concavity,abl);
		end	

	    disp('Writing ARC files')
		GRIDobj2ascii(ChiGrid,'chigrid.txt');

		switch output
		case true
			varargout{1}=ChiGrid;
		end

	case 'chi'

		if abl & strcmp(blm,'elevation')
			[S]=SetBaseLevel(DEM,FD,A,S,blm,'min_elevation',me);
		elseif abl & strcmp(blm,'drain_area');
			[S]=SetBaseLevel(DEM,FD,A,S,blm,'max_drainage_area',ma);
		elseif abl
			[S]=SetBaseLevel(DEM,FD,A,S,blm);
		end

	    disp('Calculating chi map');
		[ChiMap]=MakeChiMap(DEM,FD,A,S,ref_concavity);

	    disp('Calculating chi grid');
		if abl & strcmp(blm,'elevation')
		    [ChiGrid]=MakeChiGrid(DEM,FD,A,ref_concavity,abl,blm,me);
		elseif abl & strcmp(blm,'drain_area');
		    [ChiGrid]=MakeChiGrid(DEM,FD,A,ref_concavity,abl,blm,ma);
		elseif abl
		    [ChiGrid]=MakeChiGrid(DEM,FD,A,ref_concavity,abl,blm);
		else
			[ChiGrid]=MakeChiGrid(DEM,FD,A,ref_concavity,abl);
		end	

		disp('Writing ARC files')
		GRIDobj2ascii(ChiMap,'chimap.txt');
		GRIDobj2ascii(ChiGrid,'chigrid.txt');

		switch output
		case true
			varargout{1}=ChiMap;
			varargout{2}=ChiGrid;
		end

	case 'all'

		disp('Calculating channel steepness')
		[G,Gmin] = CompGrad(DEM,FD,A,min_ksn,ref_concavity);

		if abl & strcmp(blm,'elevation')
			[S]=SetBaseLevel(DEM,FD,A,S,blm,'min_elevation',me);
		elseif abl & strcmp(blm,'drain_area');
			[S]=SetBaseLevel(DEM,FD,A,S,blm,'max_drainage_area',ma);
		elseif abl
			[S]=SetBaseLevel(DEM,FD,A,S,blm);
		end	

		ksn=G./(A.*(A.cellsize^2)).^(-ref_concavity);
		KSN=STREAMobj2mapstruct(S,'seglength',smooth_distance,'attributes',...
			{'ksn' ksn @mean 'uparea' (A.*(A.cellsize^2)) @mean 'gradient' G @mean 'gmin' Gmin @mode});

	    disp('Calculating chi map');
		[ChiMap]=MakeChiMap(DEM,FD,A,S,ref_concavity);

	    disp('Calculating chi grid');
		if abl & strcmp(blm,'elevation')
		    [ChiGrid]=MakeChiGrid(DEM,FD,A,ref_concavity,abl,blm,me);
		elseif abl & strcmp(blm,'drain_area');
		    [ChiGrid]=MakeChiGrid(DEM,FD,A,ref_concavity,abl,blm,ma);
		elseif abl
		    [ChiGrid]=MakeChiGrid(DEM,FD,A,ref_concavity,abl,blm);
		else
			[ChiGrid]=MakeChiGrid(DEM,FD,A,ref_concavity,abl);
		end	


		disp('Writing ARC files')
		shapewrite(KSN,'ksn.shp');
		GRIDobj2ascii(ChiMap,'chimap.txt');
		GRIDobj2ascii(ChiGrid,'chigrid.txt');

		switch output
		case true
			KSNG=GRIDobj(DEM);
			for ii=1:numel(KSN)
				ix=coord2ind(DEM,KSN(ii).X,KSN(ii).Y);
				KSNG.Z(ix)=KSN(ii).ksn;
			end
			varargout{1}=KSNG;
			varargout{2}=ChiMap;
			varargout{3}=ChiGrid;
		end
	end

% Main Function End
end

function [ChiOBJ]=MakeChiMap(DEM,FD,A,S,ref_concavity);

	DA=A.*(A.cellsize^2);
	C=chitransform(S,DA,'mn',ref_concavity,'a0',1);

	% Make Empty GRIDobj
	ChiOBJ=GRIDobj(DEM);
	ChiOBJ.Z(ChiOBJ.Z==0)=NaN;

	% Populate Grid
	ChiOBJ.Z(S.IXgrid)=C;

end

function [ChiOBJ]=MakeChiGrid(DEM,FD,A,ref_concavity,abl,varargin)

	% Generate dense stream network
	S=STREAMobj(FD,'minarea',0);

	if abl & numel(varargin)==1
		blm=varargin{1};
	elseif abl & numel(varargin)==2
		blm=varargin{1};
		con=varargin{2};
	end

	if abl & strcmp(blm,'elevation')
		[S]=SetBaseLevel(DEM,FD,A,S,blm,'min_elevation',con);
	elseif abl & strcmp(blm,'drain_area');
		[S]=SetBaseLevel(DEM,FD,A,S,blm,'max_drainage_area',con);
	elseif abl
		[S]=SetBaseLevel(DEM,FD,A,S,blm);
	end

	DA=A.*(A.cellsize^2);
	C=chitransform(S,DA,'mn',ref_concavity,'a0',1);

	% Make Empty GRIDobj
	ChiOBJ=GRIDobj(A);
	ChiOBJ.Z(ChiOBJ.Z==0)=NaN;

	% Populate Grid
	ChiOBJ.Z(S.IXgrid)=C;
end

function [G,GminIX] = CompGrad(DEM,FD,A,min_ksn,ref_concavity)

	% Convert Area in Pixels to Drainage Area in Sq Meters
	DA=A.*(A.cellsize^2);

	% Find Max Drainage Area Exponent
	maxDA=max(max(DA));
	ex=ceil(log10(maxDA));

	% Build Range
	bins=logspace(4,ex,ex-4+1); % Min area is 1e4
	num_bins=numel(bins);

	% Build Masks and Gradients
	for ii=1:num_bins
		da=bins(ii);
		% Build Indices
		if ii==1
			IX=DA<=da;
		else
			IX=DA>bins(ii-1) & DA<=da;
		end

		% Calculate and Impose Min Gradient
		minG=min_ksn*(da^-ref_concavity);
		Gtemp=gradient8(imposemin(FD,DEM,minG));
		GIX=Gtemp.*IX;

		% Layer Composite
		if ii==1
			G=GIX;
		else
			G=G+GIX;
		end
	end

	% Build map of where minimum gradient has been imposed
	Gn=gradient8(DEM);
	mIX=Gn.Z~=G.Z;
	[X,Y]=getcoordinates(DEM);
	GminIX=GRIDobj(X,Y,double(mIX));
end
