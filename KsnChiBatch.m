function [varargout]=KsnChiBatch(DEM,FD,A,S,product,varargin)
	%
	% Usage:
	%	KsnChiBatch(DEM,FD,A,S,product);
	%	KsnChiBatch(DEM,FD,A,S,product,'name',value...);
	%	[outputs]=KsnChiBatch(DEM,FD,A,S,product,'output',true,...);
	%
	% Description:
	% 	Function to produce channel steepness, chi maps or chi grids for all channels within a DEM
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
	%	smooth_distance [1000] - distance in map units over which to smooth ksn measures when converting to shapefile
	% 	ref_concavity [0.50] - reference concavity (as a positive value) for calculating ksn
	% 	output [false]- switch to either output matlab files to the workspace (true) or to not only save the specified files
	%		without any workspace output (false). The number of outputs will depend on the 'product' input.
	%			'ksn' - [KSNG,ksn_ms] where KSNG is a GRIDobj with ksn values along the stream network and ksn_ms is the mapstructure suitable for creating 
	%					a shapefile
	%			'ksngrid' - [KSNgrid] where KSNgrid is a GRIDobj of interpolated ksn values
	%			'chimap' - [ChiMap] where ChiMap is a GRIDobj with chi values along the stream network
	%			'chigrid' - [ChiGrid] where ChiGrid is a GRIDobj with chi values across the entire grid
	%			'chi' - [ChiMap,ChiGrid]
	%			'all' - [KSNG,ksn_ms,KSNGrid,ChiMap,ChiGrid]
	%	ksn_method [quick] - switch between method to calculate ksn values, options are 'quick', 'trunk', or 'trib', the 'trib' method takes 3-4 times longer 
	%		than the 'quick' method. In most cases, the 'quick' method works well, but if values near tributary junctions are important, then 'trib'
	%		may be better as this calculates ksn values for individual channel segments individually. The 'trunk' option calculates steepness values
	%		of large streams independently (streams considered as trunks are controlled by the stream order value supplied to 'min_order'). The 'trunk' option
	%		may be of use if you notice anomaoloulsy high channel steepness values on main trunk streams that can result because of the way values are reach
	%		averaged.
	%	min_order [4] - minimum stream order for a stream to be considered a trunk stream, only used if 'ksn_method' is set to 'trunk'
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
	%	radius [5000] - radius of circular, moving area over which to average ksn values to produced an interpolated ksn grid if product is set to 'ksngrid' or 'all'
	%
	% Notes:
	%	Please be aware that the production of the ksngrid and/or chigrid can be time consuming, so be patient...
	%
	% Example:
	%	KsnChiBatch(DEM,FD,A,S,'ksn');
	%	[KSN,ChiMap,ChiGrid]=KsnChiBatch(DEM,FD,A,S,'output',true,'theta_ref',0.55);
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 06/18/18 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Parse Inputs
	p = inputParser;         
	p.FunctionName = 'KsnChiBatch';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'FD', @(x) isa(x,'FLOWobj'));
	addRequired(p,'A', @(x) isa(x,'GRIDobj'));
	addRequired(p,'S',@(x) isa(x,'STREAMobj'));
	addRequired(p,'product',@(x) ischar(validatestring(x,{'ksn','ksngrid','chimap','chigrid','chi','all'})));

	addParameter(p,'file_name_prefix','batch',@(x) ischar(x));
	addParameter(p,'smooth_distance',1000,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'ref_concavity',0.50,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'output',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'ksn_method','quick',@(x) ischar(validatestring(x,{'quick','trunk','trib'})));
	addParameter(p,'min_order',4,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'outlet_level_method',[],@(x) ischar(validatestring(x,{'elevation','max_out_elevation'})));
	addParameter(p,'min_elevation',[],@(x) isnumeric(x));
	addParameter(p,'conditioned_DEM',[],@(x) isa(x,'GRIDobj'));
	addParameter(p,'interp_value',0.1,@(x) isnumeric(x) && x>=0 && x<=1);
	addParameter(p,'complete_networks_only',true,@(x) isscalar(x) && islogical(x));
	addParameter(p,'radius',5000,@(x) isscalar(x) && isnumeric(x));

	parse(p,DEM,FD,A,S,product,varargin{:});
	DEM=p.Results.DEM;
	FD=p.Results.FD;
	A=p.Results.A;
	S=p.Results.S;
	product=p.Results.product;

	file_name_prefix=p.Results.file_name_prefix;
	segment_length=p.Results.smooth_distance;
	theta_ref=p.Results.ref_concavity;
	output=p.Results.output;
	ksn_method=p.Results.ksn_method;
	min_order=p.Results.min_order;
	blm=p.Results.outlet_level_method;
	me=p.Results.min_elevation;
	iv=p.Results.interp_value;
	DEMc=p.Results.conditioned_DEM;
	cno=p.Results.complete_networks_only;
	radius=p.Results.radius;

	% Check that cut off values have been provided if base level
	if strcmp(blm,'max_out_elevation') & (strcmp(product,'chigrid') | strcmp(product,'ksngrid') | strcmp(product,'all'))
		warning('"max_out_elevation" is not a valid choice for a continuous gridded product, ignoring this input')
	elseif strcmp(blm,'elevation') & isempty(me)
		error('Selected outlet level adjust method "elevation" requires that you provide an input for parameter "min_elevation"');
	end

	switch product
	case 'ksn'

		disp('Calculating channel steepness')
		
		if strcmp(blm,'elevation')
			[S]=DTSetOutlet(DEM,FD,A,S,blm,'complete_networks_only',cno,'min_elevation',me);
		elseif strcmp(blm,'max_out_elevation')
			[S]=DTSetOutlet(DEM,FD,A,S,blm,'complete_networks_only',cno);
		elseif isempty(blm) & cno
			[S]=DTSetOutlet(DEM,FD,A,S,'complete_only','complete_networks_only',cno);
		end

		% Hydrologically condition dem
		if isempty(DEMc)
			zc=mincosthydrocon(S,DEM,'interp',iv);
			DEMc=GRIDobj(DEM);
			DEMc.Z(DEMc.Z==0)=NaN;
			DEMc.Z(S.IXgrid)=zc;
		end
		
		switch ksn_method
		case 'quick'
			[ksn_ms]=KSN_Quick(DEM,DEMc,A,S,theta_ref,segment_length);
		case 'trunk'
			[ksn_ms]=KSN_Trunk(DEM,DEMc,A,S,theta_ref,segment_length,min_order);
		case 'trib'
			[ksn_ms]=KSN_Trib(DEM,DEMc,FD,A,S,theta_ref,segment_length);
		end

		disp('Writing ARC files')
		out_file=[file_name_prefix '_ksn.shp'];
		shapewrite(ksn_ms,out_file);

		switch output
		case true
			KSNG=GRIDobj(DEM);
			KSNG.Z(:,:)=NaN;
			for ii=1:numel(ksn_ms)
				ix=coord2ind(DEM,ksn_ms(ii).X,ksn_ms(ii).Y);
				KSNG.Z(ix)=ksn_ms(ii).ksn;
			end
			varargout{1}=KSNG;
			varargout{2}=ksn_ms;
		end

	case 'ksngrid'

		if strcmp(blm,'elevation')
			IDX=DEM<me;
			DEM.Z(IDX.Z)=NaN;
		else
			IDX=GRIDobj(DEM,'logical');
		end

		if isempty(DEMc)
			zc=mincosthydrocon(S,DEM,'interp',iv);
			DEMc=GRIDobj(DEM);
			DEMc.Z(DEMc.Z==0)=NaN;
			DEMc.Z(S.IXgrid)=zc;
		else
			DEMc.Z(IDX.Z)=NaN;
		end

		if strcmp(blm,'elevation')
			[S]=DTSetOutlet(DEM,FD,A,S,blm,'complete_networks_only',cno,'min_elevation',me);
		elseif strcmp(blm,'max_out_elevation')
			[S]=DTSetOutlet(DEM,FD,A,S,blm,'complete_networks_only',cno);
		elseif isempty(blm) & cno
			[S]=DTSetOutlet(DEM,FD,A,S,'complete_only','complete_networks_only',cno);
		end
		
		switch ksn_method
		case 'quick'
			[ksn_ms]=KSN_Quick(DEM,DEMc,A,S,theta_ref,segment_length);
		case 'trunk'
			[ksn_ms]=KSN_Trunk(DEM,DEMc,A,S,theta_ref,segment_length,min_order);
		case 'trib'
			[ksn_ms]=KSN_Trib(DEM,DEMc,FD,A,S,theta_ref,segment_length);
		end

		[KSNGrid]=KsnAvg(DEM,ksn_ms,radius);

		disp('Writing ARC files')
		out_file_ksng=[file_name_prefix '_ksngrid.txt'];
		GRIDobj2ascii(KSNGrid,out_file_ksng);

		switch output
		case true
			varargout{1}=KSNGrid;
		end

	case 'chimap'

		if strcmp(blm,'elevation')
			[S]=DTSetOutlet(DEM,FD,A,S,blm,'complete_networks_only',cno,'min_elevation',me);
		elseif strcmp(blm,'max_out_elevation')
			[S]=DTSetOutlet(DEM,FD,A,S,blm,'complete_networks_only',cno);
		elseif isempty(blm) & cno
			[S]=DTSetOutlet(DEM,FD,A,S,'complete_only','complete_networks_only',cno);
		end

	    disp('Calculating chi map');
		[ChiMap]=MakeChiMap(DEM,FD,A,S,theta_ref);

		disp('Writing ARC files')
		out_file_cm=[file_name_prefix '_chimap.txt'];
		GRIDobj2ascii(ChiMap,out_file_cm);

		switch output
		case true
			varargout{1}=ChiMap;
		end

	case 'chigrid'

	    disp('Calculating chi grid');
		[ChiGrid]=MakeChiGrid(DEM,FD,'theta_ref',theta_ref,'complete_networks_only',cno,'min_elevation',me);

	    disp('Writing ARC files')
		out_file=[file_name_prefix '_chigrid.txt'];
		GRIDobj2ascii(ChiGrid,out_file);

		switch output
		case true
			varargout{1}=ChiGrid;
		end

	case 'chi'

		if strcmp(blm,'elevation')
			[S]=DTSetOutlet(DEM,FD,A,S,blm,'complete_networks_only',cno,'min_elevation',me);
		elseif strcmp(blm,'max_out_elevation')
			[S]=DTSetOutlet(DEM,FD,A,S,blm,'complete_networks_only',cno);
		elseif isempty(blm) & cno
			[S]=DTSetOutlet(DEM,FD,A,S,'complete_only','complete_networks_only',cno);
		end

	    disp('Calculating chi map');
		[ChiMap]=MakeChiMap(DEM,FD,A,S,theta_ref);

	    disp('Calculating chi grid');
		[ChiGrid]=MakeChiGrid(DEM,FD,'theta_ref',theta_ref,'complete_networks_only',cno,'min_elevation',me);

		disp('Writing ARC files')
		out_file_cg=[file_name_prefix '_chigrid.txt'];
		GRIDobj2ascii(ChiGrid,out_file_cg);
		out_file_cm=[file_name_prefix '_chimap.txt'];
		GRIDobj2ascii(ChiMap,out_file_cm);

		switch output
		case true
			varargout{1}=ChiMap;
			varargout{2}=ChiGrid;
		end

	case 'all'

		if strcmp(blm,'elevation')
			IDX=DEM<me;
			DEM.Z(IDX.Z)=NaN;
		else
			IDX=GRIDobj(DEM,'logical');
		end

		disp('Calculating channel steepness')
		if strcmp(blm,'elevation')
			[S]=DTSetOutlet(DEM,FD,A,S,blm,'complete_networks_only',cno,'min_elevation',me);
		elseif strcmp(blm,'max_out_elevation')
			[S]=DTSetOutlet(DEM,FD,A,S,blm,'complete_networks_only',cno);
		elseif isempty(blm) & cno
			[S]=DTSetOutlet(DEM,FD,A,S,'complete_only','complete_networks_only',cno);
		end

		if isempty(DEMc)
			zc=mincosthydrocon(S,DEM,'interp',iv);
			DEMc=GRIDobj(DEM);
			DEMc.Z(DEMc.Z==0)=NaN;
			DEMc.Z(S.IXgrid)=zc;
		else
			DEMc.Z(IDX.Z)=NaN;
		end
		
		switch ksn_method
		case 'quick'
			[ksn_ms]=KSN_Quick(DEM,DEMc,A,S,theta_ref,segment_length);
		case 'trunk'
			[ksn_ms]=KSN_Trunk(DEM,DEMc,A,S,theta_ref,segment_length,min_order);
		case 'trib'
			[ksn_ms]=KSN_Trib(DEM,DEMc,FD,A,S,theta_ref,segment_length);
		end

		disp('Calculating interpolated ksn grid')
		[KSNGrid]=KsnAvg(DEM,ksn_ms,radius);

	    disp('Calculating chi map');
		[ChiMap]=MakeChiMap(DEM,FD,A,S,theta_ref);

	    disp('Calculating chi grid');
		[ChiGrid]=MakeChiGrid(DEM,FD,'theta_ref',theta_ref,'complete_networks_only',cno,'min_elevation',me);

		disp('Writing ARC files')
		out_file_ksn=[file_name_prefix '_ksn.shp'];
		shapewrite(ksn_ms,out_file_ksn);
		out_file_ksng=[file_name_prefix '_ksngrid.txt'];
		GRIDobj2ascii(KSNGrid,out_file_ksng);
		out_file_cg=[file_name_prefix '_chigrid.txt'];
		GRIDobj2ascii(ChiGrid,out_file_cg);
		out_file_cm=[file_name_prefix '_chimap.txt'];
		GRIDobj2ascii(ChiMap,out_file_cm);

		switch output
		case true
			KSNG=GRIDobj(DEM);
			KSNG.Z(:,:)=NaN;
			for ii=1:numel(ksn_ms)
				ix=coord2ind(DEM,ksn_ms(ii).X,ksn_ms(ii).Y);
				KSNG.Z(ix)=ksn_ms(ii).ksn;
			end
			varargout{1}=KSNG;
			varargout{2}=ksn_ms;
			varargout{3}=KSNGrid;
			varargout{4}=ChiMap;
			varargout{5}=ChiGrid;
		end
	end

% Main Function End
end

function [SC]=DTSetOutlet(DEM,FD,A,S,method,varargin)
	% Clone of Divide Tools 'SetOutlet' function.

	% Parse Inputs
	p = inputParser;         
	p.FunctionName = 'DTSetOutlet';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'FD', @(x) isa(x,'FLOWobj'));
	addRequired(p,'A', @(x) isa(x,'GRIDobj'));
	addRequired(p,'S',@(x) isa(x,'STREAMobj'));
	addRequired(p,'method',@(x) ischar(validatestring(x,{'elevation','max_out_elevation','complete_only'})));

	addParameter(p,'complete_networks_only',true,@(x) islogical(x) & isscalar(x));
	addParameter(p,'min_elevation',[],@(x) isnumeric(x));


	parse(p,DEM,FD,A,S,method,varargin{:});
	DEM=p.Results.DEM;
	FD=p.Results.FD;
	A=p.Results.A;
	S=p.Results.S;
	method=p.Results.method;

	cno=p.Results.complete_networks_only;
	me=p.Results.min_elevation;

	if ~cno & strcmp(method,'complete_networks_only')
		error('Cannot set method to complete_only and set complete_network_only to false');
	end

	if cno & ~strcmp(method,'complete_networks_only')
		S=removeedgeeffects(S,FD,DEM);
	end

	%% Initiate graphical picker if no values for either min drainage area or min elevation are provided
	if strcmp(method,'elevation') & isempty(me)
		f1=figure(1);
		set(f1,'Units','normalized','Position',[0.1 0.1 0.75 0.75]);
		hold on
		imageschs(DEM,DEM);
		title('Zoom to desired view and press enter')
		pause()
		title('Pick point on DEM to select base level elevation');
		[x,y]=ginput(1);
		hold off
		close(f1)
		el_ix=coord2ind(DEM,x,y);
		me=DEM.Z(el_ix);
		disp(['Selected elevation is : ' num2str(me) ' m']);
	end

	%% Main switch between methods
	switch method
	case 'elevation'
		st_el=getnal(S,DEM);
		idx=st_el>=me;
		IX=S.IXgrid(idx);
		W=GRIDobj(DEM);
		W.Z(IX)=1;
		W.Z=logical(W.Z);
		SC=STREAMobj(FD,W);
		% Check to see if all outlets meet the condition
		coix=streampoi(SC,'outlets','ix');
		coel=DEM.Z(coix);
		max_coel=max(coel);
		if sum(coel>me)~=0 & ~cno
			warning(['One or more stream outlets are above the provided elevation, maximum outlet elevation is ' num2str(max_coel)]);
		elseif sum(coel>me)~=0 & cno
			[xo,yo]=getoutline(DEM,true);
			% Control for incosistent output of getoutline
			sz=size(xo);
			if sz(1)==1 & sz(2)>1
				[oxy]=[xo' yo'];
			elseif sz(2)==1 & sz(1)>1
				[oxy]=[xo yo];
			end
			[coxy]=streampoi(SC,'outlets','xy');
			idx2=coel>me & ismember(coxy,oxy,'rows'); % Find streams with outlets greater than min elevation AND along boundary of DEM
			coix(idx2)=[];
			W=GRIDobj(DEM);
			W.Z(coix)=1;
			W.Z=logical(W.Z);
			SC=modify(SC,'upstreamto',W);			
		end
	case 'max_out_elevation'
		coix=streampoi(S,'outlets','ix');
		coel=DEM.Z(coix);
		max_coel=max(coel);
		st_el=getnal(S,DEM);
		idx=st_el>=max_coel;
		IX=S.IXgrid(idx);
		W=GRIDobj(DEM);
		W.Z(IX)=1;
		W.Z=logical(W.Z);
		SC=STREAMobj(FD,W);
	case 'complete_only'
		SC=removeedgeeffects(S,FD,DEM);
	end
end

function [ChiOBJ]=MakeChiMap(DEM,FD,A,S,theta_ref);

	C=chitransform(S,A,'mn',theta_ref,'a0',1);

	% Make Empty GRIDobj
	ChiOBJ=GRIDobj(DEM);
	ChiOBJ.Z(ChiOBJ.Z==0)=NaN;

	% Populate Grid
	ChiOBJ.Z(S.IXgrid)=C;

end

function [ChiOBJ]=MakeChiGrid(DEM,FD,varargin)
	% Clone of DivideTools 'ChiGrid' function with some options disabled.

	% Parse Inputs
	p = inputParser;         
	p.FunctionName = 'ChiGrid';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'FD', @(x) isa(x,'FLOWobj'));

	addParameter(p,'theta_ref',0.5,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'chi_ref_area',1,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'complete_networks_only',true,@(x) isscalar(x) && islogical(x));
	addParameter(p,'min_elevation',[],@(x) isnumeric(x));

	parse(p,DEM,FD,varargin{:});
	DEM=p.Results.DEM;
	FD=p.Results.FD;

	mn=p.Results.theta_ref;
	a0=p.Results.chi_ref_area;
	me=p.Results.min_elevation;
	cno=p.Results.complete_networks_only;

	if isempty(me)
		abl=false;
	else
		abl=true;
	end

	if cno && ~abl
		% Find nodes influenced by edge (From Topotoolbox blog)
		IXE = GRIDobj(DEM,'logical');
		IXE.Z(:,:) = true;
		IXE.Z(2:end-1,2:end-1) = false;
		IXE = influencemap(FD,IXE);
		% Rest is mine
		% Find drainage basins and all outlets
		[DB,oixi]=drainagebasins(FD);
		% Find where these share pixels other than the edge
		db=DB.Z; db=db(3:end-2,3:end-2); % Move 2 pixels in to be conservative
		ixe=IXE.Z; ixe=ixe(3:end-2,3:end-2); % Move 2 pixels in to be conservative
		dbL=db(ixe);
		% Compile list of drainage basins that are influenced by edge pixels
		dbL=unique(dbL);
		% Index list of outlets based on this
		idxi=ismember(DB.Z(oixi),dbL);
		oixi(idxi)=[];
		% Remove drainage basins based on this
		W=dependencemap(FD,oixi);
		% DEM.Z(~mask.Z)=NaN;
		% % Extract info from FLOWobj		
		% W=~isnan(DEM);
		I=W.Z(FD.ix);
		ix=double(FD.ix(I));
		ixc=double(FD.ixc(I));
		% Recalculate indices
		IX        = zeros(FD.size,'uint32');
		IX(W.Z)   = 1:nnz(W.Z);
		ix      = double(IX(ix));
		ixc     = double(IX(ixc));
		I          = ixc == 0;
		ix(I)    = [];
		ixc(I)   = [];
	elseif cno && abl
		% Recalculate flow directions after removing portions below min elevation
		IX=DEM>me;
		DEM.Z(IX.Z==false)=NaN;
		FD=FLOWobj(DEM,'preprocess','carve');
		% Find nodes influenced by edge (From Topotoolbox blog)
		IXE = GRIDobj(DEM,'logical');
		IXE.Z(:,:) = true;
		IXE.Z(2:end-1,2:end-1) = false;
		IXE = influencemap(FD,IXE);
		% Rest is mine
		% Find drainage basins and all outlets
		[DB,oixi]=drainagebasins(FD);
		% Find where these share pixels other than the edge
		db=DB.Z; db=db(3:end-2,3:end-2); % Move 2 pixels in to be conservative
		ixe=IXE.Z; ixe=ixe(3:end-2,3:end-2); % Move 2 pixels in to be conservative
		dbL=db(ixe);
		% Compile list of drainage basins that are influenced by edge pixels
		dbL=unique(dbL);
		% Index list of outlets based on this
		idxi=ismember(DB.Z(oixi),dbL);
		oixi(idxi)=[];
		% Find offending drainage basins and recalculate indicies
		W=dependencemap(FD,oixi);
		I=W.Z(FD.ix);
		ix=double(FD.ix(I));
		ixc=double(FD.ixc(I));
		% Recalculate indices
		IX        = zeros(FD.size,'uint32');
		IX(W.Z)   = 1:nnz(W.Z);
		ix      = double(IX(ix));
		ixc     = double(IX(ixc));
		I          = ixc == 0;
		ix(I)    = [];
		ixc(I)   = [];
	elseif ~cno && ~abl
		% Extract info from FLOWobj	
		W=~isnan(DEM);
		I=W.Z(FD.ix);
		ix=double(FD.ix(I));
		ixc=double(FD.ixc(I));
		% Recalculate indices
		IX        = zeros(FD.size,'uint32');
		IX(W.Z)   = 1:nnz(W.Z);
		ix      = double(IX(ix));
		ixc     = double(IX(ixc));
		I          = ixc == 0;
		ix(I)    = [];
		ixc(I)   = [];
	elseif ~cno && abl
		W=DEM>me;
		I=W.Z(FD.ix);
		ix=double(FD.ix(I));
		ixc=double(FD.ixc(I));
		% Recalculate indices
		IX        = zeros(FD.size,'uint32');
		IX(W.Z)   = 1:nnz(W.Z);
		ix      = double(IX(ix));
		ixc     = double(IX(ixc));
		I          = ixc == 0;
		ix(I)    = [];
		ixc(I)   = [];
	end

	% Generate coordinate list
	IXgrid=find(W.Z);
	[x,y]=ind2coord(DEM,IXgrid);

	% Distance between two nodes throughout grid
	d = nan(size(x));
	dedge = sqrt((x(ixc)-x(ix)).^2 + (y(ixc)-y(ix)).^2);
	d(:) = 0;
	d(ix) = dedge;

	% Cumulative trapezoidal numerical integration of draiange area grid
	DA=flowacc(FD).*(DEM.cellsize^2);
	da=DA.Z(IXgrid);
	c = zeros(size(da));
	da = ((a0) ./da).^mn;
	for r = numel(ix):-1:1;
	    c(ix(r)) = c(ixc(r)) + (da(ixc(r))+(da(ix(r))-da(ixc(r)))/2)*d(ix(r));
	end

	% Generate and populate full chi grid
	ChiOBJ=GRIDobj(DEM);
	ChiOBJ.Z(ChiOBJ.Z==0)=NaN;
	ChiOBJ.Z(IXgrid)=c;
end

function [ksn_ms]=KSN_Quick(DEM,DEMc,A,S,theta_ref,segment_length)
	g=gradient(S,DEMc);
	G=GRIDobj(DEM);
	G.Z(S.IXgrid)=g;

	Z_RES=DEMc-DEM;

	ksn=G./(A.*(A.cellsize^2)).^(-theta_ref);

	SD=GRIDobj(DEM);
	SD.Z(S.IXgrid)=S.distance;
	
	ksn_ms=STREAMobj2mapstruct(S,'seglength',segment_length,'attributes',...
		{'ksn' ksn @mean 'uparea' (A.*(A.cellsize^2)) @mean 'gradient' G @mean 'cut_fill' Z_RES @mean...
		'min_dist' SD @min 'max_dist' SD @max});

	seg_dist=[ksn_ms.max_dist]-[ksn_ms.min_dist];
	distcell=num2cell(seg_dist');
	[ksn_ms(1:end).seg_dist]=distcell{:};
	ksn_ms=rmfield(ksn_ms,{'min_dist','max_dist'});
end

function [ksn_ms]=KSN_Trunk(DEM,DEMc,A,S,theta_ref,segment_length,min_order)

	order_exp=['>=' num2str(min_order)];

    Smax=modify(S,'streamorder',order_exp);
	Smin=modify(S,'rmnodes',Smax);

	g=gradient(S,DEMc);
	G=GRIDobj(DEM);
	G.Z(S.IXgrid)=g;

	Z_RES=DEMc-DEM;

	ksn=G./(A.*(A.cellsize^2)).^(-theta_ref);

	SDmax=GRIDobj(DEM);
	SDmin=GRIDobj(DEM);
	SDmax.Z(Smax.IXgrid)=Smax.distance;
	SDmin.Z(Smin.IXgrid)=Smin.distance;

	ksn_ms_min=STREAMobj2mapstruct(Smin,'seglength',segment_length,'attributes',...
		{'ksn' ksn @mean 'uparea' (A.*(A.cellsize^2)) @mean 'gradient' G @mean 'cut_fill' Z_RES @mean...
		'min_dist' SDmin @min 'max_dist' SDmin @max});

	ksn_ms_max=STREAMobj2mapstruct(Smax,'seglength',segment_length,'attributes',...
		{'ksn' ksn @mean 'uparea' (A.*(A.cellsize^2)) @mean 'gradient' G @mean 'cut_fill' Z_RES @mean...
		'min_dist' SDmax @min 'max_dist' SDmax @max});

	ksn_ms=vertcat(ksn_ms_min,ksn_ms_max);
	seg_dist=[ksn_ms.max_dist]-[ksn_ms.min_dist];
	distcell=num2cell(seg_dist');
	[ksn_ms(1:end).seg_dist]=distcell{:};
	ksn_ms=rmfield(ksn_ms,{'min_dist','max_dist'});
end

function [ksn_ms]=KSN_Trib(DEM,DEMc,FD,A,S,theta_ref,segment_length)

	% Define non-intersecting segments
	w1=waitbar(0,'Finding network segments');
	[as]=networksegment_slim(DEM,FD,S);
	seg_bnd_ix=as.ix;
	% Precompute values or extract values needed for later
	waitbar(1/4,w1,'Calculating hydrologically conditioned stream elevations');
	z=getnal(S,DEMc);
	zu=getnal(S,DEM);
	z_res=z-zu;
	waitbar(2/4,w1,'Calculating chi values');
	g=gradient(S,DEMc);
	c=chitransform(S,A,'a0',1,'mn',theta_ref);
	d=S.distance;
	da=getnal(S,A.*(A.cellsize^2));
	ixgrid=S.IXgrid;
	waitbar(3/4,w1,'Extracting node ordered list');
	% Extract ordered list of stream indices and find breaks between streams
	s_node_list=S.orderednanlist;
	streams_ix=find(isnan(s_node_list));
	streams_ix=vertcat(1,streams_ix);
	waitbar(1,w1,'Pre computations completed');
	close(w1)
	% Generate empty node attribute list for ksn values
	ksn_nal=zeros(size(d));
	% Begin main loop through channels
	num_streams=numel(streams_ix)-1;
	w1=waitbar(0,'Calculating k_{sn} values - 0% Done');
	seg_count=1;
	for ii=1:num_streams
		% Extract node list for stream of interest
		if ii==1
			snlOI=s_node_list(streams_ix(ii):streams_ix(ii+1)-1);
		else
			snlOI=s_node_list(streams_ix(ii)+1:streams_ix(ii+1)-1);
		end

		% Determine which segments are within this stream
		[~,~,dn]=intersect(snlOI,seg_bnd_ix(:,1));
		[~,~,up]=intersect(snlOI,seg_bnd_ix(:,2));
		seg_ix=intersect(up,dn);

		num_segs=numel(seg_ix);
		dn_up=seg_bnd_ix(seg_ix,:);
		for jj=1:num_segs
			% Find positions within node list
			dnix=find(snlOI==dn_up(jj,1));
			upix=find(snlOI==dn_up(jj,2));
			% Extract segment indices of desired segment
			seg_ix_oi=snlOI(upix:dnix);
			% Extract flow distances and normalize
			dOI=d(seg_ix_oi);
			dnOI=dOI-min(dOI);
			num_bins=ceil(max(dnOI)/segment_length);
			bin_edges=[0:segment_length:num_bins*segment_length];
			% Loop through bins
			for kk=1:num_bins
				idx=dnOI>bin_edges(kk) & dnOI<=bin_edges(kk+1);
				bin_ix=seg_ix_oi(idx);
				cOI=c(bin_ix);
				zOI=z(bin_ix);
					if numel(cOI)>2
						[ksn_val,r2]=Chi_Z_Spline(cOI,zOI);
						ksn_nal(bin_ix)=ksn_val;

						% Build mapstructure
						ksn_ms(seg_count).Geometry='Line';
						ksm_ms(seg_count).BoundingBox=[min(S.x(bin_ix)),min(S.y(bin_ix));max(S.x(bin_ix)),max(S.y(bin_ix))];
						ksn_ms(seg_count).X=S.x(bin_ix);
						ksn_ms(seg_count).Y=S.y(bin_ix);
						ksn_ms(seg_count).ksn=ksn_val;
						ksn_ms(seg_count).uparea=mean(da(bin_ix));
						ksn_ms(seg_count).gradient=mean(g(bin_ix));
						ksn_ms(seg_count).cut_fill=mean(z_res(bin_ix));
						ksn_ms(seg_count).seg_dist=max(S.distance(bin_ix))-min(S.distance(bin_ix));
						ksn_ms(seg_count).chi_r2=r2;
						
						seg_count=seg_count+1;
					end
			end
		end
	perc_of_total=round((ii/num_streams)*1000)/10;
	if rem(perc_of_total,1)==0
		waitbar((ii/num_streams),w1,['Calculating k_{sn} values - ' num2str(perc_of_total) '% Done']);
	end
	
	end
	close(w1);
end

function seg = networksegment_slim(DEM,FD,S)
	% Slimmed down version of 'networksegment' from main TopoToolbox library that also removes zero and single node length segments

	%% Identify channel heads, confluences, b-confluences and outlets
	Vhead = streampoi(S,'channelheads','logical');  ihead=find(Vhead==1);  IXhead=S.IXgrid(ihead);
	Vconf = streampoi(S,'confluences','logical');   iconf=find(Vconf==1);  IXconf=S.IXgrid(iconf);
	Vout = streampoi(S,'outlets','logical');        iout=find(Vout==1);    IXout=S.IXgrid(iout);
	Vbconf = streampoi(S,'bconfluences','logical'); ibconf=find(Vbconf==1);IXbconf=S.IXgrid(ibconf);

	%% Identify basins associated to b-confluences and outlets
	DB   = drainagebasins(FD,vertcat(IXbconf,IXout));DBhead=DB.Z(IXhead); DBbconf=DB.Z(IXbconf); DBconf=DB.Z(IXconf); DBout=DB.Z(IXout);

	%% Compute flowdistance
	D = flowdistance(FD);

	%% Identify river segments
	% links between channel heads and b-confluences
	[~,ind11,ind12]=intersect(DBbconf,DBhead);
	% links between confluences and b-confluences
	[~,ind21,ind22]=intersect(DBbconf,DBconf);
	% links between channel heads and outlets
	[~,ind31,ind32]=intersect(DBout,DBhead);
	% links between channel heads and outlets
	[~,ind41,ind42]=intersect(DBout,DBconf);
	% Connecting links into segments
	IX(:,1) = [ IXbconf(ind11)' IXbconf(ind21)' IXout(ind31)'  IXout(ind41)'  ];   ix(:,1)= [ ibconf(ind11)' ibconf(ind21)' iout(ind31)'  iout(ind41)'  ];
	IX(:,2) = [ IXhead(ind12)'  IXconf(ind22)'  IXhead(ind32)' IXconf(ind42)' ];   ix(:,2)= [ ihead(ind12)'  iconf(ind22)'  ihead(ind32)' iconf(ind42)' ];

	% Compute segment flow length
	flength=double(abs(D.Z(IX(:,1))-D.Z(IX(:,2))));

	% Remove zero and one node length elements
	idx=flength>=2*DEM.cellsize;
	seg.IX=IX(idx,:);
	seg.ix=ix(idx,:);
	seg.flength=flength(idx);

	% Number of segments
	seg.n=numel(IX(:,1));
end

function [KSN,R2] = Chi_Z_Spline(c,z)

	% Resample chi-elevation relationship using cubic spline interpolation
	[~,minIX]=min(c);
	zb=z(minIX);
	chiF=c-min(c);
	zabsF=z-min(z);
	chiS=linspace(0,max(chiF),numel(chiF)).';
	zS=spline(chiF,zabsF,chiS);

	% Calculate ksn via slope
	KSN= chiS\(zS); % mn not needed because a0 is fixed to 1

	% Calculate R^2
	z_pred=chiF.*KSN;
	sstot=sum((zabsF-mean(zabsF)).^2);
	ssres=sum((zabsF-z_pred).^2);
	R2=1-(ssres/sstot);

end

function [KSNGrid] = KsnAvg(DEM,ksn_ms,radius)

	% Calculate radius
	radiuspx = ceil(radius/DEM.cellsize);

	% Record mask of current NaNs
	MASK=isnan(DEM.Z);

	% Make grid with values along channels
	KSNGrid=GRIDobj(DEM);
	KSNGrid.Z(:,:)=NaN;
	for ii=1:numel(ksn_ms)
		ix=coord2ind(DEM,ksn_ms(ii).X,ksn_ms(ii).Y);
		KSNGrid.Z(ix)=ksn_ms(ii).ksn;
	end

	% Local mean based on radius
	ISNAN=isnan(KSNGrid.Z);
    [~,L] = bwdist(~ISNAN,'e');
    ksng = KSNGrid.Z(L);           
    FLT   = fspecial('disk',radiuspx);
    ksng   = imfilter(ksng,FLT,'symmetric','same','conv');

    % Set original NaN cells back to NaN
    ksng(MASK)=NaN;

    % Output
    KSNGrid.Z=ksng;
end
