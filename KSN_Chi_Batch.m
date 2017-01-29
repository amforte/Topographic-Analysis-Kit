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
	%	file_name_prefix ['batch'] - prefix for outputs, will append the type of output, i.e. 'ksn', 'chimap', etc
	% 	segment_length [1000] - length of segments in map units for smoothing ksn values, equivalent to smoothing in Profiler
	% 	theta_ref [0.45] - reference concavity (as a positive value) for calculating ksn
	% 	output [false]- switch to either output matlab files to the workspace (true) or to not only save the specified files
	%		without any workspace output (false)
	%	ksn_method [quick] - switch between method to calculate ksn values, options are 'quick' and 'trib', the 'trib' method takes 3-4 times longer 
	%		than the 'quick' method. In most cases, the 'quick' method works well, but if values near tributary junctions are important, then 'trib'
	%		may be better as this calculates ksn values for individual channel segments individually
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
	%	[KSN,ChiMap,ChiGrid]=KSN_Chi_Batch(DEM,FD,A,S,'output',true,'theta_ref',0.55);
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

	addParamValue(p,'file_name_prefix','batch',@(x) ischar(x));
	addParamValue(p,'segment_length',1000,@(x) isscalar(x) && isnumeric(x));
	addParamValue(p,'theta_ref',0.45,@(x) isscalar(x) && isnumeric(x));
	addParamValue(p,'output',false,@(x) isscalar(x) && islogical(x));
	addParamValue(p,'ksn_method','quick',@(x) ischar(validatestring(x,{'quick','trib'})));
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

	file_name_prefix=p.Results.file_name_prefix;
	segment_length=p.Results.segment_length;
	theta_ref=p.Results.theta_ref;
	output=p.Results.output;
	ksn_method=p.Results.ksn_method;
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
		
		if abl & strcmp(blm,'elevation')
			[S]=SetBaseLevel(DEM,FD,A,S,blm,'min_elevation',me);
		elseif abl & strcmp(blm,'drain_area');
			[S]=SetBaseLevel(DEM,FD,A,S,blm,'max_drainage_area',ma);
		elseif abl
			[S]=SetBaseLevel(DEM,FD,A,S,blm);
		end	
		
		switch ksn_method
		case 'quick'
			[KSN]=KSN_Quick(DEM,A,S,theta_ref,segment_length);
		case 'trib'
			[KSN]=KSN_Trib(DEM,FD,A,S,theta_ref,segment_length);
		end

		disp('Writing ARC files')
		out_file=[file_name_prefix '_ksn.shp'];
		shapewrite(KSN,out_file);

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
		[ChiMap]=MakeChiMap(DEM,FD,A,S,theta_ref);

		disp('Writing ARC files')
		GRIDobj2ascii(ChiMap,'chimap.txt');

		switch output
		case true
			varargout{1}=ChiMap;
		end

	case 'chigrid'

	    disp('Calculating chi grid');
		if abl & strcmp(blm,'elevation')
		    [ChiGrid]=MakeChiGrid(DEM,FD,A,theta_ref,abl,blm,me);
		elseif abl & strcmp(blm,'drain_area');
		    [ChiGrid]=MakeChiGrid(DEM,FD,A,theta_ref,abl,blm,ma);
		elseif abl
		    [ChiGrid]=MakeChiGrid(DEM,FD,A,theta_ref,abl,blm);
		else
			[ChiGrid]=MakeChiGrid(DEM,FD,A,theta_ref,abl);
		end	

	    disp('Writing ARC files')
		out_file=[file_name_prefix '_chigrid.txt'];
		GRIDobj2ascii(ChiGrid,out_file);

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
		[ChiMap]=MakeChiMap(DEM,FD,A,S,theta_ref);

	    disp('Calculating chi grid');
		if abl & strcmp(blm,'elevation')
		    [ChiGrid]=MakeChiGrid(DEM,FD,A,theta_ref,abl,blm,me);
		elseif abl & strcmp(blm,'drain_area');
		    [ChiGrid]=MakeChiGrid(DEM,FD,A,theta_ref,abl,blm,ma);
		elseif abl
		    [ChiGrid]=MakeChiGrid(DEM,FD,A,theta_ref,abl,blm);
		else
			[ChiGrid]=MakeChiGrid(DEM,FD,A,theta_ref,abl);
		end	

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

		disp('Calculating channel steepness')
		if abl & strcmp(blm,'elevation')
			[S]=SetBaseLevel(DEM,FD,A,S,blm,'min_elevation',me);
		elseif abl & strcmp(blm,'drain_area');
			[S]=SetBaseLevel(DEM,FD,A,S,blm,'max_drainage_area',ma);
		elseif abl
			[S]=SetBaseLevel(DEM,FD,A,S,blm);
		end	

		switch ksn_method
		case 'quick'
			[KSN]=KSN_Quick(DEM,A,S,theta_ref,segment_length);
		case 'trib'
			[KSN]=KSN_Trib(DEM,FD,A,S,theta_ref,segment_length);
		end

	    disp('Calculating chi map');
		[ChiMap]=MakeChiMap(DEM,FD,A,S,theta_ref);

	    disp('Calculating chi grid');
		if abl & strcmp(blm,'elevation')
		    [ChiGrid]=MakeChiGrid(DEM,FD,A,theta_ref,abl,blm,me);
		elseif abl & strcmp(blm,'drain_area');
		    [ChiGrid]=MakeChiGrid(DEM,FD,A,theta_ref,abl,blm,ma);
		elseif abl
		    [ChiGrid]=MakeChiGrid(DEM,FD,A,theta_ref,abl,blm);
		else
			[ChiGrid]=MakeChiGrid(DEM,FD,A,theta_ref,abl);
		end	


		disp('Writing ARC files')
		out_file_ksn=[file_name_prefix '_ksn.shp'];
		shapewrite(KSN,out_file_ksn);
		out_file_cg=[file_name_prefix '_chigrid.txt'];
		GRIDobj2ascii(ChiGrid,out_file_cg);
		out_file_cm=[file_name_prefix '_chimap.txt'];
		GRIDobj2ascii(ChiMap,out_file_cm);

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

function [ChiOBJ]=MakeChiMap(DEM,FD,A,S,theta_ref);

	DA=A.*(A.cellsize^2);
	C=chitransform(S,DA,'mn',theta_ref,'a0',1);

	% Make Empty GRIDobj
	ChiOBJ=GRIDobj(DEM);
	ChiOBJ.Z(ChiOBJ.Z==0)=NaN;

	% Populate Grid
	ChiOBJ.Z(S.IXgrid)=C;

end

function [ChiOBJ]=MakeChiGrid(DEM,FD,A,theta_ref,abl,varargin)

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
	C=chitransform(S,DA,'mn',theta_ref,'a0',1);

	% Make Empty GRIDobj
	ChiOBJ=GRIDobj(A);
	ChiOBJ.Z(ChiOBJ.Z==0)=NaN;

	% Populate Grid
	ChiOBJ.Z(S.IXgrid)=C;
end

function [ksn_ms]=KSN_Quick(DEM,A,S,theta_ref,segment_length)

	zc=mincosthydrocon(S,DEM,'interp',0.1);
	DEMc=GRIDobj(DEM);
	DEMc.Z(DEMc.Z==0)=NaN;
	DEMc.Z(S.IXgrid)=zc;
	G=gradient8(DEMc);
	Z_RES=DEMc-DEM;

	ksn=G./(A.*(A.cellsize^2)).^(-theta_ref);
	
	ksn_ms=STREAMobj2mapstruct(S,'seglength',segment_length,'attributes',...
		{'ksn' ksn @mean 'uparea' (A.*(A.cellsize^2)) @mean 'gradient' G @mean 'cut_fill' Z_RES @mean});
end

function [ksn_ms]=KSN_Trib(DEM,FD,A,S,theta_ref,segment_length)

	% Define non-intersecting segments
	w1=waitbar(0,'Finding network segments');
	[as]=networksegment_slim(DEM,FD,S);
	seg_bnd_ix=as.ix;
	% Precompute values or extract values needed for later
	waitbar(1/4,w1,'Calculating hydrologically conditioned stream elevations');
	z=mincosthydrocon(S,DEM,'interp',0.1);
	zu=getnal(S,DEM);
	z_res=z-zu;
	waitbar(2/4,w1,'Calculating chi values');
	c=chitransform(S,A.*(A.cellsize^2),'a0',1,'mn',theta_ref);
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
						[ksn_val]=Chi_Z_Spline(cOI,zOI);
						ksn_nal(bin_ix)=ksn_val;

						% Build mapstructure
						ksn_ms(seg_count).Geometry='Line';
						ksn_ms(seg_count).X=S.x(bin_ix);
						ksn_ms(seg_count).Y=S.y(bin_ix);
						ksn_ms(seg_count).ksn=ksn_val;
						ksn_ms(seg_count).cut_fill=mean(z_res(bin_ix));
						ksn_ms(seg_count).area=mean(da(bin_ix));
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

function [KSN] = Chi_Z_Spline(c,z)

	% Resample chi-elevation relationship using cubic spline interpolation
	[~,minIX]=min(c);
	zb=z(minIX);
	chiF=c-min(c);
	zabsF=z-min(z);
	chiS=linspace(0,max(chiF),numel(chiF)).';
	zS=spline(chiF,zabsF,chiS);

	%Calculate beta
    BETA = chiS\(zS);

	KSN= BETA; %Beta.*a0^mn - if a0 set to 1, not needed
end