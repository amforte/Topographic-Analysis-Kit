function [S,zpOUT,inOUT]=ProjectedIncision(DEM,A,S,Sc,OUT,varargin)
	% Usage:
	%	[S,zpOUT,inOUT]=ProjectionIncision(DEM,A,S,Sc,OUT);
	%	[S,zpOUT,inOUT]=ProjectionIncision(DEM,A,S,Sc,OUT,'name',value);
	%
	% Description:
	% 	Function for generating maps of projected incision throughout a network based on results of stream
	%	projections from SegmentProjector. Code uses the steepness of the fit portion of the stream network
	%	and the chi values of the entire network to find the projected elevations of the original stream 
	%	network and amounts of incision implied by this, i.e. if you fit portions of streams in SegmentProjector 
	%	that you interpret as recording a former low relief landscape now uplifted and incised, this code
	%	first calculates what the stream elevations would be if no incision (but surface uplift) had occurred.
	%	This value (as a node attributed list) is output in the 'zpOUT' array. The code then subtracts this
	%	projected elevation from the current elevation of the network to estimate the amount of incision this would
	%	imply, which is stored in the 'inOUT' node attributed list array. If you projected multiple streams 
	%	in SegmentProjector (that are stored in the OUT cell array) then the code will calculate projected
	%	elevations and implied incision for each projected stream separately and then find the mean, standard
	%	deviation, minimum, and maximum values of the projected elevations and incision for the network in question.
	%	Additionally, if the Sc input has multiple outlets, the projected elevation and incision values
	%	will only be calculated for portions of the stream network S that are 1) upstream of the outlets in Sc
	%	and 2) connected to channels used to project. E.g. if Sc has two outlets, defining two connected stream
	%	networks we'll call network A and B, and you projected two streams in network A and three streams in 
	%	network B, then the portions of the resulting zpOUT and inOUT values that correspond to nodes in network 
	%	A will only be calculated using the projections from the two streams in network A. This is done in case
	%	there is spatial variability in the amount of incision.
	%
	%	Note that negative incision values indicate that the elevations of these portions of the modern stream 
	%	network are above the projected elevations. It is also important to note that this function explicitly 
	%	assumes that no change in drainage area / network topology has occurred during surface uplift.
	% 
	% Required Inputs:
	%	DEM - Digital Elevation as a GRIDobj, assumes unconditioned DEM (e.g. DEMoc from ProcessRiverBasins)
	%	A - Flow accumulation GRIDobj
	%	S - Full STREAMobj associated with the provided DEM and A grids
	%	Sc - Subset STREAMobj provided to SegmentProjector
	%	OUT - Cell array output from SegmentProjector.
	%
	% Optional Inputs:
	%	shape_name ['bsn'] - prefix on the output shapefiles names
	%	display_figure [true] - logical flag to display a figure showing the values of mean, max, min, and
	%		standard deviation of calculated incision
	%	exlcude_streams [] - optional input if you wish to exclude any of the projected streams in the OUT
	%		cell array from the calculation. Provide a list of stream numbers (e.g. if you wish to exclude
	%		the 2nd and 15th stream that you projected, you would give [2 15] as the input to this parameter)
	%	conditioned_DEM [] - option to provide a hydrologically conditioned DEM for use in this function (do 
	%		not provide a conditoned DEM for the main required DEM input!) which will be used for caculating 
	%		incision. See 'ConditionDEM' function for options for making a hydrological conditioned DEM. If 
	%		no input is provided the code defaults to using the mincosthydrocon function
	%	interp_value [0.1] - value (between 0 and 1) used for interpolation parameter in mincosthydrocon (not 
	%		used if user provides a conditioned DEM)
	%
	% Outputs:
	%	S - STREAMobj that is version of the input 'S' upstream of outlets of network in 'Sc' 
	%	zpOUT - n x 4 array of the mean, standard deviation, minimum, and maximum of the projected elevations
	%		throughout the stream network output with S
	%	inOUT - n x 4 array of the mean, standard deviation, minimum, and maximum of the incision values 
	%		throughout the stream network output with S
	%
	% 	Code also saves two shapefiles, '*_Pnts_Used.shp' and '*_Proj_Incision.shp'. '*_Pnts_Used.shp' is
	%		a point shapefile that records the channel heads of the projected stream used to calculate the 
	%		incision amounts. '*_Proj_Incision.shp' is stream shapefile containing the same outputs as 
	%		zpOUT and inOUT in the attributes
	%
	% Examples:
	%	[S,zpOUT,inOUT]=ProjectedIncision(DEM,A,S,Sc,OUT);
	%	[S,zpOUT,inOUT]=ProjectedIncision(DEM,A,S,Sc,OUT,'exclude_streams',[2 15 40]);
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 04/02/19 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	p = inputParser;
	p.FunctionName = 'ProjectedIncision';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'A',@(x) isa(x,'GRIDobj'));
	addRequired(p,'S',@(x) isa(x,'STREAMobj'));
	addRequired(p,'Sc',@(x) isa(x,'STREAMobj'));
	addRequired(p,'OUT',@(x) iscell(x));

	addParameter(p,'shape_name','bsn',@(x) ischar(x));
	addParameter(p,'display_figure',true,@(x) isscalar(x) && islogical(x));
	addParameter(p,'exclude_streams',[],@(x) isnumeric(x));
	addParameter(p,'conditioned_DEM',[],@(x) isa(x,'GRIDobj'));
	addParameter(p,'interp_value',0.1,@(x) isnumeric(x) && x>=0 && x<=1);

	parse(p,DEM,A,S,Sc,OUT,varargin{:});
	DEM=p.Results.DEM;
	A=p.Results.A;
	S=p.Results.S;
	Sc=p.Results.Sc;
	OUT=p.Results.OUT;

	shape_name=p.Results.shape_name;
	display_figure=p.Results.display_figure;
	exclude_streams=p.Results.exclude_streams;
	DEMc=p.Results.conditioned_DEM;
	iv=p.Results.interp_value;

	% Condition DEM if none is provided
	if isempty(DEMc)
		zc=mincosthydrocon(S,DEM,'interp',iv);
		DEMc=GRIDobj(DEM);
		DEMc.Z(DEMc.Z==0)=NaN;
		DEMc.Z(S.IXgrid)=zc;
	end

	% Filter streams if an input is provided to exclude_streams
	ref_num=[1:size(OUT,2)];
	if ~isempty(exclude_streams)
		idx=logical(ones(1,size(OUT,2)));
		idx(exclude_streams)=false;
		OUT=OUT(:,idx);
		ref_num=ref_num(idx);
	end

	% Find outlets of selected streams
	outix=streampoi(Sc,'outlets','ix');
	% Modify full stream network as necessary
	S=modify(S,'upstreamto',outix);
	% Make label grid for identification
	[L,nc]=conncomps(S);
	LG=GRIDobj(DEM);
	LG.Z(S.IXgrid)=L;

	% Generate empty arrays and structures
	num_proj=size(OUT,2);
	num_nodes=numel(S.x);
	zpM=zeros(num_nodes,num_proj);
	inM=zeros(num_nodes,num_proj);
	ch=struct;

	for ii=1:num_proj
		% Extract values from projected dataset
		chx=OUT{1,ii}(:,1);
		chy=OUT{1,ii}(:,2);
		c=OUT{2,ii}(:,5);
		mn=OUT{2,ii}(:,6); mn=mn(1);
		zp=OUT{2,ii}(:,8);

		% Determine which network ths stream belongs to
		chix=coord2ind(DEM,chx,chy);
		loi=LG.Z(chix);

		% Remove any nans
		idx=~isnan(c);
		c=c(idx); zp=zp(idx);

		% Find projected outlet
		outzp=min(zp);

		% Find ksn of projected stream network
		ksn=c\(zp-outzp);

		% Use modern chi and projected ksn to calculate
		% elevation of continous stream network at projected elevation
		cnal=chitransform(S,A,'a0',1,'mn',mn);
		zpnal=cnal.*ksn;
		zpnal=zpnal+outzp;

		% Use above to caclulate amount of incision
		inc=zpnal-getnal(S,DEMc);

		% Set values not in network of interest to NaN
		lidx=L~=loi;
		zpnal(lidx)=NaN;
		inc(lidx)=NaN;

		% Package output
		zpM(:,ii)=zpnal;
		inM(:,ii)=inc;

		% Build list of channel heads used in calculation
		ch(ii,1).Geometry='Point';
		ch(ii,1).X=double(chx);
		ch(ii,1).Y=double(chy);
		ch(ii,1).rivID=double(ref_num(ii));
	end

	% Generate nals
	mean_in=mean(inM,2,'omitnan');
	std_in=std(inM,0,2,'omitnan');
	min_in=min(inM,[],2,'omitnan');
	max_in=max(inM,[],2,'omitnan');

	mean_zp=mean(zpM,2,'omitnan');
	std_zp=std(zpM,0,2,'omitnan');
	min_zp=min(zpM,[],2,'omitnan');
	max_zp=max(zpM,[],2,'omitnan');

	zpOUT=[mean_zp std_zp min_zp max_zp];
	inOUT=[mean_in std_in min_in max_in];


	if display_figure
		f1=figure(1);
		set(f1,'unit','normalized','position',[0.1 0.1 0.8 0.8]);

		subplot(2,2,1);
		hold on
		[RGB]=imageschs(DEM,DEM,'colormap','gray');
		[~,R]=GRIDobj2im(DEM);
		imshow(flipud(RGB),R);
		axis xy
		plotc(S,mean_in);
		colorbar;
		title('Mean Incision');
		hold off

		subplot(2,2,2);
		hold on
		imshow(flipud(RGB),R);
		axis xy
		plotc(S,std_in);
		colorbar;
		title('StDev Incision');
		hold off

		subplot(2,2,3);
		hold on 
		imshow(flipud(RGB),R);
		axis xy
		plotc(S,min_in);
		colorbar;
		title('Min Incision');
		hold off

		subplot(2,2,4);
		hold on 
		imshow(flipud(RGB),R);
		axis xy
		plotc(S,max_in);
		colorbar;
		title('Max Incision');
		hold off
	end	

	ms=STREAMobj2mapstruct(S,'seglength',DEM.cellsize*3,'attributes',{'mean_inc' mean_in @mean 'std_inc' std_in @mean...
		'min_inc' min_in @mean 'max_inc' max_in @mean 'mean_zp' mean_zp @mean 'std_zp' std_zp @mean 'min_zp' min_zp @mean...
		'max_zp' max_zp @mean});

	pnts_name=[shape_name '_Pnts_Used.shp'];
	strm_name=[shape_name '_Proj_Incision.shp'];

	shapewrite(ch,pnts_name);
	shapewrite(ms,strm_name);

end