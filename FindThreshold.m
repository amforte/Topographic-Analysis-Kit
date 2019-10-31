function [Sn,thresh_list,xd_list]=FindThreshold(DEM,FD,A,S,num_streams,varargin)
	%
	% Usage:
	%	[Sn,thresh_list,xd_list]=FindThreshold(DEM,FD,A,S,num_streams);
	%	[Sn,thresh_list,xd_list]=FindThreshold(DEM,FD,A,S,num_streams,'name',value,...);
	% 
	% Description:
	% 	Function to interactively select an appopriate threshold area for a given stream
	%	network. Function will either have you iterate through a number of single streams,
	%	controlled by the number passed to 'num_streams', extracted from the drainage divide or 
	%	all streams within the provided drainage network if you provide 'all' to 'num_streams'. 
	%	If 'num_streams' is numeric, then the function will use the average of the user selected 
	%	minimum threshold areas to define a new stream network. If 'num_streams' is set to 'all', the
	%	function will use the user selected minimum threshold areas to define a new stream network for
	%	each individual stream, i.e. the minimum threshold area will be different for each stream base
	%	on your selections (though you can force a similar average threshold behavior with the 
	%	'remake_network' parameter).  You can use either chi-elevation or slope-area plots (the default),
	%	both plots will be displayed regardless of choice, to visually select where channels begin.  Function
	%	also outputs the lists of selected threshold areas and distance from channel head to divide. There is
	%	now also an 'auto' option which will attempt to find the colluvial-fluvial transition using the
	%	ischange function on a normalized chi-elevation plot. Using the 'auto' option requires MATLAB 2017a 
	%	or greater. The auto function may be over shoot the threshold area, especially in areas where there 
	%	is a noticeable lower relief portion of a stream, e.g. a stream draining the edge of a plateau, 
	%	one that has recently captured drainage area, etc. You can try to control for this by setting 
	%	the optional 'max_threshold' parameter, which will serve as a ceiling on the auto detected threshold.
	%
	% Required Inputs:
	%	DEM - Digital Elevation as a GRIDobj, assumes unconditioned DEM (e.g. DEMoc from ProcessRiverBasins 
	%		or output from MakeStreams)
	%	FD - Flow direction as FLOWobj
	%	A - Flow accumulation GRIDobj
	%	S - Stream network as STREAMobj (doesn't really matter what the threshold area used to calculate 
	%		this was originally)
	%	num_streams - Number of stream profiles to view and select threshold areas. If you wish to manually 
	%		select threshold areas for all streams in the provided network, provide 'all' instead of a number.
	%		If you want to instead try an auto detection method for the colluvial to fluvial transtion, provide
	%		'auto'.
	%
	% Optional Inputs;
	%	ref_concavity [0.50] - refrence concavity used to generate the chi-elevation plot. Also used if 
	%		'num_streams' is set to 'auto' as the threshold area is determined based on changes in chi-elevation
	%		relationships
	%	pick_method ['slope_area']- type of plot you wish to choose the threshold area on, valid options are:
	%		'chi' - Choose threshold areas on a chi elevation plot
	%		'slope_area' - Choose threshold areas on slope-area plot
	%	remake_network [false] - logical flag used if 'num_streams' is 'auto' or 'all'. If remake_network is false 
	%		(the default), then each channel head will have a threshold area based on the user selection or the 
	%		auto detection. If remake_network is true, then the behavior will be similar to if a numeric value is
	%		provided to 'num_streams', where the mean threshold area of all streams will be used to regenerate the
	%		network. If remake_network is false, this means the total number of channel heads of the output 'Sn' 
	%		STREAMobj will equal that of the input 'S' and thus the drainage density will basically remain unchanged.
	%		If remake_network is true then there will likely be different numbers of channel heads between Sn and S 
	%		(unless the average threshold area happens to equals the minimum threshold area used to generate S).
	%	max_threshold [] - a maximum threshold area (in squared map units) to not be exceeded by any stream when the
	%		'num_streams' is set to 'auto'.
	%
	% Outputs:
	%	Sn - New version of the STREAMobj using the mean threshold area to define streams or the user selected
	%		minimum threshold areas for each stream if 'num_streams' is set to 'all' or 'auto'
	%	thresh_list - list of chosen threshold areas
	%	xd_list - list of chosen distances from channel head to divide
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 02/09/19 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'FitThreshold';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'FD',@(x) isa(x,'FLOWobj'));
	addRequired(p,'A',@(x) isa(x,'GRIDobj'));
	addRequired(p,'S',@(x) isa(x,'STREAMobj'));
	addRequired(p,'num_streams',@(x) isnumeric(x) & isscalar(x) || ischar(validatestring(x,{'all','auto'})));

	addParameter(p,'pick_method','slope_area',@(x) ischar(validatestring(x,{'chi','slope_area'})));
	addParameter(p,'ref_concavity',0.50,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'remake_network',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'max_threshold',[],@(x) isscalar(x) && isnumeric(x) || isempty(x));

	parse(p,DEM,FD,A,S,num_streams,varargin{:});
	DEM=p.Results.DEM;
	FD=p.Results.FD;
	A=p.Results.A;
	S=p.Results.S;
	num_streams=p.Results.num_streams;

	pick_method=p.Results.pick_method;
	ref_theta=p.Results.ref_concavity;
	remake_network=p.Results.remake_network;
	mt=p.Results.max_threshold;

	% Find channel heads and flow distances
	chix=streampoi(S,'channelheads','ix');
	FLUS=flowdistance(FD);
	FLDS=flowdistance(FD,'downstream');
	DA=A.*(DEM.cellsize^2);

	% Set bin size based on DEM resolution
	bin_size=15*DEM.cellsize;

	% Parse type of operation
	if isnumeric(num_streams)
		op=1;
	elseif strcmp(num_streams,'all')
		op=2;
	elseif strcmp(num_streams,'auto')
		op=3;
	end

	if verLessThan('matlab','9.3') & op==3;
		warning('Using the "auto" version of this function requires MATLAB 2017a or later, switching to "all" streams');
		op=2;
	end

	switch op
	case 1

		% Determine number of channels and compare that to total number of channels
		num_ch=numel(chix);
		if num_streams>num_ch;
			num_streams=num_ch;
			if isdeployed
				warndlg('Number of streams to fit was greater than number of channel heads, using all streams')
			else
				warning('Number of streams to fit was greater than number of channel heads, using all streams')
			end
		end

		% Sort channels by length, fitting proceeds from largest to smallest channels
		fl=FLUS.Z(chix);
		[fl,six]=sort(fl,'descend');
		chix=chix(six);

		for ii=1:num_streams
			chOI=chix(ii);

			UP=dependencemap(FD,chOI);
			FLDSt=FLUS.*UP;

			[~,ix]=max(FLDSt);

			IX=influencemap(FD,ix);

			St=STREAMobj(FD,IX);
			z=mincosthydrocon(St,DEM,'interp',0.1);

			C=chiplot(St,z,A,'a0',1,'mn',ref_theta,'plot',false);
			[bs,ba,bc,bd,aa,ag,ac]=sa(DEM,St,A,C.chi,bin_size);

			f1=figure(1);
			set(f1,'Units','normalized','Position',[0.5 0.1 0.45 0.8],'renderer','painters');
			clf

			colormap(jet);

			switch pick_method
			case 'chi'

				ax2=subplot(2,1,2);
				hold on 
				scatter(aa,ag,5,ac,'+');
				scatter(ba,bs,20,bc,'filled','MarkerEdgeColor','k');
				xlabel('Log Drainage Area');
				ylabel('Log Gradient');
				caxis([0 max(C.chi)]);
				set(ax2,'YScale','log','XScale','log','XDir','reverse');
		        if ~verLessThan('matlab','9.5')
		            disableDefaultInteractivity(ax2);
		        end 
				hold off

				ax1=subplot(2,1,1);
				hold on
				plot(C.chi,C.elev,'-k');
				scatter(C.chi,C.elev,10,C.chi,'filled');
				xlabel('\chi');
				ylabel('Elevation (m)');
				title(['Choose hillslope to channel transition : ' num2str(num_streams-ii) ' streams remaining to pick']);
				caxis([0 max(C.chi)]);
				ax1.XColor='Red';
				ax1.YColor='Red';
		        if ~verLessThan('matlab','9.5')
		            disableDefaultInteractivity(ax1);
		        end 				
				hold off

				% Find user selected threshold area
				[c,~]=ginput(1);
				[~,cix]=min(abs(C.chi-c),[],'omitnan');
				a=C.area(cix);

				% Find xd
				cx=C.x(cix);
				cy=C.y(cix);
				ccix=coord2ind(DEM,cx,cy);
				xd=FLDS.Z(ccix);

			case 'slope_area'

				ax2=subplot(2,1,2);
				hold on
				plot(C.chi,C.elev,'-k');
				scatter(C.chi,C.elev,10,C.chi,'filled');
				xlabel('\chi');
				ylabel('Elevation (m)');
				caxis([0 max(C.chi)]);
		        if ~verLessThan('matlab','9.5')
		            disableDefaultInteractivity(ax2);
		        end 
				hold off

				ax1=subplot(2,1,1);
				hold on 
				scatter(aa,ag,5,ac,'+');
				scatter(ba,bs,20,bc,'filled','MarkerEdgeColor','k');
				xlabel('Log Drainage Area');
				ylabel('Log Gradient');
				title(['Choose hillslope to channel transition : ' num2str(num_streams-ii) ' streams remaining to pick']);
				caxis([0 max(C.chi)]);
				set(ax1,'YScale','log','XScale','log','XDir','reverse');
				ax1.XColor='Red';
				ax1.YColor='Red';

		        if ~verLessThan('matlab','9.5')
		            disableDefaultInteractivity(ax1);
		        end 
				hold off

				% Find user selected threshold area
				[a,~]=ginput(1);

				% Find xd
				[~,cix]=min(abs(C.area-a),[],'omitnan');
				cx=C.x(cix);
				cy=C.y(cix);
				ccix=coord2ind(DEM,cx,cy);
				xd=FLDS.Z(ccix);

			end
			close(f1);
			thresh_list(ii,1)=a;
			xd_list(ii,1)=xd;
		end

		mean_thresh=mean(thresh_list);
		mean_xd=mean(xd_list);

		f1=figure(1);
		set(f1,'Units','normalized','Position',[0.5 0.1 0.4 0.8],'renderer','painters');
		sbplt1=subplot(2,1,1);
		hold on
		edges=logspace(log10(min(thresh_list)),log10(max(thresh_list)),10);
		[N,e]=histcounts(thresh_list,edges);
		histogram(thresh_list,edges);
		plot([mean_thresh,mean_thresh],[0,max(N)],'-k','LineWidth',2);
		xlabel('Picked Threshold Areas (m^{2})');
		a_str=sprintf('%0.1e',mean_thresh);
		title(['Mean threshold area = ' a_str]);
        if ~verLessThan('matlab','9.5')
            disableDefaultInteractivity(sbplt1);
        end 
		hold off

		sbplt2=subplot(2,1,2);
		hold on
		[N,~]=histcounts(xd_list,10);
		histogram(xd_list,10);
		plot([mean_xd,mean_xd],[0,max(N)],'-k','LineWidth',2);
		xlabel('Mean distance from channel head to divide (m)');
		title(['Mean xd = ' num2str(round(mean_xd,1))]);
        if ~verLessThan('matlab','9.5')
            disableDefaultInteractivity(sbplt2);
        end 		
		hold off

		Sn=STREAMobj(FD,'minarea',mean_thresh,'unit','mapunits');

	case 2

		% Sort channels by length, fitting proceeds from largest to smallest channels
		fl=FLUS.Z(chix);
		[fl,six]=sort(fl,'descend');
		chix=chix(six);

		% Generate empty outputs
		xd_list=zeros(numel(chix),1);
		thresh_list=zeros(numel(chix),1);
		ix_list=cell(numel(chix),1);

		for ii=1:numel(chix)
			chOI=chix(ii);

			UP=dependencemap(FD,chOI);
			FLDSt=FLUS.*UP;

			[~,ix]=max(FLDSt);

			IX=influencemap(FD,ix);

			St=STREAMobj(FD,IX);
			z=mincosthydrocon(St,DEM,'interp',0.1);

			C=chiplot(St,z,A,'a0',1,'mn',ref_theta,'plot',false);
			[bs,ba,bc,bd,aa,ag,ac]=sa(DEM,St,A,C.chi,bin_size);

			f1=figure(1);
			set(f1,'Units','normalized','Position',[0.5 0.1 0.45 0.8],'renderer','painters');
			clf

			colormap(jet);

			switch pick_method
			case 'chi'

				ax2=subplot(2,1,2);
				hold on 
				scatter(aa,ag,5,ac,'+');
				scatter(ba,bs,20,bc,'filled','MarkerEdgeColor','k');
				xlabel('Log Drainage Area');
				ylabel('Log Gradient');
				caxis([0 max(C.chi)]);
				set(ax2,'YScale','log','XScale','log','XDir','reverse');
		        if ~verLessThan('matlab','9.5')
		            disableDefaultInteractivity(ax2);
		        end 
				hold off

				ax1=subplot(2,1,1);
				hold on
				plot(C.chi,C.elev,'-k');
				scatter(C.chi,C.elev,10,C.chi,'filled');
				xlabel('\chi');
				ylabel('Elevation (m)');
				title(['Choose hillslope to channel transition : ' num2str(numel(chix)-ii) ' streams remaining to pick']);
				caxis([0 max(C.chi)]);
				ax1.XColor='Red';
				ax1.YColor='Red';
		        if ~verLessThan('matlab','9.5')
		            disableDefaultInteractivity(ax1);
		        end 				
				hold off

				% Find user selected threshold area
				[c,~]=ginput(1);
				[~,cix]=min(abs(C.chi-c),[],'omitnan');
				a=C.area(cix);

				% Find xd
				cx=C.x(cix);
				cy=C.y(cix);
				ccix=coord2ind(DEM,cx,cy);
				xd=FLDS.Z(ccix);

				% Find ix for all streams segments
				allx=C.x(C.area>=a);
				ally=C.y(C.area>=a);
				ix_list{ii,1}=coord2ind(DEM,allx,ally);

			case 'slope_area'

				ax2=subplot(2,1,2);
				hold on
				plot(C.chi,C.elev,'-k');
				scatter(C.chi,C.elev,10,C.chi,'filled');
				xlabel('\chi');
				ylabel('Elevation (m)');
				caxis([0 max(C.chi)]);
		        if ~verLessThan('matlab','9.5')
		            disableDefaultInteractivity(ax2);
		        end 
				hold off

				ax1=subplot(2,1,1);
				hold on 
				scatter(aa,ag,5,ac,'+');
				scatter(ba,bs,20,bc,'filled','MarkerEdgeColor','k');
				xlabel('Log Drainage Area');
				ylabel('Log Gradient');
				title(['Choose hillslope to channel transition : ' num2str(numel(chix)-ii) ' streams remaining to pick']);
				caxis([0 max(C.chi)]);
				set(ax1,'YScale','log','XScale','log','XDir','reverse');
				ax1.XColor='Red';
				ax1.YColor='Red';
		        if ~verLessThan('matlab','9.5')
		            disableDefaultInteractivity(ax1);
		        end 
				hold off

				% Find user selected threshold area
				[a,~]=ginput(1);

				% Find xd
				[~,cix]=min(abs(C.area-a),[],'omitnan');
				cx=C.x(cix);
				cy=C.y(cix);
				ccix=coord2ind(DEM,cx,cy);
				xd=FLDS.Z(ccix);

				% Find ix for all streams segments
				allx=C.x(C.area>=a);
				ally=C.y(C.area>=a);
				ix_list{ii,1}=coord2ind(DEM,allx,ally);

			end
			close(f1);
			thresh_list(ii,1)=a;
			xd_list(ii,1)=xd;
		end

		f1=figure(1);
		set(f1,'Units','normalized','Position',[0.5 0.1 0.4 0.8],'renderer','painters');
		sbplt1=subplot(2,1,1);
		hold on
		edges=logspace(log10(min(thresh_list)),log10(max(thresh_list)),10);
		[N,e]=histcounts(thresh_list,edges);
		histogram(thresh_list,edges);
		xlabel('Picked Threshold Areas (m^{2})');
        if ~verLessThan('matlab','9.5')
            disableDefaultInteractivity(sbplt1);
        end 
		hold off

		sbplt2=subplot(2,1,2);
		hold on
		[N,~]=histcounts(xd_list,10);
		histogram(xd_list,10);
		xlabel('Mean distance from channel head to divide (m)');
        if ~verLessThan('matlab','9.5')
            disableDefaultInteractivity(sbplt2);
        end 
		hold off

		if remake_network
			mean_thresh=mean(thresh_list);
			Sn=STREAMobj(FD,'minarea',mean_thresh,'unit','mapunits');
		else
			% Collapse ix list and create logical raster
			ix_list=vertcat(ix_list{:});
			ix_list=unique(ix_list);
			W=GRIDobj(DEM,'logical');
			W.Z(ix_list)=true;

			Sn=STREAMobj(FD,W);
		end

	case 3
		% Generate empty outputs
		xd_list=zeros(numel(chix),1);
		thresh_list=zeros(numel(chix),1);
		ix_list=cell(numel(chix),1);

		w1=waitbar(0,'Finding stream thresholds based on chi-elevation inflection...');
		for ii=1:numel(chix)
			chOI=chix(ii);

			UP=dependencemap(FD,chOI);
			FLDSt=FLUS.*UP;

			[~,ix]=max(FLDSt);

			IX=influencemap(FD,ix);

			St=STREAMobj(FD,IX);
			z=mincosthydrocon(St,DEM,'interp',0.1);
			c=chitransform(St,A,'mn',ref_theta,'a0',1);
			da=getnal(St,DA);

			% ischange expects an ascending ordered list
			[cs,six]=sort(c);
			zs=z(six);
			das=da(six);
			
			% Normalize chi and elevation so consistent threshold can be used
			zn=zs-min(zs);
			zn=zn./max(zn);
			cn=cs./max(cs);

			TF=ischange(zn,'linear','SamplePoints',cn,'Threshold',0.05);
			chi_ix=find(TF,1,'last');

			if isempty(chi_ix)
				% Catch for if no change point is detected
				ta=DA.Z(chOI);
			else 
				ta=das(chi_ix);
				if ~isempty(mt) && ta>mt
					ta=mt;
				end
			end

			thresh_list(ii,1)=ta;
			[~,loc_ix]=min(abs(da-ta));

			% Find xd
			cx=St.x(loc_ix);
			cy=St.y(loc_ix);
			ccix=coord2ind(DEM,cx,cy);
			xd_list(ii,1)=FLDS.Z(ccix);

			% Find ix for all streams segments
			ix_list{ii,1}=St.IXgrid(da>=ta);			
			waitbar(ii/numel(chix));
		end
		close(w1);

		f1=figure(1);
		clf
		set(f1,'Units','normalized','Position',[0.5 0.1 0.4 0.8],'renderer','painters');
		sbplt1=subplot(2,1,1);
		hold on
		edges=logspace(log10(min(thresh_list)),log10(max(thresh_list)),10);
		[N,e]=histcounts(thresh_list,edges);
		histogram(thresh_list,edges);
		set(sbplt1,'XScale','log');
		xlabel('Picked Threshold Areas (m^{2})');
        if ~verLessThan('matlab','9.5')
            disableDefaultInteractivity(sbplt1);
        end 
		hold off

		sbplt2=subplot(2,1,2);
		hold on
		[N,~]=histcounts(xd_list,50);
		histogram(xd_list,50);
		xlabel('Mean distance from channel head to divide (m)');
        if ~verLessThan('matlab','9.5')
            disableDefaultInteractivity(sbplt2);
        end 
		hold off

		if remake_network
			mean_thresh=mean(thresh_list);
			Sn=STREAMobj(FD,'minarea',mean_thresh,'unit','mapunits');
		else
			% Collapse ix list and create logical raster
			ix_list=vertcat(ix_list{:});
			ix_list=unique(ix_list);
			W=GRIDobj(DEM,'logical');
			W.Z(ix_list)=true;

			Sn=STREAMobj(FD,W);
		end
	end
end

function [bs,ba,bc,bd,a,g,C]=sa(DEM,S,A,C,bin_size)
	% Modified slope area function that uses the smooth length to
	%	to determine the number of bins and uses those same bins
	%	to find mean values of chi and distance for plotting
	%	purposes

	minX=min(S.distance);
	maxX=max(S.distance);
	b=[minX:bin_size:maxX+bin_size];

	numbins=round(max([numel(b) numel(S.IXgrid)/10]));

	an=getnal(S,A.*A.cellsize^2);
	z=getnal(S,DEM);
	gn=gradient(S,z,'unit','tangent','method','robust','drop',20);
	gn=smooth(gn,3);

	% Run through STREAMobj2XY so chi and everything else are same size
	[~,~,a,g,d]=STREAMobj2XY(S,an,gn,S.distance);
	% Remove NaNs
	a(isnan(a))=[];
	g(isnan(g))=[];
	d(isnan(d))=[];
	C(isnan(C))=[];

	mina=min(a);
	maxa=max(a);

    edges = logspace(log10(mina-0.1),log10(maxa+1),numbins+1);
    try
    	% histc is deprecated
    	[ix]=discretize(a,edges);
    catch
	    [~,ix] = histc(a,edges);
	end

	ba=accumarray(ix,a,[numbins 1],@median,nan);
	bs=accumarray(ix,g,[numbins 1],@(x) mean(x(~isnan(x))),nan);
	bd=accumarray(ix,d,[numbins 1],@mean,nan);
	bc=accumarray(ix,C,[numbins 1],@mean,nan);

	% Filter negatives
	idx=bs>=0 & ba>=0 & bc>=0 & bd>=0;
	bs=bs(idx);
	ba=ba(idx);
	bc=bc(idx);
	bd=bd(idx);

	idx=a>=0 & g>=0 & d>=0 & C>=0;
	a=a(idx);
	g=g(idx);
	d=d(idx);
	C=C(idx);
end