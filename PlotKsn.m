function PlotKsn(DEM,FD,ksn,varargin)
	%
	% Usage:
	%	PlotKsn(DEM,FD,ksn);
	%	PlotKsn(DEM,FD,ksn,bnd_list);
	%	PlotKsn(DEM,FD,ksn,'knicks.shp');
	%
	% Description:
	% 	Function to plot a map of normalized channel steepness on a hillshade colored
	% 	by elevation.
	%
	% Required Inputs:
	%	DEM -  Digital Elevation as a GRIDobj used to produce the provided ksn data
	%	FD - Flow routing as a FLOWobj used to proudce the provided ksn data
	%	ksn - ksn data either as a shapefile (as ouput from KsnProfiler, ProcessRiverBasins
	%		KsnChiBatch) or a mapstructure (as output from ProcessRiverBasins or KsnChiBatch)
	% 
	% Optional Inputs:
	% 	Can provide knickpoint locations as either an array (as output from FindBasinKnicks
	%		or KsnProfiler 'bnd_list') or shapefile (as output by FindBasinKnicks or KsnProfiler)
	%
	% Examples:
	%	PlotKsn(DEMoc,FDc,MSNc); % Plot ksn map of a basin from ProcessRiverBasins
	%	PlotKsn(DEM,FD,'ksn.shp'); % Plot ksn map from a shapefile
	%	PlotKsn(DEM,FD,'ksn.shp','ksn_knicks.shp'); % Include knickpoint locations from KsnProfiler
	%	PlotKsn(DEMoc,FDc,MSNc,KnickPoints); %Include knickpoints output from FindBasinKnicks
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 06/18/18 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	if ischar(ksn) & logical(regexp(ksn,regexptranslate('wildcard','*.shp')))
		ksn=shaperead(ksn);
	elseif isstruct(ksn)
		ksn=ksn;
	else
		error('Input to "ksn" not recognized as a shapefile or a mapstructure');
	end
		
	num_seg=numel(ksn);

	sx=cell(num_seg,1);
	sy=cell(num_seg,1);
	sk=cell(num_seg,1);
	for ii=1:num_seg
		sx{ii,1}=ksn(ii,1).X(:);
		sy{ii,1}=ksn(ii,1).Y(:);
		sk{ii,1}=ones(numel(sx{ii,1}),1)*ksn(ii,1).ksn;
	end

	sx=vertcat(sx{:});
	sy=vertcat(sy{:});
	sk=vertcat(sk{:});

	ix=coord2ind(DEM,sx,sy);
	idx=isnan(ix);

	ix(idx)=[];
	sk(idx)=[];

	W=GRIDobj(DEM,'logical');
	W.Z(ix)=true;
	S=STREAMobj(FD,W);

	[~,loc,~]=unique(ix);
	sk=sk(loc);

	f1=figure(1);
	set(f1,'Visible','off');

	[RGB]=imageschs(DEM,DEM,'colormap','gray');
	[~,R]=GRIDobj2im(DEM);

	imshow(flipud(RGB),R);
	axis xy
	hold on
	colormap(KsnColor(20));
	plotc(S,sk);
	caxis([0 max(sk)]);
	c1=colorbar;
	ylabel(c1,'Normalized Channel Steepness')
	if numel(varargin)==1
		knks=varargin{1};
		if ischar(knks) & logical(regexp(knks,regexptranslate('wildcard','*.shp')))
			knk=shaperead(knks);
			knkx=[knk.X];
			knky=[knk.Y];
			scatter(knkx,knky,100,'w','p','filled','MarkerEdgeColor','k');
		else
			scatter(knks(:,1),knk(:,2),100,'w','p','filled','MarkerEdgeColor','k');
		end
	end

	hold off
	set(f1,'Visible','on','Units','normalized','Position',[0.05 0.1 0.8 0.8],'renderer','painters');
end