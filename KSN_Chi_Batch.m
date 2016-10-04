function [varargout]=KSN_Chi_Batch(DEM,FD,A,threshold_area,smooth_distance,min_ksn,ref_concavity,product,output)
	% Function to produce channel steepness, chi maps or chi grids for all channels within a DEM
	% 
	% Reqiured Inputs:
	% DEM - DEM Grid Object (assumes unconditioned DEM)
	% FD - FLOW object
	% A - GRID object of flow accumulations
	% threshold_area - area in sqaure of map units above which is considered a stream
	% smooth_distance - smoothing distance in map units for smoothing ksn values, equivalent to smoothing in Profiler
	% min_ksn - minimum ksn for calculating minimum gradients at different drainage areas 
	% ref_concavity - reference concavity (as a positive value) for calculating ksn
	% product - switch to determine which products to produce
	%		'ksn' - ksn map as a shapefile
	%		'chimap' - ascii file with chi calculated in channel networks
	%		'chigrid' - ascii file with chi calculate at all points in a grid
	%		'chi' - results for both chimap and chigrid
	%		'all' - ksn, chimap, and chigrids
	% output - switch to either output matlab files to the workspace (true) or to not only save the specified files
	%		without any workspace output (false)
	%
	% Notes:
	%	Please be aware that the production of the chigrid can be extremely time consuming, so be patient...
	%
	% Example:
	%	KSN_Chi_Batch(DEM,FD,A,1e6,1000,20,0.45,'chi',false);
	%	[KSN,ChiMap,ChiGrid]=KSN_Chi_Batch(DEM,FD,A,1e6,1000,20,0.45,'all',true);
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Last Revised Spring 2016 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	switch product
	case 'ksn'

		disp('Calculating channel steepness')
		[G,Gmin] = CompGrad(DEM,FD,A,min_ksn,ref_concavity);

		DEM_res=DEM.cellsize;
		min_area=floor(threshold_area/(DEM_res*DEM_res));
		S=STREAMobj(FD,'minarea',min_area);		

		ksn=G./(A.*(A.cellsize^2)).^(-ref_concavity);
		KSN=STREAMobj2mapstruct(S,'seglength',smooth_distance,'attributes',...
			{'ksn' ksn @mean 'uparea' (A.*(A.cellsize^2)) @mean 'gradient' G @mean 'gmin' Gmin @mode});

		disp('Writing ARC files')
		shapewrite(KSN,'ksn.shp');

		switch output
		case true
			varargout{1}=KSN;
		end


	case 'chimap'

	    disp('Calculating chi map');
		[ChiMap]=MakeChiMap(DEM,FD,A,threshold_area,ref_concavity);

		disp('Writing ARC files')
		GRIDobj2ascii(ChiMap,'chimap.txt');

		switch output
		case true
			varargout{1}=ChiMap;
		end

	case 'chigrid'

	    disp('Calculating chi grid');
	    [ChiGrid]=MakeChiGrid(FD,A,ref_concavity);

	    disp('Writing ARC files')
		GRIDobj2ascii(ChiGrid,'chigrid.txt');

		switch output
		case true
			varargout{1}=ChiGrid;
		end

	case 'chi'

	    disp('Calculating chi map');
		[ChiMap]=MakeChiMap(DEM,FD,A,threshold_area,ref_concavity);

	    disp('Calculating chi grid');
	    [ChiGrid]=MakeChiGrid(FD,A,ref_concavity);

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

		DEM_res=DEM.cellsize;
		min_area=floor(threshold_area/(DEM_res*DEM_res));
		S=STREAMobj(FD,'minarea',min_area);		

		ksn=G./(A.*(A.cellsize^2)).^(-ref_concavity);
		KSN=STREAMobj2mapstruct(S,'seglength',smooth_distance,'attributes',...
			{'ksn' ksn @mean 'uparea' (A.*(A.cellsize^2)) @mean 'gradient' G @mean 'gmin' Gmin @mode});

	    disp('Calculating chi map');
		[ChiMap]=MakeChiMap(DEM,FD,A,threshold_area,ref_concavity);

	    disp('Calculating chi grid');
	    [ChiGrid]=MakeChiGrid(FD,A,ref_concavity);


		disp('Writing ARC files')
		shapewrite(KSN,'ksn.shp');
		GRIDobj2ascii(ChiMap,'chimap.txt');
		GRIDobj2ascii(ChiGrid,'chigrid.txt');

		switch output
		case true
			varargout{1}=KSN;
			varargout{2}=ChiMap;
			varargout{3}=ChiGrid;
		end
	end

% Main Function End
end

function [ChiOBJ]=MakeChiMap(DEM,FD,A,cutoff_drainage_area,ref_concavity);

	% Function to produce a chi map for all channels within a DEM

	S=STREAMobj(FD,'minarea',cutoff_drainage_area,'unit','mapunits');
	IX=streampoi(S,'outlets','ix');

	num_outs=numel(IX);

	ChiCell=cell(num_outs,1);
	for ii=1:num_outs
		ix=IX(ii);
		Sc=STREAMobj(FD,'outlets',ix,'minarea',cutoff_drainage_area,'unit','mapunits');
		Sc=klargestconncomps(Sc,1); % Make double sure there is only one outlet
		ChiOI=chiplot(Sc,DEM,A,'a0',1,'mn',ref_concavity,'plot',false);

		% Grab x, y, and chi
		cx=ChiOI.x;
		cy=ChiOI.y;
		cc=ChiOI.chi;

		% Strip out NaNs
		idx=~isnan(cx);
		cx=cx(idx); cy=cy(idx); cc=cc(idx);

		cix=coord2ind(DEM,cx,cy);

		ChiCell{ii}=[cix cc];
	end

	% Concatenate indexes
	Chi=vertcat(ChiCell{:});
	chi_ix=Chi(:,1);
	chi_val=Chi(:,2);

	% Make Empty GRIDobj
	ChiOBJ=GRIDobj(DEM);
	idx=ChiOBJ.Z==0;
	ChiOBJ.Z(idx)=NaN;

	% Populate Grid
	ChiOBJ.Z(chi_ix)=chi_val;

end

function [ChiGrid]=MakeChiGrid(FD,A,ref_concavity)

	% Generate dense stream network
	S=STREAMobj(FD,'minarea',0);
	O=streampoi(S,'outlets','ix');
	a0=1;

	ixL=cell(numel(O),1);
	cL=cell(numel(O),1);
	for ii=1:numel(O);
		So=STREAMobj(FD,'minarea',0,'outlets',O(ii));
		a=double(a0./(A.Z(So.IXgrid)*(A.cellsize.^2)));
		x=So.distance;
		chi = netcumtrapz(x,a.^ref_concavity,So.ix,So.ixc);
		[xx,yy,c]=STREAMobj2XY(So,chi);
		idx=~isnan(xx);
		xx=xx(idx); yy=yy(idx); c=c(idx);
		ix=coord2ind(A,xx,yy);
		ixL{ii}=ix;
		cL{ii}=c;
	end

	cix=vertcat(ixL{:});
	chi=vertcat(cL{:});

	ChiGrid=GRIDobj(A);
	ChiGrid.Z(cix)=chi;

end

function z = netcumtrapz(x,y,ix,ixc)
% cumtrapz along upward direction in a directed tree network
% function from TopoToolbox chiplot by Wolfgang Schwanghart
	z = zeros(size(x));
	for lp = numel(ix):-1:1;
	    z(ix(lp)) = z(ixc(lp)) + (y(ixc(lp))+(y(ix(lp))-y(ixc(lp)))/2) *(abs(x(ixc(lp))-x(ix(lp))));
	end
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