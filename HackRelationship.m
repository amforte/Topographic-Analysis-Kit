function [C,h,drain_area,outs,globC,globh]=HackRelationship(DEM,FD,A,S,varargin)
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% THIS FUNCTION IS UNDER ACTIVE DEVELOPMENT %%%
	%% FUNCTION IS NOT FULLY DOCUMENTED IN MANUAL %%%
	%%%%%% CHANGES AND MALFUNCTIONS ARE LIKELY %%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% Usage: 
	%	[C,h,drain_area,outs,globC,globh]=HackRelationship(DEM,FD,A,S);
	%	[C,h,drain_area,outs,globC,globh]=HackRelationship(DEM,FD,A,S,'name',value,...);
	%
	% Description:
	%	Function to find the Hack relationships for all watersheds in a given landscape. Assumes
	%	that the provided STREAMobj terminates at the outlets of the watersheds for which you
	% 	wish to calculate Hack coeffecients and exponents.
	%
	% Required Inputs:
	%	DEM - Digital Elevation as a GRIDobj
	%	FD - Flow direction as FLOWobj
	%	A - Flow accumulation GRIDobj
	%	S - Stream network as STREAMobj 
	%
	% Optional Inputs:
	%	method ['trunks'] - optional parameter for controlling which data is used to 
	%		fit the Hack relationship. Options are 'trunks' (default), 'streams', or
	%		'grids'. If 'trunks', only values in trunk streams will be used. If 'streams'
	%		only values in stream networks (as defined by the input STREAMobj) will be
	%		used. If 'grids', all pixels in watersheds upstream of outlets of the provided
	%		STREAMobj will be used.
	%	relation ['original'] - optional parameter to control the form of the Hack relationship.
	%		The 'original' option will fit the original relationship as presented by Hack, e.g.
	%		length = C * drainage area ^ h. The 'inverse' will instead fit the form used by 
	%		Whipple & Tucker, 1999, drainage area = C * length ^ h.
	%	include_hillslope [false] - logical flag to either include or not include (default) the
	%		portions of the channels between the channelheads in the provided STREAMobj and the 
	%		drainage divide. Regardless of the value of this parameter, distances within the channel
	%		will be relative to the divide unless changed with 'measure_from' parameter, but if left false, 
	%		low drainage area values physically above channelhead will not be included in fits. This parameter 
	%		is ignored if 'method' is set to 'grids'.
	%	measure_from ['divide'] - optional parameter to control how distances ares measured, either measured from 
	%		the drainage divides when set to 'divide' or from the channelheads (as defined in the provided
	%		STREAMobj) if set to 'channelheads'. If 'method' is set to 'grids', then this parameter will be set to 'divide'
	%		regardless of user input.
	%	draw_fig [false] - logical flag to display a plot of the coefficient (C) and exponent (h) as
	%		a function of drainage area for all basins
	%
	% Outputs:
	%	C - n x 1 array of Hack coefficients for each watershed, there will be as many values of C as there
	%		are outlets in the provided STREAMobj.
	%	h - n x 1 array of Hack exponents for each watershed
	%	drainage_area - n x 1 array of total drainage areas (in squared map units) of each watershed
	%	outs - n x 1 array of indices of outlet locations, same as result of streampoi(S,'outlets','ix')
	%	globC - Hack coefficient of fit of all relevant points in the landscape. The points included in the
	%		this global fit will depend on the value of 'method'
	%	globh - Hack exponent of fit of all relevant points in the landscape
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 05/02/19 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	% Parse Inputs
	p=inputParser;
	p.FunctionName='HackRelationship';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'FD',@(x) isa(x,'FLOWobj'));
	addRequired(p,'A', @(x) isa(x,'GRIDobj'));
	addRequired(p,'S',@(x) isa(x,'STREAMobj'));

	addParameter(p,'method','trunks',@(x) ischar(validatestring(x,{'trunks','streams','grids'})));
	addParameter(p,'relation','original',@(x) ischar(validatestring(x,{'original','inverse'})));
	addParameter(p,'include_hillslope',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'measure_from',@(x) ischar(validatestring(x,{'divide','channelheads'})));
	addParameter(p,'draw_fig',false,@(x) isscalar(x) && islogical(x));

	parse(p,DEM,FD,A,S,varargin{:});
	DEM=p.Results.DEM;
	FD=p.Results.FD;
	A=p.Results.A;
	S=p.Results.S;

	method=p.Results.method;
	relation=p.Results.relation;
	include_hillslope=p.Results.include_hillslope;
	measure_from=p.Results.measure_from;
	draw_fig=p.Results.draw_fig;

	% Precalculate needed parameters
	outs=streampoi(S,'outlets','ix');
	FLDS=flowdistance(FD,'downstream');
	DA=A.*A.cellsize^2;
	DB=drainagebasins(FD,outs);
	num_basins=numel(outs);

	% Extract node attributed lists if method
	% is either 'trunks' or 'streams'
	switch method
	case 'trunks'
		S=trunk(S);
		% Recalculate streams if hillslope flag is thrown
		if include_hillslope
			FLUS=flowdistance(FD);
			chix=streampoi(S,'channelheads','ix');
			ix=zeros(numel(chix),1);
			for ii=1:numel(chix)
				chOI=chix(ii);

				UP=dependencemap(FD,chOI);
				FLUSt=FLUS.*UP;

				[~,ix(ii,1)]=max(FLUSt);
			end
			IX=influencemap(FD,ix);
			S=STREAMobj(FD,IX);
		end
		% Grab nal
		danal=getnal(S,DA);
		dnal=getnal(S,FLDS);
		dbnal=getnal(S,DB);
	case 'streams'
		% Recalculate streams if hillslope flag is thrown
		if include_hillslope
			FLUS=flowdistance(FD);
			chix=streampoi(S,'channelheads','ix');
			ix=zeros(numel(chix),1);
			for ii=1:numel(chix)
				chOI=chix(ii);

				UP=dependencemap(FD,chOI);
				FLUSt=FLUS.*UP;

				[~,ix(ii,1)]=max(FLUSt);
			end
			IX=influencemap(FD,ix);
			S=STREAMobj(FD,IX);
		end
		% Grab nal
		danal=getnal(S,DA);
		dnal=getnal(S,FLDS);
		dbnal=getnal(S,DB);
	end

	for ii=1:num_basins
		switch method
		case 'grids'
			IDX=DB==ii;
			da=DA.Z(IDX.Z);
			l=double(FLDS.Z(IDX.Z));
		case {'streams','trunks'}
			IDX=dbnal==ii;
			da=danal(IDX);
			l=double(dnal(IDX));
		end

		switch relation
		case 'original'
			f=fit(da,l,'power1');
		case 'inverse'
			nzidx=l>0;

			f=fit(l(nzidx),da(nzidx),'power1');
		end

		cf=coeffvalues(f);
		C(ii,1)=cf(1);
		h(ii,1)=cf(2);
		drain_area(ii,1)=max(da);
	end

	% Global fit
	switch relation
	case 'original'
		switch method
		case 'grids'
			globfit=fit(DA.Z(:),double(FLDS.Z(:)),'power1');
		case {'streams','trunks'}
			globfit=fit(danal,double(dnal),'power1');
		end
	case 'inverse'
		switch method
		case 'grids'
			flds=double(FLDS.Z(:));
			da=DA.Z(:);
			nzidx=flds>0;

			globfit=fit(flds(nzidx),da(nzidx),'power1');
		case {'streams','trunks'}
			nzidx=dnal>0;
			globfit=fit(double(dnal(nzidx)),danal(nzidx),'power1');
		end	
	end

	globcf=coeffvalues(globfit);
	globC=globcf(1);
	globh=globcf(2);

	if draw_fig
		f1=figure(1);
		clf

		subplot(2,1,1);
		hold on 
		scatter(drain_area,C,20,'k','filled');
		set(gca,'XScale','log');
		plot(xlim,[globC globC],'--b');
		plot(xlim,[mean(C) mean(C)],'--r');
		xlabel('Drainage Area');
		ylabel('C parameter');
		legend('Individual Basin Fits','All Pixel Fit','Mean of Basins','location','best');
		hold off

		subplot(2,1,2);
		hold on 
		scatter(drain_area,h,20,'k','filled');
		set(gca,'XScale','log');
		plot(xlim,[globh globh],'--b');
		plot(xlim,[mean(h) mean(h)],'--r');
		xlabel('Drainage Area');
		ylabel('h parameter');
		hold off
	end


end