function [clat,clon,crad,cpx,cpy]=BestFitSmallCircle(x,y,proj)
	% Usage:
	%	[clat,clon,crad,cpx,cpy]=BestFitSmallCircle(x,y,proj)
	%
	% Description:
	%	Function finds a best fit small circle based on a series of x,y coordinates.
	%
	% Required Inputs:
	%	x - nx1 array of x coordinates of points to fit
	%	y - nx1 array of y coordinates of points to fit
	%	proj - projection of input x,y coordinates, (e.g. data stored in DEM.georef.mstruct
	%	or DEM.georef for older versions of TopoToolbox). If an empty array is provided to 
	%	the 'proj' argument, it's assumed that the x and y coordinates are longitude and 
	%	latitude respectively. In this case, the outputs 'cpx' and 'cpy' will also be in
	%	longitude and latitude.
	%
	% Outputs:
	%	clat - latitude of center of small circle
	%	clon - longitude of center of small circle
	%	crad - radius (in degrees) of small circle
	%	cpx - x coordinates of circle perimeter (in projected coord)
	%	cpy - y coordinates of circle perimeter (in projected coord)
	%
	% Examples:
	%	[clat,clon,crad,cpx,cpy]=BestFitSmallCircle(x,y,DEM.georef.mstruct);
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 05/25/21 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Remove any NaNs
	idx=~isnan(x) | ~isnan(y);
	x=x(idx);
	y=y(idx);

	% Convert input x y points to lat lon
	if ~isempty(proj)
		[lat,lon]=projinv(proj,x,y);
	else
		lon=x; lat=y;
	end

	% Find a reasonable start point
	mlat=mean(lat);
	mlon=mean(lon);
	d=hypot(mlat-lat,mlon-lon);
	md=max(d);

	% Do minimization
	model=@bfc;
	x0=[mlat mlon md];
	est=fminsearch(model,x0);

	% Extract values
	clat=est(1);
	clon=est(2);
	crad=est(3);

	% Calculate circle perimeter in projected coordinates
	[circ_lat,circ_lon]=scircle1(clat,clon,crad);
	if ~isempty(proj)
		[cpx,cpy]=projfwd(proj,circ_lat,circ_lon);
	else
		cpx=circ_lon; cpy=circ_lat;
	end

	% Minimization function
	function [ssq]=bfc(params)
		p1=params(1);
		p2=params(2);
		p3=params(3);

		% Calculate circle
		[lat1,lon1]=scircle1(p1,p2,p3);

		% Find distance to nearest point on circle to each input point
		for ii=1:numel(lat)
			mind(ii,1)=min(hypot(lat1-lat(ii),lon1-lon(ii)));
		end

		ssq=sum(mind.^2);
	end
end