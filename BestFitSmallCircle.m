function [clat,clon,crad,cpx,cpy]=BestFitSmallCircle(x,y,proj,varargin)
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
	% Optional Inputs:
	%	trim_circle [false] - logical flag to trim circle to the approximate extent of the provided x y points,
	%		performance of this will improve if you increase the num_points beyond the default of 100.
	%	num_points [100] - number of points to generate along the small circle that is extracted
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

    p = inputParser;
    p.FunctionName = 'BestFitSmallCircle';
    addRequired(p,'x',@(x) isnumeric(x));
    addRequired(p,'y',@(x) isnumeric(x));
    addRequired(p,'proj',@(x) isstruct(x));

    addParameter(p,'trim_circle',false,@(x) islogical(x) && isscalar(x));
    addParameter(p,'num_points',100,@(x) isscalar(x) && isnumeric(x));

    parse(p,x,y,proj,varargin{:});
    x=p.Results.x;
    y=p.Results.y;
    proj=p.Results.proj;

    trim_circle=p.Results.trim_circle;
    num_points=p.Results.num_points;

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
	[circ_lat,circ_lon]=scircle1(clat,clon,crad,[],[],[],num_points);
	if ~isempty(proj)
		[cpx,cpy]=projfwd(proj,circ_lat,circ_lon);
	else
		cpx=circ_lon; cpy=circ_lat;
	end

	% Trim circle to extent 
	if trim_circle
		x0=x(1); y0=y(1);
		x1=x(end); y1=y(end);

		d0=hypot((x0-cpx),(y0-cpy));
		d1=hypot((x1-cpx),(y1-cpy));

		[~,ix0]=min(d0);
		[~,ix1]=min(d1);

		if ix0<ix1
			cpx=cpx(ix0:ix1);
			cpy=cpy(ix0:ix1);
		elseif ix0>ix1
			cpx=cpx(ix1:ix0);
			cpy=cpy(ix1:ix0);
		end
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