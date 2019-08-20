function [ds,db,dab]=ProjectSmallCircleOntoSwath(SW,x,y,cx,cy)
	% 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% THIS FUNCTION IS UNDER ACTIVE DEVELOPMENT %%%
	%% FUNCTION IS NOT FULLY DOCUMENTED IN MANUAL %%%
	%%%%%% CHANGES AND MALFUNCTIONS ARE LIKELY %%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%	
	% Usage:
	%	[ds,db,dab]=ProjectSmallCircleOntoSwath(SW,x,y,cx,cy,proj);
	%
	% Description:
	%	Function projects provided data along a family of small circles with a 
	%	common center point.
	%
	% Required Inputs:
	%	SW - SWATHobj
	%	x - nx1 array of x data to project 
	%	y - nx1 array of y data to project
	%	cx - x coordinate of center of small circles
	%	cy - y coordinate of center of small circles
	%	
	% Outputs:
	%	ds - nx1 array of distances along provided swath of projected data. Points with
	%		NaN as distances indicate that projected position of the point does not lie along
	%		the swath line
	%	db - nx1 array of distances (in map units) from the baseline of the swath. Distances
	%		will be signed.
	%	dab - nx1 array of angular distances (in radians) from the baseline along the small 
	%		circle of projection. Distances will be signed.
	%
	% Examples:
	%	[ds,db,dab]=ProjectSmallCircleOntoSwath(SW,x,y,cx,cy,DEM.georef.mstruct);
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 04/02/19 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Extract parameters from SWATHobj
	xypoints=SW.xy;
	swdist=SW.distx;
	proj=SW.georef.mstruct;

	% Calculate small circles in projected coordinates
	[dlat,dlon]=projinv(proj,x,y);
	[clat,clon]=projinv(proj,cx,cy);
	[sclat,sclon]=scircle2(clat,clon,dlat,dlon,[],[],500);
	[scx,scy]=projfwd(proj,sclat,sclon);

	% Find distance on swath
	d=zeros(size(xypoints,1),numel(x));
	w1=waitbar(0,'Finding interesctions along small circles...');
	for ii=1:size(xypoints,1)
		d(ii,:)=min(hypot(xypoints(ii,1)-scx,xypoints(ii,2)-scy));
		waitbar(ii/size(xypoints,1));
	end
	close(w1);
	[~,ix]=min(d,[],1);
	ds=swdist(ix);

	% Find arc length between swath center line and point
	xsw=xypoints(ix,1);
	ysw=xypoints(ix,2);
	[tsw,rsw]=cart2pol(xsw-cx,ysw-cy);
	[tp,rp]=cart2pol(x-cx,y-cy);
	dab=tsw-tp;
	db=dab.*rp;


	% Check for points outside of swath if the ix was either the beginning or the end
	% and set points outside swath to NaN
	if any(ix==1) | any(ix==size(xypoints,1))
		[lat1,lon1]=projinv(proj,xypoints(end,1),xypoints(end,2));
		[lat0,lon0]=projinv(proj,xypoints(1,1),xypoints(1,2));
		[eslat,eslon]=scircle2(clat,clon,vertcat(lat1,lat0),vertcat(lon1,lon0));
		[esx,esy]=projfwd(proj,eslat,eslon);
		% Define complex polygon
		esx=vertcat(esx(:,1),NaN,flipud(esx(:,2)));
		esy=vertcat(esy(:,1),NaN,flipud(esy(:,2)));
		[inp,onp]=inpolygon(x,y,esx,esy);
		in = inp | onp;
		ds(~in)=NaN;
		db(~in)=NaN;
	end
end