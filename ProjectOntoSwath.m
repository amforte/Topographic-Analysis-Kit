function [ds,db]=ProjectOntoSwath(SW,x,y,data_width,varargin)
	%
	% Usage:
	%	[ds,db]=ProjectOntoSwath(SW,x,y,data_width);
	%	[ds,db]=ProjectOntoSwath(SW,x,y,data_width,'name',value);	
	%
	% Description:
	% 	Function projects points on a SWATHobj (SW) and finds distance along (ds) and 
	%	from center line (db) of the SWATHobj, used in 'MakeCombinedSwath'
	%
	% Required Inputs:
	%	SW - SWATHobj
	%	x - xdata to project
	%	y - ydata to project
	%	data_width - width (in map units) from centerline of swath to include
	%
	% Optional Inputs:
	%	signed [false] - flag to control whether distance from the baseline
	%		"db" is signed, if false, all values will be positive. If true,
	%		values will be signed where positive values indicate distances to
	%		the right of the swath line (in the direction in which the swath
	%		nodes are defined) and negative values indicate distances to the 
	%		left of the swath line.
	%	include_concave_bend_regions [true] - flag to control whether
	%		the code will assign distances for points that lie
	%		within the triangular region within the concave bends
	%		 of a swath.
	%
	% Outputs:
	%	ds - distance along swath of each x-y point
	%	db - distance from center line in map units of each x-y point
	% 
	% Notes:
	% 	Values for points that do not project onto swath are set to NaN.
	%	Swaths with tight bends may produce unexpected results, especially when
	%		calculating signed distances.
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 09/19/20 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'ClassifyKnicks';
	addRequired(p,'SW',@(x) isa(x,'SWATHobj'));
	addRequired(p,'x',@(x) isnumeric(x));
	addRequired(p,'y',@(x) isnumeric(x));
	addRequired(p,'data_width',@(x) isnumeric(x) && isscalar(x));

	addParameter(p,'signed',false,@(x) islogical(x) && isscalar(x));
	addParameter(p,'include_concave_bend_regions',true,@(x) islogical(x) && isscalar(x));

	parse(p,SW,x,y,data_width,varargin{:});
	SW=p.Results.SW;
	x=p.Results.x;
	y=p.Results.y;
	data_width=p.Results.data_width;

	signed=p.Results.signed;
	in_cn_bnds=p.Results.include_concave_bend_regions;

	% Extract parameters from SWATHobj
	xypoints=SW.xy;
	swdist=SW.distx;
	xy0=SW.xy0;

	try
		for kk=1:numel(SW.xy0(:,1))
			[~,bend_ix(kk,1)]=min(pdist2(SW.xy,SW.xy0(kk,:)));
		end
	catch 
		% Less efficient method that doesn't require Statistics & Machine Learning Toolbox
		for kk=1:numel(SW.xy0(:,1))
			d=zeros(numel(SW.xy(:,1)),1);
			for ll=1:numel(SW.xy(:,1))
				d(ll)=hypot(SW.xy(ll,1)-SW.xy0(kk,1),SW.xy(ll,2)-SW.xy0(kk,2));
			end
			[~,bend_ix(kk,1)]=min(d);
		end
	end

	% Extract bend points, find number of segments, and find bend distances
	swxB=xy0(:,1); swyB=xy0(:,2);
	num_segs=numel(swxB)-1;

	if num_segs>1
		kk=1;
		while kk<=num_segs
			x0=swxB(kk,1);
			y0=swyB(kk,1);
			x1=swxB(kk+1,1);
			y1=swyB(kk+1,1);
			xx=x1-x0;
			yy=y1-y0;
			dist_to_bend(kk,1)=sqrt((xx^2)+(yy^2));
			kk=kk+1;
		end
		dist_to_bend=vertcat(0,dist_to_bend);
		bends=cumsum(dist_to_bend);
	else
		bends=[0;max(swdist)];
		dist_to_bend=bends;
	end

	% All points along swath line
	swx=xypoints(:,1); swy=xypoints(:,2);

	dist_in_swath=zeros(numel(x),num_segs);
	dist_from_base=zeros(numel(x),num_segs);

	for ii=1:num_segs
		% Find start and stop of segment
		swx0=swxB(ii); swx1=swxB(ii+1);
		swy0=swyB(ii); swy1=swyB(ii+1);

		% Find distances within segment
		x_dist=swx-swx0;
		y_dist=swy-swy0;
		seg_dist=sqrt((x_dist.^2)+(y_dist.^2));
		if ii>1 & ii<num_segs
			seg_dist([1:bend_ix(ii) bend_ix(ii+1)+1:end])=NaN;
		elseif ii==1
			seg_dist(bend_ix(ii+1)+1:end)=NaN;
		else
			seg_dist(1:bend_ix(ii))=NaN;
		end

		% Find distances of points from end node
		xn_dist=x-swx1;
		yn_dist=y-swy1;
		n_dist=sqrt((xn_dist.^2)+(yn_dist.^2));

		% Find angle of segment
		sw_angle=-1*atan((swy1-swy0)/(swx1-swx0));

		% Rotate all points in swath and dataset
		[n_swx,n_swy]=RotCoord(swx,swy,sw_angle,swx0,swy0);
		[n_x,n_y]=RotCoord(x,y,sw_angle,swx0,swy0);

		% Loop through all the dataset points
		num_points=numel(n_x);
		pd_in_seg=zeros(num_points,1);
		DistFromBaseLine=zeros(num_points,1);
		for jj=1:num_points
			xoi=n_x(jj);

			% Find the dataset point distance in the segment distance
			% Controlling for points beyond the end of the segment
			if abs(xoi)<=max(abs(n_swx))
				distances=n_swx-xoi;
				[~,I]=min(abs(distances)); 
				pd_in_seg(jj)=seg_dist(I);
				% Find distance from the baseline
				BaseLine=n_swy(I);
				if signed
					DistFromBaseLine(jj)=BaseLine-n_y(jj);
				else
					DistFromBaseLine(jj)=abs(BaseLine-n_y(jj));
				end

			else 
				pd_in_seg(jj)=NaN;
				DistFromBaseLine(jj)=NaN;
			end

			if in_cn_bnds
				% Control for points that are within bend overlaps
				if isnan(pd_in_seg(jj)) & n_dist(jj)<=data_width & ii~=num_segs
					pd_in_seg(jj)=seg_dist(bend_ix(ii+1));
					if signed
						% Determine sign
						[n_swx1,n_swy1]=RotCoord(swx1,swy1,sw_angle,swx0,swy0);
						if n_swy1>n_y(jj)
							DistFromBaseLine(jj)=n_dist(jj);
						else 
							DistFromBaseLine(jj)=-1*n_dist(jj);
						end
					else
						DistFromBaseLine(jj)=n_dist(jj);
					end
				end
			end 

		end

		% Index for points that project to segment
		idx=isnan(pd_in_seg) | pd_in_seg<=0;	
		pd_in_seg(idx)=NaN;
		DistFromBaseLine(idx)=NaN;

		dist_in_swath(:,ii)=pd_in_seg(:)+bends(ii);
		dist_from_base(:,ii)=DistFromBaseLine(:);
	end

	% Find distances of points that are closest to their segment line
	% [db,c]=min(abs(dist_from_base),[],2,'omitnan');
	[~,c]=min(abs(dist_from_base),[],2,'omitnan');

	r=[1:numel(c)]; r=r(:);
	ix=sub2ind(size(dist_from_base),r,c);
	ds=dist_in_swath(ix);

	db=dist_from_base(ix);

	% Set any points greater than the total swath distance to NaN;
	idx=single(ds)>=max(swdist);
	db(idx)=NaN;
	ds(idx)=NaN;
end

function [n_x,n_y]=RotCoord(x,y,theta,x0,y0)
	n_x=(x-x0).*cos(theta)-(y-y0).*sin(theta);
	n_y=(x-x0).*sin(theta)+(y-y0).*cos(theta);
end
