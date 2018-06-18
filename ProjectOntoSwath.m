function [ds,db]=ProjectOntoSwath(SW,x,y,data_width)
	%
	% Usage:
	%	[ds,db]=ProjectOntoSwath(SW,x,y,data_width);
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
	% Outputs:
	%	ds - distance along swath of each x-y point
	%	db - distance from center line in map units of each x-y point
	% 
	% Values for points that do not project onto swath are set to NaN
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 06/18/18 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Extract parameters from SWATHobj
	xypoints=SW.xy;
	swdist=SW.distx;
	xy0=SW.xy0;

	for kk=1:numel(SW.xy0(:,1))
		[~,bend_ix(kk,1)]=min(pdist2(SW.xy,SW.xy0(kk,:)));
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
				DistFromBaseLine(jj)=abs(BaseLine-n_y(jj));
			else 
				pd_in_seg(jj)=NaN;
				DistFromBaseLine(jj)=NaN;
			end

			% Control for points that are within bend overlaps
			if isnan(pd_in_seg(jj)) & n_dist(jj)<=data_width & ii~=num_segs
				pd_in_seg(jj)=seg_dist(bend_ix(ii+1));
				DistFromBaseLine(jj)=n_dist(jj);
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
	[db,c]=nanmin(dist_from_base,[],2);
	r=[1:numel(c)]; r=r(:);
	ix=sub2ind(size(dist_from_base),r,c);
	ds=dist_in_swath(ix);

	% Set any points greater than the total swath distance to NaN;
	idx=single(ds)>=max(swdist);
	db(idx)=NaN;
	ds(idx)=NaN;
end

function [n_x,n_y]=RotCoord(x,y,theta,x0,y0)
	n_x=(x-x0).*cos(theta)-(y-y0).*sin(theta);
	n_y=(x-x0).*sin(theta)+(y-y0).*cos(theta);
end
