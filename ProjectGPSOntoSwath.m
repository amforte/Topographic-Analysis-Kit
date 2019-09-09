function [ds,db,mag,unc,nc0,ec0]=ProjectGPSOntoSwath(SW,x,y,data_width,nc,ec,nu,eu)
	% Usage:
	%	[ds,db,mag,unc,nc0,ec0]=ProjectGPSOntoSwath(SW,x,y,data_width,nc,ec,nu,eu);
	%
	% Description:
	%	This function is a special version of ProjectOntoSwath for GPS velocity data that 
	%	also calculates the magnitude (north and east positive) and uncertainty of provided 
	%	GPS vectors in the direction of the swath. 
	%
	% Required Inputs:
	%	SW - SWATHobj onto which you want to project the GPS data
	%	x - nx1 array of x coordinates of GPS stations to project
	%	y - nx1 array of y coordinates of GPS stations to project
	%	data_width - width (in map units) you wish to sample, measured from swath baseline
	%	nc - north component of velocity
	%	ec - east component of velocity
	%	nu - uncertainty of the north velocity
	%	eu - uncertainty of the east velocity
	%
	% Outputs:
	% 	ds - distance along swath of station point
	% 	db - distance from baseline of swath of station point
	% 	mag - magnitude of vector along swath line, sign indicates direction of projected
	%		vector with respect to the direction of the swath. Positive values indicate that 
	%		projected vector is pointed in the direction that the swath is drawn (i.e. in the
	%		direction that the swath distance increases), negavie values indicate the projected
	%		vector points in the opposite direction
	% 	unc - uncertainty along swath line (slice of error ellipse)
	% 	nc0 - north component of projected vector (suitable for use matlab quiver plot)
	% 	ec0 - east component of projected vector (suitable for use matlab quiver plot)
	%
	% Examples:
	%	[ds,db,mag,unc,nc0,ec0]=ProjectGPSOntoSwath(SW,x,y,10000,nc,ec,nu,eu);
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 04/02/19 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

	mag_on_swath=zeros(numel(x),num_segs);
	unc_on_swath=zeros(numel(x),num_segs);
	nc0_on_swath=zeros(numel(x),num_segs);
	ec0_on_swath=zeros(numel(x),num_segs);

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

		% Do GPS vector calcs
		[mag_seg,unc_seg,ec0_seg,nc0_seg]=AngVecSw(swx0,swy0,swx1,swy1,nc,ec,nu,eu);
		% Filter gps and load out
		mag_seg(idx)=NaN;
		unc_seg(idx)=NaN;
		ec0_seg(idx)=NaN;
		nc0_seg(idx)=NaN;
		mag_on_swath(:,ii)=mag_seg;
		unc_on_swath(:,ii)=unc_seg;
		ec0_on_swath(:,ii)=ec0_seg;
		nc0_on_swath(:,ii)=nc0_seg;		

		dist_in_swath(:,ii)=pd_in_seg(:)+bends(ii);
		dist_from_base(:,ii)=DistFromBaseLine(:);
	end

	% Find distances of points that are closest to their segment line
	[db,c]=min(dist_from_base,[],2,'omitnan');
	r=[1:numel(c)]; r=r(:);
	ix=sub2ind(size(dist_from_base),r,c);
	ds=dist_in_swath(ix);
	mag=mag_on_swath(ix);
	unc=unc_on_swath(ix);
	ec0=ec0_on_swath(ix);
	nc0=nc0_on_swath(ix);

	% Set any points greater than the total swath distance to NaN;
	idx=single(ds)>=max(swdist);
	db(idx)=NaN;
	ds(idx)=NaN;
	mag(idx)=NaN;
	unc(idx)=NaN;
	ec0(idx)=NaN;
	nc0(idx)=NaN;
end

function [n_x,n_y]=RotCoord(x,y,theta,x0,y0)
	n_x=(x-x0).*cos(theta)-(y-y0).*sin(theta);
	n_y=(x-x0).*sin(theta)+(y-y0).*cos(theta);
end

function [proj_mag,unc,ec0,nc0]=AngVecSw(swx0,swy0,swx1,swy1,nc,ec,nu,eu)
	% Orient swath segment in N quadrant
	if swx0>swx1 & swy0>=swy1
		OX=swx1; OY=swy1; HX=swx0; HY=swy0;
	elseif swx0<swx1 & swy0<=swy1
		OX=swx0; OY=swy0; HX=swx1; HY=swy1;
	elseif swx0>=swx1 & swy0<swy1
		OX=swx0; OY=swy0; HX=swx1; HY=swy1;
	elseif swx0<=swx1 & swy0>swy1
		OX=swx1; OY=swy1; HX=swx0; HY=swy0;
	end

	% Find directional angles of swath and gps vector
	[swT,~]=cart2pol(HX-OX,HY-OY);
	[pT,~]=cart2pol(ec,nc);

	% Find magnitude of vector projected onto swath
	theta=swT-pT;
	mag=hypot(ec,nc);
	proj_mag=mag.*cos(theta);

	% Decompose projected vector into north and east components 
	[ec0,nc0]=pol2cart(swT,proj_mag);

	% Find uncertainty by finding radius of error ellipse at the angle of
	% the swath line
	unc=ellipserad(eu,nu,swT);

end

function [r]=ellipserad(a,b,theta);
	% a major axis
	% b minor axis
	% angle with respect to major axis

	r=(a.*b)./sqrt((b.*cos(theta)).^2 + (a.*sin(theta)).^2);
end
		