function [n_x,n_y]=ReFrameRot(SWATH,x,y)
	% Function reprojects coordinates points into a new coordinate system with 
	% a line between the endpoints of a swath profile being the new x-axis 
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Last Revised Fall 2015 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Find Angle
	xypoints=cell2mat(SWATH.xy);
	start_p=xypoints(1,:);
	stop_p=xypoints(length(xypoints),:);

	x_dist=stop_p(:,1)-start_p(:,1);
	y_dist=stop_p(:,2)-start_p(:,2);

	sw_angle=-1*atan(y_dist/x_dist);

	start_e=start_p(1,1);
	start_n=start_p(1,2);

	% Do rotation
	e=x; n=y;
	e_p=(e-start_e).*cos(sw_angle)-(n-start_n).*sin(sw_angle);
	n_p=(e-start_e).*sin(sw_angle)+(n-start_n).*cos(sw_angle);
	n_x=e_p; n_y=n_p;
end