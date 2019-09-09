function [cross_points]=CrossSwath(SW,ln)
	% Usage:
	%	[cross_points]=CrossSwath(SW,ln);
	%
	% Description:
	%	Function finds intersection of a provided line with the boundaries of 
	%	a SWATHobj
	%
	% Required Inputs:
	%	SW - SWATHobj to find the interesctions with
	%	ln - nx2 array of x,y points which define the line
	%
	% Output:
	%	cross_points - n x 7 array with a row for each separate portion of of the supplied line 
	%		that crosses the swath and columns which correspond to x position of first cross point, 
	%		y position of first cross point, nearest swath distance to first cross point, mean swath 
	%		distance for the cross section, x position of second cross point, y position of second 
	%		cross point, and nearest swath distance to second cross point.
	%
	% Examples:
	%	[cross_points]=CrossSwath(SW,ln);
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 04/02/19 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	% Find outlines of SWATH box
	try
		[verts]=SwathPolygon(SW,SW.width);
	catch
		bx=[SW.X(1,:) SW.X(:,end)' fliplr(SW.X(end,:)) fliplr(SW.X(:,1)')];
		by=[SW.Y(1,:) SW.Y(:,end)' fliplr(SW.Y(end,:)) fliplr(SW.Y(:,1)')];
		verts=horzcat(bx',by');
	end

	% Find cumulate runnind distances of line
	lxd=diff(ln(:,1)); lyd=diff(ln(:,2));
	dst=hypot(lxd,lyd);
	dst=cumsum(dst);
	dst=vertcat(0,dst);

	% Find portions of line within swath box
	[lx,ly]=interpline(ln(:,1),ln(:,2),dst,SW.dx/2);
	[in]=inpolygon(lx,ly,verts(:,1),verts(:,2));
	[lbl,nc]=bwlabel(in);

	% Process disconnected segments separately
	cross_points=zeros(nc,7);
	for ii=1:nc
		idx=lbl==ii;
		inx=lx(idx);
		iny=ly(idx);

		cross1=[inx(1) iny(1)];
		cross2=[inx(end) iny(end)];
		crossM=[mean(inx) mean(iny)];

		[~,ix1]=min(hypot(SW.xy(:,1)-cross1(1),SW.xy(:,2)-cross1(2)));
		[~,ix2]=min(hypot(SW.xy(:,1)-cross2(1),SW.xy(:,2)-cross2(2)));
		[~,ix3]=min(hypot(SW.xy(:,1)-crossM(1),SW.xy(:,2)-crossM(2)));		

		cross1=[cross1 SW.distx(ix1)];
		cross2=[cross2 SW.distx(ix2)];

		cross_points(ii,:)=[cross1 SW.distx(ix3) cross2];
	end
end

function [verts]=SwathPolygon(SW,w);

	cx=SW.xy(:,1);
	cy=SW.xy(:,2);

	dx=diff(cx);
	dy=diff(cy);

	w=w/2;

	[sw_angle,~]=cart2pol(dx,dy);
	sw_angle=vertcat(sw_angle,sw_angle(end));
	[px,py]=pol2cart(sw_angle+(pi/2),w);
	[mx,my]=pol2cart(sw_angle-(pi/2),w);

	swx=[cx+px cx+mx];
	swy=[cy+py cy+my];

	warning off
	P=polyshape(vertcat(swx(:,1),flipud(swx(:,2))),vertcat(swy(:,1),flipud(swy(:,2))));
	warning on
	P=rmholes(P);
	verts=P.Vertices;
	verts=vertcat(verts,verts(1,:));

end