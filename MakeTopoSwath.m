function [SW,SwathMat,xypoints,bends]=MakeTopoSwath(DEM,points,width,varargin)
	% Wrapper around TopoToolbox SWATHobj functionality
	%
	% Required Inputs:
	% 	DEM - DEM Grid Object with which to make topo swath
	% 	points - n x 2 matrix containing x,y points for swath, minimum are two points (start and end points).
	%		First row contains starting point and proceeds down rows, additional points besides a start and end are
	%		treated as bends in the swath. Coordinates for points must be in the same coordinate system as DEM and must
	%		lie within the DEM (cannot be coordinates on the very edge of the DEM)
	% 	width - width of swath in map units
	%
	% Optional Inputs:
	% 	sample [] - resampling distance along swath in map units, if no input is provided, code will use the cellsize of the DEM 
	%		which results in no resampling.
	% 	smooth [0] - smoothing distance, width of filter in map units over which to smooth values, default (0) results in no smoothing
	%	vex [10] - vertical exaggeration for displaying plot.
	% 	plot_figure [false] - logical flag to plot result. 
	%
	% Outputs:
	% 	SW - TopoToolbox Swath object, contains various information as a structure. Can plot path and box of swath with plot(SW) and
	%		plot version of swath profile with plotdz(SW);
	% 	SwathMat - n x 4 matrix containing distance along the swath, min elevation, mean elevation, max elevation
	% 	xypoints - n x 2 matrix containing x,y points of each swath sample point, along swath center line
	% 	bends - distances along swath of any bends, 0 if no bends
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Last Revised Spring 2018 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'MakeTopoSwath';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'points',@(x) isnumeric(x) && size(x,1)>=2 && size(x,2)==2);
	addRequired(p,'width',@(x) isscalar(x) && isnumeric(x));

	addParamValue(p,'sample',[],@(x) isscalar(x) && isnumeric(x));
	addParamValue(p,'smooth',0,@(x) isscalar(x) && isnumeric(x));
	addParamValue(p,'vex',10,@(x) isscalar(x) && isnumeric(x));
	addParamValue(p,'plot_figure',false,@(x) isscalar(x) && islogical(x));

	parse(p,DEM,points,width,varargin{:});
	DEM=p.Results.DEM;
	points=p.Results.points;
	wdth=p.Results.width;

	sample=p.Results.sample;
	smth=p.Results.smooth;
	vex=p.Results.vex;
	plot_figure=p.Results.plot_figure;

	if isempty(sample)
		sample=DEM.cellsize;
	end

	% Find Bend Points in Swath
	num_points=size(points,1);

	if num_points>2
		kk=1;
		while kk<num_points-1
			bx=points(kk,1);
			by=points(kk,2);
			ex=points(kk+1,1);
			ey=points(kk+1,2);
			xx=ex-bx;
			yy=ey-by;
			dist_to_bend(kk)=sqrt((xx^2)+(yy^2));
			kk=kk+1;
		end
		bends=cumsum(dist_to_bend);
	else
		bends=0;
	end

	% Make Swath 
	% Deal with changes in versions of TopoToolbox
	try
		SW=SWATHobj(DEM,points,'width',wdth,'dx',sample,'smooth',smth); % Older versions
	catch
		SW=SWATHobj(DEM,points(:,1),points(:,2),'width',wdth,'dx',sample,'smooth',smth); % Newer versions
	end

	% Extract useful values from swath object
	try
		elevs=cell2mat(SW.Z); % Old
	catch 
		elevs=SW.Z; % New
	end

	mean_elevs=nanmean(elevs);
	min_elevs=nanmin(elevs);
	max_elevs=nanmax(elevs);
	try
		xypoints=cell2mat(SW.xy); % Old
		swdist=cell2mat(SW.distx);
	catch
		xypoints=SW.xy; % New
		swdist=SW.distx;
	end

	SwathMat=[swdist min_elevs.' mean_elevs.' max_elevs.'];

	if plot_figure
		f1=figure(1);
		clf 
		set(f1,'Units','inches','Position',[1.0 1.5 16 4],'renderer','painters','PaperSize',[16 4],'PaperOrientation','portrait','PaperPositionMode','auto');

		hold on
		plot(swdist,min_elevs,'-b');
		plot(swdist,max_elevs,'-b');
		plot(swdist,mean_elevs,'-k');

		daspect([vex 1 1])

		yl=ylim;
		for jj=1:numel(bends)
			plot([bends(jj),bends(jj)],yl,'-k');
		end

		xlabel('Distance along swath (m)');
		ylabel('Elevation (m)');
		xlim([0 max(swdist)]);
		hold off
	end


end




