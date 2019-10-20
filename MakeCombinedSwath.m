function [SW,SwathMat,xypoints,outData]=MakeCombinedSwath(DEM,points,width,data_type,data,data_width,varargin)
	%
	% Usage:
	%	[SW,SwathMat,xypoints,outData]=MakeCombinedSwath(DEM,points,width,data_type,data,data_width);
	%	[SW,SwathMat,xypoints,outData]=MakeCombinedSwath(DEM,points,width,data_type,data,data_width,'name',value,...);
	%
	% Description:
	% 	Function to plot various additional data onto a swath profile.
	%
	% Required Inputs:
	% 	DEM - DEM Grid Object with which to make topo swath
	% 	points - n x 2 matrix containing x,y points for swath, minimum are two points (start and end points).
	%		First row contains starting point and proceeds down rows, additional points besides a start and end are
	%		treated as bends in the swath. Coordinates for points must be in the same coordinate system as DEM and must
	%		lie within the DEM (cannot be coordinates on the very edge of the DEM).If you provide an empty array this
	%		will invoke the SWATHobj behavior to display an image of the DEM on which to choose points.
	% 	width - width of swath in map units
	% 	data_type - the type of additional data you are providing to plot along with the swath, supported inputs are:
	%		'points3' - generic point dataset, expects a n x 3 matrix with values of x, y, and z
	%		'points4' - generic point dataset, expects a n x 4 matrix with values of x, y, z, and extra value. Dots will 
	%					be colored by this extra value
	%		'points5' - generic point dataset, expects a n x 5 matrix with values of x, y, z, and two extra values. Dots will
	%					colored by the first extra value (column 4) and scaled by the second extra value (column 5).
	%		'eqs' - earthquakes, expects a n x 4 matrix with x, y, depth, and magnitude. Points will be scaled by magnitude 
	%					and colored by distance from swath line. Expects depth to be positive.
	%		'gps' - gps velocity vectors, expects a n x 6 matrix, with x, y, north component, east component, north uncertainty,
	%					and east uncertainty. See 'ProjectGPSOntoSwath' for additional details
	%		'STREAMobj' - will project portions of selected stream profiles (as points) onto a swath. Expects a STREAMobj 
	%					that was generated from the provided DEM.
	%		'ksn_chandata' - will plot swath through ksn values, expects a chandata file as output from old Profiler51 code 
	%					(just in case you have some sitting around)
	%		'ksn_batch' - will plot swath through ksn values, expects the map structure output from 'KsnChiBatch' function 
	%					(i.e. run 'KsnChiBatch' with product set to 'ksn' and 'output' set to true, and provide the second 
	%					output here as 'data', i.e. run like [~,data]=KsnChiBatch(DEM,FD,A,S,'ksn'); and provide data here) 
	%		'ksn_profiler' - will plot swath through ksn values, expects the 'knl' output from 'KsnProfiler' function
	%		'basin_stats' - will plot swath through selected mean basin values as calculated from 'ProcessRiverBasins', 
	%					expects output from 'CompileBasinStats' and requires an entry to optional input 'basin_value' and  
	%					accepts optional input to 'basin_scale'. Will place point for basin at mean elevation and projected  
	%					location of the basin centroid, will color by value provided to 'basin_value' and will optionall scale  
	%					the point by the value provided to 'basin_scale'
	%		'basin_knicks' - will plot swath through knickpoints as chosen by 'FindBasinKnicks'. For 'data' provide name of folder 
	%					(or file path) where to find knickpoint files saved as a result of running 'FindBasinKnicks' on a series of  
	%					basins selected from 'ProcessRiverBasins'
	%	data - input data, form varies depending on choice of data_type
	% 	data_width - width in map units of swath through provided data. Values greater than data_width/2 from the center line 
	%					of the toposwath will not be plotted
	%
	% Optional Inputs:
	%	small_circ_center [] - option to provide a 1 x 2 array that contains the x and y coordinate of a small circle center to use
	%				to project data onto the swath, using the function 'ProjectSmallCircleOntoSwath'.
	%	dist_type ['mapunits'] - option to control how the 'data_width' is interepreted. Options are 'mapunits' or 'angle' with the
	%				default being 'mapunits'. The 'angle' option is only valid if an entry is provided to 'small_circ_center' to initiate
	%				projection along small circles. 
	% 	sample [] - resampling distance along topographic swath in map units, if no input is provided, code will use the cellsize 
	%				of the DEM which results in no resampling.
	% 	smooth [0] - smoothing distance, width of filter in map units over which to smooth values, default (0) results in no smoothing
	%	vex [10] - vertical exaggeration for the topographic swath. Note that because matlabs controls on physical axis dimensions are
	%				problematic, the vertical exaggeration controls don't work on plots that have two panels (e.g. 'ksn_batch', 'ksn_profiler',
	%				'ksn_chandata', and 'eqs')
	%	basin_value [] - required for option 'basin_stats', name (as it appears in the provided table provided to 'data') of the value 
	%				you wish to color points by
	%	basin_scale [] - optional input for option 'basin_stats', name (as it appears in the provided table provided to 'data') of the 
	%				value you wish to scale points by
	%	plot_map [true] - logical flag to plot a map displaying the location of the topographic swath and the additional data included 
	%				in the swath (red dots) and those not (white dots) based on the provided data_width parameter.
	%	cmap ['parula'] - valid name of colormap (e.g. 'jet') or a nx3 colormap array to use to color points.
	%	save_figure [false] - logical flag to save the swath figure as a pdf
	%
	% Outputs:
	%	SW - TopoToolbox Swath object, contains various information as a structure. Can plot path and box of swath with plot(SW) and
	%		plot version of swath profile with plotdz(SW);
	% 	SwathMat - n x 4 matrix containing distance along the swath, min elevation, mean elevation, max elevation
	% 	xypoints - n x 2 matrix containing x,y points of each swath sample point, along swath center line
	% 	outData - data for plotting the swath through the provided data, distances that area 'NaN' indicate those data do not
	%			fall on the swath line provided. Form of output depends on data_type:
	%		'points3' - distances, elevation, distance from base line, x coordinate, y coordinate
	%		'points4' - distances, elevation, value, distance from base line, x coordinate, y coordinate
	%		'eqs' - distances, depth, magnitude, distance from base line, x coordinate, y coordinate
	%		'STREAMobj' - distances, elevation, distance from base line, x coordinate, y coordinate
	%		'ksn_chandata' - distances, elevation, ksn, distance from base line, x coordinate, y coordinate
	%		'ksn_batch' - distances, ksn, distance from base line, x coordinate, y coordinate
	%		'ksn_profiler' - distances, ksn, distance from base line, x coordinate, y coordinate
	%		'basin_stats' - distances, mean basin elevation, 'basin_value', 'basin_scale' (if provided), distance from base line, 
	%						x coordinate, y coordinate
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 06/18/18 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'MakeCombinedSwath';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'points',@(x) isempty(x) || isnumeric(x) & size(x,1)>=2 && size(x,2)==2);
	addRequired(p,'width',@(x) isscalar(x) && isnumeric(x));
	addRequired(p,'data_type',@(x) ischar(validatestring(x,{'points3','points4','points5','eqs','gps','STREAMobj','ksn_chandata','ksn_batch','ksn_profiler','basin_stats','basin_knicks'})));
	addRequired(p,'data');
	addRequired(p,'data_width',@(x) isnumeric(x) && isscalar(x));

	addParameter(p,'file_name_prefix','Combined',@(x) ischar(x));
	addParameter(p,'small_circ_center',[],@(x) isnumeric(x) && numel(x)==2 || isempty(x));
	addParameter(p,'dist_type','mapdist',@(x) ischar(validatestring(x,{'mapdist','angle'})));	
	addParameter(p,'sample',[],@(x) isscalar(x) && isnumeric(x) || isempty(x));
	addParameter(p,'smooth',0,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'vex',10,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'basin_value',[],@(x) ischar(x));
	addParameter(p,'basin_scale',[],@(x) ischar(x) || isempty(x));
	addParameter(p,'plot_map',true,@(x) isscalar(x) && islogical(x));
	addParameter(p,'cmap','parula',@(x) ischar(x) || isnumeric(x) & size(x,2)==3);
	addParameter(p,'save_figure',false,@(x) isscalar(x) && islogical(x));


	parse(p,DEM,points,width,data_type,data,data_width,varargin{:});
	DEM=p.Results.DEM;
	points=p.Results.points;
	wdth=p.Results.width;
	data_type=p.Results.data_type;
	data=p.Results.data;
	data_width=p.Results.data_width;

	fnp=p.Results.file_name_prefix;
	small_circ_center=p.Results.small_circ_center;
	dist_type=p.Results.dist_type;
	sample=p.Results.sample;
	smth=p.Results.smooth;
	vex=p.Results.vex;
	bv=p.Results.basin_value;
	bs=p.Results.basin_scale;
	plot_map=p.Results.plot_map;
	cmap=p.Results.cmap;
	save_figure=p.Results.save_figure;

	if isempty(sample)
		sample=DEM.cellsize;
	end

	if isempty(small_circ_center)
		proj_flag=1;
	else
		proj_flag=2;
		cx=small_circ_center(1);
		cy=small_circ_center(2);
	end

	% Produce topo swath and associated datasets
	if isempty(points)
		SWt=SWATHobj(DEM);
		points=SWt.xy0;
		fig=gcf;
		close(fig);
	end
	
	[SW,SwathMat,xypoints,bends]=MakeTopoSwath(DEM,points,wdth,'sample',sample,'smooth',smth,'make_shape',false);
	swdist=SwathMat(:,1);
	min_elevs=SwathMat(:,2);
	mean_elevs=SwathMat(:,3);
	max_elevs=SwathMat(:,4);

	% Get extents of DEM
	[demx,demy]=getoutline(DEM,true);	

	% Set colormap
	colormap(cmap);

	% Perform different procedures depending on the type of data provided

	switch data_type
		case 'points3'
			x_coord=data(:,1);
			y_coord=data(:,2);
			z=data(:,3);

			% Remove any points beyond the extent of the provided DEM
			[demix]=inpolygon(x_coord,y_coord,demx,demy);	
			x_coord=x_coord(demix); y_coord=y_coord(demix); z=z(demix);

			switch proj_flag
			case 1
				[ds,db]=ProjectOntoSwath(SW,x_coord,y_coord,data_width);
				outData=[ds z db x_coord y_coord];
				idx=outData(:,3)<=(data_width/2) & ~isnan(ds);
			case 2
				[ds,db,dab]=ProjectSmallCircleOntoSwath(SW,x_coord,y_coord,cx,cy);
				switch dist_type
				case 'mapdist'
					outData=[ds z db x_coord y_coord];
					idx=abs(outData(:,3))<=(data_width/2) & ~isnan(ds);
				case 'angle'
					outData=[ds z dab x_coord y_coord];
					idx=abs(outData(:,3))<=(data_width/2) & ~isnan(ds);
				end
			end

			% Plot Swath
			f1=figure(1);
			clf 
			set(f1,'Units','normalized','Position',[0.05 0.1 0.8 0.4],'renderer','painters');

			hold on
			plot(swdist,min_elevs,'-k');
			plot(swdist,max_elevs,'-k');
			plot(swdist,mean_elevs,'-k','LineWidth',2);

			daspect([vex 1 1])

			yl=ylim;
			for jj=1:numel(bends)
				plot([bends(jj),bends(jj)],yl,'-k');
			end

			scatter(outData(idx,1),outData(idx,2),30,'k','filled');

			xlabel('Distance along swath (m)');
			ylabel('Elevation (m)');
			xlim([0 max(swdist)]);
			hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
		case 'points4'
			x_coord=data(:,1);
			y_coord=data(:,2);
			z=data(:,3);
			col=data(:,4);

			% Remove any points beyond the extent of the provided DEM
			[demix]=inpolygon(x_coord,y_coord,demx,demy);	
			x_coord=x_coord(demix); y_coord=y_coord(demix); z=z(demix); col=col(demix);

			switch proj_flag
			case 1
				[ds,db]=ProjectOntoSwath(SW,x_coord,y_coord,data_width);
				outData=[ds z col db x_coord y_coord];
				idx=outData(:,4)<=(data_width/2) & ~isnan(ds);
			case 2
				[ds,db,dab]=ProjectSmallCircleOntoSwath(SW,x_coord,y_coord,cx,cy);
				switch dist_type
				case 'mapdist'
					outData=[ds z col db x_coord y_coord];
					idx=abs(outData(:,4))<=(data_width/2) & ~isnan(ds);
				case 'angle'
					outData=[ds z col dab x_coord y_coord];
					idx=abs(outData(:,4))<=(data_width/2) & ~isnan(ds);
				end
			end			

			% Plot Swath
			f1=figure(1);
			clf 
			set(f1,'Units','normalized','Position',[0.05 0.1 0.8 0.4],'renderer','painters');

			hold on
			plot(swdist,min_elevs,'-k');
			plot(swdist,max_elevs,'-k');
			plot(swdist,mean_elevs,'-k','LineWidth',2);

			daspect([vex 1 1])

			yl=ylim;
			for jj=1:numel(bends)
				plot([bends(jj),bends(jj)],yl,'-k');
			end

			scatter(outData(idx,1),outData(idx,2),30,outData(idx,3),'filled');

			c1=colorbar;
			xlabel(c1,'User Provided Value')

			xlabel('Distance along swath (m)');
			ylabel('Elevation (m)');
			xlim([0 max(swdist)]);
			hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
		case 'points5'
			x_coord=data(:,1);
			y_coord=data(:,2);
			z=data(:,3);
			col=data(:,4)
			scle=data(:,5);

			% Remove any points beyond the extent of the provided DEM
			[demix]=inpolygon(x_coord,y_coord,demx,demy);	
			x_coord=x_coord(demix); y_coord=y_coord(demix); z=z(demix); col=col(demix); scle=scle(demix);

			switch proj_flag
			case 1
				[ds,db]=ProjectOntoSwath(SW,x_coord,y_coord,data_width);
				outData=[ds z col scle db x_coord y_coord];
				idx=outData(:,5)<=(data_width/2) & ~isnan(ds);
			case 2
				[ds,db,dab]=ProjectSmallCircleOntoSwath(SW,x_coord,y_coord,cx,cy);
				switch dist_type
				case 'mapdist'
					outData=[ds z col scle db x_coord y_coord];
					idx=abs(outData(:,5))<=(data_width/2) & ~isnan(ds);
				case 'angle'
					outData=[ds z col scle dab x_coord y_coord];
					idx=abs(outData(:,5))<=(data_width/2) & ~isnan(ds);
				end
			end

			% Plot Swath
			f1=figure(1);
			clf 
			set(f1,'Units','normalized','Position',[0.05 0.1 0.8 0.4],'renderer','painters');

			hold on
			plot(swdist,min_elevs,'-k');
			plot(swdist,max_elevs,'-k');
			plot(swdist,mean_elevs,'-k','LineWidth',2);

			daspect([vex 1 1])

			% Scale size vector
			sz_val=outData(idx,4);
			sz=(sz_val/max(sz_val))*100;
			% Create size legend
			sz_sizes=linspace(min(sz),max(sz),5);
			sz_val_scale=(sz_sizes/100)*max(sz_val);
			for ii=1:5
				sz_leg(ii)=plot(0,0,'ko','MarkerSize',sqrt(sz_sizes(ii)),'MarkerFaceColor','k','LineStyle','none');
				set(sz_leg(ii),'visible','off');
				leg_ent{ii}=num2str(sz_val_scale(ii));
			end	

			yl=ylim;
			for jj=1:numel(bends)
				plot([bends(jj),bends(jj)],yl,'-k');
			end

			scatter(outData(idx,1),outData(idx,2),sz,outData(idx,3),'filled');
			c1=colorbar('southoutside');
			xlabel(c1,'User Provided Value 1')
			legend(sz_leg,leg_ent);
			xlabel('Distance along swath (m)');
			ylabel('Elevation (m)');
			xlim([0 max(swdist)]);
			hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
		case 'eqs'
			x_coord=data(:,1);
			y_coord=data(:,2);
			depth=data(:,3);
			magnitude=data(:,4);

			% Remove any points beyond the extent of the provided DEM
			[demix]=inpolygon(x_coord,y_coord,demx,demy);	
			x_coord=x_coord(demix); y_coord=y_coord(demix); depth=depth(demix); magnitude=magnitude(demix);

			switch proj_flag
			case 1
				[ds,db]=ProjectOntoSwath(SW,x_coord,y_coord,data_width);
				outData=[ds depth magnitude db x_coord y_coord];
				idx=outData(:,4)<=(data_width/2) & ~isnan(ds);
			case 2
				[ds,db,dab]=ProjectSmallCircleOntoSwath(SW,x_coord,y_coord,cx,cy);
				switch dist_type
				case 'mapdist'
					outData=[ds depth magnitude db x_coord y_coord];
					idx=abs(outData(:,4))<=(data_width/2) & ~isnan(ds);
				case 'angle'
					outData=[ds depth magnitude dab x_coord y_coord];
					idx=abs(outData(:,4))<=(data_width/2) & ~isnan(ds);
				end
			end				

			% Plot Swath
			f1=figure(1);
			clf 
			set(f1,'Units','normalized','Position',[0.05 0.1 0.8 0.4],'renderer','painters');

			ax1=subplot(2,1,1);
			hold on
			plot(swdist,min_elevs,'-k');
			plot(swdist,max_elevs,'-k');
			plot(swdist,mean_elevs,'-k','LineWidth',2);

			yl=ylim;
			for jj=1:numel(bends)
				plot([bends(jj),bends(jj)],yl,'-k');
			end
			ylabel('Elevation (m)');
			xlim([0 max(swdist)]);
			hold off

			ax2=subplot(2,1,2);
			hold on
			% Scale size vector
			sz_val=outData(idx,3);
			sz=(sz_val/max(sz_val))*100;
			% Create size legend
			sz_sizes=linspace(min(sz),max(sz),5);
			sz_val_scale=(sz_sizes/100)*max(sz_val);
			for ii=1:5
				sz_leg(ii)=plot(0,0,'ko','MarkerSize',sqrt(sz_sizes(ii)),'MarkerFaceColor','k','LineStyle','none');
				set(sz_leg(ii),'visible','off');
				leg_ent{ii}=num2str(sz_val_scale(ii));
			end		

			scatter(outData(idx,1),outData(idx,2),sz,outData(idx,4),'filled');
			xlabel('Distance along swath (m)');
			ylabel('Depth (km')
			legend(sz_leg,leg_ent);
			xlim([0 max(swdist)]);
			c1=colorbar(ax2,'southoutside');
			xlabel(c1,'Distance from Swath Line')
			hold off

			set(ax2,'YDir','reverse');
			linkaxes([ax1,ax2],'x')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
		case 'gps'
			x_coord=data(:,1);
			y_coord=data(:,2);
			nc=data(:,3);
			ec=data(:,4);
			nu=data(:,5);
			eu=data(:,6);

			% Remove any points beyond the extent of the provided DEM
			[demix]=inpolygon(x_coord,y_coord,demx,demy);	
			x_coord=x_coord(demix); y_coord=y_coord(demix); nc=nc(demix); ec=ec(demix); nu=nu(demix); eu=eu(demix);

			switch proj_flag
			case 1
				[ds,db,mag,unc,nc0,ec0]=ProjectGPSOntoSwath(SW,x_coord,y_coord,data_width,nc,ec,nu,eu);
				outData=[ds mag unc nc0 ec0 db x_coord y_coord];
				idx=outData(:,6)<=(data_width/2) & ~isnan(ds);
			case 2
				[~,~,mag,unc,nc0,ec0]=ProjectGPSOntoSwath(SW,x_coord,y_coord,data_width,nc,ec,nu,eu);
				[ds,db,dab]=ProjectSmallCircleOntoSwath(SW,x_coord,y_coord,cx,cy);
				switch dist_type
				case 'mapdist'
					outData=[ds mag unc nc0 ec0 db x_coord y_coord];
					idx=abs(outData(:,6))<=(data_width/2) & ~isnan(ds);
				case 'angle'
					outData=[ds mag unc nc0 ec0 dab x_coord y_coord];
					idx=abs(outData(:,6))<=(data_width/2) & ~isnan(ds);
				end
			end

			% Plot Swath
			f1=figure(1);
			clf 
			set(f1,'Units','normalized','Position',[0.05 0.1 0.8 0.4],'renderer','painters');

			ax1=subplot(2,1,1);
			hold on
			plot(swdist,min_elevs,'-k');
			plot(swdist,max_elevs,'-k');
			plot(swdist,mean_elevs,'-k','LineWidth',2);

			yl=ylim;
			for jj=1:numel(bends)
				plot([bends(jj),bends(jj)],yl,'-k');
			end
			ylabel('Elevation (m)');
			xlim([0 max(swdist)]);
			hold off

			ax2=subplot(2,1,2);
			hold on	
			e1=errorbar(outData(idx,1),outData(idx,2),outData(idx,3),'.');
			e1.CapSize=0;
			e1.Color='k';
			scatter(outData(idx,1),outData(idx,2),30,outData(idx,6),'filled');
			xlabel('Distance along swath (m)');
			ylabel('Velocity in Swath Line (mm/yr)');
			xlim([0 max(swdist)]);
			c1=colorbar(ax2,'southoutside');
			xlabel(c1,'Distance from Swath Line')
			hold off

			linkaxes([ax1,ax2],'x')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		case 'STREAMobj'
			x_coord=data.x;
			y_coord=data.y;
			z=getnal(data,DEM);

			% Remove any points beyond the extent of the provided DEM
			[demix]=inpolygon(x_coord,y_coord,demx,demy);	
			x_coord=x_coord(demix); y_coord=y_coord(demix); z=z(demix);

			switch proj_flag
			case 1
				[ds,db]=ProjectOntoSwath(SW,x_coord,y_coord,data_width);
				outData=[ds z db x_coord y_coord];
				idx=outData(:,3)<=(data_width/2) & ~isnan(ds);
			case 2
				[ds,db,dab]=ProjectSmallCircleOntoSwath(SW,x_coord,y_coord,cx,cy);
				switch dist_type
				case 'mapdist'
					outData=[ds z db x_coord y_coord];
					idx=abs(outData(:,3))<=(data_width/2) & ~isnan(ds);
				case 'angle'
					outData=[ds z dab x_coord y_coord];
					idx=abs(outData(:,3))<=(data_width/2) & ~isnan(ds);
				end
			end

			% Plot Swath
			f1=figure(1);
			clf 
			set(f1,'Units','normalized','Position',[0.05 0.1 0.8 0.4],'renderer','painters');

			hold on
			plot(swdist,min_elevs,'-k');
			plot(swdist,max_elevs,'-k');
			plot(swdist,mean_elevs,'-k','LineWidth',2);

			daspect([vex 1 1])

			yl=ylim;
			for jj=1:numel(bends)
				plot([bends(jj),bends(jj)],yl,'-k');
			end

			scatter(outData(idx,1),outData(idx,2),5,'k','filled');

			xlabel('Distance along swath (m)');
			ylabel('Elevation (m)');
			xlim([0 max(swdist)]);
			hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		case 'ksn_chandata'
			x_coord=data(:,9);
			y_coord=data(:,10);
			ksn=data(:,8);
			elev=data(:,4);

			% Remove any points beyond the extent of the provided DEM
			[demix]=inpolygon(x_coord,y_coord,demx,demy);	
			x_coord=x_coord(demix); y_coord=y_coord(demix); ksn=ksn(demix); elev=elev(demix);

			switch proj_flag
			case 1
				[ds,db]=ProjectOntoSwath(SW,x_coord,y_coord,data_width);
				outData=[ds elev ksn db x_coord y_coord];
				idx=outData(:,4)<=(data_width/2) & ~isnan(ds);
			case 2
				[ds,db,dab]=ProjectSmallCircleOntoSwath(SW,x_coord,y_coord,cx,cy);
				switch dist_type
				case 'mapdist'
					outData=[ds elev ksn db x_coord y_coord];
					idx=abs(outData(:,4))<=(data_width/2) & ~isnan(ds);
				case 'angle'
					outData=[ds elev ksn dab x_coord y_coord];
					idx=abs(outData(:,4))<=(data_width/2) & ~isnan(ds);
				end
			end			

			% Plot Swath
			f1=figure(1);
			clf 
			set(f1,'Units','normalized','Position',[0.05 0.1 0.8 0.4],'renderer','painters');

			subplot(2,1,1);
			hold on

			plot(swdist,min_elevs,'-k');
			plot(swdist,max_elevs,'-k');
			plot(swdist,mean_elevs,'-k','LineWidth',2);

			yl=ylim;
			for jj=1:numel(bends)
				plot([bends(jj),bends(jj)],yl,'-k');
			end

			xlabel('Distance along swath (m)');
			ylabel('Elevation (m)');
			xlim([0 max(swdist)]);
			hold off

			subplot(2,1,2);
			hold on
			scatter(outData(idx,1),outData(idx,3),30,'k','filled');
			xlabel('Distance along swath (m)');
			ylabel('KSN');
			xlim([0 max(swdist)]);
			hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		case 'basin_stats'
			x_coord=data.center_x;
			y_coord=data.center_y;
			z=data.mean_el;
			col=data.(bv);
			if ~isempty(bs)
				scl=data.(bs);
				if ~isnumeric(scl)
					error('Value to scale points by must be numeric')
				end
			end

			if ~isnumeric(col);
				error('Value to color points by must be numeric')
			end

			% Remove any points beyond the extent of the provided DEM
			[demix]=inpolygon(x_coord,y_coord,demx,demy);	
			x_coord=x_coord(demix); y_coord=y_coord(demix); z=z(demix); col=col(demix);
			if ~isempty(bs)
				scl=scl(demix);
			end

			switch proj_flag
			case 1
				% Transform Data
				[ds,db]=ProjectOntoSwath(SW,x_coord,y_coord,data_width);

				% Assemble outData
				if isempty(bs)
					outData=[ds z col db x_coord y_coord];
				else 
					outData=[ds z col scl db x_coord y_coord];	
				end			

				% Filter based on provided data width
				if isempty(bs)
					idx=outData(:,4)<=(data_width/2) & ~isnan(ds);
				else
					idx=outData(:,5)<=(data_width/2) & ~isnan(ds);
				end
			case 2
				[ds,db,dab]=ProjectSmallCircleOntoSwath(SW,x_coord,y_coord,cx,cy);

				switch dist_type
				case 'mapdist'
					% Assemble outData
					if isempty(bs)
						outData=[ds z col db x_coord y_coord];
					else 
						outData=[ds z col scl db x_coord y_coord];	
					end			

					% Filter based on provided data width
					if isempty(bs)
						idx=abs(outData(:,4))<=(data_width/2) & ~isnan(ds);
					else
						idx=abs(outData(:,5))<=(data_width/2) & ~isnan(ds);
					end
				case 'angle'
					% Assemble outData
					if isempty(bs)
						outData=[ds z col dab x_coord y_coord];
					else 
						outData=[ds z col scl dab x_coord y_coord];	
					end			

					% Filter based on provided data width
					if isempty(bs)
						idx=abs(outData(:,4))<=(data_width/2) & ~isnan(ds);
					else
						idx=abs(outData(:,5))<=(data_width/2) & ~isnan(ds);
					end
				end
			end			


			% Plot Swath
			f1=figure(1);
			clf 
			set(f1,'Units','normalized','Position',[0.05 0.1 0.8 0.4],'renderer','painters');

			hold on
			plot(swdist,min_elevs,'-k');
			plot(swdist,max_elevs,'-k');
			plot(swdist,mean_elevs,'-k','LineWidth',2);

			daspect([vex 1 1])

			yl=ylim;
			for jj=1:numel(bends)
				plot([bends(jj),bends(jj)],yl,'-k');
			end

			if isempty(bs)
				scatter(outData(idx,1),outData(idx,2),30,outData(idx,3),'filled');
			else
				% Scale size vector
				sz_val=outData(idx,4);
				sz=(sz_val/max(sz_val))*100;
				% Create size legend
				sz_sizes=linspace(min(sz),max(sz),5);
				sz_val_scale=(sz_sizes/100)*max(sz_val);
				for ii=1:5
					sz_leg(ii)=plot(0,0,'ko','MarkerSize',sqrt(sz_sizes(ii)),'MarkerFaceColor','k','LineStyle','none');
					set(sz_leg(ii),'visible','off');
					leg_ent{ii}=num2str(sz_val_scale(ii));
				end
				scatter(outData(idx,1),outData(idx,2),sz,outData(idx,3),'filled');
				legend(sz_leg,leg_ent);
				bs_n=strrep(bs,'_',' ');
				title(['Points scaled by ' bs_n]);
			end

			c1=colorbar;
			bv_n=strrep(bv,'_',' ');
			ylabel(c1,bv_n);

			xlabel('Distance along swath (m)');
			ylabel('Elevation (m)');
			xlim([0 max(swdist)]);
			hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		case 'ksn_batch'
			% Find number of segments in provided KSN map structure	
			numSegs=numel(data);
			% Loop through stream segments to extract x,y,ksn
			streamData=zeros(numSegs,3);
			for kk=1:numSegs
				xx=mean(data(kk,1).X);
				yy=mean(data(kk,1).Y);
				ksn=mean(data(kk,1).ksn);
				streamData(kk,:)=[xx yy ksn];
			end	

			x_coord=streamData(:,1);
			y_coord=streamData(:,2);
			ksn=streamData(:,3);

			% Remove any points beyond the extent of the provided DEM
			[demix]=inpolygon(x_coord,y_coord,demx,demy);	
			x_coord=x_coord(demix); y_coord=y_coord(demix); ksn=ksn(demix);

			switch proj_flag
			case 1
				[ds,db]=ProjectOntoSwath(SW,x_coord,y_coord,data_width);
				outData=[ds ksn db x_coord y_coord];
				idx=outData(:,3)<=(data_width/2) & ~isnan(ds);
			case 2
				[ds,db,dab]=ProjectSmallCircleOntoSwath(SW,x_coord,y_coord,cx,cy);
				switch dist_type
				case 'mapdist'
					outData=[ds ksn db x_coord y_coord];
					idx=abs(outData(:,3))<=(data_width/2) & ~isnan(ds);
				case 'angle'
					outData=[ds ksn dab x_coord y_coord];
					idx=abs(outData(:,3))<=(data_width/2) & ~isnan(ds);
				end
			end			

			% Plot Swath
			f1=figure(1);
			clf 
			set(f1,'Units','normalized','Position',[0.05 0.1 0.8 0.4],'renderer','painters');

			subplot(2,1,1);
			hold on

			plot(swdist,min_elevs,'-k');
			plot(swdist,max_elevs,'-k');
			plot(swdist,mean_elevs,'-k','LineWidth',2);

			yl=ylim;
			for jj=1:numel(bends)
				plot([bends(jj),bends(jj)],yl,'-k');
			end

			xlabel('Distance along swath (m)');
			ylabel('Elevation (m)');
			xlim([0 max(swdist)]);
			hold off

			subplot(2,1,2);
			hold on
			scatter(outData(idx,1),outData(idx,2),20,'k','filled');
			xlabel('Distance along swath (m)');
			ylabel('KSN');
			xlim([0 max(swdist)]);
			hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%			
		case 'ksn_profiler'
			x_coord=data(:,1);
			y_coord=data(:,2);
			ksn=data(:,4);

			% Remove any points beyond the extent of the provided DEM
			[demix]=inpolygon(x_coord,y_coord,demx,demy);	
			x_coord=x_coord(demix); y_coord=y_coord(demix); ksn=ksn(demix);			

			switch proj_flag
			case 1
				[ds,db]=ProjectOntoSwath(SW,x_coord,y_coord,data_width);
				outData=[ds ksn db x_coord y_coord];
				idx=outData(:,3)<=(data_width/2) & ~isnan(ds);
			case 2
				[ds,db,dab]=ProjectSmallCircleOntoSwath(SW,x_coord,y_coord,cx,cy);
				switch dist_type
				case 'mapdist'
					outData=[ds ksn db x_coord y_coord];
					idx=abs(outData(:,3))<=(data_width/2) & ~isnan(ds);
				case 'angle'
					outData=[ds ksn dab x_coord y_coord];
					idx=abs(outData(:,3))<=(data_width/2) & ~isnan(ds);
				end
			end				

			% Plot Swath
			f1=figure(1);
			clf 
			set(f1,'Units','normalized','Position',[0.05 0.1 0.8 0.4],'renderer','painters');

			subplot(2,1,1);
			hold on

			plot(swdist,min_elevs,'-k');
			plot(swdist,max_elevs,'-k');
			plot(swdist,mean_elevs,'-k','LineWidth',2);

			yl=ylim;
			for jj=1:numel(bends)
				plot([bends(jj),bends(jj)],yl,'-k');
			end

			xlabel('Distance along swath (m)');
			ylabel('Elevation (m)');
			xlim([0 max(swdist)]);
			hold off

			subplot(2,1,2);
			hold on
			scatter(outData(idx,1),outData(idx,2),20,'k','filled');
			xlabel('Distance along swath (m)');
			ylabel('KSN');
			xlim([0 max(swdist)]);
			hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		case 'basin_knicks'

			if ~isdir(data)
				error('For "basin_knicks" you must provide a valid directory to "data" that contains "Knicks_*.mat" files')
			end

			current=pwd;
			cd(data);

			fileList=dir('Knicks_*.mat');
			if isempty(fileList)
				error('For "basin_knicks" you must provide a valid directory to "data" that contains "Knicks_*.mat" files')
			end
			knps=cell(numel(fileList),1);
			for jj=1:numel(fileList)
				load(fileList(jj,1).name);
				knps{jj}=[KnickTable.x_coord KnickTable.y_coord KnickTable.elevation];
			end

			knps=vertcat(knps{:});

			x_coord=knps(:,1);
			y_coord=knps(:,2);
			z=knps(:,3);

			% Remove any points beyond the extent of the provided DEM
			[demix]=inpolygon(x_coord,y_coord,demx,demy);	
			x_coord=x_coord(demix); y_coord=y_coord(demix); z=z(demix);

			switch proj_flag
			case 1
				[ds,db]=ProjectOntoSwath(SW,x_coord,y_coord,data_width);
				outData=[ds z db x_coord y_coord];
				idx=outData(:,3)<=(data_width/2) & ~isnan(ds);
			case 2
				[ds,db,dab]=ProjectSmallCircleOntoSwath(SW,x_coord,y_coord,cx,cy);
				switch dist_type
				case 'mapdist'
					outData=[ds z db x_coord y_coord];
					idx=abs(outData(:,3))<=(data_width/2) & ~isnan(ds);
				case 'angle'
					outData=[ds z dab x_coord y_coord];
					idx=abs(outData(:,3))<=(data_width/2) & ~isnan(ds);
				end
			end				

			% Plot Swath
			f1=figure(1);
			clf 
			set(f1,'Units','normalized','Position',[0.05 0.1 0.8 0.4],'renderer','painters');

			hold on
			plot(swdist,min_elevs,'-k');
			plot(swdist,max_elevs,'-k');
			plot(swdist,mean_elevs,'-k','LineWidth',2);

			daspect([vex 1 1])

			yl=ylim;
			for jj=1:numel(bends)
				plot([bends(jj),bends(jj)],yl,'-k');
			end

			scatter(outData(idx,1),outData(idx,2),20,'k','filled');

			xlabel('Distance along swath (m)');
			ylabel('Elevation (m)');
			xlim([0 max(swdist)]);
			hold off
	end

	% Make output shape
	ms=struct;
	ms(1,1).Geometry='Line';
	ms(1,1).X=SW.xy0(:,1);
	ms(1,1).Y=SW.xy0(:,2);
	ms(1,1).Type='Center';

	if ~verLessThan('matlab','9.4')
		Tverts=SwathPolygon(SW,wdth);
		ms(2,1).Geometry='Line';
		ms(2,1).X=Tverts(:,1);
		ms(2,1).Y=Tverts(:,2);
		ms(2,1).Type='TopoWdth';

		Dverts=SwathPolygon(SW,data_width);
		ms(3,1).Geometry='Line';
		ms(3,1).X=Dverts(:,1);
		ms(3,1).Y=Dverts(:,2);
		ms(3,1).Type='DataWdth';
	end

	shapewrite(ms,'SwathBounds.shp');

	if plot_map
		f2=figure(2);
		set(f2,'Units','normalized','Position',[0.05 0.1 0.6 0.6]);
		hold on
		imageschs(DEM,DEM,'colormap','gray');
		plot(SW.xy0(:,1),SW.xy0(:,2),'-g','LineWidth',0.5);
		if ~verLessThan('matlab','9.4')
			plot(Tverts(:,1),Tverts(:,2),'-g','LineWidth',0.5);
			plot(Dverts(:,1),Dverts(:,2),'-r','LineWidth',0.5);	
		end	
		scatter(x_coord(idx),y_coord(idx),20,'r','filled');
		% To avoid situations with extremely large number of points
		if ~any([strcmp(data_type,'STREAMobj') strcmp(data_type,'ksn_profiler') strcmp(data_type,'ksn_batch')])
			scatter(x_coord(~idx),y_coord(~idx),20,'w','filled');
		end
		hold off
	end

	if save_figure
		orient(f1,'Landscape')
		print(f1,'-dpdf','-bestfit',[fnp '_Swath.pdf']);
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









