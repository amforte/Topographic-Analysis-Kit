function [SW,SwathMat,xypoints,outData]=MakeCombinedSwath(DEM,points,width,data_type,data,data_width,varargin)
	% Function requires TopoToolbox to run
	%
	% Required Inputs:
	% 	DEM - DEM Grid Object with which to make topo swath
	% 	points - n x 2 matrix containing x,y points for swath, minimum are two points (start and end points).
	%		First row contains starting point and proceeds down rows, additional points besides a start and end are
	%		treated as bends in the swath. Coordinates for points must be in the same coordinate system as DEM and must
	%		lie within the DEM (cannot be coordinates on the very edge of the DEM)
	% 	width - width of swath in map units
	% 	data_type - the type of additional data you are providing to plot along with the swath, supported inputs are:
	%		'points3' - generic point dataset, expects a n x 3 matrix with values of x, y, and z
	%		'points4' - generic point dataset, expects a n x 4 matrix with values of x, y, z, and extra value. Dots will be sclaed by this extra value
	%		'eqs' - earthquakes, expects a n x 4 matrix with x, y, depth, and magnitude. Points will be scaled by magnitude and colored by distance from swath line
	%		'STREAMobj' - will project portions of selected stream profiles (as points) onto a swath. Expects a STREAMobj that was generated from the provided DEM.
	%		'ksn_chandata' - will plot swath through ksn values, expects a chandata file as output from old Profiler51 code
	%		'ksn_batch' - will plot swath through ksn values, expects the map structure output from 'KsnChiBatch' function (i.e. run 'KsnChiBatch' with product set 
	%			to 'ksn' and 'output' set to true, and provide the second output here as 'data')
	%		'ksn_profiler' - will plot swath through ksn values, expects the 'knl' output from 'KsnProfiler' function.
	%		'basin_stats' - will plot swath through selected mean basin values as calculated from 'ProcessRiverBasins', expects output from 'CompileBasinStats' and 
	%			requires an entry to optional input 'basin_value' and accepts optional input to 'basin_scale'. Will place point for basin at mean elevation and projected
	%			location of the basin centroid, will color by value provided to 'basin_value' and will optionall scale the point by the value provided to 'basin_scale'
	%	data - input data, form varies depending on choice of data_type
	% 	data_width - width in map units of swath through provided data. Values greater than data_width/2 from the center line of the toposwath will be removed
	%
	% Optional Inputs:
	% 	sample [] - resampling distance along topographic swath in map units, if no input is provided, code will use the cellsize of the DEM 
	%		which results in no resampling.
	% 	smooth [0] - smoothing distance, width of filter in map units over which to smooth values, default (0) results in no smoothing
	%	vex [10] - vertical exaggeration for displaying plot.
	%	basin_value [] - required for option 'basin_stats', name (as it appears in the provided table provided to 'data') of the value you wish to color points by
	%	basin_scale [] - optional input for option 'basin_stats', name (as it appears in the provided table provided to 'data') of the value you wish to scale points by
	%
	% Outputs:
	% SW - TopoToolbox Swath object, contains various information as a structure. Can plot path and box of swath with plot(SW) and
	%	plot version of swath profile with plotdz(SW);
	% SwathMat - n x 4 matrix containing distance along the swath, min elevation, mean elevation, max elevation
	% xypoints - n x 2 matrix containing x,y points of each swath sample point, along swath center line
	% outData - data for plotting the swath through the provided data, form of output depends on data_type:
	%		'points3' - distances, elev, distance from base line, x coordinate, y coordinate, transformed x coordinate, transformed y coordinate
	%		'points4' - distances, elev, value, distance from base line, x coordinate, y coordinate, transformed x coordinate, transformed y coordinate
	%		'eqs' - distances, depth, magnitude, distance from base line, x coordinate, y coordinate, transformed x coordinate, transformed y coordinate
	%		'ksn_chandata' - distances, elev, ksn, distance from base line, x coordinate, y coordinate, transformed x coordinate, transformed y coordinate
	%		'ksn_TT' - distances, ksn, distance from base line, x coordinate, y coordinate, transformed x coordinate, transformed y coordinate
	%		'chi_TT' - distances, chi, distance from base line, x coordinate, y coordinate, transformed x coordinate, transformed y coordinate
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Last Revised Fall 2015 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'MakeCombinedSwath';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'points',@(x) isnumeric(x) && size(x,1)>=2 && size(x,2)==2);
	addRequired(p,'width',@(x) isscalar(x) && isnumeric(x));
	addRequired(p,'data_type',@(x) ischar(validatestring(x,{'points3','points4','eqs','STREAMobj','ksn_chandata','ksn_batch','ksn_profiler','basin_stats'})));
	addRequired(p,'data');
	addRequired(p,'data_width',@(x) isnumeric(x) && isscalar(x));

	addParamValue(p,'sample',[],@(x) isscalar(x) && isnumeric(x));
	addParamValue(p,'smooth',0,@(x) isscalar(x) && isnumeric(x));
	addParamValue(p,'vex',10,@(x) isscalar(x) && isnumeric(x));
	addParamValue(p,'basin_value',[],@(x) ischar(x));
	addParamValue(p,'basin_scale',[],@(x) ischar(x));


	parse(p,DEM,points,width,data_type,data,data_width,varargin{:});
	DEM=p.Results.DEM;
	points=p.Results.points;
	wdth=p.Results.width;
	data_type=p.Results.data_type;
	data=p.Results.data;
	data_width=p.Results.data_width;

	sample=p.Results.sample;
	smth=p.Results.smooth;
	vex=p.Results.vex;
	bv=p.Results.basin_value;
	bs=p.Results.basin_scale;

	if isempty(sample)
		sample=DEM.cellsize;
	end

	% Produce topo swath and associated datasets
	[SW,SwathMat,xypoints,bends]=MakeTopoSwath(DEM,points,wdth,'sample',sample,'smooth',smth);

	sw_y=xypoints(:,2);
	sw_x=xypoints(:,1);
	[nsw_x,nsw_y]=ReFrameRot(SW,sw_x,sw_y);
	swdist=SwathMat(:,1);
	min_elevs=SwathMat(:,2);
	mean_elevs=SwathMat(:,3);
	max_elevs=SwathMat(:,4);

	% Perform different procedures depending on the type of data provided

	switch data_type
		case 'points3'
			x_coord=data(:,1);
			y_coord=data(:,2);
			z=data(:,3);

			% Transform data
			[n_x,n_y]=ReFrameRot(SW,x_coord,y_coord);

			% Find distances along swath profile of each point
			num_points=numel(x_coord);
			point_distances=zeros(num_points,1);
			DistFromBaseLine=zeros(num_points,1);
			for ii=1:num_points
				xoi=n_x(ii);

				distances=nsw_x-xoi;
				[~,I]=min(abs(distances));
				point_distances(ii)=swdist(I);

				baseLine=nsw_y(I);
				DistFromBaseLine(ii)=abs(baseLine-n_y(ii));
			end
			outData=[point_distances z DistFromBaseLine x_coord y_coord n_x n_y];

			% Filter out any data points outside of provided width
			idx=DistFromBaseLine<=(data_width/2) & point_distances>0 & point_distances<max(swdist);
			outData=outData(idx,:);

			% Plot Swath
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

			scatter(outData(:,1),outData(:,2),20,'k','filled');

			xlabel('Distance along swath (m)');
			ylabel('Elevation (m)');
			xlim([0 max(swdist)]);
			hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
		case 'points4'
			x_coord=data(:,1);
			y_coord=data(:,2);
			z=data(:,3);
			scle=data(:,4);

			% Transform data
			[n_x,n_y]=ReFrameRot(SW,x_coord,y_coord);

			% Find distances along swath profile of each point
			num_points=numel(x_coord);
			point_distances=zeros(num_points,1);
			DistFromBaseLine=zeros(num_points,1);
			for ii=1:num_points
				xoi=n_x(ii);

				distances=nsw_x-xoi;
				[~,I]=min(abs(distances));
				point_distances(ii)=swdist(I);

				baseLine=nsw_y(I);
				DistFromBaseLine(ii)=abs(baseLine-n_y(ii));
			end
			outData=[point_distances z scle DistFromBaseLine x_coord y_coord n_x n_y];

			% Filter out any data points outside of provided width
			idx=DistFromBaseLine<=(data_width/2) & point_distances>0 & point_distances<max(swdist);
			outData=outData(idx,:);

			% Plot Swath
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

			scatter(outData(:,1),outData(:,2),outData(:,3),'k','filled');

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

			% Transform data
			[n_x,n_y]=ReFrameRot(SW,x_coord,y_coord);

			% Find distances along swath profile of each point
			num_points=numel(x_coord);
			point_distances=zeros(num_points,1);
			DistFromBaseLine=zeros(num_points,1);
			for ii=1:num_points
				xoi=n_x(ii);

				distances=nsw_x-xoi;
				[~,I]=min(abs(distances));
				point_distances(ii)=swdist(I);

				baseLine=nsw_y(I);
				DistFromBaseLine(ii)=abs(baseLine-n_y(ii));
			end
			outData=[point_distances depth magnitude DistFromBaseLine x_coord y_coord n_x n_y];

			% Filter out any data points outside of provided width
			idx=DistFromBaseLine<=(data_width/2) & point_distances>0 & point_distances<max(swdist);
			outData=outData(idx,:);

			% Plot Swath
			f1=figure(1);
			clf 
			set(f1,'Units','inches','Position',[1.0 1.5 16 6],'renderer','painters','PaperSize',[16 6],'PaperOrientation','portrait','PaperPositionMode','auto');

			hold on
			plot(swdist,min_elevs,'-b');
			plot(swdist,max_elevs,'-b');
			plot(swdist,mean_elevs,'-k');

			daspect([vex 1 1])

			yl=ylim;
			for jj=1:numel(bends)
				plot([bends(jj),bends(jj)],yl,'-k');
			end

			scatter(outData(:,1),outData(:,2),outData(:,3),outData(:,4),'filled');

			xlabel('Distance along swath (m)');
			ylabel('Elevation (m)');
			xlim([0 max(swdist)]);
			hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		case 'STREAMobj'
			x_coord=data.x;
			y_coord=data.y;
			z=getnal(data,DEM);

			% Transform data
			[n_x,n_y]=ReFrameRot(SW,x_coord,y_coord);

			% Find distances along swath profile of each point
			num_points=numel(x_coord);
			point_distances=zeros(num_points,1);
			DistFromBaseLine=zeros(num_points,1);
			for ii=1:num_points
				xoi=n_x(ii);

				distances=nsw_x-xoi;
				[~,I]=min(abs(distances));
				point_distances(ii)=swdist(I);

				baseLine=nsw_y(I);
				DistFromBaseLine(ii)=abs(baseLine-n_y(ii));
			end
			outData=[point_distances z DistFromBaseLine x_coord y_coord n_x n_y];

			% Filter out any data points outside of provided width
			idx=DistFromBaseLine<=(data_width/2) & point_distances>0 & point_distances<max(swdist);
			outData=outData(idx,:);

			% Plot Swath
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

			scatter(outData(:,1),outData(:,2),5,'k','filled');

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

			% Transform data
			[n_x,n_y]=ReFrameRot(SW,x_coord,y_coord);

			% Find distances along swath profile of each point
			num_points=numel(x_coord);
			point_distances=zeros(num_points,1);
			DistFromBaseLine=zeros(num_points,1);
			for ii=1:num_points
				xoi=n_x(ii);

				distances=nsw_x-xoi;
				[~,I]=min(abs(distances));
				point_distances(ii)=swdist(I);

				baseLine=nsw_y(I);
				DistFromBaseLine(ii)=abs(baseLine-n_y(ii));
			end
			outData=[point_distances elev ksn DistFromBaseLine x_coord y_coord n_x n_y];

			% Filter out any data points outside of provided width
			idx=DistFromBaseLine<=(data_width/2) & point_distances>0 & point_distances<max(swdist);
			outData=outData(idx,:);

			% Plot Swath
			f1=figure(1);
			clf 
			set(f1,'Units','inches','Position',[1.0 1.5 16 8],'renderer','painters','PaperSize',[16 8],'PaperOrientation','portrait','PaperPositionMode','auto');

			subplot(2,1,1);
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

			subplot(2,1,2);
			hold on
			scatter(outData(:,1),outData(:,3),20,'k','filled');
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
			end

			% Transform data
			[n_x,n_y]=ReFrameRot(SW,x_coord,y_coord);

			% Find distances along swath profile of each point
			num_points=numel(x_coord);
			point_distances=zeros(num_points,1);
			DistFromBaseLine=zeros(num_points,1);
			for ii=1:num_points
				xoi=n_x(ii);

				distances=nsw_x-xoi;
				[~,I]=min(abs(distances));
				point_distances(ii)=swdist(I);

				baseLine=nsw_y(I);
				DistFromBaseLine(ii)=abs(baseLine-n_y(ii));
			end

			if isempty(bs)
				outData=[point_distances z col DistFromBaseLine x_coord y_coord n_x n_y];
			else 
				outData=[point_distances z col scl DistFromBaseLine x_coord y_coord n_x n_y];	
			end			

			% Filter out any data points outside of provided width
			idx=DistFromBaseLine<=(data_width/2) & point_distances>0 & point_distances<max(swdist);
			outData=outData(idx,:);

			% Plot Swath
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

			if isempty(bs)
				scatter(outData(:,1),outData(:,2),20,outData(:,3),'filled');
			else
				scatter(outData(:,1),outData(:,2),outData(:,4),outData(:,3),'filled');
			end

			c1=colorbar;
			ylabel(c1,bv);

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

			% Transform data
			[n_x,n_y]=ReFrameRot(SW,x_coord,y_coord);

			% Find distances along swath profile of each point
			num_points=numel(x_coord);
			point_distances=zeros(num_points,1);
			DistFromBaseLine=zeros(num_points,1);
			for ii=1:num_points
				xoi=n_x(ii);

				distances=nsw_x-xoi;
				[~,I]=min(abs(distances));
				point_distances(ii)=swdist(I);

				baseLine=nsw_y(I);
				DistFromBaseLine(ii)=abs(baseLine-n_y(ii));
			end
			outData=[point_distances ksn DistFromBaseLine x_coord y_coord n_x n_y];

			% Filter out any data points outside of provided width
			idx=DistFromBaseLine<=(data_width/2) & point_distances>0 & point_distances<max(swdist);
			outData=outData(idx,:);

			% Plot Swath
			f1=figure(1);
			clf 
			set(f1,'Units','inches','Position',[1.0 1.5 16 8],'renderer','painters','PaperSize',[16 8],'PaperOrientation','portrait','PaperPositionMode','auto');

			subplot(2,1,1);
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

			subplot(2,1,2);
			hold on
			scatter(outData(:,1),outData(:,2),20,'k','filled');
			xlabel('Distance along swath (m)');
			ylabel('KSN');
			xlim([0 max(swdist)]);
			hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%			
		case 'ksn_profiler'
			x_coord=data(:,1);
			y_coord=data(:,2);
			ksn=data(:,4);

			% Transform data
			[n_x,n_y]=ReFrameRot(SW,x_coord,y_coord);

			% Find distances along swath profile of each point
			num_points=numel(x_coord);
			point_distances=zeros(num_points,1);
			DistFromBaseLine=zeros(num_points,1);
			for ii=1:num_points
				xoi=n_x(ii);

				distances=nsw_x-xoi;
				[~,I]=min(abs(distances));
				point_distances(ii)=swdist(I);

				baseLine=nsw_y(I);
				DistFromBaseLine(ii)=abs(baseLine-n_y(ii));
			end
			outData=[point_distances ksn DistFromBaseLine x_coord y_coord n_x n_y];

			% Filter out any data points outside of provided width
			idx=DistFromBaseLine<=(data_width/2) & point_distances>0 & point_distances<max(swdist);
			outData=outData(idx,:);

			% Plot Swath
			f1=figure(1);
			clf 
			set(f1,'Units','inches','Position',[1.0 1.5 16 8],'renderer','painters','PaperSize',[16 8],'PaperOrientation','portrait','PaperPositionMode','auto');

			subplot(2,1,1);
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

			subplot(2,1,2);
			hold on
			scatter(outData(:,1),outData(:,2),20,'k','filled');
			xlabel('Distance along swath (m)');
			ylabel('KSN');
			xlim([0 max(swdist)]);
			hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%				
		case 'points3_err'
			x_coord=data(:,1);
			y_coord=data(:,2);
			z=data(:,3);
			err=data(:,4);

			% Transform data
			[n_x,n_y]=ReFrameRot(SW,x_coord,y_coord);

			% Find distances along swath profile of each point
			num_points=numel(x_coord);
			point_distances=zeros(num_points,1);
			DistFromBaseLine=zeros(num_points,1);
			for ii=1:num_points
				xoi=n_x(ii);

				distances=nsw_x-xoi;
				[~,I]=min(abs(distances));
				point_distances(ii)=swdist(I);

				baseLine=nsw_y(I);
				DistFromBaseLine(ii)=abs(baseLine-n_y(ii));
			end
			outData=[point_distances z err DistFromBaseLine x_coord y_coord n_x n_y];

			% Filter out any data points outside of provided width
			idx=DistFromBaseLine<=(data_width/2) & point_distances>0 & point_distances<max(swdist);
			outData=outData(idx,:);

			% Plot Swath
			f1=figure(1);
			clf 
			set(f1,'Units','inches','Position',[1.0 5.5 16 4],'renderer','painters','PaperSize',[16 4],'PaperOrientation','portrait','PaperPositionMode','auto');

			hold on
			plot(swdist,min_elevs,'-b');
			plot(swdist,max_elevs,'-b');
			plot(swdist,mean_elevs,'-k');

			yl=ylim;
			for jj=1:numel(bends)
				plot([bends(jj),bends(jj)],yl,'-k');
			end

			xlabel('Distance along swath (m)');
			ylabel('Elevation (m)');
			xlim([0 max(swdist)]);
			hold off

			% Plot Swath
			f2=figure(2);
			clf 
			set(f2,'Units','inches','Position',[1.0 1.5 16 4],'renderer','painters','PaperSize',[16 4],'PaperOrientation','portrait','PaperPositionMode','auto');
			hold on
			errbar(outData(:,1),outData(:,2),outData(:,3),'-r')
			scatter(outData(:,1),outData(:,2),20,'k','filled');

			xlabel('Distance along swath (m)');
			ylabel('Elevation (m)');
			xlim([0 max(swdist)]);
			hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	end

% Function End
end









