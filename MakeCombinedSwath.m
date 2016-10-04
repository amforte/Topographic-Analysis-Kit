function [SW,SwathMat,xypoints,outData]=MakeCombinedSwath(DEM,points,width,sample,smooth,data_type,data,data_width)
	% Function requires TopoToolbox to run
	%
	% Required Inputs:
	% DEM - DEM Grid Object with which to make topo swath
	% points - n x 2 matrix containing x,y points for swath, minimum are two points (start and end points).
	%	First row contains starting point and proceeds down rows, additional points besides a start and end are
	%	treated as bends in the swath. Coordinates for points must be in the same coordinate system as DEM and must
	%	lie within the DEM (cannot be coordinates on the very edge of the DEM)
	% width - width of swath in map units
	% sample - resampling distance along swath in map units, provide DEM.cellsize for no resampling, where DEM is the name of the 
	%	DEM GRIDobj.
	% smooth - smoothing distance, width of filter in map units over which to smooth values, provie 0 if no smoothing is desired
	% data_type - the type of additional data you are providing to plot along with the swath, supported inputs are:
	%		'points3' - generic point dataset, expects a n x 3 matrix with values of x, y, and z
	%		'points4' - generic point dataset, expects a n x 4 matrix with values of x, y, z, and extra value. Dots will be sclaed by this extra value
	%		'eqs' - earthquakes, expects a n x 4 matrix with x, y, depth, and magnitude. Points will be scaled by magnitude and colored by distance from swath line
	%		'ksn_chandata' - will plot swath through ksn values, expects a chandata file from Profiler
	%		'ksn_TT' - will plot swath through ksn values, expects output from KSN_TT function
	%		'chi_TT' - will plot swath through chi values, expects output of TopoToolbox chiplot function
	% data_width - width in map units of swath through provided data. Values greater than data_width/2 from the center line of the toposwath will be removed
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

	% Produce topo swath and associated datasets
	[SW,SwathMat,xypoints,bends]=MakeTopoSwath(DEM,points,width,sample,smooth,'false');

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
			idx=DistFromBaseLine<=(data_width/2);
			outData=outData(idx,:);

			% Plot Swath
			f1=figure(1);
			clf 
			set(f1,'Units','inches','Position',[1.0 1.5 16 4],'renderer','painters','PaperSize',[16 4],'PaperOrientation','portrait','PaperPositionMode','auto');

			hold on
			plot(swdist,min_elevs,'-b');
			plot(swdist,max_elevs,'-b');
			plot(swdist,mean_elevs,'-k');

			for jj=1:numel(bends)
				vline(bends(jj),'k','bend');
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
			scale=data(:,4);

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
			outData=[point_distances z scale DistFromBaseLine x_coord y_coord n_x n_y];

			% Filter out any data points outside of provided width
			idx=DistFromBaseLine<=(data_width/2);
			outData=outData(idx,:);

			% Plot Swath
			f1=figure(1);
			clf 
			set(f1,'Units','inches','Position',[1.0 1.5 16 4],'renderer','painters','PaperSize',[16 4],'PaperOrientation','portrait','PaperPositionMode','auto');

			hold on
			plot(swdist,min_elevs,'-b');
			plot(swdist,max_elevs,'-b');
			plot(swdist,mean_elevs,'-k');

			for jj=1:numel(bends)
				vline(bends(jj),'k','bend');
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
			idx=DistFromBaseLine<=(data_width/2);
			outData=outData(idx,:);

			% Plot Swath
			f1=figure(1);
			clf 
			set(f1,'Units','inches','Position',[1.0 1.5 16 6],'renderer','painters','PaperSize',[16 6],'PaperOrientation','portrait','PaperPositionMode','auto');

			hold on
			plot(swdist,min_elevs,'-b');
			plot(swdist,max_elevs,'-b');
			plot(swdist,mean_elevs,'-k');

			for jj=1:numel(bends)
				vline(bends(jj),'k','bend');
			end

			scatter(outData(:,1),outData(:,2),outData(:,3),outData(:,4),'filled');

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
			idx=DistFromBaseLine<=(data_width/2);
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

			for jj=1:numel(bends)
				vline(bends(jj),'k','bend');
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
		case 'ksn_TT'
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
			idx=DistFromBaseLine<=(data_width/2);
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

			for jj=1:numel(bends)
				vline(bends(jj),'k','bend');
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
		case 'chi_TT'
			% Extract values of interest
			x_coord=data.x;
			y_coord=data.y;
			chi=data.chi;

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
			outData=[point_distances chi DistFromBaseLine x_coord y_coord n_x n_y];

			% Filter out any data points outside of provided width
			idx=DistFromBaseLine<=(data_width/2);
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

			for jj=1:numel(bends)
				vline(bends(jj),'k','bend');
			end

			xlabel('Distance along swath (m)');
			ylabel('Elevation (m)');
			xlim([0 max(swdist)]);
			hold off

			subplot(2,1,2);
			hold on
			scatter(outData(:,1),outData(:,2),20,'k','filled');
			xlabel('Distance along swath (m)');
			ylabel('Chi');
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
			idx=DistFromBaseLine<=(data_width/2);
			outData=outData(idx,:);

			% Plot Swath
			f1=figure(1);
			clf 
			set(f1,'Units','inches','Position',[1.0 5.5 16 4],'renderer','painters','PaperSize',[16 4],'PaperOrientation','portrait','PaperPositionMode','auto');

			hold on
			plot(swdist,min_elevs,'-b');
			plot(swdist,max_elevs,'-b');
			plot(swdist,mean_elevs,'-k');

			for jj=1:numel(bends)
				vline(bends(jj),'k','bend');
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









