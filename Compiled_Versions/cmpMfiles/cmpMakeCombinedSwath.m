function cmpMakeCombinedSwath(wdir,MatFile,points,width,data_type,data,data_width,varargin)
	% Description:
	% 	Function to plot various additional data onto a swath profile.
	%
	% Required Inputs:
	%	wdir - full path of working directory
	%	MatFile - Full path of matfile output from either 'cmpMakeStreams' or the name of a single basin mat file from 'cmpProcessRiverBasins'
	% 	points - name of text file containing n x 2 matrix of x,y points for swath, minimum are two points (start and end points).
	%		First row contains starting point and proceeds down rows, additional points besides a start and end are
	%		treated as bends in the swath. Coordinates for points must be in the same coordinate system as DEM and must
	%		lie within the DEM (cannot be coordinates on the very edge of the DEM).
	% 	width - width of swath in map units
	% 	data_type - the type of additional data you are providing to plot along with the swath, supported inputs are:
	%		'points3' - generic point dataset, expects a n x 3 matrix with values of x, y, and z stored in a text file
	%		'points4' - generic point dataset, expects a n x 4 matrix with values of x, y, z, and extra value stored in a
	%					text file. Dots will be colored by this extra value
	%		'points5' - generic point dataset, expects a n x 5 matrix with values of x, y, z, and two extra values stored in
	%					a text file. Dots will colored by the first extra value (column 4) and scaled by the second extra value 
	%					(column 5).
	%		'eqs' - earthquakes, expects a n x 4 matrix with x, y, depth, and magnitude stored in a text file. Points will be 
	%					scaled by magnitude and colored by distance from swath line. Expects depth to be positive.
	%		'STREAMobj' - will project portions of selected stream profiles (as points) onto a swath. Expects a matfile containing 
	%					a STREAMobj that was generated from the provided DEM (can be the same input as MatFile, but you must provide
	%					it again)
	%		'ksn_chandata' - will plot swath through ksn values, expects full path of a *chandata.mat file as output from old Profiler51 code 
	%					(just in case you have some sitting around)
	%		'ksn_shape' - will plot swath through ksn values, expects the shapefile output from 'cmpKsnChiBatch' or 'cmpKsnProfiler' function 
	%		'basin_stats' - will plot swath through selected mean basin values as calculated from 'cmpProcessRiverBasins', 
	%					expects matfile output from 'cmpCompileBasinStats' and requires an entry to optional input 'basin_value' and  
	%					accepts optional input to 'basin_scale'. Will place point for basin at mean elevation and projected  
	%					location of the basin centroid, will color by value provided to 'basin_value' and will optionall scale  
	%					the point by the value provided to 'basin_scale'
	%		'basin_knicks' - will plot swath through knickpoints as chosen by 'cmpFindBasinKnicks'. For 'data' provide name of folder within working directory
	%					to find knickpoint files saved as a result of running 'FindBasinKnicks' on a series of basins selected from 'ProcessRiverBasins'
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
	%	cmap ['parula'] - valid name of colormap (e.g. 'jet') 
	%	save_figure [false] - logical flag to save the swath figure as a pdf
	%
	% Outputs:
	% 	SwathArray.txt - n x 6 array containing x coordinate, y coordinates, distance along the swath, min elevation, mean elevation, max elevation
	% 	SwathBends.txt - distances along swath of any bends, 0 if no bends
	%	SwathBounds.shp - polyline shapefile showing outline of swath for both topo and data and center line of swath
	% 	SwathProjectedData.txt - data for plotting the swath through the provided data, distances that area 'NaN' indicate those data do not
	%			fall on the swath line provided. Form of output depends on data_type:
	%		'points3' - distances, elevation, distance from base line, x coordinate, y coordinate
	%		'points4' - distances, elevation, value, distance from base line, x coordinate, y coordinate
	%		'eqs' - distances, depth, magnitude, distance from base line, x coordinate, y coordinate
	%		'STREAMobj' - distances, elevation, distance from base line, x coordinate, y coordinate
	%		'ksn_chandata' - distances, elevation, ksn, distance from base line, x coordinate, y coordinate
	%		'ksn_shape' - distances, ksn, distance from base line, x coordinate, y coordinate
	%		'basin_stats' - distances, mean basin elevation, 'basin_value', 'basin_scale' (if provided), distance from base line, 
	%						x coordinate, y coordinate
	%
    % Examples if running for the command line, minus OS specific way of calling main TAK function:
    %   MakeCombinedSwath /path/to/wdir Topo.mat points.txt 10000 points3 data_points.txt 20000
    %   MakeCombinedSwath /path/to/wdir Topo.mat points.txt 10000 basin_stats BasinTable.mat 20000 basin_value mean_ksn
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 06/18/18 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if isdeployed
		if ~isempty(varargin)
			varargin=varargin{1};
		end
	end

	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'cmpMakeCombinedSwath';
	addRequired(p,'wdir',@(x) ischar(x));
	addRequired(p,'MatFile',@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.mat'))));
	addRequired(p,'points',@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.txt'))));
	addRequired(p,'width',@(x) isscalar(x) && isnumeric(x));
	addRequired(p,'data_type',@(x) ischar(validatestring(x,{'points3','points4','points5','eqs','STREAMobj','ksn_chandata','ksn_shape','basin_stats','basin_knicks'})));
	addRequired(p,'data');
	addRequired(p,'data_width',@(x) isnumeric(x) && isscalar(x));

	addParameter(p,'small_circ_center',[],@(x) isnumeric(x) && numel(x)==2);
	addParameter(p,'dist_type','mapdist',@(x) ischar(validatestring(x,{'mapdist','angle'})));	
	addParameter(p,'sample',[],@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'smooth',0,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'vex',10,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'basin_value',[],@(x) ischar(x));
	addParameter(p,'basin_scale',[],@(x) ischar(x));
	addParameter(p,'plot_map',true,@(x) isscalar(x) && islogical(x));
	addParameter(p,'cmap','parula',@(x) ischar(x) || isnumeric(x) & size(x,2)==3);
	addParameter(p,'save_figure',false,@(x) isscalar(x) && islogical(x));


	parse(p,wdir,MatFile,points,width,data_type,data,data_width,varargin{:});
	wdir=p.Results.wdir;
	MatFile=p.Results.MatFile;
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
	plot_map=p.Results.plot_map;
	cmap=p.Results.cmap;
	save_figure=p.Results.save_figure;

	% Determine the type of input
	MatFile=fullfile(wdir,MatFile);
	D=whos('-file',MatFile);
	VL=cell(numel(D),1);
	for ii=1:numel(D);
		VL{ii}=D(ii,1).name;
	end

	if any(strcmp(VL,'DEM')) & any(strcmp(VL,'FD')) & any(strcmp(VL,'A')) & any(strcmp(VL,'S'))
		load(MatFile,'DEM');
	elseif any(strcmp(VL,'DEMcc')) & any(strcmp(VL,'FDc')) & any(strcmp(VL,'Ac')) & any(strcmp(VL,'Sc'))
		load(MatFile,'DEMoc','Sc');
		DEM=DEMoc;
	end

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

	if strcmp(data_type,'basin_stats') & isempty(bv)
		error('You must provide an argument to "basin_value" when plotting "basin_stats"')
	end

	if ~isnumeric(points)
		pts=readtable(fullfile(wdir,points));
		clear points;
		points(:,1)=pts.(1);
		points(:,2)=pts.(2);
	end	

	% Produce topo swath and associated datasets
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
			data=readtable(fullfile(wdir,data));
			x_coord=data.(1);
			y_coord=data.(2);
			z=data.(3);

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

			OT=array2table(outData,'VariableNames',{'distance_along_swath','elevation','distance_from_center','x_coord','y_coord'});
			writetable(OT,fullfile(wdir,'SwathProjectedData.txt'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
		case 'points4'
			data=readtable(fullfile(wdir,data));
			x_coord=data.(1);
			y_coord=data.(2);
			z=data.(3);
			col=data.(4);

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

			OT=array2table(outData,'VariableNames',{'distance_along_swath','elevation','data1','distance_from_center','x_coord','y_coord'});
			writetable(OT,fullfile(wdir,'SwathProjectedData.txt'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
		case 'points5'
			data=readtable(fullfile(wdir,data));
			x_coord=data.(1);
			y_coord=data.(2);
			z=data.(3);
			col=data.(4)
			scle=data.(5);

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

			OT=array2table(outData,'VariableNames',{'distance_along_swath','elevation','data1','data2','distance_from_center','x_coord','y_coord'});
			writetable(OT,fullfile(wdir,'SwathProjectedData.txt'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
		case 'eqs'
			data=readtable(fullfile(wdir,data));
			x_coord=data.(1);
			y_coord=data.(2);
			depth=data.(3);
			magnitude=data.(4);

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


			OT=array2table(outData,'VariableNames',{'distance_along_swath','depth','magntiude','distance_from_center','x_coord','y_coord'});
			writetable(OT,fullfile(wdir,'SwathProjectedData.txt'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		case 'STREAMobj'
			data=fullfile(wdir,data);
			SD=whos('-file',data);
			VS=cell(numel(SD),1);
			for ii=1:numel(SD)
				VS{ii}=SD(ii,1).name;
			end

			if any(strcmp(VS,'Sn'))
				load(data,'Sn');
				S=Sn;
			elseif any(strcmp(VS,'Sc'))
				load(data,'Sc');
				S=Sc;
			elseif any(strcmp(VS,'S'));
				load(data,'S');
			end

			x_coord=S.x;
			y_coord=S.y;
			z=getnal(S,DEM);

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

			OT=array2table(outData,'VariableNames',{'distance_along_swath','elevation','distance_from_center','x_coord','y_coord'});
			writetable(OT,fullfile(wdir,'SwathProjectedData.txt'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		case 'ksn_chandata'
			load(fullfile(wdir,data),'chandata');
			x_coord=chandata(:,9);
			y_coord=chandata(:,10);
			ksn=chandata(:,8);
			elev=chandata(:,4);

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

			OT=array2table(outData,'VariableNames',{'distance_along_swath','elevation','ksn','distance_from_center','x_coord','y_coord'});
			writetable(OT,fullfile(wdir,'SwathProjectedData.txt'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		case 'basin_stats'
			load(fullfile(wdir,data),'T');
			x_coord=T.center_x;
			y_coord=T.center_y;
			z=T.mean_el;
			col=T.(bv);
			if ~isempty(bs)
				scl=T.(bs);
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

			if isempty(bs)
				OT=array2table(outData,'VariableNames',{'distance_along_swath','elevation',bv,'distance_from_center','x_coord','y_coord'});
				writetable(OT,fullfile(wdir,'SwathProjectedData.txt'));
			else
				OT=array2table(outData,'VariableNames',{'distance_along_swath','elevation',bv,bs,'distance_from_center','x_coord','y_coord'});
				writetable(OT,fullfile(wdir,'SwathProjectedData.txt'));
			end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		case 'ksn_shape'
			data=shaperead(fullfile(wdir,data));
			% Find number of segments in provided KSN map structure	
			numSegs=numel(data);
			% Loop through stream segments to extract x,y,ksn
			streamData=zeros(numSegs,3);
			for kk=1:numSegs
				xx=nanmean(data(kk,1).X);
				yy=nanmean(data(kk,1).Y);
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

			OT=array2table(outData,'VariableNames',{'distance_along_swath','ksn','distance_from_center','x_coord','y_coord'});
			writetable(OT,fullfile(wdir,'SwathProjectedData.txt'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%			
		case 'basin_knicks'

			fileList=dir(fullfile(wdir,data,'Knicks_*.mat'));
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

			OT=array2table(outData,'VariableNames',{'distance_along_swath','elevation','distance_from_center','x_coord','y_coord'});
			writetable(OT,fullfile(wdir,'SwathProjectedData.txt'));
	end

	% Make output shape
	ms=struct;
	ms(1,1).Geometry='Line';
	ms(1,1).X=SW.xy0(:,1);
	ms(1,1).Y=SW.xy0(:,2);
	ms(1,1).Type='Center';

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

	shapewrite(ms,fullfile(wdir,'SwathBounds.shp'));

	% Write out other info as table
	ST=table;
	ST.x_coord=xypoints(:,1);
	ST.y_coord=xypoints(:,2);
	ST.distance=SwathMat(:,1);
	ST.min_z=SwathMat(:,2);
	ST.mean_z=SwathMat(:,3);
	ST.max_z=SwathMat(:,4);
	writetable(ST,fullfile(wdir,'SwathArray.txt'));

	BT=table;
	BT.bends_distance=bends;
	writetable(BT,fullfile(wdir,'SwathBends.txt'));

	if plot_map
		f2=figure(2);
		set(f2,'Units','normalized','Position',[0.05 0.1 0.6 0.6]);
		hold on
		imageschs(DEM,DEM,'colormap','gray');
		plot(SW.xy0(:,1),SW.xy0(:,2),'-g','LineWidth',0.5);
		plot(Tverts(:,1),Tverts(:,2),'-g','LineWidth',0.5);
		plot(Dverts(:,1),Dverts(:,2),'-r','LineWidth',0.5);		
		scatter(x_coord(idx),y_coord(idx),20,'r','filled');
		scatter(x_coord(~idx),y_coord(~idx),20,'w','filled');
		hold off
	end

	if save_figure
		orient(f1,'Landscape')
		print(f1,'-dpdf','-bestfit','Swath.pdf');
		close(f1);
	end

	msgbox('Close all figures to complete code execution');

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









