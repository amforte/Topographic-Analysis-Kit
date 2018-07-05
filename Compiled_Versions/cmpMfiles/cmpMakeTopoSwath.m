function cmpMakeTopoSwath(wdir,MatFile,points,width,varargin)
	% Description:
	% 	Wrapper around TopoToolbox SWATHobj functionality
	%
	% Required Inputs:
	%	wdir - full path of working directory
	%	MatFile - Full path of matfile output from either 'cmpMakeStreams' or the name of a single basin mat file from 'cmpProcessRiverBasins'
	% 	points - name of text file containing n x 2 matrix of x,y points for swath, minimum are two points (start and end points).
	%		First row contains starting point and proceeds down rows, additional points besides a start and end are
	%		treated as bends in the swath. Coordinates for points must be in the same coordinate system as DEM and must
	%		lie within the DEM (cannot be coordinates on the very edge of the DEM).
	% 	width - width of swath in map units
	%
	% Optional Inputs:
	% 	sample [] - resampling distance along swath in map units, if no input is provided, code will use the cellsize of the DEM 
	%		which results in no resampling.
	% 	smooth [0] - smoothing distance, width of filter in map units over which to smooth values, default (0) results in no smoothing
	%	vex [10] - vertical exaggeration for displaying plot.
	%	plot_as_points [false] - logical flag to switch plot type to distributions of points
	%	plot_as_heatmap [false] - logical flag to switch plot type to a heat map
	%	save_figure [false] - logical flag to save the swath figure as a pdf (this will also set 'plot_figure' to true)
	%
	% Outputs:
	% 	SwathArray.txt - n x 6 array containing x coordinate, y coordinates, distance along the swath, min elevation, mean elevation, max elevation
	% 	SwathBends.txt - distances along swath of any bends, 0 if no bends
	%	SwathBounds.shp - polyline shapefile showing outline of swath and center line of swath
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 07/02/18 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if isdeployed
		if ~isempty(varargin)
			varargin=varargin{1};
		end
	end

	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'cmpMakeTopoSwath';
	addRequired(p,'wdir',@(x) ischar(x));
	addRequired(p,'MatFile',@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.mat'))));
	addRequired(p,'points',@(x)  ~isempty(regexp(x,regexptranslate('wildcard','*.txt'))));
	addRequired(p,'width',@(x) isscalar(x) && isnumeric(x));

	addParameter(p,'sample',[],@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'smooth',0,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'vex',10,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'plot_as_points',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'plot_as_heatmap',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'save_figure',false,@(x) isscalar(x) && islogical(x));

	parse(p,wdir,MatFile,points,width,varargin{:});
	wdir=p.Results.wdir;
	MatFile=p.Results.MatFile;
	points=p.Results.points;
	wdth=p.Results.width;

	sample=p.Results.sample;
	smth=p.Results.smooth;
	vex=p.Results.vex;
	plot_figure=p.Results.plot_figure;
	plot_as_points=p.Results.plot_as_points;
	plot_as_heatmap=p.Results.plot_as_heatmap;
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
		load(MatFile,'DEMoc','FDc','Ac','Sc','DEMcc');
		DEM=DEMoc;
	end

	if isempty(sample)
		sample=DEM.cellsize;
	end

	if plot_as_points & plot_as_heatmap
		error('Please only set one of "plot_as_points" and "plot_as_heatmap" to true');
	end

	if save_figure
		plot_figure=true;
	end

	if ~isnumeric(points)
		pts=readtable(fullfile(wdir,points));
		clear points;
		points(:,1)=pts.(1);
		points(:,2)=pts.(2);
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

	f1=figure(1);
	clf 
	set(f1,'Units','normalized','Position',[0.05 0.1 0.8 0.4],'renderer','painters');

	hold on

	if plot_as_points
		for ii=1:size(elevs,1)
			scatter(swdist,elevs(ii,:),1,'k','.');
		end

	elseif plot_as_heatmap
		el_range=linspace(min(min_elevs)-1,max(max_elevs)+1,101);
		el_range_p=linspace(min(min_elevs)-1,max(max_elevs)+1,100);
		C=zeros(100,numel(swdist));

		cmap = jet(256);
		cmap(1,:) = 1;
		colormap(cmap);

		for ii=1:numel(swdist)
			[N,~]=histcounts(SW.Z(:,ii),el_range);
			N=N';
			mi=min_elevs(ii);
			ma=max_elevs(ii);
			idx=el_range_p>ma | el_range_p<mi;
			N(idx)=-1;
			C(:,ii)=N;
		end

		imagesc(swdist,el_range_p,C);
		plot(swdist,min_elevs,'-k');
		plot(swdist,max_elevs,'-k');		
	else
		xx=vertcat(swdist,flipud(swdist));
		yy=horzcat(min_elevs,fliplr(max_elevs));
		patch(xx,yy,[0.8 0.8 0.8]);

		plot(swdist,min_elevs,'-k');
		plot(swdist,max_elevs,'-k');
		plot(swdist,mean_elevs,'-k','LineWidth',2);
	end

	daspect([vex 1 1])

	yl=ylim;
	for jj=1:numel(bends)
		plot([bends(jj),bends(jj)],yl,'-k');
	end

	xlabel(['Distance along swath (m) : VEX = ' num2str(vex)]);
	ylabel('Elevation (m)');
	xlim([0 max(swdist)]);
	hold off


	if save_figure
		orient(f1,'Landscape')
		if plot_as_points
			print(f1,'-dtiff',fullfile(wdir,'Swath.tif'));
		else
			print(f1,'-dpdf','-bestfit',fullfile(wdir,'Swath.pdf'));
		end
		close(f1);
	end

	% Make output shape
	ms=struct;
	ms(1,1).Geometry='Line';
	ms(1,1).X=SW.xy0(:,1);
	ms(1,1).Y=SW.xy0(:,2);
	ms(1,1).Type='Center';

	verts=SwathPolygon(SW,wdth);
	ms(2,1).Geometry='Line';
	ms(2,1).X=verts(:,1);
	ms(2,1).Y=verts(:,2);
	ms(2,1).Type='TopoWdth';

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


