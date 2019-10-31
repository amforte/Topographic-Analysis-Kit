function [SWcell,points]=MakeSerialSwath(DEM,points,divisions,sw_length,varargin)
	% Usage:
	% [SWcell,points]=MakeSerialSwath(DEM,points,divisions,sw_length);
	% [SWcell,points]=MakeSerialSwath(DEM,points,divisions,sw_length,'name',value);
	%
	% Description:
	%	Function to create a series of swath profiles perpendicular to a provided line
	%
	% Required Inputs:
	% 	DEM - DEM Grid Object with which to make topo swath
	% 	points - n x 2 matrix containing x,y points for the line along which series of perpendicular swaths will be generated,
	%		minimum are two points (start and end points). First row contains starting point and proceeds down rows, additional 
	%		points besides a start and end are treated as bends in the line. Coordinates for points must be in the same coordinate 
	%		system as DEM and must lie within the DEM (cannot be coordinates on the very edge of the DEM). If you provide an empty 
	%		array this will invoke a display of the DEM that you can use to manually define the line.
	%	divisions - scalar that controls the number of swaths that will be generated. How this parameter is interpreted depends on
	%		value of the optional 'div_type' parameter. If 'div_type' is 'number' then the value provided to 'divisions' will be the 
	%		total number of swaths created and the width of these swaths will equal the length of the line (defined by 'points')
	%		divided by the number of swaths. If 'div_type' is 'width' then the value provided to 'divisions' will be the width of each 
	%		swath in map units. The code will start producing swaths from the beginning of the line defined by points and will proceed 
	%		until it can not produce a swath that is 'width' wide (e.g. if the total length of the line defined by points is 11000 meters
	%	 	and 2000 is provided to 'divisions' with 'div_type' set to 'width' then the result will be 5, 2000 meter wide swaths and the 
	%		last 1000 meter length of the line defined by 'points' will be ignored).
	%	sw_length - length of individual swaths (perpendicular to the line defined by 'points')
	%	
	% Optional Inputs:
	%	div_type ['number'] - controls how the 'divisions' parameter is interpreted. Viable entries are 'number' or 'width', see
	%		above the 'division' parameter for more info.
	%	alignment ['center'] - controls how the serial swaths are drawn relative to the line defined by 'points'. Vialble entries are
	%		'center', 'left', 'right', and 'between'. Behavior is as follows:
	%			'center' - Swaths will be drawn so that their centers intersect the line defined by 'points'. 
	%			'right' - Swaths will be drawn to the righthand side of the line. This is based on the order in which points are defined, 
	%				e.g. if the xy coordinates provided to 'points' define an east-west oriented line and the order of the points go from west 
	%				to east then a 'right' alignment will draw swaths south of the line.
	%			'left' - Swaths will be drawn to the lefthand side of the line. This is based on the order in which points are defined, 
	%				e.g. if the xy coordinates provided to 'points' define an east-west oriented line and the order of the points go from west 
	%				to east then a 'left' alignment will draw swaths north of the line.
	%			'between' - Swaths will be drawn between two lines, one provided to the required 'points' input and one provided to the optional
	%				'points2' input. If 'alignment' is set to 'between' and you provide an input to 'points' but none to 'points2', this will
	%				generate an error. If 'alignment' is set to 'between' and you leave 'points' empty, then you can graphically select both 
	%				lines between which the serial swaths will be drawn. The code uses the shorter of the two provided lines to determine the
	%				number of swaths (and the widths of those swaths if 'div_type' is set to 'width'). Note that setting 'alignment' to 'between'
	%				means the input to the provided 'sw_length' parameter will be ignored.
	%	'points2' [] - option to provide a n x matrix containing x,y points for defining the other half of the area within which to draw swaths
	%		if 'alignment' is set to 'between'.
	%	add_grids [] - option to provide a cell array of additional grids to be sampled in the swath. The expected input is a nx2 cell array,
	%		where the first column is a GRIDobj and the second column is a string identifying what this grid is (will be used as the label on the
	%		resultant swaths).
	% 	sample [cellsize of DEM] - resampling distance along swath in map units, if no input is provided, code will use the cellsize of the DEM 
	%		which results in no resampling. 
	% 	smooth [0] - smoothing distance, width of filter in map units over which to smooth values, default (0) results in no smoothing
	%	plot_map [true] - logical to turn on plotting of map showing location of all swaths
	%	plot_individual [false] - logical to turn on plotting of individual maps and swaths (i.e. there will be one figure per swath)
	%
	% Outputs:
	%	SWcell - Cell array containing TopoToolbox SWATHobj's resultant from running the function. If only a DEM is provided, then the
	%		cell array will have one column and rows correspond to individual swaths. If additional grids are provided, then additional
	%		columns correspond to these additional grids.
	%	points - n x 2 array defining the points for the reference line used to draw the swaths. Will be identical to the 'points' input
	%		if provided, this is provided as an output in case the user doesn't intially provide an input to 'points' and uses the graphical
	%		selection option. Note that if 'alignment' is set to 'between' then the output of 'points' will be a cell array with two columns
	%		containing the n x 2 arrays that define the two lines between which the serial swaths were generated.
	%
	% Examples:
	%	[SWcell,points]=MakeSerialSwath(DEM,points,100,1000);
	%	[SWcell,points]=MakeSerialSwath(DEM,points,1000,1000,'div_type','width');
	%	[SWcell,points]=MakeSerialSwath(DEM,points,100,1000,'alignment','center');	
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 04/02/19 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'MakeSerialSwath';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'points',@(x) isempty(x) || isnumeric(x) & size(x,1)>=2 && size(x,2)==2);
	addRequired(p,'divisions',@(x) isscalar(x) && isnumeric(x));
	addRequired(p,'sw_length',@(x) isscalar(x) && isnumeric(x));

	addParameter(p,'points2',[],@(x) isnumeric(x) & size(x,1)>=2 && size(x,2)==2 || isempty(x));
	addParameter(p,'div_type','number',@(x) ischar(validatestring(x,{'number','width'})));
	addParameter(p,'alignment','center',@(x) ischar(validatestring(x,{'center','right','left','between'})));
	addParameter(p,'add_grids',[],@(x) isa(x,'cell') && size(x,2)==2 || isempty(x));
	addParameter(p,'sample',[],@(x) isscalar(x) && isnumeric(x) || isempty(x));
	addParameter(p,'smooth',0,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'plot_map',true,@(x) isscalar(x) && islogical(x));
	addParameter(p,'plot_individual',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'make_shape',true,@(x) islogical(x)); % Hidden parameter to fix problems with generating shapefiles here for compiled version
	addParameter(p,'out_dir',[],@(x) isdir(x));

	parse(p,DEM,points,divisions,sw_length,varargin{:});
	DEM=p.Results.DEM;
	points=p.Results.points;
	divisions=p.Results.divisions;
	sw_length=p.Results.sw_length;

	points2=p.Results.points2;
	div_type=p.Results.div_type;
	alignment=p.Results.alignment;
	AG=p.Results.add_grids;
	sample=p.Results.sample;
	smth=p.Results.smooth;
	plot_map=p.Results.plot_map;
	plot_individual=p.Results.plot_individual;
	make_shape=p.Results.make_shape;
	out_dir=p.Results.out_dir;

	if isempty(sample)
		sample=DEM.cellsize;
	end

	if isempty(out_dir)
		out_dir=pwd;
	end	

	if ~isempty(points) & strcmp(alignment,'between') & isempty(points2)
		error('If "alignment" is set to "between" then an entry must be provided for "points2"')
	end

	% Initiate graphics to choose points if 'points' is empty
	if isempty(points) & ~strcmp(alignment,'between')
		f1=figure(1);
		set(f1,'Units','normalized','Position',[0.05 0.1 0.9 0.9]);
		clf
        imagesc(DEM)
    	if ~verLessThan('matlab','9.5')
        disableDefaultInteractivity(gca);
	    end 
        title('Draw line along which to generate serial swaths (double-click to end)')
        [points] = getline;
        close(f1);
    elseif isempty(points) & strcmp(alignment,'between')
		f1=figure(1);
		set(f1,'Units','normalized','Position',[0.05 0.1 0.9 0.9]);
		clf
		hold on
        imagesc(DEM)
		if ~verLessThan('matlab','9.5')
	        disableDefaultInteractivity(gca);
	    end         
        title('Draw first bounding line for serial swaths (double-click to end)')
        [points1] = getline;
        plot(points1(:,1),points1(:,2),'-r');
        title('Draw second bounding line for serial swaths (double-click to end)')
        [points2] = getline;
        hold off
        close(f1);

        points={points1,points2};
    elseif ~isempty(points) & ~isempty(points2) & strcmp(alignment,'between')
    	points1=points;
    	points={points1,points2};
    end

	% Get extents of DEM
	[demx,demy]=getoutline(DEM,true);	

    % Discriminate between 'between' and other options
    if strcmp(alignment,'between');
    	num_bounds=2;

    	% Check if either bound is outside DEM
    	in_bnd1=inpolygon(points1(:,1),points1(:,2),demx,demy);
    	in_bnd2=inpolygon(points2(:,1),points2(:,2),demx,demy);

    	if nnz(in_bnd1)~=numel(in_bnd1)
    		error('Portion of first boundary for serial swaths lies outside the DEM');
    	elseif nnz(in_bnd2)~=numel(in_bnd2)
    		error('Portion of second boundary for serial swaths lies outside the DEM');
    	end

    else
    	num_bounds=1;
    end

    switch num_bounds
    case 1 % If alignment is left, right, or center

		% Find Total Distance Along Swath
		xd=diff(points(:,1)); yd=diff(points(:,2));
		dst=sqrt((xd.^2)+(yd.^2));
		dst=cumsum(dst);
		dst=vertcat(0,dst);
		tot_dst=dst(end);

		% Determine type of divisions
		if strcmp(div_type,'number')
			number=divisions;
			trim_flag=false;
			if tot_dst/number < DEM.cellsize*2
				warning('Width of resultant swaths are near that of a single grid cell, did you input a width for swaths instead of a number of swaths?')
			end
		elseif strcmp(div_type,'width')
			rmndr=rem(tot_dst,divisions);
			number=(tot_dst-rmndr)/divisions;
			if rmndr~=0
				trim_flag=true;
				tot_dst=tot_dst-rmndr;
			else
				trim_flag=false
			end
		end

		% Find swath edges and midpoints in distance space
		edg=linspace(0,tot_dst,number+1);
		wdth=unique(round(diff(edg)));
		mp=edg(2:end)-(wdth/2);

		% Find swath edges and midpoints in cartesian coordinates
		[nx,ny,nd]=interpline(points(:,1),points(:,2),dst,1);

		if trim_flag
			nx=nx(nd<=tot_dst);
			ny=ny(nd<=tot_dst);
			nd=nd(nd<=tot_dst);
		end

		edg_x=zeros(numel(edg),1);
		edg_y=zeros(numel(edg),1);
		for ii=1:numel(edg)
			[~,ix]=min(abs(edg(ii)-nd));
			edg_x(ii)=nx(ix);
			edg_y(ii)=ny(ix);
		end

		mp_x=zeros(numel(mp),1);
		mp_y=zeros(numel(mp),1);
		for ii=1:numel(mp)
			[~,ix]=min(abs(mp(ii)-nd));
			mp_x(ii)=nx(ix);
			mp_y(ii)=ny(ix);
		end

		% Calculate swath start and end points
		dx=diff(edg_x); dy=diff(edg_y);
		[ang_alng,~]=cart2pol(dx,dy);

		switch alignment
		case 'center'
			[px,py]=pol2cart(ang_alng+(pi/2),sw_length/2);
			[mx,my]=pol2cart(ang_alng-(pi/2),sw_length/2);

			start_x=mp_x+px; start_y=mp_y+py;
			stop_x=mp_x+mx; stop_y=mp_y+my;
		case 'left'
			[px,py]=pol2cart(ang_alng+(pi/2),sw_length);

			start_x=mp_x+px; start_y=mp_y+py;
			stop_x=mp_x; stop_y=mp_y;
		case 'right'
			[mx,my]=pol2cart(ang_alng-(pi/2),sw_length);

			start_x=mp_x; start_y=mp_y;
			stop_x=mp_x+mx; stop_y=mp_y+my;
		end


		% Check if any of the calculated swath points lie outside the dem
		start_in=inpolygon(start_x,start_y,demx,demy);
		stop_in=inpolygon(stop_x,stop_y,demx,demy);

		if nnz(start_in)~=numel(start_in)
			error('Some of the serial swath start points lie outside the DEM, adjust the length of the swaths or the alignment');
		elseif nnz(stop_in)~=numel(stop_in)
			error('Some of the serial swath stop points lie outside the DEM, adjust the length of the swaths or the alignment');
		end

	case 2 % If alignment is between

		% Find Total Distance Along Edges of Swath Zone
		xd1=diff(points1(:,1)); yd1=diff(points1(:,2));
		dst1=sqrt((xd1.^2)+(yd1.^2));
		dst1=cumsum(dst1);
		dst1=vertcat(0,dst1);
		tot_dst1=dst1(end);

		xd2=diff(points2(:,1)); yd2=diff(points2(:,2));
		dst2=sqrt((xd2.^2)+(yd2.^2));
		dst2=cumsum(dst2);
		dst2=vertcat(0,dst2);
		tot_dst2=dst2(end);	

		% Determine type of divisions
		[min_dst,min_ix]=min([tot_dst1 tot_dst2]);
		[max_dst]=max([tot_dst1 tot_dst2]);

		if strcmp(div_type,'number')
			number=divisions;
			if min_dst/number < DEM.cellsize*2
				warning('Width of resultant swaths are near that of a single grid cell, did you input a width for swaths instead of a number of swaths?')
			end

			% Find swath edges and midpoints in distance space
			edg_min=linspace(0,min_dst,number+1);
			wdth=unique(round(diff(edg_min)));
			mp_min=edg_min(2:end)-(wdth/2);

			edg_max=linspace(0,max_dst,number+1);
			wdth_max=unique(round(diff(edg_max)));
			mp_max=edg_max(2:end)-(wdth_max/2);	

			% Find swath edges and midpoints in cartesian coordinates
			[nx1,ny1,nd1]=interpline(points1(:,1),points1(:,2),dst1,1);
			[nx2,ny2,nd2]=interpline(points2(:,1),points2(:,2),dst2,1);

			if min_ix==1
				start_x=zeros(numel(mp_min),1);
				start_y=zeros(numel(mp_min),1);
				for ii=1:numel(mp_min)
					[~,ix]=min(abs(mp_min(ii)-nd1));
					start_x(ii)=nx1(ix);
					start_y(ii)=ny1(ix);
				end

				stop_x=zeros(numel(mp_max),1);
				stop_y=zeros(numel(mp_max),1);
				for ii=1:numel(mp_max)
					[~,ix]=min(abs(mp_max(ii)-nd2));
					stop_x(ii)=nx2(ix);
					stop_y(ii)=ny2(ix);
				end
			else
				start_x=zeros(numel(mp_max),1);
				start_y=zeros(numel(mp_max),1);
				for ii=1:numel(mp_max)
					[~,ix]=min(abs(mp_max(ii)-nd1));
					start_x(ii)=nx1(ix);
					start_y(ii)=ny1(ix);
				end

				stop_x=zeros(numel(mp_min),1);
				stop_y=zeros(numel(mp_min),1);
				for ii=1:numel(mp_min)
					[~,ix]=min(abs(mp_min(ii)-nd2));
					stop_x(ii)=nx2(ix);
					stop_y(ii)=ny2(ix);
				end
			end

		elseif strcmp(div_type,'width')
			rmndr=rem(min_dst,divisions);
			number=(min_dst-rmndr)/divisions;
			if rmndr~=0
				trim_flag=true;
				min_dst=min_dst-rmndr;
			else
				trim_flag=false
			end

			% Find swath edges and midpoints in distance space
			edg_min=linspace(0,min_dst,number+1);
			wdth=unique(round(diff(edg_min)));
			mp_min=edg_min(2:end)-(wdth/2);

			edg_max=linspace(0,max_dst,number+1);
			wdth_max=unique(round(diff(edg_max)));
			mp_max=edg_max(2:end)-(wdth_max/2);	

			% Find swath edges and midpoints in cartesian coordinates
			[nx1,ny1,nd1]=interpline(points1(:,1),points1(:,2),dst1,1);
			[nx2,ny2,nd2]=interpline(points2(:,1),points2(:,2),dst2,1);

			if min_ix==1
				if trim_flag
					nx1=nx1(nd1<=min_dst);
					ny1=ny1(nd1<=min_dst);
					nd1=nd1(nd1<=min_dst);
				end

				start_x=zeros(numel(mp_min),1);
				start_y=zeros(numel(mp_min),1);
				for ii=1:numel(mp_min)
					[~,ix]=min(abs(mp_min(ii)-nd1));
					start_x(ii)=nx1(ix);
					start_y(ii)=ny1(ix);
				end

				stop_x=zeros(numel(mp_max),1);
				stop_y=zeros(numel(mp_max),1);
				for ii=1:numel(mp_max)
					[~,ix]=min(abs(mp_max(ii)-nd2));
					stop_x(ii)=nx2(ix);
					stop_y(ii)=ny2(ix);
				end
			else
				if trim_flag
					nx2=nx2(nd2<=min_dst);
					ny2=ny2(nd2<=min_dst);
					nd2=nd2(nd2<=min_dst);
				end

				start_x=zeros(numel(mp_max),1);
				start_y=zeros(numel(mp_max),1);
				for ii=1:numel(mp_max)
					[~,ix]=min(abs(mp_max(ii)-nd1));
					start_x(ii)=nx1(ix);
					start_y(ii)=ny1(ix);
				end

				stop_x=zeros(numel(mp_min),1);
				stop_y=zeros(numel(mp_min),1);
				for ii=1:numel(mp_min)
					[~,ix]=min(abs(mp_min(ii)-nd2));
					stop_x(ii)=nx2(ix);
					stop_y(ii)=ny2(ix);
				end
			end
		end
	end

	% Calculate Swaths
	if ~isempty(AG)
		num_grids=size(AG,1);
	else
		num_grids=0;
	end

	SWcell=cell(number,1+num_grids);
	PLOTcell=cell(number,1+num_grids);
	if ~verLessThan('matlab','9.4');
		VERTcell=cell(number,1);
	end

	w1=waitbar(0,'Generating swaths');
	for ii=1:number
		SWcell{ii,1}=SWATHobj(DEM,vertcat(start_x(ii),stop_x(ii)),vertcat(start_y(ii),stop_y(ii)),'width',wdth,'dx',sample,'smooth',smth);
		if ~isempty(AG)
			for jj=1:num_grids
				AGoi=AG{jj,1};
				SWcell{ii,jj+1}=SWATHobj(AGoi,vertcat(start_x(ii),stop_x(ii)),vertcat(start_y(ii),stop_y(ii)),'width',wdth);
			end
		end

		PLOTcell{ii,1}=[SWcell{ii,1}.distx min(SWcell{ii,1}.Z,[],'omitnan').' mean(SWcell{ii,1}.Z,'omitnan').' max(SWcell{ii,1}.Z,[],'omitnan').'];
		if ~isempty(AG)
			for jj=1:num_grids
				PLOTcell{ii,jj+1}=[SWcell{ii,jj+1}.distx min(SWcell{ii,jj+1}.Z,[],'omitnan').' mean(SWcell{ii,jj+1}.Z,'omitnan').' max(SWcell{ii,jj+1}.Z,[],'omitnan').'];
			end
		end

		if ~verLessThan('matlab','9.4');
			VERTcell{ii}=SwathPolygon(SWcell{ii,1},wdth);
		end
		waitbar(ii/number,w1);
	end
	close(w1);

	if plot_map
		f1=figure(1);
		clf
		set(f1,'Units','normalized','Position',[0.05 0.1 0.6 0.9]);
		cmap=jet(number);


		switch num_bounds
		case 1

			sbplt1=subplot(3,1,1);
			hold on
			imagesc(DEM);
			for ii=1:number
				if ~verLessThan('matlab','9.4');
					plot(VERTcell{ii}(:,1),VERTcell{ii}(:,2),'Color',cmap(ii,:),'LineWidth',2);
				else
					plot(SWcell{ii,1}.xy0(:,1),SWcell{ii,1}.xy0(:,2),'Color',cmap(ii,:),'LineWidth',2);
				end
			end
			plot(points(:,1),points(:,2),'-k','LineWidth',2);
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(sbplt1);
		    end 
			hold off

			sbplt2=subplot(3,1,2);
			hold on 
			for ii=1:number
				plot(PLOTcell{ii,1}(:,1),PLOTcell{ii,1}(:,3),'LineWidth',2,'Color',cmap(ii,:));
				plot(PLOTcell{ii,1}(:,1),PLOTcell{ii,1}(:,2),'LineWidth',0.5,'Color',cmap(ii,:));
				plot(PLOTcell{ii,1}(:,1),PLOTcell{ii,1}(:,4),'LineWidth',0.5,'Color',cmap(ii,:));
			end
			title('All Elevation Swaths')
			xlabel('Distance Along Swaths (m)');
			ylabel('Elevation Along Swaths (m)');
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(sbplt2);
		    end 			
			hold off

			% Control for slight variability in number of points sampled along x direction
			for ii=1:number
				num_dist(ii,1)=numel(PLOTcell{ii,1}(:,3));
			end
			mnd=min(num_dist);

			% Assemble arrays
			for ii=1:number
				mean_comp(:,ii)=PLOTcell{ii,1}(1:mnd,3);
				max_comp(:,ii)=PLOTcell{ii,1}(1:mnd,4);
				min_comp(:,ii)=PLOTcell{ii,1}(1:mnd,2);
			end

			sbplt3=subplot(3,1,3);
			hold on
			plt(1)=plot(PLOTcell{1,1}(1:mnd,1),mean(mean_comp,2),'LineWidth',2,'Color','k');
			plt(2)=plot(PLOTcell{1,1}(1:mnd,1),max(mean_comp,[],2),'LineWidth',0.5,'Color','k');
			plot(PLOTcell{1,1}(1:mnd,1),min(mean_comp,[],2),'LineWidth',0.5,'Color','k');
			plt(3)=plot(PLOTcell{1,1}(1:mnd,1),mean(mean_comp,2)+std(mean_comp,0,2),'--','LineWidth',0.5,'Color','k');
			plot(PLOTcell{1,1}(1:mnd,1),mean(mean_comp,2)-std(mean_comp,0,2),'--','LineWidth',0.5,'Color','k');
			plt(4)=plot(PLOTcell{1,1}(1:mnd,1),max(max_comp,[],2),':','LineWidth',0.5,'Color','k');
			plot(PLOTcell{1,1}(1:mnd,1),min(min_comp,[],2),':','LineWidth',0.5,'Color','k');
			title('Mean Topography Along Swath')
			legend(plt,{'Mean of Mean','Min and Max of Mean','Std of Mean','All Extremes'},'location','best');
			xlabel('Distance Along Swaths (m)');
			ylabel('Elevation Along Swaths (m)');
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(sbplt3);
		    end 			
			hold off

		case 2

			sbplt1=subplot(2,1,1);
			hold on
			imagesc(DEM);
			for ii=1:number
				if ~verLessThan('matlab','9.4');
					plot(VERTcell{ii}(:,1),VERTcell{ii}(:,2),'Color',cmap(ii,:),'LineWidth',2);
				else
					plot(SWcell{ii,1}.xy0(:,1),SWcell{ii,1}.xy0(:,2),'Color',cmap(ii,:),'LineWidth',2);
				end
			end
			plot(points1(:,1),points1(:,2),'-k','LineWidth',2);
			plot(points2(:,1),points2(:,2),'-k','LineWidth',2);
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(sbplt1);
		    end 			
			hold off

			sbplt2=subplot(2,1,2);
			hold on 
			for ii=1:number
				plot(PLOTcell{ii,1}(:,1),PLOTcell{ii,1}(:,3),'LineWidth',2,'Color',cmap(ii,:));
				plot(PLOTcell{ii,1}(:,1),PLOTcell{ii,1}(:,2),'LineWidth',0.5,'Color',cmap(ii,:));
				plot(PLOTcell{ii,1}(:,1),PLOTcell{ii,1}(:,4),'LineWidth',0.5,'Color',cmap(ii,:));
			end
			title('All Elevation Swaths')
			xlabel('Distance Along Swaths (m)');
			ylabel('Elevation Along Swaths (m)');
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(sbplt2);
		    end 			
			hold off

		end

	end

	if plot_individual
		disp('Generating individual plots...')
		for ii=1:number
			f=figure(ii+1);
			clf
			set(f,'Units','normalized','Position',[0.05 0.1 0.6 0.6]);

			subplot(2+num_grids,1,1)
			hold on
			imagesc(DEM);
			if ~verLessThan('matlab','9.4')
				plot(VERTcell{ii}(:,1),VERTcell{ii}(:,2),'Color','k','LineWidth',2);
			else
				plot(SWcell{ii,1}.xy0(:,1),SWcell{ii,1}.xy0(:,2),'Color','k','LineWidth',2);
			end

			switch num_bounds
			case 1
				plot(points(:,1),points(:,2),'-k','LineWidth',2);
			case 2
				plot(points1(:,1),points1(:,2),'-k','LineWidth',2);
				plot(points2(:,1),points2(:,2),'-k','LineWidth',2);
			end
			title(['Swath Number ' num2str(ii)]);
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(gca);
		    end 			
			hold off

			subplot(2+num_grids,1,2)
			hold on 
			plot(PLOTcell{ii,1}(:,1),PLOTcell{ii,1}(:,3),'LineWidth',2,'Color','k');
			plot(PLOTcell{ii,1}(:,1),PLOTcell{ii,1}(:,2),'LineWidth',0.5,'Color','k');
			plot(PLOTcell{ii,1}(:,1),PLOTcell{ii,1}(:,4),'LineWidth',0.5,'Color','k');
			xlabel('Distance Along Swaths (m)');
			ylabel('Elevation Along Swaths (m)');
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(gca);
		    end 			
			hold off

			if ~isempty(AG)
				for jj=1:num_grids
					subplot(2+num_grids,1,2+jj)
					hold on 
					plot(PLOTcell{ii,1+jj}(:,1),PLOTcell{ii,1+jj}(:,3),'LineWidth',2,'Color','k');
					plot(PLOTcell{ii,1+jj}(:,1),PLOTcell{ii,1+jj}(:,2),'LineWidth',0.5,'Color','k');
					plot(PLOTcell{ii,1+jj}(:,1),PLOTcell{ii,1+jj}(:,4),'LineWidth',0.5,'Color','k');
					xlabel('Distance Along Swaths (m)');
					ylabel(AG{jj,2});
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(gca);
				    end 
					hold off
				end
			end
		end	
	end


	if make_shape
		% Make output shape
		ms=struct;
		switch num_bounds
		case 1
			ms(1,1).Geometry='Line';
			ms(1,1).X=points(:,1);
			ms(1,1).Y=points(:,2);
			ms(1,1).Type='SerialBase';

			for ii=1:number
				if ~verLessThan('matlab','9.4')	
					ms(ii+1,1).Geometry='Line';
					ms(ii+1,1).X=VERTcell{ii}(:,1);
					ms(ii+1,1).Y=VERTcell{ii}(:,2);
					ms(ii+1,1).Type=['Swth' num2str(ii) 'Box'];
				else
					ms(ii+1,1).Geometry='Line';
					ms(ii+1,1).X=SWcell{ii,1}.xy0(:,1);
					ms(ii+1,1).Y=SWcell{ii,1}.xy0(:,2);
					ms(ii+1,1).Type=['Swth' num2str(ii) 'Cntr'];
				end
			end
		case 2
			ms(1,1).Geometry='Line';
			ms(1,1).X=points1(:,1);
			ms(1,1).Y=points1(:,2);
			ms(1,1).Type='SerialBound1';

			ms(2,1).Geometry='Line';
			ms(2,1).X=points2(:,1);
			ms(2,1).Y=points2(:,2);
			ms(2,1).Type='SerialBound2';
			for ii=1:number
				if ~verLessThan('matlab','9.4')	
					ms(ii+2,1).Geometry='Line';
					ms(ii+2,1).X=VERTcell{ii}(:,1);
					ms(ii+2,1).Y=VERTcell{ii}(:,2);
					ms(ii+2,1).Type=['Swth' num2str(ii) 'Box'];
				else
					ms(ii+2,1).Geometry='Line';
					ms(ii+2,1).X=SWcell{ii,1}.xy0(:,1);
					ms(ii+2,1).Y=SWcell{ii,1}.xy0(:,2);
					ms(ii+2,1).Type=['Swth' num2str(ii) 'Cntr'];
				end
			end
		end
		shapewrite(ms,fullfile(out_dir,'SerialSwathBounds.shp'));
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








