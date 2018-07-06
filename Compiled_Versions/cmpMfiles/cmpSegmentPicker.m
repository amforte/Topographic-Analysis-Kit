function cmpSegmentPicker(wdir,MatFile,basin_num,varargin)
	% 	Function to select a segment of a stream network from the top of the stream, and plot the long profile 
	% 	and chi-Z relationship of that segment,also outputs the extraced portion of the stream network and chi structure 
	% 	(out of 'chiplot'). Allows user to iteratively select different parts of the stream network and display. 
	% 	Keeps running dataset of all the streams you pick and accept.
	%
	% Required Inputs:
	%	wdir - full path of working directory
	%	MatFile - Name of matfile output from either 'cmpMakeStreams' or the name of a single basin mat file from 'cmpProcessRiverBasins'
	%	basin_num - basin number from ProcessRiverBasins for output name or other identifying number for the set of streams you will pick
	%
	% Optional Inputs
	%	conditioned_DEM [] - option to provide full path of a hydrologically conditioned DEM for use in this function, expects the mat file as saved by 
	%		'cmpConditionDEM'. See 'cmpConditionDEM' function for options for making a hydrological conditioned DEM. If no input is provided the code defaults 
	%		to using the mincosthydrocon function.
	%	new_stream_net [] - option to provide name of a matfile containing a new stream network as as output from another function (e.g. cmpFindThreshold) to use
	%		instead of the stream network saved in the MatFile provided to the function. This new stream network must have been generated from the
	%		DEM stored in the provided MatFile
	%	direction ['down'] - expects either 'up' or 'down', default is 'down', if 'up' assumes individual selections are points above
	%		which you wish to extract and view stream profiles (i.e. a pour point), if 'down' assumes individual
	%		selections are channel heads if specific streams you wish to extract and view stream profiles. 
	%	method ['new_picks'] - expects either 'new_picks' or 'prev_picks', default is 'new_picks' if no input is provided. If 'prev_picks' is
	%			 given, the user must also supply an input for the 'picks' input (see below)
	%	plot_style ['refresh'] - expects either 'refresh' or 'keep', default is 'refresh' if no input is provided. If 'refresh' is given, the plots reset
	%			after each new stream pick, but if 'keep' is given, all selected streams remain on both the map (as thick red lines) and the
	%			chi-z/longitudinal profile/slope-area plots.
	%	plot_type ['vector'] - expects either 'vector' or 'grid', default is 'vector'. Controls whether all streams are drawn as individual lines ('vector') or if
	%		   the stream network is plotted as a grid and downsampled ('grid'). The 'grid' option is much faster with large datasets, 
	%			but can result in inaccurate choices. The 'vector' option is easier to see, but can be very slow to load and interact with.
	%	complete_networks_only [false] - if true, the code will filter out portions of the stream network that are incomplete prior to choosing
	%			streams
	%	picks - expects name of a textfile containing a n x 3 matrix with columns as x coordinates, y coordinates, and an identifying number OR the name of a point shapefile 
	%			with a single value column of identifying numbers. Will interpret this input as a list of channel heads if 'direction' is 'down' and a list of channel outlets 
	%			if 'direction' is 'up'.
	%	ref_concavity [0.50] - reference concavity for calculating Chi-Z, default is 0.50
	%	min_elev [] - minimum elevation below which the code stops extracting channel information, only used if 'direction'
	%			   is 'down'
	%	max_area [] - maximum drainage area above which the code stops extracting channel information, only used if 'direction'
	%			   is 'down'
	%	recalc [false] - only valid if either min_elev or max_area are specified. If recalc is false (default) then extraction of 
	%			 streams stops downstream of the condition specified in either min_elev or max_area, but chi is not recalculated 
	%			 and distances will remain tied to the original stream (i.e. distances from the outlet will be relative to the outlet
	%			 of the stream if it continued to the edge of the DEM, not where it stops extracting the stream profile). If recalc is true, 
	%	     	 then chi and distance are recalculated (i.e. the outlet as determined by the min_elev or max_area condition
	%			 will have a chi value of zero and a distance from mouth value of zero).
	%	threshold_area [1e6] - used to redraw downsampled stream network if 'plot_type' is set to 'grid'
	%	interp_value [0.1] - value (between 0 and 1) used for interpolation parameter in mincosthydrocon (not used if user provides a conditioned DEM)
	%	bin_size [500] - bin size (in map units) for binning slope area data.
	%
	% Outputs:
	%	Saves an output called 'PickedSegements_*.mat' with the provided basin number containing these results:
	%		StreamSgmnts - Cell array of selected stream segments as STREAMobj
	%		ChiSgmnts - Cell array of selected chi structures 
	%		SlpAreaSgmnts - Cell array of slope area data
	%		Sc - Single STREAMobj containing all the streams chosen.
	%		and if 'down' is selected:
	%			Heads - nx3 matrix of channel heads you picked with x cooord, y coordinate, and pick number as the columns
	%		and if 'up' is selected:
	%			Outlets - nx3 matrix of outlets you picked with x cooord, y coordinate, and pick number as the columns (valid input to 'ProcessRiverBasins'
	%			as 'river_mouths' parameter)
	%	Also saves a shapefile of the selected stream network
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
	p.FunctionName = 'cmpSegmentPicker';
	addRequired(p,'wdir',@(x) ischar(x));
	addRequired(p,'MatFile',@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.mat'))));
	addRequired(p,'basin_num',@(x) isnumeric(x));

	addParameter(p,'direction','down',@(x) ischar(validatestring(x,{'down','up'})));
	addParameter(p,'method','new_picks',@(x) ischar(validatestring(x,{'new_picks','prev_picks'})));
	addParameter(p,'plot_type','vector',@(x) ischar(validatestring(x,{'vector','grid'})));	
	addParameter(p,'plot_style','refresh',@(x) ischar(validatestring(x,{'refresh','keep'})));
	addParameter(p,'complete_networks_only',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'ref_concavity',0.50,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'min_elev',[],@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'max_area',[],@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'recalc',false,@(x) isscalar(x));
	addParameter(p,'picks',[],@(x) (isnumeric(x) && size(x,2)==3) | ischar(x));
	addParameter(p,'threshold_area',1e6,@(x) isnumeric(x));
	addParameter(p,'conditioned_DEM',[],@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.mat'))));
	addParameter(p,'new_stream_net',[],@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.mat'))));
	addParameter(p,'interp_value',0.1,@(x) isnumeric(x) && x>=0 && x<=1);
	addParameter(p,'bin_size',500,@(x) isscalar(x) && isnumeric(x));

	parse(p,wdir,MatFile,basin_num,varargin{:});
	wdir=p.Results.wdir;
	MatFile=p.Results.MatFile;
	basin_num=p.Results.basin_num;

	cno=p.Results.complete_networks_only;
	direction=p.Results.direction;
	method=p.Results.method;
	theta_ref=p.Results.ref_concavity;
	plot_type=p.Results.plot_type;
	plot_style=p.Results.plot_style;
	threshold_area=p.Results.threshold_area;
	points=p.Results.picks;
	iv=p.Results.interp_value;
	DEMc=p.Results.conditioned_DEM;
	nsn=p.Results.new_stream_net;
	bin_size=p.Results.bin_size;


	% Determine the type of input
	MatFile=fullfile(wdir,MatFile);
	D=whos('-file',MatFile);
	VL=cell(numel(D),1);
	for ii=1:numel(D);
		VL{ii}=D(ii,1).name;
	end

	if any(strcmp(VL,'DEM')) & any(strcmp(VL,'FD')) & any(strcmp(VL,'A')) & any(strcmp(VL,'S'))
		load(MatFile,'DEM','FD','A','S');
	elseif any(strcmp(VL,'DEMoc')) & any(strcmp(VL,'FDc')) & any(strcmp(VL,'Ac')) & any(strcmp(VL,'Sc'))
		load(MatFile,'DEMoc','FDc','Ac','Sc');
		DEM=DEMoc;
		FD=FDc;
		A=Ac;
		S=Sc;
	end	

	if ~isempty(nsn)
		nsn=fullfile(wdir,nsn);
		SD=whos('-file',nsn);
		VS=cell(numel(SD),1);
		for ii=1:numel(SD)
			VS{ii}=SD(ii,1).name;
		end

		if strcmp(VS,'Sn')
			load(nsn,'Sn');
			S=Sn;
		elseif strcmp(VS,'Sc')
			load(nsn,'Sc');
			S=Sc;
		elseif strcmp(VS,'S');
			load(nsn,'S');
		end

		if ~validatealignment(S,DEM)
			error('Supplied new stream net does not appear to have been generated from the same DEM supplied to the function');
		end
	end


	% Catch errors
	if strcmp(direction,'up') && ~isempty(p.Results.min_elev)
		warning('Input for minimum elevation is ignored when picking streams up from a pour point')
	elseif strcmp(direction,'up') && ~isempty(p.Results.max_area)
		warning('Input for maximum drainage area is ignored when picking streams up from a pour point')
	elseif ~isempty(p.Results.min_elev) && ~isempty(p.Results.max_area)
		error('Cannot specify both a minimum elevation and a maximum area, please provide one input only')
	elseif strcmp(method,'prev_picks') && isempty(p.Results.picks)
		error('If you choose the previous picks method you must provide a list of outlets or channel heads') 
	end

	% Parse different inputs
	if strcmp(direction,'down') & strcmp(plot_style,'keep');
		plot_switch='down_keep';
	elseif strcmp(direction,'down') & strcmp(plot_style,'refresh');
		plot_switch='down_ref';
	elseif strcmp(direction,'up') & strcmp(plot_style,'keep');
		plot_switch='up_keep';
	elseif strcmp(direction,'up') & strcmp(plot_style,'refresh');
		plot_switch='up_ref';
	end

    if strcmp(method,'prev_picks') & isempty(points)
        error('Please provide a m x 3 array of points or the valid path to a shapefile of points');
    elseif strcmp(method,'prev_picks') & ~isempty(regexp(points,regexptranslate('wildcard','*.shp')))
        disp('Reading shapefile')
    	points=fullfile(wdir,points);
        pt_shp=shaperead(points);
        fn=fieldnames(pt_shp);
        pt=horzcat([pt_shp.(fn{4})]',[pt_shp.X]',[pt_shp.Y]');
    elseif strcmp(method,'prev_picks') & ~isempty(regexp(points,regexptranslate('wildcard','*.txt')))
    	points=fullfile(wdir,points);
    	pT=readtable(points);
    	pt=[pT.(1) pT.(2) pT.(3)];
    end

    % Remove edges if flag is thrown
    if cno
    	S=removeedgeeffects(S,FD,DEM);
    end

	% Hydrologically condition dem
	if isempty(DEMc)
		zc=mincosthydrocon(S,DEM,'interp',iv);
		DEMc=GRIDobj(DEM);
		DEMc.Z(DEMc.Z==0)=NaN;
		DEMc.Z(S.IXgrid)=zc;
	else
		DEMc=fullfile(wdir,DEMc);
		load(DEMc,'DEMc');
		if ~validatealignment(DEMc,DEM)
			error('Provided conditioned DEM does not appear to align with the original DEM')
		end
	end

	switch plot_type
	case 'grid'
		DA=A.*(A.cellsize^2);
		DA.Z(DA.Z<threshold_area)=0;
		LA=log10(DA);
	end

	colcol=colorcube(25);
	% Strip out greys and whites
	colcol=colcol(1:20,:);
	
	switch method
	case 'new_picks'
		switch plot_switch
		% Extract points downstream from a channel head selection
		case 'down_keep'

			str1='N';
			str2='Y';

			ii=1;
			close all
			f1=figure(1);
			set(f1,'Units','normalized','Position',[0.05 0.1 0.45 0.8],'renderer','painters');
			clf
			switch plot_type
			case 'grid'
				hold on
				imageschs(DEM,LA,'colormap','parula','colorbar',false);
				hold off
			case 'vector'
				hold on
				imageschs(DEM,DEM,'colormap','parula','colorbar',false);
				plot(S,'-w');
				hold off
			end

			while strcmpi(str2,'Y')         
				while strcmpi(str1,'N')    	
					
					% Reset short circuit switch
					short_circ=0;

					figure(f1)
					hold on
					title('Choose point near channel head of interest')
					hold off
					[x,y]=ginput(1);
					pOI=[x y];

					% Find nearest channel head
					[ch]=streampoi(S,'channelheads','xy');
					distance=sqrt(sum(bsxfun(@minus, ch, pOI).^2,2));
					chOI=ch(distance==min(distance),:);

					% Build logical raster
					ix=coord2ind(DEM,chOI(:,1),chOI(:,2));
					IX=GRIDobj(DEM);
					IX.Z(ix)=1;
					[ixmat,X,Y]=GRIDobj2mat(IX);
					ixmat=logical(ixmat);
					IX=GRIDobj(X,Y,ixmat);

					% Extract stream from channel head to outlet
					Sn_t=modify(S,'downstreamto',IX);

					% Check if additional constraints have been specified
					if ~isempty(p.Results.max_area)
						AR=A.*(A.cellsize^2);
						ar=getnal(Sn_t,AR);

						sx=Sn_t.x; sy=Sn_t.y;
						d=Sn_t.distance;
						[d_s,d_ix]=sort(d,'ascend');
						ar_s=ar(d_ix);
						sx=sx(d_ix); sy=sy(d_ix);

						ma=p.Results.max_area;

						ix2=find(ar_s>=ma,1,'last');
						if isempty(ix2)
							warning('Input maximum drainage area is too large, extracting full stream')
							Sn=Sn_t;
							short_circ=1;
						elseif ix2==numel(ar)
							error('Input maximum drainage area is too small, no portion of the stream selected')
						else
							xn=sx(ix2);
							yn=sy(ix2);

							ix2=coord2ind(DEM,xn,yn);
							IX2=GRIDobj(DEM);
							IX2.Z(ix2)=1;
							[ix2mat,X,Y]=GRIDobj2mat(IX2);
							ix2mat=logical(ix2mat);
							IX2=GRIDobj(X,Y,ix2mat);

							Sn=modify(Sn_t,'upstreamto',IX2);
						end

					elseif ~isempty(p.Results.min_elev);
						el=getnal(Sn_t,DEMc);

						sx=Sn_t.x; sy=Sn_t.y;
						d=Sn_t.distance;
						[d_s,d_ix]=sort(d,'ascend');
						el_s=el(d_ix);
						sx=sx(d_ix); sy=sy(d_ix);

						me=p.Results.min_elev;

						ix2=find(el_s>=me,1,'first');
						if ix2==1
							warning('Input minimum elevation is too low, extracting full stream')
							Sn=Sn_t;
							short_circ=1;
						elseif isempty(ix2)
							error('Input minimum elevation is too high, no portion of the stream selected')
						else
							xn=sx(ix2);
							yn=sy(ix2);

							ix2=coord2ind(DEM,xn,yn);
							IX2=GRIDobj(DEM);
							IX2.Z(ix2)=1;
							[ix2mat,X,Y]=GRIDobj2mat(IX2);
							ix2mat=logical(ix2mat);
							IX2=GRIDobj(X,Y,ix2mat);

							Sn=modify(Sn_t,'upstreamto',IX2);
						end
					else
						Sn=Sn_t;
					end

					hold on
					SP=plot(Sn,'-r','LineWidth',2);
					hold off

		            qa=questdlg('Is this the stream segment you wanted?','Stream Selection','No','Yes','Yes');

		            switch qa
		            case 'Yes'
		                str1 = 'Y';
		            case 'No'
		                str1 = 'N';
		                delete(SP);
		            end
				end

				if isempty(p.Results.min_elev) && isempty(p.Results.max_area) 
					C=chiplot(Sn,DEMc,A,'a0',1,'mn',theta_ref,'plot',false);
				elseif ~isempty(p.Results.min_elev) | ~isempty(p.Results.max_area) && p.Results.recalc
					C=chiplot(Sn,DEMc,A,'a0',1,'mn',theta_ref,'plot',false);
				elseif short_circ==1;
					C=chiplot(Sn,DEMc,A,'a0',1,'mn',theta_ref,'plot',false);
				elseif ~isempty(p.Results.min_elev) | ~isempty(p.Results.max_area) && p.Results.recalc==0
					C_t=chiplot(Sn_t,DEMc,A,'a0',1,'mn',theta_ref,'plot',false);
					C_n=chiplot(Sn,DEMc,A,'a0',1,'mn',theta_ref,'plot',false);

					% Find extents
					txyz=[C_t.x C_t.y C_t.elev];
					nxyz=[C_n.x C_n.y C_n.elev];
					ix3=ismember(txyz,nxyz,'rows');
					% Remake chi structure
					C=struct;
					C.mn=C_n.mn;
					C.beta=C_n.beta;
					C.betase=C_n.betase;
					C.a0=C_n.a0;
					C.ks=C_n.ks;
					C.R2=C_n.R2;
					C.chi=C_t.chi(ix3);
					C.x=C_t.x(ix3); C.y=C_t.y(ix3);
					C.elev=C_t.elev(ix3); C.elevbl=C_t.elevbl(ix3);
					C.distance=C_t.distance(ix3); C.pred=C_t.pred(ix3);
					C.area=C_t.area(ix3); C.res=C_t.res(ix3);
				end

				f2=figure(2);
				set(f2,'Units','normalized','Position',[0.5 0.1 0.45 0.8],'renderer','painters');

				subplot(3,1,1);
				hold on
				plot(C.chi,C.elev,'Color',colcol(mod(ii,20)+1,:));
				xlabel('\chi')
				ylabel('Elevation (m)')
				title('\chi - Z')
				hold off

				subplot(3,1,2);
				hold on
				if isempty(p.Results.min_elev) && isempty(p.Results.max_area) 
					plotdz(Sn,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
					plotdz(Sn,DEMc,'dunit','km','Color',colcol(mod(ii,20)+1,:));
				elseif ~isempty(p.Results.min_elev) | ~isempty(p.Results.max_area) && p.Results.recalc
					plotdz(Sn,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
					plotdz(Sn,DEMc,'dunit','km','Color',colcol(mod(ii,20)+1,:));
				elseif short_circ==1;
					plotdz(Sn,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
					plotdz(Sn,DEMc,'dunit','km','Color',colcol(mod(ii,20)+1,:));
				elseif ~isempty(p.Results.min_elev) | ~isempty(p.Results.max_area) && p.Results.recalc==0
					Cu=chiplot(Sn_t,DEM,A,'a0',1,'mn',theta_ref,'plot',false);
					plot((Cu.distance(ix3))./1000,Cu.elev(ix3),'Color',[0.5 0.5 0.5]);
					% plot(C.distance./1000,C.elev,'-k');
					plot(C.distance./1000,C.elev,'Color',colcol(mod(ii,20)+1,:));
				end
				xlabel('Distance from Mouth (km)')
				ylabel('Elevation (m)')
				% legend('Unconditioned DEM','Conditioned DEM','location','best');
				title('Long Profile')
				hold off

				saax=subplot(3,1,3);
				hold on
				if isempty(p.Results.min_elev) && isempty(p.Results.max_area) 
					[bs,ba,aa,ag]=sa(DEMc,Sn,A,bin_size);
				elseif ~isempty(p.Results.min_elev) | ~isempty(p.Results.max_area) && p.Results.recalc								
					[bs,ba,aa,ag]=sa(DEMc,Sn,A,bin_size);
				elseif short_circ==1;
					[bs,ba,aa,ag]=sa(DEMc,Sn,A,bin_size);
				elseif ~isempty(p.Results.min_elev) | ~isempty(p.Results.max_area) && p.Results.recalc==0
					[bs,ba,aa,ag]=sa(DEMc,Sn_t,A,bin_size);
				end

				scatter(aa,ag,5,[0.5 0.5 0.5],'+');
				scatter(ba,bs,20,colcol(mod(ii,20)+1,:),'filled');	
				set(saax,'Xscale','log','Yscale','log','XDir','reverse');
				xlabel('Log Drainage Area');
				ylabel('Log Slope');	
				hold off

				StreamSgmnts{ii}=Sn;
				ChiSgmnts{ii}=C;
				SlpAreaSgmnts{ii,1}=[bs ba];
				SlpAreaSgmnts{ii,2}=[ag aa];				
				Heads(ii,1)=chOI(:,1);
				Heads(ii,2)=chOI(:,2);
				Heads(ii,3)=ii;

				ii=ii+1;

	            qa2=questdlg('Continue picking streams?','Stream Selection','No','Yes','Yes');
	            switch qa2
	            case 'Yes'
	            	str2 = 'Y';
	                str1 = 'N';
	            case 'No'
	                str2 = 'N';
	            end				
			end

		% Extract segements upstream from a pour point selection
		case 'up_keep'

			str1='N';
			str2='Y';

			ii=1;
	
			f1=figure(1);
			set(f1,'Units','normalized','Position',[0.05 0.1 0.45 0.8],'renderer','painters');
			clf
			switch plot_type
			case 'grid'
				hold on
				imageschs(DEM,LA,'colormap','parula','colorbar',false);
				hold off
			case 'vector'
				hold on
				imageschs(DEM,DEM,'colormap','parula','colorbar',false);
				plot(S,'-w');
				hold off
			end

			while strcmpi(str2,'Y')       
				while strcmpi(str1,'N')   	

					figure(f1)

					hold on
					title('Choose point above to which calculate Chi-Z')
					hold off

					[x,y]=ginput(1);

					% Build logical raster
					[xn,yn]=snap2stream(S,x,y);
					ix=coord2ind(DEM,xn,yn);
					IX=GRIDobj(DEM);
					IX.Z(ix)=1;
					[ixmat,X,Y]=GRIDobj2mat(IX);
					ixmat=logical(ixmat);
					IX=GRIDobj(X,Y,ixmat);

					Sn=modify(S,'upstreamto',IX);

					hold on
					SP=plot(Sn,'-r','LineWidth',2);
					hold off

		            qa=questdlg('Is this the stream segment you wanted?','Stream Selection','No','Yes','Yes');
		            switch qa
		            case 'Yes'
		                str1 = 'Y';
		            case 'No'
		                str1 = 'N';
		                delete(SP);
		            end
				end

				C=chiplot(Sn,DEMc,A,'a0',1,'mn',theta_ref,'plot',false);

				f2=figure(2);
				set(f2,'Units','normalized','Position',[0.5 0.1 0.45 0.8],'renderer','painters');

				subplot(3,1,1);
				hold on
				plot(C.chi,C.elev,'Color',colcol(mod(ii,20)+1,:));
				xlabel('\chi')
				ylabel('Elevation (m)')
				title('\chi - Z')
				hold off

				subplot(3,1,2);
				hold on
				plotdz(Sn,DEMc,'dunit','km','Color',colcol(mod(ii,20)+1,:));
				xlabel('Distance from Mouth (km)')
				ylabel('Elevation (m)')
				title('Long Profile')
				hold off

				saax=subplot(3,1,3);
				hold on
				[bs,ba,aa,ag]=sa(DEMc,Sn,A,bin_size);
				scatter(aa,ag,5,[0.5 0.5 0.5],'+');
				scatter(ba,bs,20,colcol(mod(ii,20)+1,:),'filled');
				set(saax,'Xscale','log','Yscale','log','XDir','reverse');
				xlabel('Log Drainage Area');
				ylabel('Log Slope');	
				hold off

				StreamSgmnts{ii}=Sn;
				ChiSgmnts{ii}=C;
				SlpAreaSgmnts{ii,1}=[bs ba];
				SlpAreaSgmnts{ii,2}=[ag aa];	
				Outlets(ii,1)=xn;
				Outlets(ii,2)=yn;
				Outlets(ii,3)=ii;

				ii=ii+1;


	            qa2=questdlg('Continue picking streams?','Stream Selection','No','Yes','Yes');
	            switch qa2
	            case 'Yes'
	            	str2 = 'Y';
	                str1 = 'N';
	            case 'No'
	                str2 = 'N';
	            end
			end

		case 'down_ref'

			str1='N';
			str2='Y';

			ii=1;

			f1=figure(1);
			set(f1,'Units','normalized','Position',[0.05 0.1 0.45 0.8],'renderer','painters');
			clf
			switch plot_type
			case 'grid'
				hold on
				imageschs(DEM,LA,'colormap','parula','colorbar',false);
				hold off
			case 'vector'
				hold on
				imageschs(DEM,DEM,'colormap','parula','colorbar',false);
				plot(S,'-w');
				hold off
			end

			while strcmpi(str2,'Y')     
				while strcmpi(str1,'N')    	
					
					% Reset short circuit switch
					short_circ=0;

					figure(1)

					hold on
					title('Choose point near channel head of interest')
					hold off
					[x,y]=ginput(1);
					pOI=[x y];

					% Find nearest channel head
					[ch]=streampoi(S,'channelheads','xy');
					distance=sqrt(sum(bsxfun(@minus, ch, pOI).^2,2));
					chOI=ch(distance==min(distance),:);

					% Build logical raster
					ix=coord2ind(DEM,chOI(:,1),chOI(:,2));
					IX=GRIDobj(DEM);
					IX.Z(ix)=1;
					[ixmat,X,Y]=GRIDobj2mat(IX);
					ixmat=logical(ixmat);
					IX=GRIDobj(X,Y,ixmat);

					% Extract stream from channel head to outlet
					Sn_t=modify(S,'downstreamto',IX);

					% Check if additional constraints have been specified
					if ~isempty(p.Results.max_area)
						AR=A.*(A.cellsize^2);
						ar=getnal(Sn_t,AR);

						sx=Sn_t.x; sy=Sn_t.y;
						d=Sn_t.distance;
						[d_s,d_ix]=sort(d,'ascend');
						ar_s=ar(d_ix);
						sx=sx(d_ix); sy=sy(d_ix);

						ma=p.Results.max_area;

						ix2=find(ar_s>=ma,1,'last');
						if isempty(ix2)
							warning('Input maximum drainage area is too large, extracting full stream')
							Sn=Sn_t;
							short_circ=1;
						elseif ix2==numel(ar)
							error('Input maximum drainage area is too small, no portion of the stream selected')
						else
							xn=sx(ix2);
							yn=sy(ix2);

							ix2=coord2ind(DEM,xn,yn);
							IX2=GRIDobj(DEM);
							IX2.Z(ix2)=1;
							[ix2mat,X,Y]=GRIDobj2mat(IX2);
							ix2mat=logical(ix2mat);
							IX2=GRIDobj(X,Y,ix2mat);

							Sn=modify(Sn_t,'upstreamto',IX2);
						end

					elseif ~isempty(p.Results.min_elev);
						el=getnal(Sn_t,DEMc);

						sx=Sn_t.x; sy=Sn_t.y;
						d=Sn_t.distance;
						[d_s,d_ix]=sort(d,'ascend');
						el_s=el(d_ix);
						sx=sx(d_ix); sy=sy(d_ix);

						me=p.Results.min_elev;

						ix2=find(el_s>=me,1,'first');
						if ix2==1
							warning('Input minimum elevation is too low, extracting full stream')
							Sn=Sn_t;
							short_circ=1;
						elseif isempty(ix2)
							error('Input minimum elevation is too high, no portion of the stream selected')
						else
							xn=sx(ix2);
							yn=sy(ix2);

							ix2=coord2ind(DEM,xn,yn);
							IX2=GRIDobj(DEM);
							IX2.Z(ix2)=1;
							[ix2mat,X,Y]=GRIDobj2mat(IX2);
							ix2mat=logical(ix2mat);
							IX2=GRIDobj(X,Y,ix2mat);

							Sn=modify(Sn_t,'upstreamto',IX2);
						end
					else
						Sn=Sn_t;
					end

					hold on
					SP=plot(Sn,'-r','LineWidth',2);
					hold off


		            qa=questdlg('Is this the stream segment you wanted?','Stream Selection','No','Yes','Yes');
		            switch qa
		            case 'Yes'
		                str1 = 'Y';
		            case 'No'
		                str1 = 'N';
		                delete(SP);
		            end

				end

				if isempty(p.Results.min_elev) && isempty(p.Results.max_area) 
					C=chiplot(Sn,DEMc,A,'a0',1,'mn',theta_ref,'plot',false);
				elseif ~isempty(p.Results.min_elev) | ~isempty(p.Results.max_area) && p.Results.recalc
					C=chiplot(Sn,DEMc,A,'a0',1,'mn',theta_ref,'plot',false);
				elseif short_circ==1;
					C=chiplot(Sn,DEMc,A,'a0',1,'mn',theta_ref,'plot',false);
				elseif ~isempty(p.Results.min_elev) | ~isempty(p.Results.max_area) && p.Results.recalc==0
					C_t=chiplot(Sn_t,DEMc,A,'a0',1,'mn',theta_ref,'plot',false);
					C_n=chiplot(Sn,DEMc,A,'a0',1,'mn',theta_ref,'plot',false);

					% Find extents
					txyz=[C_t.x C_t.y C_t.elev];
					nxyz=[C_n.x C_n.y C_n.elev];
					ix3=ismember(txyz,nxyz,'rows');
					% Remake chi structure
					C=struct;
					C.mn=C_n.mn;
					C.beta=C_n.beta;
					C.betase=C_n.betase;
					C.a0=C_n.a0;
					C.ks=C_n.ks;
					C.R2=C_n.R2;
					C.chi=C_t.chi(ix3);
					C.x=C_t.x(ix3); C.y=C_t.y(ix3);
					C.elev=C_t.elev(ix3); C.elevbl=C_t.elevbl(ix3);
					C.distance=C_t.distance(ix3); C.pred=C_t.pred(ix3);
					C.area=C_t.area(ix3); C.res=C_t.res(ix3);
				end

				f2=figure(2);
				set(f2,'Units','normalized','Position',[0.5 0.1 0.45 0.8],'renderer','painters');

				subplot(3,1,1);
				hold on
				plot(C.chi,C.elev,'-k');
				xlabel('\chi')
				ylabel('Elevation (m)')
				title('\chi - Z')
				hold off

				subplot(3,1,2);
				hold on
				if isempty(p.Results.min_elev) && isempty(p.Results.max_area) 
					plotdz(Sn,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
					plotdz(Sn,DEMc,'dunit','km','Color','k');
				elseif ~isempty(p.Results.min_elev) | ~isempty(p.Results.max_area) && p.Results.recalc
					plotdz(Sn,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
					plotdz(Sn,DEMc,'dunit','km','Color','k');
				elseif short_circ==1;
					plotdz(Sn,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
					plotdz(Sn,DEMc,'dunit','km','Color','k');
				elseif ~isempty(p.Results.min_elev) | ~isempty(p.Results.max_area) && p.Results.recalc==0
					Cu=chiplot(Sn_t,DEM,A,'a0',1,'mn',theta_ref,'plot',false);
					plot((Cu.distance(ix3))./1000,Cu.elev(ix3),'Color',[0.5 0.5 0.5]);
					plot(C.distance./1000,C.elev,'-k');
				end
				xlabel('Distance from Mouth (km)')
				ylabel('Elevation (m)')
				legend('Unconditioned DEM','Conditioned DEM','location','best');
				title('Long Profile')
				hold off


				saax=subplot(3,1,3);
				hold on
				if isempty(p.Results.min_elev) && isempty(p.Results.max_area) 
					[bs,ba,aa,ag]=sa(DEMc,Sn,A,bin_size);
				elseif ~isempty(p.Results.min_elev) | ~isempty(p.Results.max_area) && p.Results.recalc								
					[bs,ba,aa,ag]=sa(DEMc,Sn,A,bin_size);
				elseif short_circ==1;
					[bs,ba,aa,ag]=sa(DEMc,Sn,A,bin_size);
				elseif ~isempty(p.Results.min_elev) | ~isempty(p.Results.max_area) && p.Results.recalc==0
					[bs,ba,aa,ag]=sa(DEMc,Sn_t,A,bin_size);
				end		
				scatter(aa,ag,5,[0.5 0.5 0.5],'+');
				scatter(ba,bs,20,'k','filled');
				set(saax,'Xscale','log','Yscale','log','XDir','reverse');
				xlabel('Log Drainage Area');
				ylabel('Log Slope');	
				hold off								

				StreamSgmnts{ii}=Sn;
				ChiSgmnts{ii}=C;
				SlpAreaSgmnts{ii,1}=[bs ba];
				SlpAreaSgmnts{ii,2}=[ag aa];	
				Heads(ii,1)=chOI(:,1);
				Heads(ii,2)=chOI(:,2);
				Heads(ii,3)=ii;

				ii=ii+1;

	            qa2=questdlg('Continue picking streams?','Stream Selection','No','Yes','Yes');
	            switch qa2
	            case 'Yes'
	            	str2 = 'Y';
	                str1 = 'N';
	                close figure 2
	            case 'No'
	                str2 = 'N';
	            end	

			end

		% Extract segements upstream from a pour point selection
		case 'up_ref'

			str1='N';
			str2='Y';

			ii=1;
			f1=figure(1);
			set(f1,'Units','normalized','Position',[0.05 0.1 0.45 0.8],'renderer','painters');
			clf
			switch plot_type
			case 'grid'
				hold on
				imageschs(DEM,LA,'colormap','parula','colorbar',false);
				hold off
			case 'vector'
				hold on
				imageschs(DEM,DEM,'colormap','parula','colorbar',false);
				plot(S,'-w');
				hold off
			end

			while strcmpi(str2,'Y')        
				while strcmpi(str1,'N')   	
					
					figure(1);

					hold on
					title('Choose point above to which calculate Chi-Z')
					hold off
					[x,y]=ginput(1);

					% Build logical raster
					[xn,yn]=snap2stream(S,x,y);
					ix=coord2ind(DEM,xn,yn);
					IX=GRIDobj(DEM);
					IX.Z(ix)=1;
					[ixmat,X,Y]=GRIDobj2mat(IX);
					ixmat=logical(ixmat);
					IX=GRIDobj(X,Y,ixmat);

					Sn=modify(S,'upstreamto',IX);

					hold on
					SP=plot(Sn,'-r','LineWidth',2);
					hold off

		            qa=questdlg('Is this the stream segment you wanted?','Stream Selection','No','Yes','Yes');
		            switch qa
		            case 'Yes'
		                str1 = 'Y';
		            case 'No'
		                str1 = 'N';
		                delete(SP);
		            end
				end

				C=chiplot(Sn,DEMc,A,'a0',1,'mn',theta_ref,'plot',false);

				f2=figure(2);
				set(f2,'Units','normalized','Position',[0.5 0.1 0.45 0.8],'renderer','painters');
				subplot(3,1,1);
				hold on
				plot(C.chi,C.elev,'-k');
				xlabel('\chi')
				ylabel('Elevation (m)')
				title('\chi - Z')
				hold off

				subplot(3,1,2);
				hold on
				plotdz(Sn,DEMc,'dunit','km','Color','k');
				xlabel('Distance from Mouth (km)')
				ylabel('Elevation (m)')
				title('Long Profile')
				hold off

				saax=subplot(3,1,3);
				hold on
				[bs,ba,aa,ag]=sa(DEMc,Sn,A,bin_size);
				scatter(aa,ag,5,[0.5 0.5 0.5],'+');
				scatter(ba,bs,20,'k','filled');
				set(saax,'Xscale','log','Yscale','log','XDir','reverse');
				xlabel('Log Drainage Area');
				ylabel('Log Slope');	
				hold off

				StreamSgmnts{ii}=Sn;
				ChiSgmnts{ii}=C;
				SlpAreaSgmnts{ii,1}=[bs ba];
				SlpAreaSgmnts{ii,2}=[ag aa];	
				Outlets(ii,1)=xn;
				Outlets(ii,2)=yn;
				Outlets(ii,3)=ii;

				ii=ii+1;

	            qa2=questdlg('Continue picking streams?','Stream Selection','No','Yes','Yes');
	            switch qa2
	            case 'Yes'
	            	str2 = 'Y';
	                str1 = 'N';
	                close figure 2
	            case 'No'
	                str2 = 'N';
	            end	

			end
		end
		
	% Previous Picks
	case 'prev_picks'
		switch direction
		case 'down'
			heads=pt;
			[num_heads,~]=size(heads);

			w1=waitbar(0,'Extracting segments');
			for ii=1:num_heads
				short_circ=0;
				x=heads(ii,1); y=heads(ii,2);
				pOI=[x y];

				% Find nearest channel head
				[ch]=streampoi(S,'channelheads','xy');
				distance=sqrt(sum(bsxfun(@minus, ch, pOI).^2,2));
				chOI=ch(distance==min(distance),:);

				% Build logical raster
				ix=coord2ind(DEM,chOI(:,1),chOI(:,2));
				IX=GRIDobj(DEM);
				IX.Z(ix)=1;
				[ixmat,X,Y]=GRIDobj2mat(IX);
				ixmat=logical(ixmat);
				IX=GRIDobj(X,Y,ixmat);

				% Extract stream from channel head to outlet
				Sn_t=modify(S,'downstreamto',IX);

				% Check if additional constraints have been specified
				if ~isempty(p.Results.max_area)
					AR=A.*(A.cellsize^2);

					sx=Sn_t.x; sy=Sn_t.y;
					d=Sn_t.distance;
					[d_s,d_ix]=sort(d,'ascend');
					ar_s=ar(d_ix);
					sx=sx(d_ix); sy=sy(d_ix);

					ix2=find(ar_s>=ma,1,'last');
					if isempty(ix2)
						warning('Input maximum drainage area is too large, extracting full stream')
						Sn=Sn_t;
						short_circ=1;
					elseif ix2==numel(ar)
						error('Input maximum drainage area is too small, no portion of the stream selected')
					else
						xn=sx(ix2);
						yn=sy(ix2);

						ix2=coord2ind(DEM,xn,yn);
						IX2=GRIDobj(DEM);
						IX2.Z(ix2)=1;
						[ix2mat,X,Y]=GRIDobj2mat(IX2);
						ix2mat=logical(ix2mat);
						IX2=GRIDobj(X,Y,ix2mat);

						Sn=modify(Sn_t,'upstreamto',IX2);
					end

				elseif ~isempty(p.Results.min_elev);
					el=getnal(Sn_t,DEMc);

					sx=Sn_t.x; sy=Sn_t.y;
					d=Sn_t.distance;
					[d_s,d_ix]=sort(d,'ascend');
					el_s=el(d_ix);
					sx=sx(d_ix); sy=sy(d_ix);

					me=p.Results.min_elev;

					ix2=find(el_s>=me,1,'first');
					if ix2==1
						warning('Input minimum elevation is too low, extracting full stream')
						Sn=Sn_t;
					elseif isempty(ix2)
						error('Input minimum elevation is too high, no portion of the stream selected')
					else
						xn=sx(ix2);
						yn=sy(ix2);

						ix2=coord2ind(DEM,xn,yn);
						IX2=GRIDobj(DEM);
						IX2.Z(ix2)=1;
						[ix2mat,X,Y]=GRIDobj2mat(IX2);
						ix2mat=logical(ix2mat);
						IX2=GRIDobj(X,Y,ix2mat);

						Sn=modify(Sn_t,'upstreamto',IX2);
					end
				else
					Sn=Sn_t;
				end

				[bs,ba,aa,ag]=sa(DEMc,Sn,A,bin_size);

				if isempty(p.Results.min_elev) && isempty(p.Results.max_area) 
					C=chiplot(Sn,DEMc,A,'a0',1,'mn',theta_ref,'plot',false);
				elseif ~isempty(p.Results.min_elev) | ~isempty(p.Results.max_area) && p.Results.recalc
					C=chiplot(Sn,DEMc,A,'a0',1,'mn',theta_ref,'plot',false);
				elseif short_circ==1;
					C=chiplot(Sn,DEMc,A,'a0',1,'mn',theta_ref,'plot',false);
				elseif ~isempty(p.Results.min_elev) | ~isempty(p.Results.max_area) && p.Results.recalc==0
					C_t=chiplot(Sn_t,DEMc,A,'a0',1,'mn',theta_ref,'plot',false);
					C_n=chiplot(Sn,DEMc,A,'a0',1,'mn',theta_ref,'plot',false);

					% Find extents
					txyz=[C_t.x C_t.y C_t.elev];
					nxyz=[C_n.x C_n.y C_n.elev];
					ix3=ismember(txyz,nxyz,'rows');

					% Remake chi structure
					C=C_t;
					C.chi=C_t.chi(ix3);
					C.x=C_t.x(ix3); C.y=C_t.y(ix3);
					C.elev=C_t.elev(ix3); C.elevbl=C_t.elevbl(ix3);
					C.distance=C_t.distance(ix3); C.pred=C_t.pred(ix3);
					C.area=C_t.area(ix3); C.res=C_t.res(ix3);
				end

				StreamSgmnts{ii}=Sn;
				ChiSgmnts{ii}=C;
				SlpAreaSgmnts{ii,1}=[bs ba];
				SlpAreaSgmnts{ii,2}=[ag aa];	
				Heads(ii,1)=chOI(:,1);
				Heads(ii,2)=chOI(:,2);
				Heads(ii,3)=ii;
				waitbar(ii/num_heads);
			end
			close(w1);

		case 'up'
			outlets=pt;
			[num_outs,~]=size(outlets);

			w1=waitbar(0,'Extracting segments');
			for ii=1:num_outs
				x=outlets(ii,1); y=outlets(ii,2);
				% Build logical raster
				[xn,yn]=snap2stream(S,x,y);
				ix=coord2ind(DEM,xn,yn);
				IX=GRIDobj(DEM);
				IX.Z(ix)=1;
				[ixmat,X,Y]=GRIDobj2mat(IX);
				ixmat=logical(ixmat);
				IX=GRIDobj(X,Y,ixmat);


				Sn=modify(S,'upstreamto',IX);
				Sn=klargestconncomps(Sn,1);
				C=chiplot(Sn,DEMc,A,'a0',1,'mn',theta_ref,'plot',false);
				[bs,ba,aa,ag]=sa(DEMc,Sn,A,bin_size);

				StreamSgmnts{ii}=Sn;
				ChiSgmnts{ii}=C;
				SlpAreaSgmnts{ii,1}=[bs ba];
				SlpAreaSgmnts{ii,2}=[ag aa];	
				Outlets(ii,1)=xn;
				Outlets(ii,2)=yn;
				Outlets(ii,3)=ii;
				waitbar(ii/num_outs);
			end
			close(w1);
		end
	end

	% Clean up and generate outputs
	num_picks=numel(StreamSgmnts);
	if num_picks==1
		Sc=StreamSgmnts{1};
	else
		Sc=StreamSgmnts{1};
		for ii=2:num_picks
			Sc=union(Sc,StreamSgmnts{ii});
		end
	end

	fileOut=fullfile(wdir,['PickedSegments_' num2str(basin_num) '.mat']);
	switch direction
	case 'up'
		save(fileOut,'StreamSgmnts','ChiSgmnts','SlpAreaSgmnts','Outlets','Sc');
	case 'down'
		save(fileOut,'StreamSgmnts','ChiSgmnts','SlpAreaSgmnts','Heads','Sc');	
	end

	MSn=STREAMobj2mapstruct(Sc);
	shapewrite(MSn,fullfile(wdir,'SelectedStreams.shp'));

% Main Function End
end

function [bs,ba,a,g]=sa(DEM,S,A,bin_size)
	% Modified slope area function that uses the smooth length to
	%	to determine the number of bins and uses those same bins
	%	to find mean values of chi and distance for plotting
	%	purposes

	minX=min(S.distance);
	maxX=max(S.distance);
	b=[minX:bin_size:maxX+bin_size];

	numbins=round(max([numel(b) numel(S.IXgrid)/10]));

	an=getnal(S,A.*A.cellsize^2);
	z=getnal(S,DEM);
	gn=gradient(S,z,'unit','tangent','method','robust','drop',20);
	gn=smooth(gn,3);

	% Run through STREAMobj2XY so chi and everything else are same size
	[~,~,a,g]=STREAMobj2XY(S,an,gn);
	% Remove NaNs
	a(isnan(a))=[];
	g(isnan(g))=[];

	mina=min(a);
	maxa=max(a);

    edges = logspace(log10(mina-0.1),log10(maxa+1),numbins+1);
    try
    	% histc is deprecated
    	[ix]=discretize(a,edges);
    catch
	    [~,ix] = histc(a,edges);
	end

	ba=accumarray(ix,a,[numbins 1],@median,nan);
	bs=accumarray(ix,g,[numbins 1],@(x) mean(x(~isnan(x))),nan);

	% Filter negatives
	idx=bs>=0 & ba>=0;
	bs=bs(idx);
	ba=ba(idx);

	idx=a>=0 & g>=0;
	a=a(idx);
	g=g(idx);
end





