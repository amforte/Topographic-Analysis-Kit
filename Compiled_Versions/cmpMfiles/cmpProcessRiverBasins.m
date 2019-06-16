function cmpProcessRiverBasins(wdir,MakeStreamsMat,river_mouths,basin_dir,varargin)
	% Function takes grid object outputs from MakeStreams script (DEM,FD,A,S), a series of x,y coordinates of river mouths,
	% 	and outputs clipped dem, stream network, various topographic metrics, and river values (ks, ksn, chi). Saves mat files
	% 	for use in other codes, if you want to produce files that are usable in an outside GIS program, be sure to set 
	% 	'write_arc_files' to true.
	
	% Required Inputs:
	% 		wdir - full path of working directory
	% 		MakeStreamsMat - name or location of the mat file saved after running 'cmpMakeStreams'
	% 		river_mouths - locations of river mouths (i.e. pour points) above which you wish to extract basins, can take one of three forms:
	% 			1) name of a text file containing a nx3 array of river mouths with x, y, and a number identifying the stream/basin of interest 
	% 				(must be same projection as DEM) saved as a '.txt'. No extra columns should be present, the code should work with or without headers as 
	% 				long as headers are restricted to a single line. The text file of Outlets saved by BasinPicker, can be used as the river mouths
	% 			2) a single value that will be interpreted as an elevation that  the code will use this to autogenerate river mouths at this elevation.
	% 			3) point shapefile with one numeric user input field (e.g. the default 'ID' field generated by ArcGIS) that will be used as the
	% 				river mouth ID (must be same projection as DEM).
	% 		basin_dir - name of folder to store basin files (if specified folder does not exist in current directory, code will create it)
	
	% Optional Inputs:
	% 		conditioned_DEM [] - option to provide a hydrologically conditioned DEM for use in this function, expects an ascii grid as saved by 'cmpConditionDEM'
	% 			See 'cmpConditionDEM' function for options for making a hydrological conditioned DEM. If no input is provided the code defaults to using the 
	% 			mincosthydrocon function.
	% 		interp_value [0.1] - value (between 0 and 1) used for interpolation parameter in mincosthydrocon (not used if user provides a conditioned DEM)
	% 		threshold_area [1e6] - minimum accumulation area to define streams in meters squared
	% 		segment_length [1000] - smoothing distance in meters for averaging along ksn, suggested value is 1000 meters
	% 		ref_concavity [0.5] - reference concavity for calculating ksn, suggested value is 0.45
	% 		ksn_method [quick] - switch between method to calculate ksn values, options are 'quick', 'trunk', or 'trib', the 'trib' method takes 3-4 times longer 
	% 			than the 'quick' method. In most cases, the 'quick' method works well, but if values near tributary junctions are important, then 'trib'
	% 			may be better as this calculates ksn values for individual channel segments individually. The 'trunk' option calculates steepness values
	% 			of large streams independently (streams considered as trunks are controlled by the stream order value supplied to 'min_order'). The 'trunk' option
	% 			may be of use if you notice anomaoloulsy high channel steepness values on main trunk streams that can result because of the way values are reach
	% 			averaged.
	% 		min_order [4] - minimum stream order for a stream to be considered a trunk stream, only used if 'ksn_method' is set to 'trunk'
	% 		write_arc_files [false] - set value to true to output a ascii's of various grids and a shapefile of the ksn, false to not output arc files
	% 		add_grids [] - option to provide the name of the mat file produced by running 'cmpPrepareAddGrids'. Use this function if you want this function 
	% 			to calculate statistics for additional grids (e.g. precipitation).
	% 		add_cat_grids [] - option provide the name of the mat file produced by running 'cmpPrepareCatAddGrids'. Use this if you want to calculate statisitcs
	% 			related to categorical data stored in a shapefile (e.g. geologic map).
	% 		resample_method ['nearest'] - method to use in the resample function on additional grids (if required). Acceptable inputs are 'nearest', 'bilinear', 
	% 			or 'bicubic'. Method 'nearest' is appropriate if you do not want the resampling to interpolate between values (e.g. if an additinal grid has specific values
	% 			that correlate to a property like rock type) and either 'bilinear' or 'bicubic' is appropriate if you want smooth variations between nodes. 
	% 		gradient_method ['arcslope'] - function used to calculate gradient, either 'arcslope' (default) or 'gradient8'. The 'arcslope' function calculates
	% 			gradient the same way as ArcGIS by fitting a plane to the 8-connected neighborhood and 'gradient8' returns the steepest descent for the same
	% 			8-connected neighborhood. 'gradient8' will generally return higher values than 'arcslope'.
	% 		calc_relief [false] - option to calculate local relief. Can provide an array of radii to use with 'relief_radii' option.
	% 		relief_radii [2500] - a 1d vector (column or row) of radii to use for calculating local relief, values must be in map units. If more than one value is provided
	% 			the function assumes you wish to calculate relief at all of these radii. Note, the local relief function is slow so providing multiple radii will
	% 			slow code performance. Saved outputs will be in a m x 2 cell array, with the columns of the cell array corresponding to the GRIDobj and the input radii.
	%		ksn_radius [5000] - radius of circular, moving area over which to average ksn values for making an interpolated ksn grid
	%
	% Notes:
	% 		-The code will perform a check of the river_mouths input to confirm that 1) there are no duplicate ID numbers (it will dump your ID numbers and create new
	% 			ID numbers if this is the case and output a text file contatining the river mouth locations with their new ID nubmers) and 2) that no provided river mouths 
	% 			are outside the boundaries of the DEM (it will remove these IDs if this the case).
	%
 	%  	Examples if running for the command line, minus OS specific way of calling main TAK function:
 	%   ProcessRiverBasins /path/to/wdir Topo.mat river_mouths.txt Basins
 	%   ProcessRiverBasins /path/to/wdir Topo.mat river_mouths.shp Basins
	%   ProcessRiverBasins /path/to/wdir Topo.mat 500 Basins
 	%   ProcessRiverBasins /path/to/wdir Topo.mat river_mouths.txt Basins add_cat_grids AddCatGrids.mat add_grids AddGrids.mat
 	%   ProcessRiverBasins /path/to/wdir Topo.mat river_mouths.shp Basins calc_relief true relief_radii [1000 2500 5000]
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
	p.FunctionName = 'ProcessRiverBasins';
	addRequired(p,'wdir',@(x) ischar(x));
	addRequired(p,'MakeStreamsMat',@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.mat'))));
	addRequired(p,'river_mouths',@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.txt')))|| isnumeric(x) & isscalar(x) || ~isempty(regexp(x,regexptranslate('wildcard','*.shp'))));
	addRequired(p,'basin_dir',@(x) ischar(x));

	addParamValue(p,'ref_concavity',0.5,@(x) isscalar(x) && isnumeric(x));
	addParamValue(p,'threshold_area',1e6,@(x) isscalar(x) && isnumeric(x));
	addParamValue(p,'segment_length',1000,@(x) isscalar(x) && isnumeric(x));
	addParamValue(p,'write_arc_files',false,@(x) isscalar(x));
	addParamValue(p,'ksn_method','quick',@(x) ischar(validatestring(x,{'quick','trunk','trib'})));
	addParamValue(p,'min_order',4,@(x) isscalar(x) && isnumeric(x));
	addParamValue(p,'add_grids',[],@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.mat'))));
	addParamValue(p,'add_cat_grids',[],@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.mat'))));
	addParamValue(p,'resample_method','nearest',@(x) ischar(validatestring(x,{'nearest','bilinear','bicubic'})));
	addParamValue(p,'gradient_method','arcslope',@(x) ischar(validatestring(x,{'arcslope','gradient8'})));
	addParamValue(p,'calc_relief',false,@(x) isscalar(x));
	addParamValue(p,'relief_radii',[2500],@(x) isnumeric(x) && size(x,2)==1 || size(x,1)==1);
	addParamValue(p,'conditioned_DEM',[],@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.mat'))));
	addParamValue(p,'interp_value',0.1,@(x) isnumeric(x) && x>=0 && x<=1);
	addParamValue(p,'ksn_radius',5000,@(x) isnumeric(x) && isscalar(x));

	parse(p,wdir,MakeStreamsMat,river_mouths,basin_dir,varargin{:});
	wdir=p.Results.wdir;
	MATfile=p.Results.MakeStreamsMat;
	river_mouths=p.Results.river_mouths;
	basin_dir=p.Results.basin_dir;

	theta_ref=p.Results.ref_concavity;
	threshold_area=p.Results.threshold_area;
	segment_length=p.Results.segment_length;
	write_arc_files=p.Results.write_arc_files;
	ksn_method=p.Results.ksn_method;
	min_order=p.Results.min_order;
	AGmat=p.Results.add_grids;
	ACGmat=p.Results.add_cat_grids;
	resample_method=p.Results.resample_method;
	gradient_method=p.Results.gradient_method;
	calc_relief=p.Results.calc_relief;
	relief_radii=p.Results.relief_radii;
	iv=p.Results.interp_value;
	DEMhc=p.Results.conditioned_DEM;
	radius=p.Results.ksn_radius;



	% Set redo_flag
	redo_flag=false;

	% Load in required data from cmpMakeStreams
	MATfile=fullfile(wdir,MATfile);
	load(MATfile,'DEM','FD','S','A');

	% Load in conditioned DEM if provided
	if ~isempty(DEMhc)
		DEMhc=fullfile(wdir,DEMhc);
		load(DEMhc,'DEMc');
		if ~validatealignment(DEMc,DEM)
			error('Provided conditioned DEM does not appear to align with the original DEM')
		end
		DEMhc=DEMc;
	end

	% Navigate into dir
	basin_path=fullfile(wdir,basin_dir);
	if ~isdir(basin_path)
		mkdir(basin_path);
	end

	% Load Additional Grids if present and perform check on dimensions and cellsize of additional grids and resample if necessary
	if ~isempty(AGmat)
		AGmat=fullfile(wdir,AGmat);
		load(AGmat,'AG');
		num_grids=size(AG,1);
		for jj=1:num_grids
			AGoi=AG{jj,1};
			if ~validatealignment(AGoi,DEM);
				disp(['Resampling ' AG{jj,2} ' GRIDobj to be the same resolution and dimensions as the input DEM by the ' resample_method ' method']);
				AG{jj,1}=resample(AGoi,DEM,resample_method);
			end
		end
	else 
		AG=[];
	end

	% Peform check on segment length
	if (DEM.cellsize*3)>segment_length
		segment_length=DEM.cellsize*3;
		warning(['Provided segment_length is incompatible with DEM resolution, segment_length reset to ' num2str(segment_length)])
	end

	% Load Additional Categorical Grids if present
	if ~isempty(ACGmat);
		ACGmat=fullfile(wdir,ACGmat);
		load(ACGmat,'ACG');
	else
		ACG=[];
	end

	if ~isempty(regexp(river_mouths,regexptranslate('wildcard','*.shp')))
		disp('Reading shapefile and snapping river mouths to stream network')
		river_mouths=fullfile(wdir,river_mouths);
		rm_ms=shaperead(river_mouths);
		rm_t=struct2table(rm_ms);
		if ~strcmp(rm_t.Geometry(1),'Point')
			error('Shapefile provided as "river_mouths" does not appear to be a point shapefile');
		end
		fn=rm_t.Properties.VariableNames;
		xi=rm_t.X;
		yi=rm_t.Y;
		riv_nums=rm_t.(fn{4});

		num_basins=numel(xi);

		% Perform check if there are duplicate river_mouth IDs and junk them if so
		if numel(riv_nums)~=numel(unique(riv_nums))
			riv_nums=[1:numel(riv_nums)]';
			warning('Duplicate values present in "river_mouths" IDs, IDs have been reassigned');
			redo_flag=true;
		end	

		% Perform check that no river mouths are outside extent of DEM and remove
		[demx,demy]=getoutline(DEM,true);
		[demix]=inpolygon(xi,yi,demx,demy);	
		xi=xi(demix); yi=yi(demix); riv_nums=riv_nums(demix);

		% Snap to streams
		[xn,yn]=snap2stream(S,xi,yi);
		RM=[xn yn riv_nums];

		num_basins=numel(xn);

		if redo_flag
			csvwrite(fullfile(wdir,'river_mouths_updated.txt'),RM);
		end

	elseif ~isempty(regexp(river_mouths,regexptranslate('wildcard','*.txt')))
		river_mouths=fullfile(wdir,river_mouths);
		river_mouths=readtable(river_mouths);

		disp('Snapping river mouths to stream network')
		xi=river_mouths.(1);
		yi=river_mouths.(2);
		riv_nums=river_mouths.(3);

		num_basins=numel(xi);

		% Perform check if there are duplicate river_mouth IDs and junk them if so
		if numel(riv_nums)~=numel(unique(riv_nums))
			riv_nums=[1:numel(riv_nums)]';
			warning('Duplicate values present in "river_mouths" IDs, IDs have been reassigned');
			redo_flag=true;
		end	

		% Perform check that no river mouths are outside extent of DEM and remove
		[demx,demy]=getoutline(DEM,true);
		[demix]=inpolygon(xi,yi,demx,demy);	
		xi=xi(demix); yi=yi(demix); riv_nums=riv_nums(demix);

		% Snap to streams
		[xn,yn]=snap2stream(S,xi,yi);
		RM=[xn yn riv_nums];

		num_basins=numel(xn);

		if redo_flag
			csvwrite(fullfile(basin_path,'river_mouths_updated.txt'),RM);
		end

	elseif isscalar(river_mouths)
		disp('Generating river mouths based on provided elevation')
		sz=getnal(S,DEM);
		ix1=S.IXgrid;
		ix1(sz<river_mouths)=[];
		W=GRIDobj(DEM,'logical');
		W.Z(ix1)=true;
		Stemp=STREAMobj(FD,W);
		oxy=streampoi(Stemp,'outlets','xy');
		num_basins=size(oxy,1);
		olist=[1:num_basins]';
		RM=[oxy olist];
		csvwrite(fullfile(basin_path,'river_mouths.txt'),RM);
	end

	% Check for zeros and replace and warn
	if nnz(RM(:,3))~=numel(RM(:,3))
		warning('Zeros present in "river_mouths" IDs, IDs for zeros have been reassigned')
		zeroIDX=RM(:,3)==0;
		maxRM=max(RM(:,3));
		numZeros=sum(zeroIDX);
		zeroIX=find(zeroIDX);
		for ii=1:numZeros
			RM(zeroIX(ii),3)=maxRM+ii;
		end
		csvwrite(fullfile(basin_path,'river_mouths_updated.txt'),RM);
	end

	w1=waitbar(0,['Working on Basin Number 1 of ' num2str(num_basins) ' total basins']);
	for ii=1:num_basins
		xx=RM(ii,1);
		yy=RM(ii,2);
		basin_num=RM(ii,3);

		RiverMouth=[xx yy basin_num];

		% Build dependence map and clip out drainage basins
		I=dependencemap(FD,xx,yy);
		DEMoc=crop(DEM,I,nan);
		FDc=crop(FD,I);
		Ac=crop(A,I,nan);

		% Calculate drainage area
		dep_map=GRIDobj2mat(I);
		num_pix=sum(sum(dep_map));
		drainage_area=(num_pix*DEMoc.cellsize*DEMoc.cellsize)/(1e6);

		% Calculate hypsometry
		[rb,eb]=hypscurve(DEMoc,100);
		hyps=[rb eb];

		% Find weighted centroid of drainage basin
		[Cx,Cy]=FindCentroid(DEMoc);
		Centroid=[Cx Cy];

		% Generate new stream map
		Sc=STREAMobj(FDc,'minarea',threshold_area,'unit','mapunits');

		% Check to make sure the stream object isn't empty, this shouldn't occur anymore unless a bad pour point was provided...
		if isempty(Sc.x)
			warning(['Input threshold drainage area is too large for basin ' num2str(basin_num) ' decreasing threshold area for this basin']);
			new_thresh=threshold_area;
			while isempty(Sc.x)
				new_thresh=new_thresh/2;
				Sc=STREAMobj(FDc,'minarea',new_thresh,'unit','mapunits');
			end
		end

		% Calculate chi and create chi map
		Cc=chitransform(Sc,Ac,'a0',1,'mn',theta_ref);
		ChiOBJc=GRIDobj(DEMoc);
		ChiOBJc.Z(Sc.IXgrid)=Cc;

		% Calculate gradient
		switch gradient_method
		case 'gradient8'
			Goc=gradient8(DEMoc);
		case 'arcslope'
			Goc=arcslope(DEMoc);
		end

		% Hydrologically Condition DEM
		if isempty(DEMhc)
			zcon=mincosthydrocon(Sc,DEMoc,'interp',iv);
		else
			DEMhcc=crop(DEMhc,I,nan);
			zcon=getnal(Sc,DEMhcc);
		end
		DEMcc=GRIDobj(DEMoc);
		DEMcc.Z(DEMcc.Z==0)=NaN;
		DEMcc.Z(Sc.IXgrid)=zcon;

		% Find best fit concavity	
		SLc=klargestconncomps(Sc,1);	
		Chic=chiplot(SLc,DEMcc,Ac,'a0',1,'plot',false);

		% Calculate ksn
		switch ksn_method
		case 'quick'
			[MSc]=KSN_Quick(DEMoc,DEMcc,Ac,Sc,Chic.mn,segment_length);
			[MSNc]=KSN_Quick(DEMoc,DEMcc,Ac,Sc,theta_ref,segment_length);
		case 'trunk'
			[MSc]=KSN_Trunk(DEMoc,DEMcc,Ac,Sc,Chic.mn,segment_length,min_order);
			[MSNc]=KSN_Trunk(DEMoc,DEMcc,Ac,Sc,theta_ref,segment_length,min_order);
		case 'trib'
			% Overide choice if very small basin as KSN_Trib will fail for small basins
			if drainage_area>2.5
				[MSc]=KSN_Trib(DEMoc,DEMcc,FDc,Ac,Sc,Chic.mn,segment_length);
				[MSNc]=KSN_Trib(DEMoc,DEMcc,FDc,Ac,Sc,theta_ref,segment_length);
			else
				[MSc]=KSN_Quick(DEMoc,DEMcc,Ac,Sc,Chic.mn,segment_length);
				[MSNc]=KSN_Quick(DEMoc,DEMcc,Ac,Sc,theta_ref,segment_length);
			end
		end

		% Calculate basin wide ksn statistics
		min_ksn=nanmin([MSNc.ksn]);
		mean_ksn=nanmean([MSNc.ksn]);
		max_ksn=nanmax([MSNc.ksn]);
		std_ksn=nanstd([MSNc.ksn]);
		se_ksn=std_ksn/sqrt(numel(MSNc)); % Standard error

		% Calculate basin wide gradient statistics
		min_grad=nanmin(Goc.Z(:));
		mean_grad=nanmean(Goc.Z(:));
		max_grad=nanmax(Goc.Z(:));
		std_grad=nanstd(Goc.Z(:));
		se_grad=std_grad/sqrt(sum(~isnan(Goc.Z(:)))); % Standard error

		% Calculate basin wide elevation statistics
		min_z=nanmin(DEMoc.Z(:));
		mean_z=nanmean(DEMoc.Z(:));
		max_z=nanmax(DEMoc.Z(:));
		std_z=nanstd(DEMoc.Z(:));
		se_z=std_z/sqrt(sum(~isnan(DEMoc.Z(:)))); % Standard error

		KSNc_stats=[mean_ksn se_ksn std_ksn min_ksn max_ksn];
		Gc_stats=double([mean_grad se_grad std_grad min_grad max_grad]);
		Zc_stats=double([mean_z se_z std_z min_z max_z]);

		% Find outlet elevation
		out_ix=coord2ind(DEMoc,xx,yy);
		out_el=double(DEMoc.Z(out_ix));

		% Save base file
		FileName=fullfile(basin_path,['Basin_' num2str(basin_num) '_Data.mat']);
		save(FileName,'RiverMouth','DEMcc','DEMoc','out_el','drainage_area','hyps','FDc','Ac','Sc','SLc','Chic','Goc','MSc','MSNc','KSNc_stats','Gc_stats','Zc_stats','Centroid','ChiOBJc','ksn_method','gradient_method','theta_ref','-v7.3');

		if strcmp(ksn_method,'trunk')
			save(FileName,'min_order','-append');
		end
		
		% Make interpolated ksn grid
		try 
			[KsnOBJc] = KsnAvg(DEMoc,MSNc,radius);
			save(FileName,'KsnOBJc','radius','-append');
		catch
			warning(['Interpolation of KSN grid failed for basin ' num2str(RiverMouth(:,3))]);
			save(FileName,'radius','-append');
		end

		% If additional grids are present, append them to the mat file
		if ~isempty(AG)
			num_grids=size(AG,1);
			AGc=cell(size(AG));
			for jj=1:num_grids
				AGcOI=crop(AG{jj,1},I,nan);
				AGc{jj,1}=AGcOI;
				AGc{jj,2}=AG{jj,2};
				mean_AGc=nanmean(AGcOI.Z(:));
				min_AGc=nanmin(AGcOI.Z(:));
				max_AGc=nanmax(AGcOI.Z(:));
				std_AGc=nanstd(AGcOI.Z(:));
				se_AGc=std_AGc/sqrt(sum(~isnan(AGcOI.Z(:))));
				AGc_stats(jj,:)=[mean_AGc se_AGc std_AGc min_AGc max_AGc];
			end
			save(FileName,'AGc','AGc_stats','-append');				
		end

		if ~isempty(ACG)
			num_grids=size(ACG,1);
			ACGc=cell(size(ACG));
			for jj=1:num_grids
				ACGcOI=crop(ACG{jj,1},I,nan);
				ACGc{jj,1}=ACGcOI;
				ACGc{jj,3}=ACG{jj,3};
				edg=ACG{jj,2}.Numbers;
				edg=edg+0.5;
				edg=vertcat(0.5,edg);
				[N,~]=histcounts(ACGcOI.Z(:),edg);
				T=ACG{jj,2};
				T.Counts=N';
				ACGc{jj,2}=T;
				ACGc_stats(jj,1)=[mode(ACGcOI.Z(:))];
			end
			save(FileName,'ACGc','ACGc_stats','-append');	
		end				

		if calc_relief
			num_rlf=numel(relief_radii);
			rlf=cell(num_rlf,2);
			rlf_stats=zeros(num_rlf,6);
			for jj=1:num_rlf
				% Calculate relief
				radOI=relief_radii(jj);
				rlf{jj,2}=radOI;
				rlfOI=localtopography(DEMoc,radOI);
				rlf{jj,1}=rlfOI;
				% Calculate stats
				mean_rlf=nanmean(rlfOI.Z(:));
				min_rlf=nanmin(rlfOI.Z(:));
				max_rlf=nanmax(rlfOI.Z(:));
				std_rlf=nanstd(rlfOI.Z(:));
				se_rlf=std_rlf/sqrt(sum(~isnan(rlfOI.Z(:))));
				rlf_stats(jj,:)=[mean_rlf se_rlf std_rlf min_rlf max_rlf radOI];
			end
			save(FileName,'rlf','rlf_stats','-append');
		end

		if write_arc_files
			% Replace NaNs in DEM with -32768
			Didx=isnan(DEMoc.Z);
			DEMoc_temp=DEMoc;
			DEMoc_temp.Z(Didx)=-32768;

			DEMFileName=fullfile(basin_path,['Basin_' num2str(basin_num) '_DEM.txt']);
			GRIDobj2ascii(DEMoc_temp,DEMFileName);
			CHIFileName=fullfile(basin_path,['Basin_' num2str(basin_num) '_CHI.txt']);
			GRIDobj2ascii(ChiOBJc,CHIFileName);
			KSNFileName=fullfile(basin_path,['Basin_' num2str(basin_num) '_KSN.shp']);
			shapewrite(MSNc,KSNFileName);

			if calc_relief
				for jj=1:num_rlf
					RLFFileName=fullfile(basin_path,['Basin_' num2str(basin_num) '_RLF_' num2str(rlf{jj,2}) '.txt']);
					GRIDobj2ascii(rlf{jj,1},RLFFileName);
				end
			end

			if ~isempty(AG);
				for jj=1:num_grids
					AGcFileName=fullfile(basin_path,['Basin_' num2str(basin_num) '_' AGc{jj,2} '.txt']);
					GRIDobj2ascii(AGc{jj,1},AGcFileName);
				end
			end

			if ~isempty(ACG);
				for jj=1:num_grids
					ACGcFileName=fullfile(basin_path,['Basin_' num2str(basin_num) '_' ACGc{jj,3} '.txt']);
					GRIDobj2ascii(ACGc{jj,1},ACGcFileName);
				end
			end
		end

		waitbar(ii/num_basins,w1,['Completed ' num2str(ii) ' of ' num2str(num_basins) ' total basins'])
	end

	close(w1)

end


function [ksn_ms]=KSN_Quick(DEM,DEMc,A,S,theta_ref,segment_length)
	g=gradient(S,DEMc);
	G=GRIDobj(DEM);
	G.Z(S.IXgrid)=g;

	Z_RES=DEMc-DEM;

	ksn=G./(A.*(A.cellsize^2)).^(-theta_ref);

	SD=GRIDobj(DEM);
	SD.Z(S.IXgrid)=S.distance;
	
	ksn_ms=STREAMobj2mapstruct(S,'seglength',segment_length,'attributes',...
		{'ksn' ksn @mean 'uparea' (A.*(A.cellsize^2)) @mean 'gradient' G @mean 'cut_fill' Z_RES @mean...
		'min_dist' SD @min 'max_dist' SD @max});

	seg_dist=[ksn_ms.max_dist]-[ksn_ms.min_dist];
	distcell=num2cell(seg_dist');
	[ksn_ms(1:end).seg_dist]=distcell{:};
	ksn_ms=rmfield(ksn_ms,{'min_dist','max_dist'});
end

function [ksn_ms]=KSN_Trunk(DEM,DEMc,A,S,theta_ref,segment_length,min_order)

	order_exp=['>=' num2str(min_order)];

    Smax=modify(S,'streamorder',order_exp);
	Smin=modify(S,'rmnodes',Smax);

	g=gradient(S,DEMc);
	G=GRIDobj(DEM);
	G.Z(S.IXgrid)=g;

	Z_RES=DEMc-DEM;

	ksn=G./(A.*(A.cellsize^2)).^(-theta_ref);

	SDmax=GRIDobj(DEM);
	SDmin=GRIDobj(DEM);
	SDmax.Z(Smax.IXgrid)=Smax.distance;
	SDmin.Z(Smin.IXgrid)=Smin.distance;

	ksn_ms_min=STREAMobj2mapstruct(Smin,'seglength',segment_length,'attributes',...
		{'ksn' ksn @mean 'uparea' (A.*(A.cellsize^2)) @mean 'gradient' G @mean 'cut_fill' Z_RES @mean...
		'min_dist' SDmin @min 'max_dist' SDmin @max});

	ksn_ms_max=STREAMobj2mapstruct(Smax,'seglength',segment_length,'attributes',...
		{'ksn' ksn @mean 'uparea' (A.*(A.cellsize^2)) @mean 'gradient' G @mean 'cut_fill' Z_RES @mean...
		'min_dist' SDmax @min 'max_dist' SDmax @max});

	ksn_ms=vertcat(ksn_ms_min,ksn_ms_max);
	seg_dist=[ksn_ms.max_dist]-[ksn_ms.min_dist];
	distcell=num2cell(seg_dist');
	[ksn_ms(1:end).seg_dist]=distcell{:};
	ksn_ms=rmfield(ksn_ms,{'min_dist','max_dist'});
end

function [ksn_ms]=KSN_Trib(DEM,DEMc,FD,A,S,theta_ref,segment_length)

	% Define non-intersecting segments
	[as]=networksegment_slim(DEM,FD,S);
	seg_bnd_ix=as.ix;
	% Precompute values or extract values needed for later
	z=getnal(S,DEMc);
	zu=getnal(S,DEM);
	z_res=z-zu;
	g=gradient(S,DEMc);
	c=chitransform(S,A,'a0',1,'mn',theta_ref);
	d=S.distance;
	da=getnal(S,A.*(A.cellsize^2));
	ixgrid=S.IXgrid;
	% Extract ordered list of stream indices and find breaks between streams
	s_node_list=S.orderednanlist;
	streams_ix=find(isnan(s_node_list));
	streams_ix=vertcat(1,streams_ix);
	% Generate empty node attribute list for ksn values
	ksn_nal=zeros(size(d));
	% Begin main loop through channels
	num_streams=numel(streams_ix)-1;
	seg_count=1;
	for ii=1:num_streams
		% Extract node list for stream of interest
		if ii==1
			snlOI=s_node_list(streams_ix(ii):streams_ix(ii+1)-1);
		else
			snlOI=s_node_list(streams_ix(ii)+1:streams_ix(ii+1)-1);
		end

		% Determine which segments are within this stream
		[~,~,dn]=intersect(snlOI,seg_bnd_ix(:,1));
		[~,~,up]=intersect(snlOI,seg_bnd_ix(:,2));
		seg_ix=intersect(up,dn);

		num_segs=numel(seg_ix);
		dn_up=seg_bnd_ix(seg_ix,:);
		for jj=1:num_segs
			% Find positions within node list
			dnix=find(snlOI==dn_up(jj,1));
			upix=find(snlOI==dn_up(jj,2));
			% Extract segment indices of desired segment
			seg_ix_oi=snlOI(upix:dnix);
			% Extract flow distances and normalize
			dOI=d(seg_ix_oi);
			dnOI=dOI-min(dOI);
			num_bins=ceil(max(dnOI)/segment_length);
			bin_edges=[0:segment_length:num_bins*segment_length];
			% Loop through bins
			for kk=1:num_bins
				idx=dnOI>bin_edges(kk) & dnOI<=bin_edges(kk+1);
				bin_ix=seg_ix_oi(idx);
				cOI=c(bin_ix);
				zOI=z(bin_ix);
					if numel(cOI)>2
						[ksn_val,r2]=Chi_Z_Spline(cOI,zOI);
						ksn_nal(bin_ix)=ksn_val;

						% Build mapstructure
						ksn_ms(seg_count).Geometry='Line';
						ksm_ms(seg_count).BoundingBox=[min(S.x(bin_ix)),min(S.y(bin_ix));max(S.x(bin_ix)),max(S.y(bin_ix))];
						ksn_ms(seg_count).X=S.x(bin_ix);
						ksn_ms(seg_count).Y=S.y(bin_ix);
						ksn_ms(seg_count).ksn=ksn_val;
						ksn_ms(seg_count).uparea=mean(da(bin_ix));
						ksn_ms(seg_count).gradient=mean(g(bin_ix));
						ksn_ms(seg_count).cut_fill=mean(z_res(bin_ix));
						ksn_ms(seg_count).seg_dist=max(S.distance(bin_ix))-min(S.distance(bin_ix));
						ksn_ms(seg_count).chi_r2=r2;
						
						seg_count=seg_count+1;
					end
			end
		end
	end
end

function seg = networksegment_slim(DEM,FD,S)
	% Slimmed down version of 'networksegment' from main TopoToolbox library that also removes zero and single node length segments

	%% Identify channel heads, confluences, b-confluences and outlets
	Vhead = streampoi(S,'channelheads','logical');  ihead=find(Vhead==1);  IXhead=S.IXgrid(ihead);
	Vconf = streampoi(S,'confluences','logical');   iconf=find(Vconf==1);  IXconf=S.IXgrid(iconf);
	Vout = streampoi(S,'outlets','logical');        iout=find(Vout==1);    IXout=S.IXgrid(iout);
	Vbconf = streampoi(S,'bconfluences','logical'); ibconf=find(Vbconf==1);IXbconf=S.IXgrid(ibconf);

	%% Identify basins associated to b-confluences and outlets
	DB   = drainagebasins(FD,vertcat(IXbconf,IXout));DBhead=DB.Z(IXhead); DBbconf=DB.Z(IXbconf); DBconf=DB.Z(IXconf); DBout=DB.Z(IXout);

	%% Compute flowdistance
	D = flowdistance(FD);

	%% Identify river segments
	% links between channel heads and b-confluences
	[~,ind11,ind12]=intersect(DBbconf,DBhead);
	% links between confluences and b-confluences
	[~,ind21,ind22]=intersect(DBbconf,DBconf);
	% links between channel heads and outlets
	[~,ind31,ind32]=intersect(DBout,DBhead);
	% links between channel heads and outlets
	[~,ind41,ind42]=intersect(DBout,DBconf);
	% Connecting links into segments
	IX(:,1) = [ IXbconf(ind11)' IXbconf(ind21)' IXout(ind31)'  IXout(ind41)'  ];   ix(:,1)= [ ibconf(ind11)' ibconf(ind21)' iout(ind31)'  iout(ind41)'  ];
	IX(:,2) = [ IXhead(ind12)'  IXconf(ind22)'  IXhead(ind32)' IXconf(ind42)' ];   ix(:,2)= [ ihead(ind12)'  iconf(ind22)'  ihead(ind32)' iconf(ind42)' ];

	% Compute segment flow length
	flength=double(abs(D.Z(IX(:,1))-D.Z(IX(:,2))));

	% Remove zero and one node length elements
	idx=flength>=2*DEM.cellsize;
	seg.IX=IX(idx,:);
	seg.ix=ix(idx,:);
	seg.flength=flength(idx);

	% Number of segments
	seg.n=numel(IX(:,1));
end

function [KSN,R2] = Chi_Z_Spline(c,z)

	% Resample chi-elevation relationship using cubic spline interpolation
	[~,minIX]=min(c);
	zb=z(minIX);
	chiF=c-min(c);
	zabsF=z-min(z);
	chiS=linspace(0,max(chiF),numel(chiF)).';
	zS=spline(chiF,zabsF,chiS);

	% Calculate ksn via slope
	KSN= chiS\(zS); % mn not needed because a0 is fixed to 1

	% Calculate R^2
	z_pred=chiF.*KSN;
	sstot=sum((zabsF-mean(zabsF)).^2);
	ssres=sum((zabsF-z_pred).^2);
	R2=1-(ssres/sstot);

end

function [KSNGrid] = KsnAvg(DEM,ksn_ms,radius)

	% Calculate radius
	radiuspx = ceil(radius/DEM.cellsize);

	% Record mask of current NaNs
	MASK=isnan(DEM.Z);

	% Make grid with values along channels
	KSNGrid=GRIDobj(DEM);
	KSNGrid.Z(:,:)=NaN;
	for ii=1:numel(ksn_ms)
		ix=coord2ind(DEM,ksn_ms(ii).X,ksn_ms(ii).Y);
		KSNGrid.Z(ix)=ksn_ms(ii).ksn;
	end

	% Local mean based on radius
	ISNAN=isnan(KSNGrid.Z);
    [~,L] = bwdist(~ISNAN,'e');
    ksng = KSNGrid.Z(L);           
    FLT   = fspecial('disk',radiuspx);
    ksng   = imfilter(ksng,FLT,'symmetric','same','conv');

    % Set original NaN cells back to NaN
    ksng(MASK)=NaN;

    % Output
    KSNGrid.Z=ksng;
end
