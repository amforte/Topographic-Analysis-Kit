function [MS]=Basin2Shape(DEM,location_of_data_files,varargin)
	%
	% Usage:
	%	[MapStruct]=Basin2Shape(DEM,location_of_data_files);
	%	[MapStruct]=Basin2Shape(DEM,location_of_data_files,'name',value,...);	
	%
	% Description:
	% 	Function to take the outputs from 'ProcessRiverBasins' and 'SubDivideBigBasins' and produce a single shapefile showing the outlines of polygons
	% 	and with commonly desired attributes from the results of 'ProcessRiverBasins' etc. See below for a full list of fields that the output shapefile
	% 	will include. If additional grids were provided to 'ProcessRiverBasins', mean and standard error values for those grids will be auto-populated in
	% 	the shapefile and the name of the fields will be the character array provided in the second column of additional grids input. This function also
	% 	allows you to input a list of additional fields you wish to include (see Optional Inputs below). If you would rather create a GRIDobj with specified
	% 	values, use 'Basin2Raster'.
	%
	% Required Inputs:
	%		DEM - GRIDobj of the DEM originally used as input for 'ProcessRiverBasins' for the basins of interest.
	% 		location_of_data_files - full path of folder which contains the mat files from 'ProcessRiverBasins'
	%
	% Optional Inputs:
	%		location_of_subbasins ['SubBasins'] - name of folder that contains subbasins of interest (if you created subbasins using
	%			"SubDivideBigBasins"), expected to be within the main Basin folder provided with "location_of_data_files". Note that if you do not provide
	%			the correct directory name for the location of the subbasins, subbasin values will not be included in the output regardless of your choice
	%			for the "include" parameter.
	%		shape_name ['basins'] - name for the shapefile to be export, must have no spaces to be a valid name for ArcGIS and should NOT include the '.shp';
	%		include ['all'] - parameter to specify which basins to include in building the shapfile. The default 'all' will include all basin mat files in the 
	%			folder you specify. Providing 'subdivided' will check to see if a given main basin was subdivided using 'SubdivideBigBasins' and then only include 
	%			the subdivided versions of that basin (i.e. the original main basin for those subbasins will not be included in the shapefile). Providing 'bigonly'
	%			will only include the original basins produced by 'ProcessRiverBasins' even if 'SubDivideBigBasins' was run. If 'SubDivideBigBasins' was never run,
	%			result of 'all' and 'bigonly' will be the same.
	%		extra_field_values [] - cell array of extra field values you wish to include. The first column in this cell array must be the river basin number
	%			(i.e. the identifying number in the third column of the RiverMouth input to ProcessRiverBasins or the number generated for the basin in
	%			SubDivideBigBasins). Only one row per river basin number is allowed and ALL river basin numbers in the basins being processed must have a value
	%			associated with them. Additional columns are interpreted as the values with which you wish to populate the extra fields. These can either be character
	%			arrays or numbers, other values will results in an error. 
	%		extra_field_names [] - a 1 x m cell array of field names, as characters (no spaces as this won't work with shapefile attributes), associated with the field values. 
	%			These must be in the same order as values are given in extra_field_values. If for example your extra_field_values cell array is 3 columns with the river number, 
	%			sample name, and erosion rate then your extra_field_names cell array should include entries for 'sample_name' and 'erosion_rate' in that order.
	%		new_concavity [] - a 1 x m array of concavity values to recalculate normalized channel steepness statistics (mean, standard error and/or standard deviation) using the
	%			provided concavities. The default method for this is very quick, but an approximation. If you are willing to wait and you want the ksn statistics at the new 
	%			concavity to be exact, change 'new_ksn_method' to 'exact'.
	%		new_ksn_method ['approximate'] - parameter to control how a new concavity is calculated if an entry is provided to 'new_concavity', options are 'approximate' (the default)
	%			and 'exact'. Setting this to exact will slow the calculation time considerably.
	%		segment_length [1000] - smoothing distance for ksn if new_concavities are provided and new_ksn_method is set to 'exact', otherwise ignored.	
	%		uncertainty ['se'] - parameter to control which measure of uncertainty is included, expects 'se' for standard error (default), 'std' for standard deviation, or 'both'
	%			to include both standard error and deviation.
	%		populate_categories [false] - logical flag to add entries that indicate the percentage of a watershed occupied by each category from a categorical grid, e.g. if you
	%			provided an entry for 'add_cat_grids' to ProcessRiverBasins that was a geologic map that had three units, 'Q', 'Mz', and 'Pz' and you set 'populate_categories' 
	%			to true there will be field names in the resulting shapefile named 'Q', 'Mz', and 'Pz' and the values stored in those columns will correspond to the percentage 
	%			of each basin covered by each unit for each basin. Setting populate_categories to true will not have any effect if no entry was provided to 'add_cat_grids' when
	%			running ProcessRiverBasins.
	%
	% Output:
	%		Outputs a mapstructure (MS) and saves a shapefile with the following default fields:
	%			river_mouth - river mouth number provided to ProcessRiverBasins
	%			drainage_area - drainage area of basin in km^2
	%			center_x - x coordinate of basin in projected coordinates
	%			center_y - y coordinate of basin in projected coordinates
	%			outlet_elevation - elevation of pour point in m
	%			mean_el - mean elevation of basin in meters
	%			max_el - maximum elevation of basin in meters
	%			mean_ksn - mean channel steepenss
	%			mean_gradient - mean gradient
	%		Either standard errors, standard deviations or both will be populated for elevation, ksn, and gradient depending on value of 'uncertainty'
	%		Mean and standard error / standard deviation / both values will be populated for any additional grids
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 06/18/18 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	

	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'Basin2Shape';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'location_of_data_files',@(x) isdir(x));

	addParameter(p,'location_of_subbasins','SubBasins',@(x) ischar(x) || isempty(x));
	addParameter(p,'shape_name','basins',@(x) ischar(x));
	addParameter(p,'include','all',@(x) ischar(validatestring(x,{'all','subdivided','bigonly'})));
	addParameter(p,'extra_field_values',[],@(x) isa(x,'cell'));
	addParameter(p,'extra_field_names',[],@(x) isa(x,'cell') & size(x,1)==1);
	addParameter(p,'new_concavity',[],@(x) isnumeric(x));
	addParameter(p,'new_ksn_method','approximate',@(x) ischar(validatestring(x,{'approximate','exact'})));
	addParameter(p,'segment_length',1000,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'uncertainty','se',@(x) ischar(validatestring(x,{'se','std','both'})));
	addParameter(p,'populate_categories',false,@(x) isscalar(x) && islogical(x))
	addParameter(p,'suppress_shape_write',false,@(x) isscalar(x) && islogical(x))

	parse(p,DEM,location_of_data_files,varargin{:});
	DEM=p.Results.DEM;
	location_of_data_files=p.Results.location_of_data_files;

	location_of_subbasins=p.Results.location_of_subbasins;
	shape_name=p.Results.shape_name;
	include=p.Results.include;
	efv=p.Results.extra_field_values;
	efn=p.Results.extra_field_names;
	new_concavity=p.Results.new_concavity;
	new_ksn_method=p.Results.new_ksn_method;
	segment_length=p.Results.segment_length;
	uncertainty=p.Results.uncertainty;
	pc=p.Results.populate_categories;
	ssw=p.Results.suppress_shape_write;

	% Deal with variability in format of locations
	[sub_head,~,~]=fileparts(location_of_subbasins);
	if isempty(sub_head)
		location_of_subbasins=[location_of_data_files filesep location_of_subbasins];
	end

	% Switch for which basins to include
	switch include
	case 'all'
		FileList1=dir([location_of_data_files filesep '*_Data.mat']);
		FileList2=dir([location_of_subbasins filesep '*_DataSubset*.mat']);
		FileList=vertcat(FileList1,FileList2);
		num_files=numel(FileList);
	case 'bigonly'
		FileList=dir([location_of_data_files filesep '*_Data.mat']);
		num_files=numel(FileList);
	case 'subdivided'
		AllFullFiles=dir([location_of_data_files filesep '*_Data.mat']);
		num_basins=numel(AllFullFiles);
		basin_nums=zeros(num_basins,1);
		for jj=1:num_basins
			fileName=AllFullFiles(jj,1).name;
			basin_nums(jj)=sscanf(fileName,'%*6s %i'); %%%
		end

		FileCell=cell(num_basins,1);
		for kk=1:num_basins
			basin_num=basin_nums(kk);
			SearchAllString=[location_of_data_files filesep '*_' num2str(basin_num) '_Data.mat'];
			SearchSubString=[location_of_subbasins filesep '*_' num2str(basin_num) '_DataSubset*.mat'];

			if numel(dir(SearchSubString))>0
				Files=dir(SearchSubString);
			else
				Files=dir(SearchAllString);
			end

			FileCell{kk}=Files;
		end
		FileList=vertcat(FileCell{:});
		num_files=numel(FileList);
	end

	% Sort by drainage area to ensure proper drawing order
	dalist=zeros(num_files,1);
	for ii=1:num_files
		FileName=[FileList(ii,1).folder filesep FileList(ii,1).name];
		load(FileName,'drainage_area');
		dalist(ii,1)=drainage_area;
	end
	[~,six]=sort(dalist,'descend');
	FileList=FileList(six);

	% Initiate Map Structure
	MS=struct;

	% Begin main loop
	w1=waitbar(0,'Building polygons');
	for ii=1:num_files;
		FileName=[FileList(ii,1).folder filesep FileList(ii,1).name];
		DB=GRIDobj(DEM);

		load(FileName,'DEMoc','RiverMouth','drainage_area','out_el','KSNc_stats','Zc_stats','Gc_stats','Centroid','hyps','Chic','DEMcc','Sc','Ac','theta_ref');

		I=~isnan(DEMoc.Z);
		[X,Y]=getcoordinates(DEMoc);
		xmat=repmat(X,numel(Y),1);
		ymat=repmat(Y,1,numel(X));

		xix=xmat(I);
		yix=ymat(I);

		ix=coord2ind(DEM,xix,yix);

		DB.Z(ix)=RiverMouth(:,3);

		[ms_temp,~,~]=GRIDobj2polygon(DB);

		% Populate default fields in output map structure
		MS(ii,1).Geometry='Polygon';
		MS(ii,1).X=ms_temp.X;
		MS(ii,1).Y=ms_temp.Y;
		MS(ii,1).ID=ii;
		MS(ii,1).river_mouth=ms_temp.gridval;
		MS(ii,1).center_x=Centroid(1);
		MS(ii,1).center_y=Centroid(2);
		MS(ii,1).drainage_area=drainage_area;
		MS(ii,1).outlet_elevation=out_el;
		MS(ii,1).mean_el=Zc_stats(1);
		MS(ii,1).max_el=Zc_stats(5);
		MS(ii,1).mean_ksn=KSNc_stats(1);
		MS(ii,1).mean_gradient=Gc_stats(1);

		switch uncertainty
		case 'se'
			MS(ii,1).se_el=Zc_stats(2);
			MS(ii,1).se_ksn=KSNc_stats(2);
			MS(ii,1).se_gradient=Gc_stats(2);
		case 'std'
			MS(ii,1).std_el=Zc_stats(3);
			MS(ii,1).std_ksn=KSNc_stats(3);
			MS(ii,1).std_gradient=Gc_stats(3);
		case 'both'
			MS(ii,1).se_el=Zc_stats(2);
			MS(ii,1).se_ksn=KSNc_stats(2);
			MS(ii,1).se_gradient=Gc_stats(2);
			MS(ii,1).std_el=Zc_stats(3);
			MS(ii,1).std_ksn=KSNc_stats(3);
			MS(ii,1).std_gradient=Gc_stats(3);
		end

		if ~isempty(new_concavity)
			load(FileName,'MSNc');
			for jj=1:numel(new_concavity)
				switch new_ksn_method
				case 'approximate'
					[mean_ksn,std_ksn,se_ksn]=ksn_convert_approx(MSNc,new_concavity(jj));
				case 'exact'
					[mean_ksn,std_ksn,se_ksn]=ksn_convert_exact(FileName,segment_length,new_concavity(jj));
				end
				ksn_cat_name=matlab.lang.makeValidName(['mean_ksn_' num2str(new_concavity(jj))]);
				MS(ii,1).(ksn_cat_name)=mean_ksn;
				switch uncertainty
				case 'se'
					ksn_cat_name_se=matlab.lang.makeValidName(['se_ksn_' num2str(new_concavity(jj))]);
					MS(ii,1).(ksn_cat_name_se)=se_ksn;
				case 'std'
					ksn_cat_name_std=matlab.lang.makeValidName(['std_ksn_' num2str(new_concavity(jj))]);
					MS(ii,1).(ksn_cat_name_std)=std_ksn;
				case 'both'
					ksn_cat_name_se=matlab.lang.makeValidName(['se_ksn_' num2str(new_concavity(jj))]);
					MS(ii,1).(ksn_cat_name_se)=se_ksn;
					ksn_cat_name_std=matlab.lang.makeValidName(['std_ksn_' num2str(new_concavity(jj))]);
					MS(ii,1).(ksn_cat_name_std)=std_ksn;					
				end
			end
		end		

		MS(ii,1).hyp_int=double(abs(trapz((hyps(:,2)-min(hyps(:,2)))/(max(hyps(:,2))-min(hyps(:,2))),hyps(:,1)/100)));
		MS(ii,1).theta=Chic.mn;

		c=chiplot(Sc,DEMcc,Ac,'a0',1,'mn',theta_ref,'plot',false);
		MS(ii,1).chi_r2=c.R2;
		
		% Determine if a georef structure exists and if so, produce lat-lon locations for sample points
		if ~isempty(DEM.georef)
			try
				% Check how projection is stored (control for older TopoToolbox methods)
				if isfield(DEM.georef,'mstruct')
					proj=DEM.georef.mstruct;
				else
					proj=DEM.georef;
				end
				% Project
				[s_lat,s_lon]=projinv(proj,RiverMouth(:,1),RiverMouth(:,2));
				% Populate lat lon coordinates
				MS(ii,1).outlet_lat=s_lat;
				MS(ii,1).outlet_lon=s_lon;
			catch
				str=sprintf('Projection is either not supported or you do not have the Mapping Toolbox,\n unable to convert river mouth coordinates to lat-lon');
				disp(str);
			end
		else
			str=sprintf('GRIDobj does not have projection information, unable to convert river mouth coordinates to lat-lon');
			disp(str);
		end

		% Check for additional grids within the process river basins output
		VarList=whos('-file',FileName);

		VarInd=find(strcmp(cellstr(char(VarList.name)),'KSNQc_stats'));
		if ~isempty(VarInd)
			load(FileName,'KSNQc_stats');

			MS(ii,1).mean_ksn_q=KSNQc_stats(:,1);

			switch uncertainty
			case 'se'
				MS(ii,1).se_ksn_q=KSNQc_stats(:,2);
			case 'std'
				MS(ii,1).std_ksn_q=KSNQc_stats(:,3);
			case 'both'
				MS(ii,1).se_ksn_q=KSNQc_stats(:,2);
				MS(ii,1).std_ksn_q=KSNQc_stats(:,3);
			end
		end


		VarInd=find(strcmp(cellstr(char(VarList.name)),'AGc'));
		if ~isempty(VarInd)
			load(FileName,'AGc','AGc_stats');
			num_grids=size(AGc,1);

			for kk=1:num_grids
				mean_prop_name=['mean_' AGc{kk,2}];		
				MS(ii,1).(mean_prop_name)=double(AGc_stats(kk,1));

				switch uncertainty
				case 'se'
					se_prop_name=['se_' AGc{kk,2}];
					MS(ii,1).(se_prop_name)=double(AGc_stats(kk,2));
				case 'std'
					std_prop_name=['std_' AGc{kk,2}];
					MS(ii,1).(std_prop_name)=double(AGc_stats(kk,3));
				case 'both'
					se_prop_name=['se_' AGc{kk,2}];
					MS(ii,1).(se_prop_name)=double(AGc_stats(kk,2));
					std_prop_name=['std_' AGc{kk,2}];
					MS(ii,1).(std_prop_name)=double(AGc_stats(kk,3));
				end
			end
		end	
	

		VarInd=find(strcmp(cellstr(char(VarList.name)),'rlf'));
		if ~isempty(VarInd)
			load(FileName,'rlf','rlf_stats');
			num_grids=size(rlf,1);

			for kk=1:num_grids
				mean_prop_name=['mean_rlf' num2str(rlf{kk,2})];
				se_prop_name=['se_rlf' num2str(rlf{kk,2})];
				MS(ii,1).(mean_prop_name)=double(rlf_stats(kk,1));
				MS(ii,1).(se_prop_name)=double(rlf_stats(kk,2));

				switch uncertainty
				case 'se'
					se_prop_name=['se_rlf' num2str(rlf{kk,2})];
					MS(ii,1).(se_prop_name)=double(rlf_stats(kk,2));
				case 'std'
					std_prop_name=['std_rlf' num2str(rlf{kk,2})];
					MS(ii,1).(std_prop_name)=double(rlf_stats(kk,3));
				case 'both'
					se_prop_name=['se_rlf' num2str(rlf{kk,2})];
					MS(ii,1).(se_prop_name)=double(rlf_stats(kk,2));
					std_prop_name=['std_rlf' num2str(rlf{kk,2})];
					MS(ii,1).(std_prop_name)=double(rlf_stats(kk,3));
				end
			end
		end		

		% Check for the presence of extra fields provided at input
		if ~isempty(efv)
			bnl=cell2mat(efv(:,1));

			ix=find(bnl==RiverMouth(:,3));
			% Check to make sure a single entry exists for each basin number
			if ~isempty(ix) & numel(ix)==1
				efvOI=efv(ix,2:end); % Strip out the basin number column
				num_efv=size(efvOI,2);

				for kk=1:num_efv
					field_name=efn{kk};
					field_value=efvOI{kk};
					% Check to see if field value is a number or string
					if ischar(field_value)
						MS(ii,1).(field_name)=field_value;
					elseif isnumeric(field_value)
						MS(ii,1).(field_name)=double(field_value);
					else
						if isdeployed
							errordlg(['Extra field value provided for ' field_name ' is neither numeric or a character'])
						end
						error(['Extra field value provided for ' field_name ' is neither numeric or a character']);
					end
				end
			elseif numel(ix)>1
				if isdeployed
					errordlg(['More than one entry was provided for extra fields for basin ' num2str(RiverMouth(:,3))])
				end
				error(['More than one entry was provided for extra fields for basin ' num2str(RiverMouth(:,3))]);
			elseif isempty(ix)
				if isdeployed
					errordlg(['No one entry was provided for extra field values for basin ' num2str(RiverMouth(:,3))])
				end
				error(['No one entry was provided for extra field values for basin ' num2str(RiverMouth(:,3))]);
			end
		end

		VarInd=find(strcmp(cellstr(char(VarList.name)),'ACGc'));
		if ~isempty(VarInd)
			load(FileName,'ACGc','ACGc_stats');
			num_grids=size(ACGc,1);

			for kk=1:num_grids
				mode_prop_name=['mode_' ACGc{kk,3}];	
				perc_prop_name=['mode_' ACGc{kk,3} '_pct'];
				ix=find(ACGc{kk,2}.Numbers==ACGc_stats(kk,1),1);	
				MS(ii,1).(mode_prop_name)=ACGc{kk,2}.Categories{ix};
				total_nodes=sum(ACGc{kk,2}.Counts);	
				MS(ii,1).(perc_prop_name)=double((ACGc{kk,2}.Counts(ix)/total_nodes)*100);

				if pc
					ACG_T=ACGc{kk,2};
					total_nodes=sum(ACG_T.Counts);
					for ll=1:numel(ACG_T.Categories)
						cat_name=matlab.lang.makeValidName(ACG_T.Categories{ll});
						MS(ii,1).(cat_name)=double((ACG_T.Counts(ll)/total_nodes)*100);
					end
				end
			end
		end

		waitbar(ii/num_files);
	end
	close(w1);

	[head_dir,~,~]=fileparts(location_of_data_files);
	if ~isempty(head_dir)
		out_shape_name=[head_dir filesep shape_name '.shp'];
	else
		out_shape_name=[shape_name '.shp'];
	end
	shapewrite(MS,out_shape_name);

end

function [mean_ksn,std_ksn,se_ksn]=ksn_convert_approx(okm,new_ref_concavity)

	g=[okm.gradient];
	a=[okm.uparea];

	ksn_calc=g./a.^-new_ref_concavity;

	mean_ksn=mean(ksn_calc,'omitnan');
	std_ksn=std(ksn_calc,'omitnan');
	se_ksn=std_ksn/sqrt(numel(ksn_calc));

end

function [mean_ksn,std_ksn,se_ksn]=ksn_convert_exact(FN,segment_length,new_ref_concavity)
	% Determine ksn method
	load(FN,'DEMoc','DEMcc','FDc','Ac','Sc','ksn_method');

	% Calculate ksn
	switch ksn_method
	case 'quick'
		[MSNc]=KSN_Quick(DEMoc,DEMcc,Ac,Sc,new_ref_concavity,segment_length);
	case 'trunk'
		[MSNc]=KSN_Trunk(DEMoc,DEMcc,Ac,Sc,new_ref_concavity,segment_length,min_order);			
	case 'trib'
		% Overide choice if very small basin as KSN_Trib will fail for small basins
		if drainage_area>2.5
			[MSNc]=KSN_Trib(DEMoc,DEMcc,FDc,Ac,Sc,new_ref_concavity,segment_length);
		else
			[MSNc]=KSN_Quick(DEMoc,DEMcc,Ac,Sc,new_ref_concavity,segment_length);
		end
	end

	% Calculate basin wide ksn statistics
	mean_ksn=mean([MSNc.ksn],'omitnan');
	std_ksn=std([MSNc.ksn],'omitnan');
	se_ksn=std_ksn/sqrt(numel(MSNc)); % Standard error
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
