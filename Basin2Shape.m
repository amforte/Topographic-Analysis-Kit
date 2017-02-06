function [MS]=Basin2Shape(DEM,location_of_data_files,varargin)
	% Function to take the outputs from 'ProcessRiverBasins' and 'SubDivideBigBasins' and produce a single shapefile showing the outlines of polygons
	% and with commonly desired attributes from the results of 'ProcessRiverBasins' etc. See below for a full list of fields that the output shapefile
	% will include. If additional grids were provided to 'ProcessRiverBasins', mean and standard error values for those grids will be auto-populated in
	% the shapefile and the name of the fields will be the character array provided in the second column of additional grids input. This function also
	% allows you to input a list of additional fields you wish to include (see Optional Inputs below). If you would rather create a GRIDobj with specified
	% values, use 'Basin2Raster'.
	%
	% Required Inputs:
	%		DEM - GRIDobj of the DEM originally used as input for 'ProcessRiverBasins' for the basins of interest.
	% 		location_of_data_files - full path of folder which contains the mat files from 'ProcessRiverBasins'
	%
	% Optional Inputs:
	%		shape_name ['basins'] - name for the shapefile to be export, must have no spaces to be a valid name for ArcGIS and should NOT include the '.shp';
	%		include ['all'] - parameter to specify which basins to include in building the shapfile. The default 'all' will include all basin mat files in the 
	%			folder you specify. The optional 'subdivided' will check to see if a given main basin was subdivided using 'SubdivideBigBasin' and then only include 
	%			the subdivided versions of that basin (i.e. the original main basin for those subbasins will not be included in the shapefile). 
	%		extra_field_values [] - cell array of extra field values you wish to include. The first column in this cell array must be the river basin number
	%			(i.e. the identifying number in the third column of the RiverMouth input to ProcessRiverBasins or the number generated for the basin in
	%			SubDivideBigBasins). Only one row per river basin number is allowed and ALL river basin numbers in the basins being processed must have a value
	%			associated with them. Additional columns are interpreted as the values with which you wish to populate the extra fields. These can either be character
	%			arrays or numbers, other values will results in an error. 
	%		extra_field_names [] - a 1 x m cell array of field names, as characters (no spaces as this won't work with shapefile attributes), associated with the field values. 
	%			These must be in the same order as values are given in extra_field_values. If for example your extra_field_values cell array is 3 columns with the river number, 
	%			sample name, and erosion rate then your extra_field_names cell array should include entries for 'sample_name' and 'erosion_rate' in that order. 
	%
	% Output:
	%		Outputs a mapstructure and saves a shapefile with the following default fields:
	%			river_mouth - river mouth number provided to ProcessRiverBasins
	%			drainage_area - drainage area of basin in km^2
	%			outlet_elevation - elevation of pour point in m
	%			mean_el - mean elevation of basin in meters
	%			max_el - maximum elevation of basin in meters
	%			se_el - standard error on basin mean elevation
	%			mean_ksn - mean channel steepenss
	%			se_ksn - standard error on channel steepness
	%			mean_gradient - mean gradient
	%			se_gradient - standard error on gradient
	%		Mean and standard error values will be populated for any additional grids
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Last Revised Winter 2017 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	

	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'Basin2Shape';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'location_of_data_files',@(x) isdir(x));

	addParamValue(p,'shape_name','basins',@(x) ischar(x));
	addParamValue(p,'include','all',@(x) ischar(validatestring(x,{'all','subdivided'})));
	addParamValue(p,'extra_field_values',[],@(x) isa(x,'cell'));
	addParamValue(p,'extra_field_names',[],@(x) isa(x,'cell') & size(x,1)==1);

	parse(p,DEM,location_of_data_files,varargin{:});
	DEM=p.Results.DEM;
	location_of_data_files=p.Results.location_of_data_files;

	shape_name=p.Results.shape_name;
	include=p.Results.include;
	efv=p.Results.extra_field_values;
	efn=p.Results.extra_field_names;

	current=pwd;
	cd(location_of_data_files);

	% Switch for which basins to include
	switch include
	case 'all'
		FileList=dir('Basin*.mat');
		num_files=numel(FileList);
	case 'subdivided'
		AllFullFiles=dir('*_Data.mat');
		num_basins=numel(AllFullFiles);
		basin_nums=zeros(num_basins,1);
		for jj=1:num_basins
			fileName=AllFullFiles(jj,1).name;
			basin_nums(jj)=sscanf(fileName,'%*6s %i'); %%%
		end

		FileCell=cell(num_basins,1);
		for kk=1:num_basins
			basin_num=basin_nums(kk);
			SearchAllString=['*_' num2str(basin_num) '_Data.mat'];
			SearchSubString=['*_' num2str(basin_num) '_DataSubset*.mat'];

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

	% Initiate Map Structure
	MS=struct;

	w1=waitbar(0,'Building polygons');
	for ii=1:num_files;
		FileName=FileList(ii,1).name;
		DB=GRIDobj(DEM);

		load(FileName,'DEMoc','RiverMouth','drainage_area','out_el','KSNc_stats','Zc_stats','Gc_stats');

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
		MS(ii,1).drainage_area=drainage_area;
		MS(ii,1).outlet_elevation=out_el;
		MS(ii,1).mean_el=Zc_stats(1);
		MS(ii,1).max_el=Zc_stats(5);
		MS(ii,1).se_el=Zc_stats(2);
		MS(ii,1).mean_ksn=KSNc_stats(1);
		MS(ii,1).se_ksn=KSNc_stats(2);
		MS(ii,1).mean_gradient=Gc_stats(1);
		MS(ii,1).se_gradient=Gc_stats(2);
		
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
		VarInd=find(strcmp(cellstr(char(VarList.name)),'AGc'));

		if ~isempty(VarInd)
			load(FileName,'AGc','AGc_stats');
			num_grids=size(AGc,1);

			for kk=1:num_grids
				mean_prop_name=['mean_' AGc{kk,2}];
				se_prop_name=['se_' AGc{kk,2}];
				MS(ii,1).(mean_prop_name)=double(AGc_stats(kk,1));
				MS(ii,1).(se_prop_name)=double(AGc_stats(kk,2));	
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
						error(['Extra field value provided for ' field_name ' is neither numeric or a character']);
					end
				end
			elseif numel(ix)>1
				error(['More than one entry was provided for extra fields for basin ' num2str(RiverMouth(:,3))]);
			elseif isempty(ix)
				error(['No one entry was provided for extra field values for basin ' num2str(RiverMouth(:,3))]);
			end
		end


		waitbar(ii/num_files);
	end
	close(w1);

	cd(current);
	out_shape_name=[shape_name '.shp'];
	shapewrite(MS,out_shape_name);
