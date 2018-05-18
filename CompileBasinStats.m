function [T]=CompileBasinStats(location_of_data_files,varargin)
	% Function to take the outputs from 'ProcessRiverBasins' and 'SubDivideBigBasins' and produce a single shapefile showing the outlines of polygons
	% 	and with commonly desired attributes from the results of 'ProcessRiverBasins' etc. See below for a full list of fields that the output shapefile
	% 	will include. If additional grids were provided to 'ProcessRiverBasins', mean and standard error values for those grids will be auto-populated in
	% 	the shapefile and the name of the fields will be the character array provided in the second column of additional grids input. This function also
	% 	allows you to input a list of additional fields you wish to include (see Optional Inputs below). If you would rather create a GRIDobj with specified
	% 	values, use 'Basin2Raster'.
	%
	% Required Inputs:
	% 		location_of_data_files - full path of folder which contains the mat files from 'ProcessRiverBasins'
	%
	% Optional Inputs:
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
	%		uncertainty ['se'] - parameter to control which measure of uncertainty is included, expects 'se' for standard error (default), 'std' for standard deviation, or 'both'
	%			to include both standard error and deviation.
	%		populate_categories [false] - logical flag to add entries that indicate the percentage of a watershed occupied by each category from a categorical grid, e.g. if you
	%			provided an entry for 'add_cat_grids' to ProcessRiverBasins that was a geologic map that had three units, 'Q', 'Mz', and 'Pz' and you set 'populate_categories' 
	%			to true there will be field names in the resulting shapefile named 'Q', 'Mz', and 'Pz' and the values stored in those columns will correspond to the percentage 
	%			of each basin covered by each unit for each basin. Setting populate_categories to true will not have any effect if no entry was provided to 'add_cat_grids' when
	%			running ProcessRiverBasins.
	%		means_by_category [] - method to calculate means of various continuous values within by categories. Requires that a categorical grid was input to ProcessRiverBasins.
	%			Expects a cell 1 x m cell array where the first entry is the name of the category to use (i.e. name for categorical grid you provided to ProcessRiverBasins) and
	%			following entries are names of grids you wish to use to find means by categories, e.g. an example array might be {'geology','ksn','rlf2500','gradient'} if you 
	%			were interested in looking for patterns in channel steepness, 2.5 km^2 relief, and gradient as a function of rock type/age. Valid inputs for the grid names are:
	%				'ksn' - uses interpolated channel steepness grid
	%				'gradient' - uses gradient grid
	%				'rlf####' - where #### is the radius you provided to ProcessRiverBasins (requires that 'calc_relief' was set to true when running ProcessRiverBasins
	%				'NAME' - where NAME is the name of an additional grid provided with the 'add_grids' option to ProcessRiverBasins
	%
	% Output:
	%		Outputs a table (T) with the following default fields:
	%			river_mouth - river mouth number provided to ProcessRiverBasins
	%			drainage_area - drainage area of basin in km^2
	%			out_x - x coordinate of basin mouth
	%			out_y - y coordinate of basin mouth
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
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Last Revised Spring 2018 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	

	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'CompileBasinStats';
	addRequired(p,'location_of_data_files',@(x) isdir(x));

	addParamValue(p,'include','all',@(x) ischar(validatestring(x,{'all','subdivided','bigonly'})));
	addParamValue(p,'extra_field_values',[],@(x) isa(x,'cell'));
	addParamValue(p,'extra_field_names',[],@(x) isa(x,'cell') && size(x,1)==1);
	addParamValue(p,'uncertainty','se',@(x) ischar(validatestring(x,{'se','std','both'})));
	addParamValue(p,'populate_categories',false,@(x) isscalar(x) && islogical(x))
	addParamValue(p,'means_by_category',[],@(x) isa(x,'cell' && size(x,1)>=2));


	parse(p,location_of_data_files,varargin{:});
	location_of_data_files=p.Results.location_of_data_files;

	include=p.Results.include;
	efv=p.Results.extra_field_values;
	efn=p.Results.extra_field_names;
	uncertainty=p.Results.uncertainty;
	pc=p.Results.populate_categories;
	mbc=p.Results.means_by_category;

	if ~isempty(mbc)
		disp('Means_by_category option has not been fully tested')
	end

	current=pwd;
	cd(location_of_data_files);

	% Switch for which basins to include
	switch include
	case 'all'
		FileList=dir('Basin*.mat');
		num_files=numel(FileList);
	case 'bigonly'
		FileList=dir('*_Data.mat');
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

	% Initiate Table
	T=table;

	w1=waitbar(0,'Compiling table');
	warning off
	for ii=1:num_files;
		FileName=FileList(ii,1).name;

		load(FileName,'DEMoc','RiverMouth','drainage_area','out_el','KSNc_stats','Zc_stats','Gc_stats','Centroid');

		% Populate default fields in Table
		T.ID(ii)=ii;
		T.river_mouth(ii)=RiverMouth(3);
		T.out_x(ii)=RiverMouth(1);
		T.out_y(ii)=RiverMouth(2);
		T.center_x(ii)=Centroid(1);
		T.center_y(ii)=Centroid(2);
		T.drainage_area(ii)=drainage_area;
		T.outlet_elevation(ii)=out_el;
		T.mean_el(ii)=Zc_stats(1);
		T.max_el(ii)=Zc_stats(5);
		T.mean_ksn(ii)=KSNc_stats(1);
		T.mean_gradient(ii)=Gc_stats(1);

		switch uncertainty
		case 'se'
			T.se_el(ii)=Zc_stats(2);
			T.se_ksn(ii)=KSNc_stats(2);
			T.se_gradient(ii)=Gc_stats(2);
		case 'std'
			T.std_el(ii)=Zc_stats(3);
			T.std_ksn(ii)=KSNc_stats(3);
			T.std_gradient(ii)=Gc_stats(3);
		case 'both'
			T.se_el(ii)=Zc_stats(2);
			T.se_ksn(ii)=KSNc_stats(2);
			T.se_gradient(ii)=Gc_stats(2);
			T.std_el(ii)=Zc_stats(3);
			T.std_ksn(ii)=KSNc_stats(3);
			T.std_gradient(ii)=Gc_stats(3);
		end

		% Check for additional grids within the process river basins output
		VarList=whos('-file',FileName);
		VarInd=find(strcmp(cellstr(char(VarList.name)),'AGc'));

		if ~isempty(VarInd)
			load(FileName,'AGc','AGc_stats');
			num_grids=size(AGc,1);

			for kk=1:num_grids
				mean_prop_name=['mean_' AGc{kk,2}];		
				T.(mean_prop_name)(ii)=double(AGc_stats(kk,1));

				switch uncertainty
				case 'se'
					se_prop_name=['se_' AGc{kk,2}];
					T.(se_prop_name)(ii)=double(AGc_stats(kk,2));
				case 'std'
					std_prop_name=['std_' AGc{kk,2}];
					T.(std_prop_name)(ii)=double(AGc_stats(kk,3));
				case 'both'
					se_prop_name=['se_' AGc{kk,2}];
					T.(se_prop_name)(ii)=double(AGc_stats(kk,2));
					std_prop_name=['std_' AGc{kk,2}];
					T.(std_prop_name)(ii)=double(AGc_stats(kk,3));
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
				T.(mean_prop_name)(ii)=double(rlf_stats(kk,1));
				T.(se_prop_name)(ii)=double(rlf_stats(kk,2));

				switch uncertainty
				case 'se'
					se_prop_name=['se_rlf' num2str(rlf{kk,2})];
					T.(se_prop_name)(ii)=double(rlf_stats(kk,2));
				case 'std'
					std_prop_name=['std_rlf' num2str(rlf{kk,2})];
					T.(std_prop_name)(ii)=double(rlf_stats(kk,3));
				case 'both'
					se_prop_name=['se_rlf' num2str(rlf{kk,2})];
					T.(se_prop_name)(ii)=double(rlf_stats(kk,2));
					std_prop_name=['std_rlf' num2str(rlf{kk,2})];
					T.(std_prop_name)(ii)=double(rlf_stats(kk,3));
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
						T.(field_name)(ii)=field_value;
					elseif isnumeric(field_value)
						T.(field_name)(ii)=double(field_value);
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

		VarInd=find(strcmp(cellstr(char(VarList.name)),'ACGc'));
		if ~isempty(VarInd)
			load(FileName,'ACGc','ACGc_stats');
			num_grids=size(ACGc,1);

			for kk=1:num_grids
				mode_prop_name=['mode_' ACGc{kk,3}];
				ix=find(ACGc{kk,2}.Numbers==ACGc_stats(kk,1),1);
				T.(mode_prop_name){ii}=ACGc{kk,2}.Categories{ix};

				if pc
					ACG_T=ACGc{kk,2};
					total_nodes=sum(ACG_T.Counts);
					for ll=1:numel(ACG_T.Categories)
						cat_name=matlab.lang.makeValidName(ACG_T.Categories{ll});
						cat_name=[ACGc{kk,3} '_perc_' cat_name];
						T.(cat_name)(ii)=double((ACG_T.Counts(ll)/total_nodes)*100);
					end
				end

				if ~isempty(mbc)
					% Partition input
					cg=mbc(1);
					dg=mbc(2:end);
					num_dg=numel(dg);
					% Find categorical grid of interest
					cix=find(ACGc(:,3),cg);
					ACG=ACGc{cix,1}; % GRID
					ACG_T=ACGc{cix,2}; %look up table
					% Iterate through categories
					for ll=1:numel(ACG_T.Categories)
						IDX=GRIDobj(ACG,'logical');
						IDX.Z=ismember(ACG.Z,ACG_T.Numbers(ll));
						cat_name=matlab.lang.makeValidName(ACG_T.Categories{ll});
						for mm=1:num_dg
							dgOI=dg{mm};
							if strcmp(dgOI,'ksn')
								load(FileName,'KsnOBJc');
								cat_nameN=['mksn_' cat_name];
								T.(cat_nameN)(ii)=nanmean(KsnOBJc.Z(IDX.Z));
							elseif strcmp(dgOI,'gradient')
								load(FileName,'Goc');
								cat_nameN=['mgrad_' cat_name];
								T.(cat_nameN)(ii)=nanmean(Goc.Z(IDX.Z));
							elseif regexp(dgOI,regexptranslate('wildcard','rlf*'))
								rlfval=str2num(strrep(dgOI,'rlf',''));
								rlfix=find(rlf(:,2)==rlfval);
								if ~isempty(rlfix)
									Rg=rlf{rlfix,1};
									cat_nameN=['mr' num2str(rlfval) '_' cat_name];
									T.(cat_nameN)(ii)=nanmean(Rg.Z(IDX.Z));	
								end								
							else 
								dgix=find(strcmp(AGc(:,2),dgOI));
								AGcOI=AGc{dgix,1};
								cat_nameN=['m' AGc{dgix,2} '_' cat_name];
								T.(cat_nameN)(ii)=nanmean(AGcOI.Z(IDX.Z));
							end
						end
					end
				end
			end
		end	


		waitbar(ii/num_files);
	end
	warning on
	close(w1);

	cd(current);

end
