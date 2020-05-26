function [T]=CompileBasinStats(location_of_data_files,varargin)
	%
	% Usage:
	%	[T]=CompileBasinStats(location_of_data_files);
	%	[T]=CompileBasinStats(location_of_data_files,'name',value,...);
	%
	% Description:
	% 	Function to take the outputs from 'ProcessRiverBasins' and 'SubDivideBigBasins' and produce a Matlab table that summarizes the results of ProcessRiverBasins
	%	and optionally SubDivideBigBasins. This table is a required input for 'BasinStatsPlots'. If additional grids were provided to 'ProcessRiverBasins', mean and 
	%	standard error values for those grids will be included in the table. This function also allows you to input a list of additional fields you wish to include 
	%	(see Optional Inputs below). There are also a variety of additional parameters / quantities that can be calculated if you provided a categorical grid
	%	to 'ProcessRiverBasins'.
	%
	% Required Inputs:
	% 		location_of_data_files - full path of folder which contains the mat files from 'ProcessRiverBasins'
	%
	% Optional Inputs:
	%		location_of_subbasins ['SubBasins'] - name of folder that contains subbasins of interest (if you created subbasins using
	%			"SubDivideBigBasins"), expected to be within the main Basin folder provided with "location_of_data_files". Note that if you do not provide
	%			the correct directory name for the location of the subbasins, subbasin values will not be included in the output regardless of your choice
	%			for the "include" parameter.
	%		include ['all'] - parameter to specify which basins to include in building the shapfile. The default 'all' will include all basin mat files in the 
	%			folder you specify. Providing 'subdivided' will check to see if a given main basin was subdivided using 'SubdivideBigBasins' and then only include 
	%			the subdivided versions of that basin (i.e. the original main basin for those subbasins will not be included in the table). Providing 'bigonly'
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
	%		dist_along_azimuth [] - option to calculate distances along a given azimuth for all basins. Expects an single numeric input, interpreted as an azimuth in degrees 
	%		filter_by_category [false] - logical flag to recalculate selected mean values based on filtering by particular categories within a categorical grid (provided to
	%			ProcessRiverBasins as 'add_cat_grids'). Requires entries to 'filter_type', 'cat_grid', and 'cat_values'. Will produce filtered values for channel steepness, gradient,
	%			and mean  elevation by default along with any additonal grids present (i.e. grids provided with 'add_grids' to ProcessRiverBasins).
	%		filter_type ['exclude'] - behavior of filter, if 'filter_by_categories' is set to true. Valid inputs are 'exclude', 'include', or 'mode'. If set to 'exclude', the filtered 
	%			means will be calculated excluding any portions of grids have the values of 'cat_values' in the 'cat_grid'. If set to 'include', filtered means will only be calculated 
	%			for portions of grids that are within specified categories. If set to 'mode', filtered means will be calculated based on the modal value of the categorical grid by basin,
	%			e.g. if the mode of basin 1 is 'grMz' and the mode of basin 2 is 'T', then the filtered mean will be calculated based on nodes that are 'grMz' in basin 1 and are 'T' in 
	%			basin 2. The idea behind this filter is if you wish to find characteristic stats for particular categories. If filter type is 'mode' then an entry for 'cat_values' is not
	%			required.
	%		cat_grid [] - name of categorical grid to use as filter, must be the same as the name provided to ProcessRiverBasins (i.e. third column in the cell array provided to
	%			'add_cat_grids').
	%		cat_values [] - 1xm cell array of categorical values of interest to use in filter. These must match valid categories in the lookup table as output from CatPoly2GRIDobj
	%			(i.e. second colmun in cell array provided to 'add_cat_grids')
	%		populate_categories [false] - logical flag to add entries that indicate the percentage of a watershed occupied by each category from a categorical grid, e.g. if you
	%			provided an entry for 'add_cat_grids' to ProcessRiverBasins that was a geologic map that had three units, 'Q', 'Mz', and 'Pz' and you set 'populate_categories' 
	%			to true there will be field names in the resulting shapefile named 'Q', 'Mz', and 'Pz' and the values stored in those columns will correspond to the percentage 
	%			of each basin covered by each unit for each basin. Setting populate_categories to true will not have any effect if no entry was provided to 'add_cat_grids' when
	%			running ProcessRiverBasins.
	%		means_by_category [] - method to calculate means of various continuous values within by categories. Requires that a categorical grid(s) was input to ProcessRiverBasins.
	%			Expects a cell 1 x m cell array where the first entry is the name of the category to use (i.e. name for categorical grid you provided to ProcessRiverBasins) and
	%			following entries are names of grids you wish to use to find means by categories, e.g. an example array might be {'geology','ksn','rlf2500','gradient'} if you 
	%			were interested in looking for patterns in channel steepness, 2.5 km^2 relief, and gradient as a function of rock type/age. Valid inputs for the grid names are:
	%				'ksn' - uses channel steepness map structure with user provided reference concavity
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
	% 
	% Examples:
	%		[T]=CompileBasinStats('/Users/You/basin_files');
	%		[T]=CompileBasinStats('/Users/You/basin_files','means_by_category',{'geology','gradient','rlf2500','rlf5000'})
	%
	%		To include recalculated channel steepness values at difference reference concavities
	%		[T]=CompileBasinStats('/Users/You/basin_files','new_concavity',[0.45 0.55 0.60]);
	%
	%		To recalculate means excluding any area of watersheds that are mapped as either 'Q' or 'Water' in the geology dataset provided to ProcessRiverBasins
	%		[T]=CompileBasinStats('/Users/You/basin_files','filter_by_categories',true,'cat_grid','geology','cat_values',{'Q','Water'},'filter_type','exclude'); 
	%
	%		To recalculate means only in the areas mapped as 'grMZ', 'grPz', or 'grpC' in the geology dataset provided to ProcessRiverBasins
	%		[T]=CompileBasinStats('/Users/You/basin_files','filter_by_categories',true,'cat_grid','geology','cat_values',{'grMz','grPz','grpC'},'filter_type','include');
	%
	% Notes
	%		-If you use 'filter_by_category' to create filtered means and uncertainites, note that the filtered value for channel steepness is calcuated using the 
	%		interpolated 'KsnOBJc', not the stream values like the the value reported in mean_ksn in the output table.
	%		
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 06/18/18 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	

	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'CompileBasinStats';
	addRequired(p,'location_of_data_files',@(x) isdir(x));

	addParameter(p,'location_of_subbasins','SubBasins',@(x) ischar(x) || isempty(x));
	addParameter(p,'include','all',@(x) ischar(validatestring(x,{'all','subdivided','bigonly'})));
	addParameter(p,'extra_field_values',[],@(x) isa(x,'cell') || isempty(x));
	addParameter(p,'extra_field_names',[],@(x) isa(x,'cell') && size(x,1)==1 || isempty(x));
	addParameter(p,'new_concavity',[],@(x) isnumeric(x));
	addParameter(p,'new_ksn_method','approximate',@(x) ischar(validatestring(x,{'approximate','exact'})));
	addParameter(p,'segment_length',1000,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'dist_along_azimuth',[],@(x) isnumeric(x) && isscalar(x) && x>=0 && x<=360 || isempty(x));
	addParameter(p,'uncertainty','se',@(x) ischar(validatestring(x,{'se','std','both'})));
	addParameter(p,'populate_categories',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'means_by_category',[],@(x) isa(x,'cell') && size(x,2)>=2 || isempty(x));
	addParameter(p,'filter_by_category',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'filter_type','exclude',@(x) ischar(validatestring(x,{'exclude','include','mode'})));
	addParameter(p,'cat_grid',[],@(x) ischar(x));
	addParameter(p,'cat_values',[],@(x) isa(x,'cell') && size(x,1)==1);


	parse(p,location_of_data_files,varargin{:});
	location_of_data_files=p.Results.location_of_data_files;

	location_of_subbasins=p.Results.location_of_subbasins;
	include=p.Results.include;
	efv=p.Results.extra_field_values;
	efn=p.Results.extra_field_names;
	new_concavity=p.Results.new_concavity;
	new_ksn_method=p.Results.new_ksn_method;
	segment_length=p.Results.segment_length;
	az=p.Results.dist_along_azimuth;
	uncertainty=p.Results.uncertainty;
	pc=p.Results.populate_categories;
	mbc=p.Results.means_by_category;
	fbc=p.Results.filter_by_category;
	ft=p.Results.filter_type;
	cgn=p.Results.cat_grid;
	cgv=p.Results.cat_values;

	% Check required entries
	if fbc && ~strcmp(ft,'mode') && isempty(cgn) | isempty(cgv)
		if isdeployed
			errordlg('For "include" or "exclude" filters, entries must be provided for both "cat_grid" and "cat_values"')
		end
		error('For "include" or "exclude" filters, entries must be provided for both "cat_grid" and "cat_values"');
	elseif fbc && strcmp(ft,'mode') && isempty(cgn)
		if isdeployed
			errordlg('For "mode" filter, entry must be provided for "cat_grid"')
		end
		error('For "mode" filter, entry must be provided for "cat_grid"');
	end

	% Deal with variability in format of locations
	if ~isempty(location_of_subbasins)
		[sub_head,~,~]=fileparts(location_of_subbasins);
		if isempty(sub_head)
			location_of_subbasins=[location_of_data_files filesep location_of_subbasins];
		end
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

	% Catch for an empty list of files
	if num_files==0
		error('The directory provided to "location_of_data_files" was a valid directory, but contained no valid basin files. Check that you have provided the correct file path.');
	end

	% Initiate Table
	T=table;

	if ~isempty(mbc)
		w1=waitbar(0,'Compiling table and calculating means by categories');
	elseif fbc
		w1=waitbar(0,'Compiling table and calculating filtered means');
	else
		w1=waitbar(0,'Compiling table');
	end

	warning off
	for ii=1:num_files;
		FileName=[FileList(ii,1).folder filesep FileList(ii,1).name];

		load(FileName,'DEMoc','RiverMouth','drainage_area','out_el','KSNc_stats','Zc_stats','Gc_stats','Centroid','hyps','Chic','DEMcc','Sc','Ac','theta_ref');

		% Populate default fields in Table
		T.ID(ii,1)=ii;
		T.river_mouth(ii,1)=RiverMouth(3);
		T.out_x(ii,1)=RiverMouth(1);
		T.out_y(ii,1)=RiverMouth(2);
		T.center_x(ii,1)=Centroid(1);
		T.center_y(ii,1)=Centroid(2);
		T.drainage_area(ii,1)=drainage_area;
		T.outlet_elevation(ii,1)=out_el;
		T.mean_el(ii,1)=Zc_stats(1);
		T.max_el(ii,1)=Zc_stats(5);
		switch uncertainty
		case 'se'
			T.se_el(ii,1)=Zc_stats(2);
		case 'std'
			T.std_el(ii,1)=Zc_stats(3);
		case 'both'
			T.se_el(ii,1)=Zc_stats(2);
			T.std_el(ii,1)=Zc_stats(3);
		end

		T.mean_ksn(ii,1)=KSNc_stats(1);
		switch uncertainty
		case 'se'
			T.se_ksn(ii,1)=KSNc_stats(2);
		case 'std'
			T.std_ksn(ii,1)=KSNc_stats(3);
		case 'both'
			T.se_ksn(ii,1)=KSNc_stats(2);
			T.std_ksn(ii,1)=KSNc_stats(3);
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
				T.(ksn_cat_name)(ii,1)=mean_ksn;
				switch uncertainty
				case 'se'
					ksn_cat_name_se=matlab.lang.makeValidName(['se_ksn_' num2str(new_concavity(jj))]);
					T.(ksn_cat_name_se)(ii,1)=se_ksn;
				case 'std'
					ksn_cat_name_std=matlab.lang.makeValidName(['std_ksn_' num2str(new_concavity(jj))]);
					T.(ksn_cat_name_std)(ii,1)=std_ksn;
				case 'both'
					ksn_cat_name_se=matlab.lang.makeValidName(['se_ksn_' num2str(new_concavity(jj))]);
					T.(ksn_cat_name_se)(ii,1)=se_ksn;
					ksn_cat_name_std=matlab.lang.makeValidName(['std_ksn_' num2str(new_concavity(jj))]);
					T.(ksn_cat_name_std)(ii,1)=std_ksn;					
				end
			end
		end

		T.mean_gradient(ii,1)=Gc_stats(1);
		switch uncertainty
		case 'se'
			T.se_gradient(ii,1)=Gc_stats(2);
		case 'std'
			T.std_gradient(ii,1)=Gc_stats(3);
		case 'both'
			T.se_gradient(ii,1)=Gc_stats(2);
			T.std_gradient(ii,1)=Gc_stats(3);
		end

		T.hypsometry{ii,1}=hyps;
		T.hyp_integral(ii,1)=abs(trapz((hyps(:,2)-min(hyps(:,2)))/(max(hyps(:,2))-min(hyps(:,2))),hyps(:,1)/100));
		T.concavity(ii,1)=Chic.mn;

		c=chiplot(Sc,DEMcc,Ac,'a0',1,'mn',theta_ref,'plot',false);
		T.chi_R_squared(ii,1)=c.R2;
		c_trunk=chiplot(trunk(Sc),DEMcc,Ac,'a0',1,'mn',theta_ref,'plot',false);
		T.chi_R_squared_trunk(ii,1)=c_trunk.R2;

		% Check for additional grids within the process river basins output
		VarList=whos('-file',FileName);
		AgInd=find(strcmp(cellstr(char(VarList.name)),'AGc'));
		RlfInd=find(strcmp(cellstr(char(VarList.name)),'rlf'));
		AcgInd=find(strcmp(cellstr(char(VarList.name)),'ACGc'));

		if ~isempty(AgInd)
			load(FileName,'AGc','AGc_stats');
			num_grids=size(AGc,1);

			for kk=1:num_grids
				mean_prop_name=['mean_' AGc{kk,2}];		
				T.(mean_prop_name)(ii,1)=double(AGc_stats(kk,1));

				switch uncertainty
				case 'se'
					se_prop_name=['se_' AGc{kk,2}];
					T.(se_prop_name)(ii,1)=double(AGc_stats(kk,2));
				case 'std'
					std_prop_name=['std_' AGc{kk,2}];
					T.(std_prop_name)(ii,1)=double(AGc_stats(kk,3));
				case 'both'
					se_prop_name=['se_' AGc{kk,2}];
					T.(se_prop_name)(ii,1)=double(AGc_stats(kk,2));
					std_prop_name=['std_' AGc{kk,2}];
					T.(std_prop_name)(ii,1)=double(AGc_stats(kk,3));
				end
			end

			ag_flag=true;
		else
			ag_flag=false;
		end		


		if ~isempty(RlfInd)
			load(FileName,'rlf','rlf_stats');
			num_grids=size(rlf,1);

			for kk=1:num_grids
				mean_prop_name=['mean_rlf' num2str(rlf{kk,2})];
				T.(mean_prop_name)(ii,1)=double(rlf_stats(kk,1));

				switch uncertainty
				case 'se'
					se_prop_name=['se_rlf' num2str(rlf{kk,2})];
					T.(se_prop_name)(ii,1)=double(rlf_stats(kk,2));
				case 'std'
					std_prop_name=['std_rlf' num2str(rlf{kk,2})];
					T.(std_prop_name)(ii,1)=double(rlf_stats(kk,3));
				case 'both'
					se_prop_name=['se_rlf' num2str(rlf{kk,2})];
					T.(se_prop_name)(ii,1)=double(rlf_stats(kk,2));
					std_prop_name=['std_rlf' num2str(rlf{kk,2})];
					T.(std_prop_name)(ii,1)=double(rlf_stats(kk,3));
				end
			end
			rlf_flag=true;
		else
			rlf_flag=false;
		end		

		% Calculate filtered values

		if fbc & ~isempty(AcgInd)
			load(FileName,'ACGc');
			% Isolate Cat Grid and lookup table of interest
			cix=find(strcmp(ACGc(:,3),cgn));
			CG=ACGc{cix,1};
			cgt=ACGc{cix,2};
			% Create filter 
			F=GRIDobj(CG,'logical');
			if strcmp(ft,'include')
				% Find entries that match values of interest
				vcix=ismember(cgt.Categories,cgv);
				vnix=cgt.Numbers(vcix);
				F.Z=ismember(CG.Z,vnix);
			elseif strcmp(ft,'exclude')
				% Find entries that match values of interest
				vcix=ismember(cgt.Categories,cgv);
				vnix=cgt.Numbers(vcix);
				F.Z=~ismember(CG.Z,vnix);
			elseif strcmp(ft,'mode')
				load(FileName,'ACGc_stats');
				F.Z=ismember(CG.Z,ACGc_stats(cix,1));
			end
			% Apply filter
			load(FileName,'DEMoc','Goc','MSNc');

			T.mean_el_f(ii,1)=mean(DEMoc.Z(F.Z),'omitnan');
			switch uncertainty
			case 'se'
				T.se_el_f(ii,1)=std(DEMoc.Z(F.Z),'omitnan')/sqrt(sum(~isnan(DEMoc.Z(F.Z))));
			case 'std'
				T.std_el_f(ii,1)=std(DEMoc.Z(F.Z),'omitnan');
			case 'both'
				T.se_el_f(ii,1)=std(DEMoc.Z(F.Z),'omitnan')/sqrt(sum(~isnan(DEMoc.Z(F.Z))));
				T.std_el_f(ii,1)=std(DEMoc.Z(F.Z),'omitnan');
			end

			T.mean_gradient_f(ii,1)=mean(Goc.Z(F.Z),'omitnan');
			switch uncertainty
			case 'se'
				T.se_gradient_f(ii,1)=std(Goc.Z(F.Z),'omitnan')/sqrt(sum(~isnan(Goc.Z(F.Z))));
			case 'std'
				T.std_gradient_f(ii,1)=std(Goc.Z(F.Z),'omitnan');
			case 'both'
				T.se_gradient_f(ii,1)=std(Goc.Z(F.Z),'omitnan')/sqrt(sum(~isnan(Goc.Z(F.Z))));
				T.std_gradient_f(ii,1)=std(Goc.Z(F.Z),'omitnan');
			end

			KSNG=GRIDobj(CG);
			KSNG.Z(:,:)=NaN;
			for kk=1:numel(MSNc)
				ix=coord2ind(CG,MSNc(kk).X,MSNc(kk).Y);
				KSNG.Z(ix)=MSNc(kk).ksn;
			end

			T.mean_ksn_f(ii,1)=mean(KSNG.Z(F.Z),'omitnan');
			switch uncertainty
			case 'se'
				T.se_ksn_f(ii,1)=std(KSNG.Z(F.Z),'omitnan')/sqrt(sum(~isnan(KSNG.Z(F.Z))));
			case 'std'
				T.std_ksn_f(ii,1)=std(KSNG.Z(F.Z),'omitnan');
			case 'both'
				T.se_ksn_f(ii,1)=std(KSNG.Z(F.Z),'omitnan')/sqrt(sum(~isnan(KSNG.Z(F.Z))));
				T.std_ksn_f(ii,1)=std(KSNG.Z(F.Z),'omitnan');
			end

			if ag_flag
				ag_grids=size(AGc,1);
				for kk=1:ag_grids
					agG=AGc{kk,1};
					mean_prop_name=['mean_' AGc{kk,2} '_f'];		
					T.(mean_prop_name)(ii,1)=mean(agG.Z(F.Z),'omitnan');

					switch uncertainty
					case 'se'
						se_prop_name=['se_' AGc{kk,2} '_f'];
						T.(se_prop_name)(ii,1)=std(agG.Z(F.Z),'omitnan')/sqrt(sum(~isnan(agG.Z(F.Z))));
					case 'std'
						std_prop_name=['std_' AGc{kk,2} '_f'];
						T.(std_prop_name)(ii,1)=std(agG.Z(F.Z),'omitnan');
					case 'both'
						se_prop_name=['se_' AGc{kk,2} '_f'];
						T.(se_prop_name)(ii,1)=std(agG.Z(F.Z),'omitnan')/sqrt(sum(~isnan(agG.Z(F.Z))));
						std_prop_name=['std_' AGc{kk,2} '_f'];
						T.(std_prop_name)(ii,1)=std(agG.Z(F.Z),'omitnan');
					end
				end
			end

			if rlf_flag
				rlf_grids=size(rlf,1);
				for kk=1:rlf_grids
					rlfG=rlf{kk,1};
					mean_prop_name=['mean_rlf' num2str(rlf{kk,2}) '_f'];
					T.(mean_prop_name)(ii,1)=mean(rlfG.Z(F.Z),'omitnan');

					switch uncertainty
					case 'se'
						se_prop_name=['se_rlf' num2str(rlf{kk,2}) '_f'];
						T.(se_prop_name)(ii,1)=std(rlfG.Z(F.Z),'omitnan')/sqrt(sum(~isnan(rlfG.Z(F.Z))));
					case 'std'
						std_prop_name=['std_rlf' num2str(rlf{kk,2}) '_f'];
						T.(std_prop_name)(ii,1)=std(rlfG.Z(F.Z),'omitnan');
					case 'both'
						se_prop_name=['se_rlf' num2str(rlf{kk,2}) '_f'];
						T.(se_prop_name)(ii,1)=std(rlfG.Z(F.Z),'omitnan')/sqrt(sum(~isnan(rlfG.Z(F.Z))));
						std_prop_name=['std_rlf' num2str(rlf{kk,2})];
						T.(std_prop_name)(ii,1)=std(rlfG.Z(F.Z),'omitnan');
					end
				end
			end
			% Generate column to record filter
			if ~strcmp(ft,'mode')
				filt_name=join(cgv);
				filt_name=filt_name{1};
				T.filter{ii,1}=[ft ' ' filt_name];
			else
				T.filter{ii,1}=[cgn ' mode'];
			end

		elseif fbc & isempty(AcgInd)
			if isdeployed
				errordlg('No Categorical Grids were provided to ProcessRiverBasins so filtered values cannot be calculated')
			end
			error('No Categorical Grids were provided to ProcessRiverBasins so filtered values cannot be calculated');
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
						T.(field_name){ii,1}=field_value;
					elseif isnumeric(field_value)
						T.(field_name)(ii,1)=double(field_value);
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

		if ~isempty(AcgInd)
			load(FileName,'ACGc','ACGc_stats');
			num_grids=size(ACGc,1);

			for kk=1:num_grids
				mode_prop_name=['mode_' ACGc{kk,3}];
				perc_prop_name=['mode_' ACGc{kk,3} '_percent'];
				ix=find(ACGc{kk,2}.Numbers==ACGc_stats(kk,1),1);
				T.(mode_prop_name){ii,1}=ACGc{kk,2}.Categories{ix};
				total_nodes=sum(ACGc{kk,2}.Counts);				
				T.(perc_prop_name)(ii,1)=double((ACGc{kk,2}.Counts(ix)/total_nodes)*100);

				if pc
					ACG_T=ACGc{kk,2};
					total_nodes=sum(ACG_T.Counts);
					for ll=1:numel(ACG_T.Categories)
						cat_name=ACG_T.Categories{ll};
						cat_name=matlab.lang.makeValidName([ACGc{kk,3} '_perc_' cat_name]);
						T.(cat_name)(ii,1)=double((ACG_T.Counts(ll)/total_nodes)*100);
					end
				end

				if ~isempty(mbc)
					warn_flag=false;
					% Partition input
					cg=mbc(1);
					dg=mbc(2:end);
					num_dg=numel(dg);
					% Find categorical grid of interest
					cix=find(strcmp(ACGc(:,3),cg));
					ACG=ACGc{cix,1}; % GRID
					ACG_T=ACGc{cix,2}; %look up table
					% Iterate through categories
					for ll=1:numel(ACG_T.Categories)
						IDX=GRIDobj(ACG,'logical');
						IDX.Z=ismember(ACG.Z,ACG_T.Numbers(ll));
						cat_name=ACG_T.Categories{ll};
						for mm=1:num_dg
							dgOI=dg{mm};
							if strcmp(dgOI,'ksn')
								load(FileName,'MSNc');
								KSNG=GRIDobj(ACG);
								KSNG.Z(:,:)=NaN;
								for oo=1:numel(MSNc)
									ix=coord2ind(ACG,MSNc(oo).X,MSNc(oo).Y);
									KSNG.Z(ix)=MSNc(oo).ksn;
								end
								cat_nameN=matlab.lang.makeValidName(['mksn_' cat_name]);
								T.(cat_nameN)(ii,1)=mean(KSNG.Z(IDX.Z),'omitnan');
							elseif strcmp(dgOI,'gradient')
								load(FileName,'Goc');
								cat_nameN=matlab.lang.makeValidName(['mgrad_' cat_name]);
								T.(cat_nameN)(ii,1)=mean(Goc.Z(IDX.Z),'omitnan');
							elseif regexp(dgOI,regexptranslate('wildcard','rlf*'))
								rlfval=str2num(strrep(dgOI,'rlf',''));
								rlfix=find(cell2mat(rlf(:,2))==rlfval);
								if ~isempty(rlfix)
									Rg=rlf{rlfix,1};
									cat_nameN=matlab.lang.makeValidName(['mr' num2str(rlfval) '_' cat_name]);
									T.(cat_nameN)(ii,1)=mean(Rg.Z(IDX.Z),'omitnan');	
								end								
							else 
								try
									dgix=find(strcmp(AGc(:,2),dgOI));
									AGcOI=AGc{dgix,1};
									cat_nameN=matlab.lang.makeValidName(['m' AGc{dgix,2} '_' cat_name]);
									T.(cat_nameN)(ii,1)=mean(AGcOI.Z(IDX.Z),'omitnan');
								catch
									warn_flag=true;
								end
							end
						end
					end
				end
			end
		end	

		T.file_path{ii,1}=FileName;

		waitbar(ii/num_files);
	end
	warning on

	if ~isempty(az)
		% Rotate by center of dataset
		x0=mean(T.center_x);
		y0=mean(T.center_y);
		% Convert provided azimuth
		azn=az-90;
		% Do Rotation
		d=(T.center_x-x0).*cosd(azn)-(T.center_y-y0).*sind(azn);
		% Normalize distance
		d=d-min(d);	
		% Add to the table
		az_name=['dist_along_' num2str(round(az))];
		T=addvars(T,d,'NewVariableNames',az_name,'After','hyp_integral');
	end

	if ~isempty(mbc)
		if warn_flag==true
			if isdeployed
				warndlg('One or more input for grid names to "means_by_category" was not recognized, table compiled without this entry')
			end
			warning('One or more input for grid names to "means_by_category" was not recognized, table compiled without this entry')
		end
	end

	close(w1);
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
