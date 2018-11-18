function cmpCompileBasinStats(wdir,location_of_data_files,varargin)
	% Description:
	% 	Function to take the outputs from 'ProcessRiverBasins' and 'SubDivideBigBasins' and produce a Matlab table that summarizes the results of ProcessRiverBasins
	%	and optionally SubDivideBigBasins. This table is a required input for 'BasinStatsPlots'. If additional grids were provided to 'ProcessRiverBasins', mean and 
	%	standard error values for those grids will be included in the table. This function also allows you to input a list of additional fields you wish to include 
	%	(see Optional Inputs below). There are also a variety of additional parameters / quantities that can be calculated if you provided a categorical grid
	%	to 'ProcessRiverBasins'.
	%
	% Required Inputs:
	%		wdir - full path of working directory
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
	%		file_name_prefix ['Basins'] - parameter to specify a file name prefix for the output tables 
	%		extra_field_values [] - name of text file of extra field values you wish to include. The first column in this file must be the river basin number
	%			(i.e. the identifying number in the third column of the RiverMouth input to ProcessRiverBasins or the number generated for the basin in
	%			SubDivideBigBasins). Only one row per river basin number is allowed and ALL river basin numbers in the basins being processed must have a value
	%			associated with them. Additional columns are interpreted as the values with which you wish to populate the extra fields. These can either be character
	%			arrays or numbers, other values will results in an error. The function will use the header names within this file to name fields in the output shapefile 
	%		new_concavity [] - a 1 x m array of concavity values to recalculate normalized channel steepness statistics (mean, standard error and/or standard deviation) using the
	%			provided concavities.
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
	%		cat_values [] - name of text file containing single row of comma separated categorical values of interest to use in filter. These must match valid categories in the field of 
	%			interest in the shapefile provided to PrepareAddCatGrids to prepare the grid name in the 'cat_grid' function
	%		populate_categories [false] - logical flag to add entries that indicate the percentage of a watershed occupied by each category from a categorical grid, e.g. if you
	%			provided an entry for 'add_cat_grids' to ProcessRiverBasins that was a geologic map that had three units, 'Q', 'Mz', and 'Pz' and you set 'populate_categories' 
	%			to true there will be field names in the resulting shapefile named 'Q', 'Mz', and 'Pz' and the values stored in those columns will correspond to the percentage 
	%			of each basin covered by each unit for each basin. Setting populate_categories to true will not have any effect if no entry was provided to 'add_cat_grids' when
	%			running ProcessRiverBasins.
	%		means_by_category [] - method to calculate means of various continuous values within by categories. Requires that a categorical grid(s) was input to ProcessRiverBasins.
	%			Expects a text file containing a single row with comma separated entries, where the first entry is the name of the category to use (i.e. name for categorical grid 
	%			you provided to PrepareAddCatGrids) following entries are names of grids you wish to use to find means by categories, e.g. an example single row in an input table 
	%			would be 'geology,ksn,rlf2500,gradient' (without the quotes) if you were interested in looking for patterns in channel steepness, 2.5 km^2 relief, and gradient as 
	%			a function of rock type/age. Valid inputs for the grid names are:
	%				ksn - uses channel steepness map structure with user provided reference concavity
	%				gradient - uses gradient grid
	%				rlf#### - where #### is the radius you provided to ProcessRiverBasins (requires that 'calc_relief' was set to true when running ProcessRiverBasins
	%				NAME - where NAME is the name of an additional grid provided with the 'add_grids' option to ProcessRiverBasins
	%
	% Output:
	%		Outputs a table as a text with the following default fields:
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
	%		Also saves a matfile for use in 'cmpBasinStatsPlots' or 'cmpMakeCombinedSwath'
	%
	%
	% Notes
	%		-If you use 'filter_by_category' to create filtered means and uncertainites, note that the filtered value for channel steepness is calcuated using the 
	%		interpolated 'KsnOBJc', not the stream values like the the value reported in mean_ksn in the output table.
	%		
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
	p.FunctionName = 'cmpCompileBasinStats';
	addRequired(p,'wdir',@(x) ischar(x));
	addRequired(p,'location_of_data_files',@(x) ischar(x));

	addParameter(p,'location_of_subbasins','SubBasins',@(x) ischar(x));
	addParameter(p,'include','all',@(x) ischar(validatestring(x,{'all','subdivided','bigonly'})));
	addParameter(p,'extra_field_values',[],@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.txt'))));
	addParameter(p,'new_concavity',[],@(x) isnumeric(x));
	addParameter(p,'dist_along_azimuth',[],@(x) isnumeric(x) && isscalar(x) && x>=0 && x<=360);
	addParameter(p,'uncertainty','se',@(x) ischar(validatestring(x,{'se','std','both'})));
	addParameter(p,'populate_categories',false,@(x) isscalar(x) && islogical(x))
	addParameter(p,'means_by_category',[],@(x) ischar(x));
	addParameter(p,'filter_by_category',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'filter_type','exclude',@(x) ischar(validatestring(x,{'exclude','include','mode'})));
	addParameter(p,'cat_grid',[],@(x) ischar(x));
	addParameter(p,'cat_values',[],@(x) ischar(x));
	addParameter(p,'file_name_prefix','Basins',@(x) ischar(x));


	parse(p,wdir,location_of_data_files,varargin{:});
	wdir=p.Results.wdir;
	location_of_data_files=p.Results.location_of_data_files;

	location_of_subbasins=p.Results.location_of_subbasins;
	include=p.Results.include;
	efv=p.Results.extra_field_values;
	new_concavity=p.Results.new_concavity;
	az=p.Results.dist_along_azimuth;
	uncertainty=p.Results.uncertainty;
	pc=p.Results.populate_categories;
	mbc=p.Results.means_by_category;
	fbc=p.Results.filter_by_category;
	ft=p.Results.filter_type;
	cgn=p.Results.cat_grid;
	cgv=p.Results.cat_values;
	fnp=p.Results.file_name_prefix;

	% Check required entries
	if fbc && ~strcmp(ft,'mode') && isempty(cgn) | isempty(cgv)
		error('For "include" or "exclude" filters, entries must be provided for both "cat_grid" and "cat_values"');
	elseif fbc && strcmp(ft,'mode') && isempty(cgn)
		error('For "mode" filter, entry must be provided for "cat_grid"');
	end

	% Switch for which basins to include
	switch include
	case 'all'
		FileList1=dir(fullfile(wdir,location_of_data_files,'*_Data.mat'));
		FileList2=dir(fullfile(wdir,location_of_data_files,location_of_subbasins, '*_DataSubset*.mat'));
		FileList=vertcat(FileList1,FileList2);
		num_files=numel(FileList);
	case 'bigonly'
		FileList=dir(fullfile(wdir,location_of_data_files,'*_Data.mat'));
		num_files=numel(FileList);
	case 'subdivided'
		AllFullFiles=dir(fullfile(wdir,location_of_data_files,'*_Data.mat'));
		num_basins=numel(AllFullFiles);
		basin_nums=zeros(num_basins,1);
		for jj=1:num_basins
			fileName=AllFullFiles(jj,1).name;
			basin_nums(jj)=sscanf(fileName,'%*6s %i'); %%%
		end

		FileCell=cell(num_basins,1);
		for kk=1:num_basins
			basin_num=basin_nums(kk);
			SearchAllString=fullfile(wdir,location_of_data_files,['*_' num2str(basin_num) '_Data.mat']);
			SearchSubString=fullfile(wdir,location_of_data_files,location_of_subbasins,['*_' num2str(basin_num) '_DataSubset*.mat']);

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

	if ~isempty(cgv)
		CGVT=readtable(fullfile(wdir,cgv));
		cgv=CGVT.Properties.VariableNames;
	end

	if ~isempty(mbc)
		MBCT=readtable(fullfile(wdir,mbc));
		mbc=MBCT.Properties.VariableNames;
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
		FileName=fullfile(FileList(ii,1).folder,FileList(ii,1).name);

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
				[mean_ksn,std_ksn,se_ksn]=ksn_convert(MSNc,new_concavity(jj));
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

			T.mean_el_f(ii,1)=nanmean(DEMoc.Z(F.Z));
			switch uncertainty
			case 'se'
				T.se_el_f(ii,1)=nanstd(DEMoc.Z(F.Z))/sqrt(sum(~isnan(DEMoc.Z(F.Z))));
			case 'std'
				T.std_el_f(ii,1)=nanstd(DEMoc.Z(F.Z));
			case 'both'
				T.se_el_f(ii,1)=nanstd(DEMoc.Z(F.Z))/sqrt(sum(~isnan(DEMoc.Z(F.Z))));
				T.std_el_f(ii,1)=nanstd(DEMoc.Z(F.Z));
			end

			T.mean_gradient_f(ii,1)=nanmean(Goc.Z(F.Z));
			switch uncertainty
			case 'se'
				T.se_gradient_f(ii,1)=nanstd(Goc.Z(F.Z))/sqrt(sum(~isnan(Goc.Z(F.Z))));
			case 'std'
				T.std_gradient_f(ii,1)=nanstd(Goc.Z(F.Z));
			case 'both'
				T.se_gradient_f(ii,1)=nanstd(Goc.Z(F.Z))/sqrt(sum(~isnan(Goc.Z(F.Z))));
				T.std_gradient_f(ii,1)=nanstd(Goc.Z(F.Z));
			end

			KSNG=GRIDobj(CG);
			KSNG.Z(:,:)=NaN;
			for kk=1:numel(MSNc)
				ix=coord2ind(CG,MSNc(kk).X,MSNc(kk).Y);
				KSNG.Z(ix)=MSNc(kk).ksn;
			end

			T.mean_ksn_f(ii,1)=nanmean(KSNG.Z(F.Z));
			switch uncertainty
			case 'se'
				T.se_ksn_f(ii,1)=nanstd(KSNG.Z(F.Z))/sqrt(sum(~isnan(KSNG.Z(F.Z))));
			case 'std'
				T.std_ksn_f(ii,1)=nanstd(KSNG.Z(F.Z));
			case 'both'
				T.se_ksn_f(ii,1)=nanstd(KSNG.Z(F.Z))/sqrt(sum(~isnan(KSNG.Z(F.Z))));
				T.std_ksn_f(ii,1)=nanstd(KSNG.Z(F.Z));
			end

			ag_grids=size(AGc,1);
			for kk=1:ag_grids
				agG=AGc{kk,1};
				mean_prop_name=['mean_' AGc{kk,2} '_f'];		
				T.(mean_prop_name)(ii,1)=nanmean(agG.Z(F.Z));

				switch uncertainty
				case 'se'
					se_prop_name=['se_' AGc{kk,2} '_f'];
					T.(se_prop_name)(ii,1)=nanstd(agG.Z(F.Z))/sqrt(sum(~isnan(agG.Z(F.Z))));
				case 'std'
					std_prop_name=['std_' AGc{kk,2} '_f'];
					T.(std_prop_name)(ii,1)=nanstd(agG.Z(F.Z));
				case 'both'
					se_prop_name=['se_' AGc{kk,2} '_f'];
					T.(se_prop_name)(ii,1)=nanstd(agG.Z(F.Z))/sqrt(sum(~isnan(agG.Z(F.Z))));
					std_prop_name=['std_' AGc{kk,2} '_f'];
					T.(std_prop_name)(ii,1)=nanstd(agG.Z(F.Z));
				end
			end

			if rlf_flag
				rlf_grids=size(rlf,1);
				for kk=1:rlf_grids
					rlfG=rlf{kk,1};
					mean_prop_name=['mean_rlf' num2str(rlf{kk,2}) '_f'];
					T.(mean_prop_name)(ii,1)=nanmean(rlfG.Z(F.Z));

					switch uncertainty
					case 'se'
						se_prop_name=['se_rlf' num2str(rlf{kk,2}) '_f'];
						T.(se_prop_name)(ii,1)=nanstd(rlfG.Z(F.Z))/sqrt(sum(~isnan(rlfG.Z(F.Z))));
					case 'std'
						std_prop_name=['std_rlf' num2str(rlf{kk,2}) '_f'];
						T.(std_prop_name)(ii,1)=nanstd(rlfG.Z(F.Z));
					case 'both'
						se_prop_name=['se_rlf' num2str(rlf{kk,2}) '_f'];
						T.(se_prop_name)(ii,1)=nanstd(rlfG.Z(F.Z))/sqrt(sum(~isnan(rlfG.Z(F.Z))));
						std_prop_name=['std_rlf' num2str(rlf{kk,2})];
						T.(std_prop_name)(ii,1)=nanstd(rlfG.Z(F.Z));
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
			error('No Categorical Grids were provided to ProcessRiverBasins so filtered values cannot be calculated');
		end

		% Check for the presence of extra fields provided at input
		if ~isempty(efv)
			efv=readtable(fullfile(wdir,efv));
			bnl=efv.(1);

			efn=efv.Properties.VariableNames;
			efn=efn{2:end};

			ix=find(bnl==RiverMouth(:,3));
			% Check to make sure a single entry exists for each basin number
			if ~isempty(ix) & numel(ix)==1
				num_efv=numel(efn);

				for kk=1:num_efv
					field_name=efn{kk};
					field_value=efv.(field_name)(ix);
					% Check to see if field value is a number or string
					if ischar(field_value)
						T.(field_name){ii,1}=field_value;
					elseif isnumeric(field_value)
						T.(field_name)(ii,1)=double(field_value);
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
								T.(cat_nameN)(ii,1)=nanmean(KSNG.Z(IDX.Z));
							elseif strcmp(dgOI,'gradient')
								load(FileName,'Goc');
								cat_nameN=matlab.lang.makeValidName(['mgrad_' cat_name]);
								T.(cat_nameN)(ii,1)=nanmean(Goc.Z(IDX.Z));
							elseif regexp(dgOI,regexptranslate('wildcard','rlf*'))
								rlfval=str2num(strrep(dgOI,'rlf',''));
								rlfix=find(cell2mat(rlf(:,2))==rlfval);
								if ~isempty(rlfix)
									Rg=rlf{rlfix,1};
									cat_nameN=matlab.lang.makeValidName(['mr' num2str(rlfval) '_' cat_name]);
									T.(cat_nameN)(ii,1)=nanmean(Rg.Z(IDX.Z));	
								end								
							else 
								try
									dgix=find(strcmp(AGc(:,2),dgOI));
									AGcOI=AGc{dgix,1};
									cat_nameN=matlab.lang.makeValidName(['m' AGc{dgix,2} '_' cat_name]);
									T.(cat_nameN)(ii,1)=nanmean(AGcOI.Z(IDX.Z));
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
			warning('One or more input for grid names to "means_by_category" was not recognized, table compiled without this entry')
		end
	end

	close(w1);

	save(fullfile(wdir,[fnp '_Table.mat']),'T','-v7.3');
	writetable(T,fullfile(wdir,[fnp '_Table.txt']));

end

function [mean_ksn,std_ksn,se_ksn]=ksn_convert(okm,new_ref_concavity)

	g=[okm.gradient];
	a=[okm.uparea];

	ksn_calc=g./a.^-new_ref_concavity);

	mean_ksn=mean(ksn_calc,'omitnan');
	std_ksn=std(ksn_calc,'omitnan');
	se_ksn=std_ksn/sqrt(numel(ksn_calc));

end
