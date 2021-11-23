function RetroAddGrids2Basins(location_of_data_files,DEM,behavior,varargin)
	%
	% Usage:
	%		RetroAddGrids2Basins(locations_of_data_files,DEM,behavior,'name','value',...);
	%
	% Description:
	%		Function operates on an existing set of basins and/or subbasins produced by 'ProcessRiverBasins' 
	%		and/or 'SubdivideBigBasins' and adds additional grids or additional categorical grids (or both) 
	%		to the basin files. You can control whether this action will overwrite exisitng or append
	%		to existing additional grids or additional categorical grids within those basins.		
	%		
	% Required Inputs:
	% 		location_of_data_files - full path of folder which contains the mat files from 'ProcessRiverBasins'	
	%		DEM - GRIDobj of original DEM provided to ProcessRiverBasins
	%		behavior - parameter to specify how to deal with if existing additional grids or additional categorical
	%			grids are found in the specified basins. Valid inputs are 'overwrite' or 'append'. The 'overwrite' option
	%			will ignore any exisitng additional grids or additional categorical grids and overwrite what is there, 
	%			where as append will add the new inputs specified here to the relevant arrays in the basin files. If you 
	%			specify append and there are no existing grids, this will not generate an error, it will just have no effect.
	%
	% Optional Inputs:
	%		add_grids [] - An additional grids input as you would provide if running ProcessRiverBasins, 
	%			see the details in the 'add_grids' parameter of 'ProcessRiverBasins' for additional info. Either 
	%			'add_grids' or 'add_cat_grids' must be provided (otherwise there is no reason to run this function.)
	%		add_cat_grids [] - An additional categorical grids input as you would provide if running 'ProcessRiverBasins',
	%			see the details in the 'add_cat_grids' parameter of 'ProcessRiverBasins' for additional info. Either 
	%			'add_grids' or 'add_cat_grids' must be provided (otherwise there is no reason to run this function.)
	%		include ['all'] - parameter to specify which basins you wish to add data to. The default 'all' 
	%			will include all basin mat files in the folder you specify. Providing 'subdivided' will check to see 
	%			if a given main basin was subdivided using 'SubdivideBigBasins' and then only process the subdivided 
	%			versions of that basin (i.e. the original main basin for those subbasins will not have the additional grids). 
	%			Providing 'bigonly' will only process the original basins produced by 'ProcessRiverBasins' even if 
	%			'SubDivideBigBasins' was run. If 'SubDivideBigBasins' was never run, result of 'all' and 'bigonly' will be the same.		
	%		location_of_subbasins ['SubBasins'] - name of folder that contains subbasins of interest (if you created 
	%			subbasins using"SubDivideBigBasins"), expected to be within the main Basin folder provided with "location_of_data_files".
	%			Note that if you do not provide the correct directory name for the location of the subbasins, subbasins
	%			will not be modified regardless of your choicefor the "include" parameter.
	%		resample_method ['nearest'] - method to use in the resample function on additional grids (if required). 
	%			Acceptable inputs are 'nearest', 'bilinear', or 'bicubic'. Method 'nearest' is appropriate if you do 
	%			not want the resampling to interpolate between values (e.g. if an additinal grid has specific values 
	%			that correlate to a property like rock type) and either 'bilinear' or 'bicubic' is appropriate if you 
	%			want smooth variations between nodes. If you are providing additional categorical grids, resampling will be
	%			performed if necessary, but will use a 'nearest' resampling method regardless of the choice here.
	% 
	% Examples:
	% 		RetroAddGrids2Basins('basins',DEM,'append','add_grids',{PRE,'precip'});
	% 		RetroAddGrids2Basins('basins',DEM,'append','add_grids',{PRE,'precip'},'include','subdivided','location_of_subbasins','my_subs');	
	%		RetroAddGrids2Basins('basins',DEM,'overwrite','add_cat_grids',{GEO,geo_table,'geology'});
	%	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 05/16/20 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'RetroAddGrids2Basins';
	addRequired(p,'location_of_data_files',@(x) isdir(x));
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'behavior',@(x) ischar(validatestring(x,{'overwrite','append'})));

	addParameter(p,'location_of_subbasins','SubBasins',@(x) ischar(x) || isempty(x));
	addParameter(p,'include','all',@(x) ischar(validatestring(x,{'all','subdivided','bigonly'})));
	addParameter(p,'add_grids',[],@(x) isa(x,'cell') && size(x,2)==2 || isempty(x));
	addParameter(p,'add_cat_grids',[],@(x) isa(x,'cell') && size(x,2)==3 || isempty(x));
	addParameter(p,'resample_method','nearest',@(x) ischar(validatestring(x,{'nearest','bilinear','bicubic'})));

	parse(p,location_of_data_files,DEM,behavior,varargin{:});
	location_of_data_files=p.Results.location_of_data_files;
	DEM=p.Results.DEM;
	behavior=p.Results.behavior;

	location_of_subbasins=p.Results.location_of_subbasins;
	include=p.Results.include;
	AG=p.Results.add_grids;
	ACG=p.Results.add_cat_grids;
	resample_method=p.Results.resample_method;

	% Check to make sure that something meaningful has been provided
	if isempty(AG) & isempty(ACG)
		error('You must provide an input to either "add_grids" or "add_cat_grids" or there is nothing for this function to do')
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


	% Perform check on dimensions and cellsize of additional grids and resample if necessary
	if ~isempty(AG)
		num_grids=size(AG,1);
		for jj=1:num_grids
			AGoi=AG{jj,1};
			if ~validatealignment(AGoi,DEM);
				disp(['Resampling ' AG{jj,2} ' GRIDobj to be the same resolution and dimensions as the input DEM by the ' resample_method ' method']);
				AG{jj,1}=resample(AGoi,DEM,resample_method);
			end
		end
	end


	% Perform check on dimensions and cellsize of additional grids and resample if necessary
	if ~isempty(ACG)
		num_grids=size(ACG,1);
		for jj=1:num_grids
			ACGoi=ACG{jj,1};
			if ~validatealignment(ACGoi,DEM);
				disp(['Resampling ' ACG{jj,3} ' GRIDobj to be the same resolution and dimensions as the input DEM by the nearest method']);
				ACG{jj,1}=resample(ACGoi,DEM,'nearest');
			end
		end
	end

	% Start main loop
	w1=waitbar(0,'Adding Grids to Basins...');
	for ii=1:num_files
		% Load in DEM
		FileName=[FileList(ii,1).folder filesep FileList(ii,1).name];
		load(FileName,'DEMoc');
		% Build mask
		[x,y]=getcoordinates(DEMoc);
		[X,Y]=meshgrid(x,y);
		ix=isnan(DEMoc.Z(:));
		xx=X(ix); yy=Y(ix);
		gix=coord2ind(DEM,xx,yy);
		I=GRIDobj(DEM,'logical');
		I.Z(gix)=true;

		% Operate on any additional grids
		if ~isempty(AG)
			switch behavior
			case 'overwrite'
				num_grids=size(AG,1);
				AGc=cell(size(AG));
				for jj=1:num_grids
					AGcOI=crop(AG{jj,1},I,nan);
					AGc{jj,1}=AGcOI;
					AGc{jj,2}=AG{jj,2};
					mean_AGc=mean(AGcOI.Z(:),'omitnan');
					min_AGc=min(AGcOI.Z(:),[],'omitnan');
					max_AGc=max(AGcOI.Z(:),[],'omitnan');
					std_AGc=std(AGcOI.Z(:),'omitnan');
					se_AGc=std_AGc/sqrt(sum(~isnan(AGcOI.Z(:))));
					AGc_stats(jj,:)=[mean_AGc se_AGc std_AGc min_AGc max_AGc];
				end
				save(FileName,'AGc','AGc_stats','-append');
			case 'append'
				vI=who('-file',FileName);
				if any(ismember(vI,'AGc'))
					load(FileName,'AGc','AGc_stats');
					num_existing=size(AGc,1);
					pos_vec=num_existing+1:num_grids+num_existing;
					num_grids=size(AG,1);
					for jj=1:num_grids
						AGcOI=crop(AG{jj,1},I,nan);
						AGc{pos_vec(jj),1}=AGcOI;
						AGc{pos_vec(jj),2}=AG{jj,2};
						mean_AGc=mean(AGcOI.Z(:),'omitnan');
						min_AGc=min(AGcOI.Z(:),[],'omitnan');
						max_AGc=max(AGcOI.Z(:),[],'omitnan');
						std_AGc=std(AGcOI.Z(:),'omitnan');
						se_AGc=std_AGc/sqrt(sum(~isnan(AGcOI.Z(:))));
						AGc_stats(pos_vec(jj),:)=[mean_AGc se_AGc std_AGc min_AGc max_AGc];
					end
					save(FileName,'AGc','AGc_stats','-append');						

				else
					num_grids=size(AG,1);
					AGc=cell(size(AG));
					for jj=1:num_grids
						AGcOI=crop(AG{jj,1},I,nan);
						AGc{jj,1}=AGcOI;
						AGc{jj,2}=AG{jj,2};
						mean_AGc=mean(AGcOI.Z(:),'omitnan');
						min_AGc=min(AGcOI.Z(:),[],'omitnan');
						max_AGc=max(AGcOI.Z(:),[],'omitnan');
						std_AGc=std(AGcOI.Z(:),'omitnan');
						se_AGc=std_AGc/sqrt(sum(~isnan(AGcOI.Z(:))));
						AGc_stats(jj,:)=[mean_AGc se_AGc std_AGc min_AGc max_AGc];
					end
					save(FileName,'AGc','AGc_stats','-append');					
				end
			end				
		end

		% Operate on any additional categorical grids
		if ~isempty(ACG)
			switch behavior
			case 'overwrite'
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
			case 'append'
				vI=who('-file',FileName);
				if any(ismember(vI,'ACGc'))
					num_existing=size(ACGc,1);
					pos_vec=num_existing+1:num_grids+num_existing;
					num_grids=size(ACG,1);
					for jj=1:num_grids
						ACGcOI=crop(ACG{jj,1},I,nan);
						ACGc{pos_vec(jj),1}=ACGcOI;
						ACGc{pos_vec(jj),3}=ACG{jj,3};
						edg=ACG{jj,2}.Numbers;
						edg=edg+0.5;
						edg=vertcat(0.5,edg);
						[N,~]=histcounts(ACGcOI.Z(:),edg);
						T=ACG{jj,2};
						T.Counts=N';
						ACGc{pos_vec(jj),2}=T;
						ACGc_stats(pos_vec(jj),1)=[mode(ACGcOI.Z(:))];
					end
					save(FileName,'ACGc','ACGc_stats','-append');
				else
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
			end
		end
		waitbar(ii/num_files);
	end
	close(w1);
end
