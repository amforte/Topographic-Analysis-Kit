function [OUT]=Basin2Raster(DEM,valueOI,location_of_data_files,varargin)
	%
	% Usage:
	%	[RASTER]=Basin2Raster(DEM,valueOI,location_of_data_files);
	%	[RASTER]=Basin2Raster(DEM,valueOI,location_of_data_files,'name',value,...);
	%
	% Description:
	% 	Function takes outputs from 'ProcessRiverBasins' function and produces a single GRIDobj with individual drainage
	% 	basins (as selected by 'ProcessRiverBasins' and 'SubDivideBigBasins') assinged various values
	%
	% Required Inputs:
	%	DEM - GRIDobj of full extent of datasets
	% 	valueOI - value to assign to basins, acceptable inputs are:
	%		'ksn' - mean ksn value of basin
	%		'gradient' - mean gradient of basin
	%		'elevation' - mean elevation of basin
	%		'relief' - mean relief of basin (must specify the radius of interest with the 'relief_radius' parameter)
	%		'chir2' - R^2 value of chi-z fit (proxy for disequilibrium)
	%		'drainage_area' - drainage area in km2 of basin
	%		'hypsometric_integral' - hypsometric integral of basin
	%		'id' - basin ID number (i.e third column RiverMouth output)
	%   	'theta' - best fit concavity resultant from the topo toolbox chiplot function 
	%		'NAME' - where name is the name provided for an extra grid (i.e. entry to second column of 'add_grid' or entry to 
	%			third column of 'add_cat_grid'), value input will be mean for additional grid names or mode for additional 
	%			categorical grid names
	%	location_of_data_files - full path of folder which contains the mat files from 'ProcessRiverBasins' as a string
	%
	% Optional Inputs:
	%	file_name_prefix ['basins'] - prefix for outputs, will automatically append the type of output, i.e. 'ksn', 'elevation', etc
	%	location_of_subbasins ['SubBasins'] - name of folder that contains subbasins of interest (if you created subbasins using
	%		"SubDivideBigBasins"), expected to be within the main Basin folder provided with "location_of_data_files".
	%	method ['subdivided'] - method used for subdividing watersheds. If you used 'ProcessRiversBasins' and then
	%		'SubDivideBigBasins' or if you only used 'ProcessRiverBasins' but did not pick any nested catchments, i.e.
	%		none of the river mouths supplied to 'ProcessRiverBasins' were within the catchment boundaries of other 
	%		watersheds for which you provided river mouths, then you should use use 'subdivided' which is the default
	%		so you do not need to specify a value for this property. If you picked nested catchments manually and then
	%		ran 'ProcessRiverBasins' you should use 'nested'.
	%	relief_radius [2500] - relief radius to use if 'valueOI' is set to 'relief'
	%		 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 01/27/19 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'Basin2Raster';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'valueOI',@(x) ischar(x));
	addRequired(p,'location_of_data_files',@(x) isdir(x));

	addParameter(p,'location_of_subbasins','SubBasins',@(x) ischar(x) || isempty(x));
	addParameter(p,'file_name_prefix','basins',@(x) ischar(x));
	addParameter(p,'method','subdivided',@(x) ischar(validatestring(x,{'subdivided','nested'})));
	addParameter(p,'relief_radius',2500,@(x) isscalar(x) && isnumeric(x));

	parse(p,DEM,valueOI,location_of_data_files,varargin{:});
	DEM=p.Results.DEM;
	valueOI=p.Results.valueOI;
	location_of_data_files=p.Results.location_of_data_files;

	location_of_subbasins=p.Results.location_of_subbasins;
	file_name_prefix=p.Results.file_name_prefix;
	method=p.Results.method;
	rr=p.Results.relief_radius;

	OUT=GRIDobj(DEM);
	OUT=OUT-32768;

	[head_dir,~,~]=fileparts(location_of_data_files);
	file_name_prefix=[head_dir filesep file_name_prefix];

	switch method
	case 'subdivided'
		%% Build File List
		% Get Basin Numbers
		AllFullFiles=dir([location_of_data_files filesep '*_Data.mat']);
		num_basins=numel(AllFullFiles);
		basin_nums=zeros(num_basins,1);
		for jj=1:num_basins
			FileName=AllFullFiles(jj,1).name;
			basin_nums(jj)=sscanf(FileName,'%*6s %i'); %%%
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

		w1=waitbar(0,'Building raster...');
		for ii=1:num_files;
			FileName=[FileList(ii,1).folder filesep FileList(ii,1).name];
			switch valueOI
			case 'ksn'
				load(FileName,'DEMoc','KSNc_stats');
				val=KSNc_stats(:,1);
			case 'gradient'
				load(FileName,'DEMoc','Gc_stats');
				val=Gc_stats(:,1);
			case 'elevation'
				load(FileName,'DEMoc','Zc_stats');
				val=Zc_stats(:,1);
			case 'chir2'
				load(FileName,'DEMoc','Sc','Ac','theta_ref');
				c=chiplot(Sc,DEMoc,Ac,'a0',1,'mn',theta_ref,'plot',false);
				val=c.R2;
			case 'drainage_area'
				load(FileName,'DEMoc','drainage_area');
				val=drainage_area;
			case 'hypsometric_integral'
				load(FileName,'DEMoc','hyps');
				val=abs(trapz((hyps(:,2)-min(hyps(:,2)))/(max(hyps(:,2))-min(hyps(:,2))),hyps(:,1)/100));
			case 'id'
				load(FileName,'DEMoc','RiverMouth');
				val=RiverMouth(:,3);
			case 'theta'
				load(FileName,'DEMoc','Chic');
				val=Chic.mn;
			case 'relief'
				VarList=whos('-file',FileName);
				RLFInd=find(strcmp(cellstr(char(VarList.name)),'rlf'));
				if isempty(RLFInd)
					error('Relief does appear to have been calculated for these basins')
				end
				load(FileName,'DEMoc','rlf','rlf_stats');
				ix=find(rlf(:,2)==rr);
				if isempty(ix)
					error('Input relief radius was not found in relief outputs, please check to make sure relief radius is correct')
				end
				val=rlf_stats(ix,1);
			otherwise
				VarList=whos('-file',FileName);
				AGInd=find(strcmp(cellstr(char(VarList.name)),'AGc'));
				ACGInd=find(strcmp(cellstr(char(VarList.name)),'ACGc'));
				if ~isempty(AGInd)
					load(FileName,'DEMoc','AGc');
					AGInd=true;
				end

				if ~isempty(ACGInd)
					load(FileName,'DEMoc','ACGc');
					ACGInd=true;
				end

				if AGInd & any(strcmp(AGc(:,2),valueOI))
					Nix=find(strcmp(AGc(:,2),valueOI));
					load(FileName,'AGc_stats');
					val=AGc_stats(Nix,1);
				elseif ACGInd & any(strcmp(ACGc(:,3),valueOI))
					Nix=find(trcmp(ACGc(:,3),valueOI));
					load(FileName,'ACGc_stats');
					val=ACGc_stats(Nix,1);
				else
					error('Name provided to "valueOI" does not match name in either additional grids or additional categorical grids');
				end
			end

			I=~isnan(DEMoc.Z);
			[X,Y]=getcoordinates(DEMoc);
			xmat=repmat(X,numel(Y),1);
			ymat=repmat(Y,1,numel(X));

			xix=xmat(I);
			yix=ymat(I);

			ix=coord2ind(DEM,xix,yix);
			val_list=ones(numel(ix),1).*val;
			OUT.Z(ix)=val_list;
			waitbar(ii/num_files);
		end
		close(w1);

	case 'nested'
		% Build list of indices
		AllFiles=dir([location_of_data_files filesep '*_Data.mat']);
		num_basins=numel(AllFiles);

		ix_cell=cell(num_basins,1);
		basin_list=zeros(num_basins,1);
		fileCell=cell(num_basins,1);
		for jj=1:num_basins
			FileName=AllFiles(jj,1).name;
			FileCell{jj}=FileName;

			load(FileName,'DEMoc');
			[x,y]=getcoordinates(DEMoc);
			xg=repmat(x,numel(y),1);
			yg=repmat(y,1,numel(x));
			xl=xg(~isnan(DEMoc.Z));
			yl=yg(~isnan(DEMoc.Z));
			ix=coord2ind(DEM,xl,yl);

			ix_cell{jj}=ix;

			basin_list(jj,1)=numel(ix);
		end

		% Sort basin size list in descending order
		[~,six]=sort(basin_list,'descend');
		% Apply sorting index to fileCell and ix_cell
		FileCell=FileCell(six);
		ix_cell=ix_cell(six);

		w1=waitbar(0,'Building raster...');
		for ii=1:num_basins
			FileName=FileCell{ii};

			switch valueOI
			case 'ksn'
				load(FileName,'DEMoc','KSNc_stats');
				val=KSNc_stats(:,1);
			case 'gradient'
				load(FileName,'DEMoc','Gc_stats');
				val=Gc_stats(:,1);
			case 'elevation'
				load(FileName,'DEMoc','Zc_stats');
				val=Zc_stats(:,1);
			case 'chir2'
				load(FileName,'DEMoc','Chic');
				val=Chic.R2;
			case 'drainage_area'
				load(FileName,'DEMoc','drainage_area');
				val=drainage_area;
			case 'hypsometric_integral'
				load(FileName,'DEMoc','hyps');
				val=abs(trapz((hyps(:,2)-min(hyps(:,2)))/(max(hyps(:,2))-min(hyps(:,2))),hyps(:,1)/100));
			case 'id'
				load(FileName,'DEMoc','RiverMouth');
				val=RiverMouth(:,3);
			case 'theta'
				load(FileName,'DEMoc','Chic');
				val=Chic.mn;
			case 'relief'
				VarList=whos('-file',FileName);
				RLFInd=find(strcmp(cellstr(char(VarList.name)),'rlf'));
				if isempty(RLFInd)
					error('Relief does appear to have been calculated for these basins')
				end
				load(FileName,'DEMoc','rlf','rlf_stats');
				ix=find(rlf(:,2)==rr);
				if isempty(ix)
					error('Input relief radius was not found in relief outputs, please check to make sure relief radius is correct')
				end
				val=rlf_stats(ix,1);
			otherwise
				VarList=whos('-file',FileName);
				AGInd=find(strcmp(cellstr(char(VarList.name)),'AGc'));
				ACGInd=find(strcmp(cellstr(char(VarList.name)),'ACGc'));
				if ~isempty(AGInd)
					load(FileName,'DEMoc','AGc');
					AGInd=true;
				end

				if ~isempty(ACGInd)
					load(FileName,'DEMoc','ACGc');
					ACGInd=true;
				end

				if AGInd & any(strcmp(AGc(:,2),valueOI))
					Nix=find(strcmp(AGc(:,2),valueOI));
					load(FileName,'AGc_stats');
					val=AGc_stats(Nix,1);
				elseif ACGInd & any(strcmp(ACGc(:,3),valueOI))
					Nix=find(trcmp(ACGc(:,3),valueOI));
					load(FileName,'ACGc_stats');
					val=ACGc_stats(Nix,1);
				else
					error('Name provided to "valueOI" does not match name in either additional grids or additional categorical grids');
				end
			end

			I=~isnan(DEMoc.Z);
			[X,Y]=getcoordinates(DEMoc);
			xmat=repmat(X,numel(Y),1);
			ymat=repmat(Y,1,numel(X));

			ix=ix_cell{ii};

			val_list=ones(numel(ix),1).*val;
			OUT.Z(ix)=val_list;
			waitbar(ii/num_files);
		end
		close(w1);
	end

	switch valueOI
	case 'ksn'
		out_file=[file_name_prefix '_ksn.txt'];
		GRIDobj2ascii(OUT,out_file);
	case 'gradient'
		out_file=[file_name_prefix '_gradient.txt'];
		GRIDobj2ascii(OUT,out_file);
	case 'elevation'
		out_file=[file_name_prefix '_elevation.txt'];
		GRIDobj2ascii(OUT,out_file);
	case 'chir2'
		out_file=[file_name_prefix '_chiR2.txt'];
		GRIDobj2ascii(OUT,out_file);
	case 'drainage_area'
		out_file=[file_name_prefix '_drainarea.txt'];
		GRIDobj2ascii(OUT,out_file);
	case 'hypsometric_integral'
		out_file=[file_name_prefix '_hyps_int.txt'];
		GRIDobj2ascii(OUT,out_file);
	case 'id'
		out_file=[file_name_prefix '_id.txt'];
		GRIDobj2ascii(OUT,out_file);
	case 'theta'
		out_file=[file_name_prefix '_theta.txt'];
		GRIDobj2ascii(OUT,out_file);
	case 'relief'
		out_file=[file_name_prefix '_relief.txt'];
		GRIDobj2ascii(OUT,out_file);		
	otherwise
		out_file=[file_name_prefix '_' valueOI '.txt'];
		GRIDobj2ascii(OUT,out_file);
	end

