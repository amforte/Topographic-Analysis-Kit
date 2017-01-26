function [OUT]=Basin2Raster(DEM,valueOI,location_of_data_files,varargin)
	% Function takes outputs from 'ProcessRiverBasins' function and produces a single GRIDobj with individual drainage
	% basins (as selected by 'ProcessRiverBasins' and 'SubDivideBigBasins') assigned various values. If you would rather
	% create a shapefile, use 'Basin2Shape'.
	%
	% Required Inputs:
	%	DEM - GRIDobj of full extent of datasets
	% 	valueOI - value to assign to basins, acceptable inputs are:
	%		'ksn' - mean ksn value of basin
	%		'gradient' - mean gradient of basin
	%		'elevation' - mean elevation of basin
	%		'chir2' - R^2 value of chi-z fit (proxy for disequilibrium)
	%		'drainage_area' - drainage area in km2 of basin
	%		'id' - basin ID number (i.e third column RiverMouth output)
	%   	'theta' - best fit concavity resultant from the topo toolbox chiplot function (will not work if method for 
	%			ProcessRiverBasin was 'segment')
	%	location_of_data_files - full path of folder which contains the mat files from 'ProcessRiverBasins' as a string
	%
	% Optional Inputs:
	%	method ['subdivided'] - method used for subdividing watersheds. If you used 'ProcessRiversBasins' and then
	%		'SubDivideBigBasins' or if you only used 'ProcessRiverBasins' but did not pick any nested catchments, i.e.
	%		none of the river mouths supplied to 'ProcessRiverBasins' were within the catchment boundaries of other 
	%		watersheds for which you provided river mouths, then you should use use 'subdivided' which is the default
	%		so you do not need to specify a value for this property. If you picked nested catchments manually and then
	%		ran 'ProcessRiverBasins' you should use 'nested'.
	%		 
	% 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Last Revised Spring 2016 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'Basin2Raster';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'valueOI',@(x) ischar(validatestring(x,{'ksn','gradient','chir2','elevation','drainage_area','theta','id'})));
	addRequired(p,'location_of_data_files',@(x) ischar(x));

	addParamValue(p,'method','subdivided',@(x) ischar(validatestring(x,{'subdivided','nested'})));

	parse(p,DEM,valueOI,location_of_data_files,varargin{:});
	DEM=p.Results.DEM;
	valueOI=p.Results.valueOI;
	location_of_data_files=p.Results.location_of_data_files;

	method=p.Results.method;

	OUT=GRIDobj(DEM);
	OUT=OUT-32768;

	current=pwd;
	cd(location_of_data_files);

	switch method
	case 'subdivided'
		%% Build File List
		% Get Basin Numbers
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
		fileList=vertcat(FileCell{:});
		num_files=numel(fileList);

		for ii=1:num_files;
			disp(['Working on ' num2str(ii) ' of ' num2str(num_files)]);
			FileName=fileList(ii,1).name;
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
			case 'id'
				load(FileName,'DEMoc','RiverMouth');
				val=RiverMouth(:,3);
			case 'theta'
				load(FileName,'DEMcc','Sc','Ac');
				C=chiplot(Sc,DEMcc,Ac,'a0',1,'plot',false);
				val=C.mn;
			end

			if strcmp(valueOI,'theta')
				I=~isnan(DEMcc.Z);
				[X,Y]=getcoordinates(DEMcc);
				xmat=repmat(X,numel(Y),1);
				ymat=repmat(Y,1,numel(X));
			else
				I=~isnan(DEMoc.Z);
				[X,Y]=getcoordinates(DEMoc);
				xmat=repmat(X,numel(Y),1);
				ymat=repmat(Y,1,numel(X));
			end

			xix=xmat(I);
			yix=ymat(I);

			ix=coord2ind(DEM,xix,yix);
			val_list=ones(numel(ix),1).*val;
			OUT.Z(ix)=val_list;
		end

	case 'nested'
		% Build list of indices
		allfiles=dir('*_Data.mat');
		num_basins=numel(allfiles);

		ix_cell=cell(num_basins,1);
		basin_list=zeros(num_basins,1);
		fileCell=cell(num_basins,1);
		for jj=1:num_basins
			fileName=allfiles(jj,1).name;
			fileCell{jj}=fileName;

			load(fileName,'DEMoc');
			[x,y]=getcoordinates(DEMoc);
			xg=repmat(x,numel(y),1);
			yg=repmat(y,1,numel(x));
			xl=xg(~isnan(DEMc.Z));
			yl=yg(~isnan(DEMc.Z));
			ix=coord2ind(DEM,xl,yl);

			ix_cell{jj}=ix;

			basin_list(jj,1)=numel(ix);
		end

		% Sort basin size list in descending order
		[~,six]=sort(basin_list,'descend');
		% Apply sorting index to fileCell and ix_cell
		fileCell=fileCell(six);
		ix_cell=ix_cell(six);

		for ii=1:num_basins
			FileName=fileCell{ii};
			disp(['Working on ' num2str(ii) ' of ' num2str(num_files)]);

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
			case 'id'
				load(FileName,'DEMoc','RiverMouth');
				val=RiverMouth(:,3);
			case 'theta'
				load(FileName,'DEMcc','Sc','Ac');
				C=chiplot(Sc,DEMc,Ac,'a0',1,'plot',false);
				val=C.mn;
			end

			if strcmp(valueOI,'theta')
				I=~isnan(DEMcc.Z);
				[X,Y]=getcoordinates(DEMcc);
				xmat=repmat(X,numel(Y),1);
				ymat=repmat(Y,1,numel(X));
			else
				I=~isnan(DEMoc.Z);
				[X,Y]=getcoordinates(DEMoc);
				xmat=repmat(X,numel(Y),1);
				ymat=repmat(Y,1,numel(X));
			end

			ix=ix_cell{ii};

			val_list=ones(numel(ix),1).*val;
			OUT.Z(ix)=val_list;
		end
	end

	cd(current);
