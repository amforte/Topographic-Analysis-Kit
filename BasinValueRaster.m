function [OUT]=BasinValueRaster(DEM,valueOI,location_of_data_files)
	% Function takes outputs from 'ProcessRiverBasins' function and produces a single GRIDobj with individual drainage
	% basins (as selected by 'ProcessRiverBasins' and 'SubDivideBigBasins') assinged various values
	%
	%Inputs:
	% DEM - GRIDobj of full extent of datasets
	% valueOI - value to assign to basins, acceptable inputs are:
	%	'ksn' - mean ksn value of basin
	%	'gradient' - mean gradient of basin
	%	'chir2' - R^2 value of chi-z fit (proxy for disequilibrium)
	%	'drainage_area' - drainage area in km2 of basin
    	%   'theta' - best fit concavity resultant from the topo toolbox chiplot function
	% location_of_data_files - full path of folder which contains the mat files from 'ProcessRiverBasins' as a string
	% 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Last Revised Spring 2016 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	OUT=GRIDobj(DEM);
	OUT=OUT-32768;

	current=pwd;
	cd(location_of_data_files);

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
			load(FileName,'DEMc','KSNc_stats');
			val=KSNc_stats(:,1);
		case 'gradient'
			load(FileName,'DEMc','Gc_stats');
			val=Gc_stats(:,1);
		case 'chir2'
			load(FileName,'DEMc','Chic');
			val=Chic.R2;
		case 'drainage_area'
			load(FileName,'DEMc','drainage_area');
			val=drainage_area;
		case 'theta'
			load(FileName,'DEMc','Sc','Ac');
			C=chiplot(Sc,DEMc,Ac,'a0',1,'plot',false);
			val=C.mn;
		end

		I=~isnan(DEMc.Z);
		[X,Y]=getcoordinates(DEMc);
		xmat=repmat(X,numel(Y),1);
		ymat=repmat(Y,1,numel(X));

		xix=xmat(I);
		yix=ymat(I);

		ix=coord2ind(DEM,xix,yix);
		val_list=ones(numel(ix),1).*val;
		OUT.Z(ix)=val_list;
	end

	cd(current);
