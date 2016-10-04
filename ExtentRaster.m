function [OUT]=ExtentRaster(DEM,location_of_data_files)
	% Function takes outputs from 'ProcessRiverBasins' function and produces a single raster with the extent of selected
	% basins assigned the basin number, suitable for producing extent shapefiles in ArcGIS.
	%
	%Inputs:
	% DEM - GRIDobj of full extent of datasets
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
		load(FileName,'DEMc','RiverMouth');
		basin_num=RiverMouth(:,3);

		I=~isnan(DEMc.Z);
		[X,Y]=getcoordinates(DEMc);
		xmat=repmat(X,numel(Y),1);
		ymat=repmat(Y,1,numel(X));

		xix=xmat(I);
		yix=ymat(I);

		ix=coord2ind(DEM,xix,yix);
		val_list=ones(numel(ix),1).*basin_num;
		OUT.Z(ix)=val_list;
	end

	cd(current);