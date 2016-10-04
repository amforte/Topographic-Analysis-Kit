function [DEM,FD,A,S]=MakeStreams(location_of_dem,threshold_area,save_output,varargin)
	% Function takes a dem and outputs the necessary base datasets for use in other TopoToolbox functions
	%
	%Inputs:
	% location_of_dem - full path of dem file as either an ascii text file (recommended) or geotiff
	% threshold_area - minimum accumulation area to define streams in meters squared
	%
	%Outputs:
	% DEM - GRIDobj of the DEM
	% FD - FLOWobj from the supplied DEM
	% A - Flow accumulation grid (GRIDobj)
	% S - STREAMobj derived from the DEM
	% save_output - set to true if you wish to save a mat file containing the DEM, FD, A, and S and a shapefile of 
	%	the stream network produced, set to false if you do not wish to save the output as a mat file. 
	%	If you set the value to true, you must follow this with a file name as a string.
	% 
	%
	% examples [DEM,FD,A,S]=MakeStreams('/Users/forte/topo/dem.tif',1e6,false);
	%		   [DEM,FD,A,S]=MakeStreams('/Users/forte/topo/dem.tif',1e6,true,'AreaFiles');	
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Last Revised Spring 2016 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	current=pwd;

	% Check for filename
	switch save_output
	case true
		if length(varargin)~=1
			error('Error: You must provide a filename as a string if you wish to save outputs');
		else
			;
		end
	end


	disp('Loading and processing DEM')
	DEM=GRIDobj(location_of_dem);

	% Clippin gout NaNs, max and min elev set arbitrarily large but below (or above) the (-)32,768 internal NaN value.
	max_elev=10000;
	min_elev=-1000;
	IDX=DEM<max_elev & DEM>min_elev;
	DEM=crop(DEM,IDX,nan);
	cd(current);

	disp('Calculating Flow Direction')
	FD=FLOWobj(DEM,'preprocess','carve','verbose',true);

	disp('Calculating Flow Accumulation')
	A=flowacc(FD);

	disp('Extracting total stream network')
	DEM_res=DEM.cellsize;
	min_area=floor(threshold_area/(DEM_res*DEM_res));
	isstream=A>min_area;
	S=STREAMobj(FD,isstream);


	switch save_output
	case true
        %Save output if switch thrown
        fileNameBase=varargin{1};
        MatFileName=[fileNameBase '.mat'];
        ShpFileName=[fileNameBase '.shp'];
		save(MatFileName,'DEM','FD','A','S','-v7.3');
		MS=STREAMobj2mapstruct(S);
		shapewrite(MS,ShpFileName);
	case false
		;
	end


end