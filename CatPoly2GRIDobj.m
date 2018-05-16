function [OUT,look_table]=CatPoly2GRIDobj(DEM,poly_shape,field)
	% Function to convert a categorical polygon shape file (e.g. a digitzed geologic map) to a GRIDobj. Can
	%	be useful for use in 'ProcessRiverBasins'
	%
	% Required Inputs:
	% 	DEM - DEM that you want the output to match 
	%	poly_shape - name or path to shapefile containing the categorical data
	%	field - field name of categorical data within the shapefile 
	%
	% Outputs:
	%	OUT - GRIDobj of the same size as DEM where values correspond to categorical data
	%		as defined in the look_table
	%	look_table - nx2 table with columns Numbers and Categories that serves as a lookup table 
	%		to convert between the numbers and the original categories.
	%
	% Example:
	%	[GEO,geo_table]=CatPoly2GRIDobj(DEM,'geologic_map.shp','rtype');
	% 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Last Revised Spring 2018 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Read in shape and covert to a table
	PS=shaperead(poly_shape);
	TS=struct2table(PS);

	% Separate out the field of interest
	Foi=TS.(field);

	% Generate unique liss and the output lookup table
	Categories=unique(Foi);
	Numbers=[1:numel(Categories)]; Numbers=Numbers';
	look_table=table(Numbers,Categories);

	% Replace categorical with number
	for ii=1:numel(PS)
		Eoi=PS(ii,1).(field);
		ix=find(strcmp(Categories,Eoi));
		PS(ii,1).replace_number=Numbers(ix);
	end

	% Run polygon2GRIDobj
	[OUT]=polygon2GRIDobj(DEM,PS,'replace_number');
end