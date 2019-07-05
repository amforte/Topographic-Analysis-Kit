function [OUT,look_table]=CatPoly2GRIDobj(DEM,poly_shape,field,varargin)
	%
	% Usage:
	%	[GRIDobj,look_up_table]=CatPoly2GRIDobj(DEM,poly_shape,field);
	%
	% Description:
	% 	Function to convert a categorical polygon shape file (e.g. a digitzed geologic map) to a GRIDobj. Can
	%	be useful for use in 'ProcessRiverBasins'
	%
	% Required Inputs:
	% 	DEM - DEM that you want the output to match 
	%	poly_shape - name or path to shapefile containing the categorical data
	%	field - field name of categorical data within the shapefile 
	%
	% Optional Input:
	%	table_in - Optional input to manually provide a table to use as the look_table. This table must contain
	%		two columns where the first column is n x 1 array of scalar integers and the second column is a n x 1
	%		cell array containing the category names. These categories must match entries in the 'field' of the
	%		provided 'poly_shape' or all of the values in the 'OUT' will be 0.
	%		 
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
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 07/05/19 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'CatPoly2GRIDobj';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'poly_shape',@(x) ischar(x));
	addRequired(p,'field',@(x) ischar(x));

	addParameter(p,'table_in',[],@(x) istable(x));

	parse(p,DEM,poly_shape,field,varargin{:});
	DEM=p.Results.DEM;
	poly_shape=p.Results.poly_shape;
	field=p.Results.field;

	table_in=p.Results.table_in;

	% Validate table if provided
	if ~isempty(table_in)
		if size(table_in,2)~=2
			error('Provided "table_in" must have two columns');
		elseif table_in.(1)(1)~=0
			error('First entry of first column of "table_in" must equal 0');
		elseif ~iscell(table_in.(2))
			error('Second column of "table_in" must be a cell array');
		elseif ~ischar(table_in.(2){1})
			error('Entries in second column of "table_in" must be characters');
		end
	end
			

	% Read in shape and covert to a table
	PS=shaperead(poly_shape);
	TS=struct2table(PS);

	% Separate out the field of interest
	Foi=TS.(field);

	if isempty(table_in)
		% Generate unique list and the output lookup table
		Categories=unique(Foi);
		Numbers=[1:numel(Categories)]; Numbers=Numbers';
		% Add an 'undef' category to deal with zeros that appear because of read errors
		Numbers=vertcat(0,Numbers);
		Categories=vertcat('undef',Categories);
		look_table=table(Numbers,Categories);
	else 
		% Load provided table into specific table format
		Categories=table_in.(2);
		Numbers=table_in.(1);
		look_table=table(Numbers,Categories);
	end

	% Replace categorical with number
	for ii=1:numel(PS)
		Eoi=PS(ii,1).(field);
		ix=find(strcmp(Categories,Eoi));
		PS(ii,1).replace_number=Numbers(ix);
	end

	% Run polygon2GRIDobj
	[OUT]=polygon2GRIDobj(DEM,PS,'replace_number');

	% Remove nonexistent category-number pairs from look_table
	pres=unique(OUT.Z(:));
	ix=ismember(look_table.Numbers,pres);
	look_table=look_table(ix,:);
end