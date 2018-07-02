function Mat2Arc(mat_file,file_prefix)
	%
	% Usage:
	%	Mat2Arc(mat_file,file_prefix);
	%
	% Description:
	% 	Function converts all valid topotoolbox files contained within a mat file
	%	to Arc compatible outputs. Specifically converts any GRIDobjs to
	%	ascii files, any STREAMobjs to shapefiles, any FLOWobjs to ArcGIS 
	%	flow direction grids saved as an ascii file, and any valid mapstructures
	%	to shapefiles.
	%
	% Input:
	%	mat_file - name or path to matfile of interest
	%	file_prefix - characters to add to the front of all output files
	%
	% Example:
	%	Mat file 'Data.mat' containing:
	%		DEM - GRIDobj
	%		FD - FLOWobj
	%		S - STREAMobj
	%		MS_KSN - Mapstructure with ksn data
	%		
	%	Mat2Arc('Data.mat','Example'); would produce:
	%		Example_DEM.txt - ascii of the DEM
	%		Example_FD.txt - ascii of the Arc flow direction grid
	%		Example_S.shp - shapefile of the streams in STREAMobj
	%		Example_MS_KSN.shp - shapefile of the stream data with ksn
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 06/18/18 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	d=whos('-file',mat_file);

	num_vars=numel(d);

	for ii=1:num_vars

		classOI=d(ii,1).class;
		varNM=d(ii,1).name;

		if strcmp(classOI,'GRIDobj')
			varOI=load(mat_file,varNM);
			out_name=[file_prefix '_' varNM '.txt'];
			GRIDobj2ascii(varOI.(varNM),out_name);
		elseif strcmp(classOI,'STREAMobj')
			varOI=load(mat_file,varNM);
			out_name=[file_prefix '_' varNM '.shp'];
			MS=STREAMobj2mapstruct(varOI.(varNM));
			shapewrite(MS,out_name);
		elseif strcmp(classOI,'FLOWobj')
			varOI=load(mat_file,varNM);
			out_name=[file_prefix '_' varNM '.txt'];
			G=FLOWobj2GRIDobj(varOI.(varNM));
			GRIDobj2ascii(G,out_name);
		elseif strcmp(classOI,'struct')
			varOI=load(mat_file,varNM);
			st=varOI.(varNM);
			fn=fieldnames(st);
			if any(strcmp(fn,'Geometry')) & any(strcmp(fn,'X')) & any(strcmp(fn,'Y'))
				out_name=[file_prefix '_' varNM '.shp'];
				shapewrite(st,out_name);
			end
		end
	end
end
		