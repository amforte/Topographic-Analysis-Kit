function cmpMat2Arc(wdir,MatFile,file_prefix)
	% Description:
	% 	Function converts all valid topotoolbox files contained within a mat file
	%	to Arc compatible outputs. Specifically converts any GRIDobjs to
	%	ascii files, any STREAMobjs to shapefiles, any FLOWobjs to ArcGIS 
	%	flow direction grids saved as an ascii file, and any valid mapstructures
	%	to shapefiles.
	%
	% Input:
	%	wdir - full path of working directory
	%	MatFile - full path to matfile of interest
	%	file_prefix - characters to add to the front of all output files
	%
   	% Examples if running for the command line, minus OS specific way of calling main TAK function:
    %   Mat2Arc /path/to/wdir Topo.mat outputs
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 07/02/18 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	MatFile=fullfile(wdir,MatFile);
	d=whos('-file',MatFile);

	num_vars=numel(d);

	for ii=1:num_vars

		classOI=d(ii,1).class;
		varNM=d(ii,1).name;

		if strcmp(classOI,'GRIDobj')
			varOI=load(mat_file,varNM);
			out_name=fullfile(wdir,[file_prefix '_' varNM '.txt']);
			GRIDobj2ascii(varOI.(varNM),out_name);
		elseif strcmp(classOI,'STREAMobj')
			varOI=load(mat_file,varNM);
			out_name=fullfile(wdir,[file_prefix '_' varNM '.shp']);
			MS=STREAMobj2mapstruct(varOI.(varNM));
			shapewrite(MS,out_name);
		elseif strcmp(classOI,'FLOWobj')
			varOI=load(mat_file,varNM);
			out_name=fullfile(wdir,[file_prefix '_' varNM '.txt']);
			G=FLOWobj2GRIDobj(varOI.(varNM));
			GRIDobj2ascii(G,out_name);
		elseif strcmp(classOI,'struct')
			varOI=load(mat_file,varNM);
			st=varOI.(varNM);
			fn=fieldnames(st);
			if any(strcmp(fn,'Geometry')) & any(strcmp(fn,'X')) & any(strcmp(fn,'Y'))
				out_name=fullfile(wdir,[file_prefix '_' varNM '.shp']);
				shapewrite(st,out_name);
			end
		end
	end
end
		