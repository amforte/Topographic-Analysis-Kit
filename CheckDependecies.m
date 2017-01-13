function CheckDependecies()

	p=ver;

	ix=find(strcmp(cellstr(char(p.Name)),'TopoToolbox'));
	if isempty(ix)
		warning('Fatal error: TopoToolbox does not appear to be installed or on your path, available from https://topotoolbox.wordpress.com/')
	end

	if str2num(p(ix).Version)<2.2
		warning('Version of TopoToolbox is pre 2.2, please download newer version, available from https://topotoolbox.wordpress.com/')
	end

	ix=find(strcmp(cellstr(char(p.Name)),'Image Processing Toolbox'));
	if isempty(ix)
		warning('Fatal error: You do not have a license for the Image Processing Toolbox, TopoToolbox will not function properly')
	end

	ix=find(strcmp(cellstr(char(p.Name)),'Mapping Toolbox'));
	if isempty(ix)
		warning('You do not have a license for the Mapping Toolbox, functions to output shapefiles will not function')
	end
	
	if exist('vline')~=2
		warning('Dependecy function "hline and vline" is not installed or on your path, available from https://www.mathworks.com/matlabcentral/fileexchange/1039-hline-and-vline')
	end

	if exist('cbrewer.m')~=2
		warning('Dependecy function "cbrewer" is not installed or on your path, available from https://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---colorbrewer-schemes-for-matlab')
	end

	if exist('errbar.m')~=2
		warning('Dependecy function "errbar" is not installed or on your path, available from https://www.mathworks.com/matlabcentral/fileexchange/50472-errbar')
	end

	if exist('real2rgb')~=2
		warning('Dependecy function "real2rgb" is not installed or on your path, available from http://www.mathworks.com/matlabcentral/fileexchange/16233-sc-powerful-image-rendering')
	end		

end