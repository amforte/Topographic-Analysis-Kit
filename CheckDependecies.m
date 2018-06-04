function CheckTAKDependecies()

	p=ver;

	warn_flag=false;

	ix=find(strcmp(cellstr(char(p.Name)),'TopoToolbox'));
	if isempty(ix)
		warning('Fatal error: TopoToolbox does not appear to be installed or on your path, available from https://github.com/wschwanghart/topotoolbox/')
	end

	if str2num(p(ix).Version)<2.3
		warning('Version of TopoToolbox is pre 2.3, please download newest version, available from https://github.com/wschwanghart/topotoolbox/')
	end

	ix=find(strcmp(cellstr(char(p.Name)),'Image Processing Toolbox'));
	if isempty(ix)
		warning('Fatal error: You do not have a license for the Image Processing Toolbox, TopoToolbox will not function properly')
		warn_flag=true;
	end

	ix=find(strcmp(cellstr(char(p.Name)),'Mapping Toolbox'));
	if isempty(ix)
		warning('Fatal error: You do not have a license for the Mapping Toolbox, TopoToolbox will not function properly')
		warn_flag=true;
	end


	ix=find(strcmp(cellstr(char(p.Name)),'Optimization Toolbox'));
	if isempty(ix)
		warning('You do not have a license for the Optimization Toolbox, some functions will not work properly')
		warn_flag=true;
	end	

	ix=find(strcmp(cellstr(char(p.Name)),'Statistics and Machine Learning Toolbox'));
	if isempty(ix)
		warning('You do not have a license for the Statistics and Machine Learning Toolbox, some functions will not work properly')
		warn_flag=true;
	end	


	ix=find(strcmp(cellstr(char(p.Name)),'MATLAB'));
	if str2num(p(ix).Version)<9.4
		warning('Topographic Analysis Kit is optimized for MATLAB version post 2018a, some functions may not work properly')
	end	

	if warn_flag
		warning('Your system is missing one or more of the toolboxes necessary for all Topographic Analysis Kit functions to run, consider using the compiled versions of the codes')
	end

end