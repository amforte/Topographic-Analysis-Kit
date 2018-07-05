function cmpConditionDEM(wdir,MatFile,method,varargin)
	% 	Wrapper around the variety of methods provided by TopoToolbox for smoothing a stream profile. With the exception of
	% 	'quantc_grid' and 'mingrad' these methods will only modify elevations along the stream network provided to the code.
	% 	See the relevant parent functions for a more in depth description of the behavior of these individual methods. Produces
	%	one figure that compares the long profile of the longest stream within the dataset using the uncondtioned and 
	%	conditioned DEM to provide a quick method of evaluating the result. These methods vary in their complexity and processing
	%	times so it is recommended you understand your choice. Using the 'mincost' method is a good starting place before
	%	exploring some of the more complicated methods. Outputs the resulting DEM as an ascii file for use in other 'cmp*' codes
	%	our outside GIS applications.
	%
	% Required Inputs:
	%	wdir - full path of working directory
	%	MatFile - Name of matfile output from 'cmpMakeStreams'
	%	method - method of conditioning, valid inputs are as follows:
	%		'mincost' - uses the 'mincosthydrocon' function, valid optional inputs are 'mc_method' and 'fillp'.
	%		'mingrad' - uses the 'imposemin' function, valid optional inputs are 'ming'. Note that providing a large minimum gradient
	%			to this code can carve the stream well below the topography.
	%		'quantc' - uses the 'quantcarve' function (for STREAMobjs), valid optional inputs are 'tau','ming', and 'split'. Requires the
	%			Optimization Toolbox and if 'split' is set to true, requires Parallel Processing Toolbox.
	%		'quantc_grid' - uses the 'quantcarve' function for (GRIDobjs), valid optional inputs are 'tau'. Requires the Optimization Toolbox.
	%			This is a computationally expensive calculation and because it operates it on the whole grid, it can take a long time and/or
	%			fail on large grids. The 'quantc' method which only operates on the stream network is significantly fasters and less prone
	%			to failure.
	%		'smooth' - uses the 'smooth' function, valid optional inputs are 'sm_method','split','stiffness','stiff_tribs', and 'positive' depending 
	%			on inputs to optional parameters may require Optimization Toolbox ('sm_method'='regularization' and 'positive'=true) and Parallel
	%			Processing Toolbox ('split'=true).
	%		'crs' - uses the 'crs' function, valid optional inputs are 'stiffness', 'tau', 'ming', 'stiff_tribs', 'knicks', and 'split'. Requires
	%			Optimization Toolbox.
	%		'crslin' - uses the 'crslin' function, valid optional inputs are 'stiffness', 'stiff_tribs', 'ming', 'imposemin', 'attachtomin', 
	%			'attachheads', 'discardflats','precisecoords'
	%			
	% Optional Inputs:
	%	file_name ['cond_DEM'] - name for output ascii file containing conditioned DEM
	%	new_stream_net [] - option to provide name of a matfile containing a new stream network as as output from another function (e.g. cmpFindThreshold) to use
	%		instead of the stream network saved in the MatFile provided to the function. This new stream network must have been generated from the
	%		DEM stored in the provided MatFile
	%	mc_method [interp] - method for 'mincost', valid inputs are 'minmax' or 'interp'
	%	fillp [0.1] - scalar value between 0 and 1 controlling the ratio of carving to filling for 'mincost'
	%	ming [0] - minimum gradient [m/m] in downslope direction, used in 'mingrad','quantc','crs','crslin'
	% 	tau [0.5] - quantile for carving, used in 'quantc', 'quantc_grid', 'crs'.
	%	split [true] - logical flag to utilized parallel processing to independently process tributaries, used in 'quantc_grid', 'smooth', and 'crs'
	%	sm_method ['regularization'] - method for 'smooth', valid inputs are 'regularization' and 'movmean'. 
	%	stiffness [10] - scalar positive value for stiffness penalty, used in 'smooth', 'crs', and 'crslin'
	%	stiff_tribs [true] - logical flag to relax the stiffness penalty at tributary junctions, used in 'smooth', 'crs', and 'crslin'
	% 	knicks [] - nx2 matrix of x and y locations of knickpoints where stiffness penalty should be relaxed, used in 'crs' and 'crslin'
	%	imposemin [false] -logical flag to preprocess DEM with imposemin during crslin
	%	attachtomin [false] - logical flag to prevent elevations from going below profile minima, used in crslin
	%	attachheads [false] - logical flag to fix the channel head elevations, used in crslin
	%	discardflats [false] - logical flag to discard flat portions of profiles, used in crslin
	%	maxcurvature [] - maximum convex curvature at any vertices along profile, used in crslin
	%	precisecoords [] - nx3 matrix with x, y, and z coordinates of points that the smoothed profile must past through, used in crslin
	%
	% Outputs:
	%		georeferenced ascii file of the processed dem
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 07/02/18 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if isdeployed
		if ~isempty(varargin)
			varargin=varargin{1};
		end
	end

	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'cmpConditionDEM';
	addRequired(p,'wdir',@(x) ischar(x));
	addRequired(p,'MatFile',@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.mat'))));
	addRequired(p,'method',@(x) ischar(validatestring(x,{'mincost','mingrad','quantc','quantc_grid','smooth','crs','crslin'})))

	addParameter(p,'file_name','cond_DEM',@(x) ischar(x));
	addParameter(p,'new_stream_net',[],@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.mat'))));	
	addParameter(p,'mc_method','interp',@(x) ischar(validatestring(x,{'minmax','interp'})));
	addParameter(p,'fillp',0.1,@(x) isscalar(x) && isnumeric(x) && x>=0 && x<=1);
	addParameter(p,'ming',0,@(x) isscalar(x) && isnumeric(x) && x>=0);
	addParameter(p,'tau',0.5,@(x) isscalar(x) && isnumeric(x) && x>=0 && x<=1);
	addParameter(p,'split',true,@(x) islogical(x))
	addParameter(p,'sm_method','regularization',@(x) ischar(validatestring(x,{'regularization','movmean'})));
	addParameter(p,'stiffness',10,@(x) isscalar(x) && isnumeric(x) && x>=0);
	addParameter(p,'stiff_tribs',true,@(x) islogical(x));
	addParameter(p,'positive',true,@(x) islogical(x));
	addParameter(p,'knicks',[],@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.txt'))));	
	addParameter(p,'imposemin',false,@(x) islogical(x));
	addParameter(p,'attachtomin',false,@(x) islogical(x));
	addParameter(p,'attachheads',false,@(x) islogical(x));
	addParameter(p,'discardflats',false,@(x) islogical(x));
	addParameter(p,'maxcurvature',[],@(x) isnumeric(x) && isscalar(x));
	addParameter(p,'precisecoords',[],@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.txt'))));	

	parse(p,wdir,MatFile,method,varargin{:});
	wdir=p.Results.wdir;
	MatFile=p.Results.MatFile;
	method=p.Results.method;

	nsn=p.Results.new_stream_net;
	file_name=p.Results.file_name;
	mc_method=p.Results.mc_method;
	fillp=p.Results.fillp;
	ming=p.Results.ming;
	tau=p.Results.tau;
	split=p.Results.split;
	sm_method=p.Results.sm_method;
	K=p.Results.stiffness;
	st=p.Results.stiff_tribs;
	po=p.Results.positive;
	knicks=p.Results.knicks;
	%CRSLIN Parameters
	p1=p.Results.imposemin;
	p2=p.Results.attachtomin;
	p3=p.Results.attachheads;
	p4=p.Results.discardflats;
	p5=p.Results.maxcurvature;
	p6=p.Results.precisecoords;

	% Load in results of matfile
	MatFile=fullfile(wdir,MatFile);
	load(MatFile,'DEM','FD','S');

	if ~isempty(nsn)
		nsn=fullfile(wdir,nsn);
		SD=whos('-file',nsn);
		VS=cell(numel(SD),1);
		for ii=1:numel(SD)
			VS{ii}=SD(ii,1).name;
		end

		if any(strcmp(VS,'Sn'))
			load(nsn,'Sn');
			S=Sn;
		elseif any(strcmp(VS,'Sc'))
			load(nsn,'Sc');
			S=Sc;
		elseif any(strcmp(VS,'S'));
			load(nsn,'S');
		end

		if ~validatealignment(S,DEM)
			error('Supplied new stream net does not appear to have been generated from the same DEM supplied to the function');
		end
	end	

	switch method
	case 'mincost'
		zc=mincosthydrocon(S,DEM,mc_method,fillp);
		DEMc=GRIDobj(DEM);
		DEMc.Z(DEMc.Z==0)=NaN;
		DEMc.Z(S.IXgrid)=zc;
	case 'mingrad'
		DEMc=imposemin(FD,DEM,ming);
	case 'quantc'
		% Split parameter expects a 1 or 2 not a logical for quantcarve
		if split
			sp=2;
		elseif ~split
			sp=1;
		end
		[zc]=quantcarve(S,DEM,tau,'split',sp,'mingradient',ming);
		DEMc=GRIDobj(DEM);
		DEMc.Z(DEMc.Z==0)=NaN;
		DEMc.Z(S.IXgrid)=zc;	
	case 'quantc_grid'
		DEMc=quantcarve(FD,DEM,tau);
	case 'smooth'
		zc=smooth(S,DEM,'method',sm_method,'split',split,'K',K,'nstribs',st,'positive',po);
		DEMc=GRIDobj(DEM);
		DEMc.Z(DEMc.Z==0)=NaN;
		DEMc.Z(S.IXgrid)=zc;
	case 'crs'
		% Split parameter expects a 0, 1, or 2
		if split
			sp=2;
		elseif ~split
			sp=0;
		end	

		if isempty(knicks)
			knicksix=[];
		else
			knicks=readtable(knicks);
			knicks=table2array(knicks);
			knicksix=coord2ind(DEM,knicks(:,1),knicks(:,2));
		end

		[zc]=crs(S,DEM,'K',K,'tau',tau,'mingradient',ming,'split',sp,'nonstifftribs',st,'knickpoints',knicksix);
		DEMc=GRIDobj(DEM);
		DEMc.Z(DEMc.Z==0)=NaN;
		DEMc.Z(S.IXgrid)=zc;		
	case 'crslin'
		knicks=readtable(knicks);
		knicks=table2array;
		[zc,~,~]=crslin(S,DEM,'K',K,'mingradient',ming,'nonstifftribs',st,'knickpoints',knicks,'imposemin',p1,'attachtomin',p2,'attachheads',p3,'maxcurvature',p4,'precisecoords',p5);
		DEMc=GRIDobj(DEM);
		DEMc.Z(DEMc.Z==0)=NaN;
		DEMc.Z(S.IXgrid)=zc;	
	end

	f1=figure(1);
	set(f1,'Units','normalized','Position',[0.05 0.1 0.45 0.8],'renderer','painters');
	clf

	SL=trunk(klargestconncomps(S,1));

	subplot(2,2,1:2)
	hold on
	plotdz(SL,DEM,'color','k');
	plotdz(SL,DEMc,'color','r');
	legend('Unconditioned DEM','Conditioned DEM','location','best');
	hold off

	subplot(2,2,3)
	hold on
	plotdz(SL,DEM-DEMc,'color','k');
	legend('Elevation Difference between Un- and Conditioned DEM','location','best');
	hold off

	subplot(2,2,4)
	hold on
	imagesc(DEM-DEMc);
	colorbar;
	title('Elevation Difference Between Un- and Conditioned DEM');
	hold off

	outFile=fullfile(wdir,'ConditionedDEM.txt');
	GRIDobj2ascii(DEMc,outFile);
	outMat=fullfile(wdir,'ConditionedDEM.mat');
	save(outMat,'DEMc');
end