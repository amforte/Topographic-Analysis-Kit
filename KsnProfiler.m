function [knl,ksn_master,bnd_list,kn_list,Sc]=KsnProfiler(DEM,FD,A,S,varargin)
	%
	% Usage:
	%	[knl,ksn_master,bnd_list,kn_list,Sc]=KsnProfiler(DEM,FD,A,S);
	%	[knl,ksn_master,bnd_list,kn_list,Sc]=KsnProfiler(DEM,FD,A,S,'name',value,...);
	%
	% Description:
	% 	Function to interactively select channel heads and define segements over which to calculate channel steepness values.
	% 	This function is designed to be similar to the operation of Profiler_51, with some improvements. Function will display map
	%	with the stream network and expects the user to select a location near a channel head of interest. The user will be then 
	%	prompted to confirm that the defined stream is the desired choice. Finally, displays of the chi-z and longitudinal profile 
	%	of the selected river will appear and the user is expected to define (with mouse clicks) any obvious segments with different 
	%	channel steepness (or concavity) on either the chi-z plot, the stream profile, or a slope-area plot (see 'pick_method' option). 
	%	When done selecting press enter/return. The user will be prompted whether they wish to continue picking streams or if they are done. 
	%	When done picking streams, the function will output five different products (see below) and produce a shapefile of the selected streams 
	%	with ksn, concavity, area, and gradient.
	%
	% Picking Bounds:
	%	The function expects that you will select bounds for defining different ksn segments on the appropriate plot (the one with red
	%	axes, based on what you provide to 'pick_method') with LEFT MOUSE CLICKS. If you click any other mouse button or press any key
	%	(other than the return/enter key) the code will recognize this as selecting a knickpoint location. This is distinguished from a 
	%	bound in that it's positions will be logged and recorded in the 'kn_list' output, but it will NOT be used to define different
	%	ksn segments. This is to provide the user the ability to 'mark' locations on the profile that you do not wish to use as bounds.
	%	Both bounds and these marked locations are provided to inputs to the companion 'ClassifyKnicks' function (where the bounds and
	%	these knickpoints / marked locations are distinguished from each other).
	%	
	% Required Inputs:
	%	DEM - Digital Elevation as a GRIDobj, assumes unconditioned DEM (e.g. DEMoc from ProcessRiverBasins or output from MakeStreams)
	%	FD - Flow direction as FLOWobj
	%	S - Stream network as STREAMobj
	%	A - Flow accumulation GRIDobj
	%
	%%%%%%%%%%%%%%%%%%
	% Optional Inputs:
	%
	%%% Restart Picking
	%	restart [] - providing an entry to this parameter allows the user to restart a run, either a run that you succesfully completed
	%			but want to restart or a run that failed part way through either because of an error or because you aborted out. While 
	%			the code is running, it will save data necessary to restart in a mat file called '*_restart.mat'. If the code succesfully
	%			completes, this '*_restart.mat' file will be deleted. DO NOT DELETE THIS FILE WHILE THE CODE IS RUNNING OR IF THE CODE FAILS
	%			AND YOU WISH TO SALVAGE THE RUN. You can also call use restart if you just wish to restart picking streams from a previously 
	%			completed run. If you run the code with an 'input_method' other than 'interactive' and the code succesfully completes (i.e
	%			you fit all the streams selected via the input method you choose and you did not stop the code early) then running with restart
	%			will not do anything. If you wish to restart, you do not need to define any of the original parameters, these are saved in the
	%			output files and will be loaded in, you only need to provide the four required inputs (see example) along with the restart parameter. 
	%			Valid inputs to restart are:
	%		'continue' - will restart the run. If used with a completed or failed 'interactive' run will repopulate the map with already picked 
	%			streams and you can continue picking. If using with a non interactive input method that either failed or you aborted, will start 
	%			on the next stream in the sequence.
	%		'skip' - only a meaningful input for a non interactive run. This will skip the next stream segment in the sequence. This would be useful
	%			if a particular stream segment causes the code to error, this way you can skip that stream in a restart without having to modifying
	%			the stream network.
	%
	%%% Main Options
	% 	input_method ['interactive'] - parameter which controls how streams of interest are supplied:
	%		'interactive' - user picks streams of interest by selecting channelheads on a map, this option will also iteratively build a 
	%			channel steepness map as the user picks more streams.
	%		'all_streams' - will use the supplied STREAMobj and iterate through all channel heads. There is an internal parameter to avoid 
	%			selecting streams that are too short to properly fit (mostly relevant if 'junction method' is set to 'check'). The default 
	%			value is ~4 * the DEM cellisze, the user can change this value by providing an input for the optional parameter 'min_channel_length',
	%			input should be in map units and greater than the default. You can use a code like 'SegmentPicker' to select portions of a STREAMobj
	%		'stream_length' - will use supplied STREAMobj and entry to 'min_length_to_extract' to iterate through all streams that are longer than
	%			the length provided to 'min_length_to_extract'. There is an internal parameter to avoid selecting streams that are too short to fit
	%			(mostly relevant if 'junction method' is set to 'check'). The default value is ~4 * the DEM cellsize, the user can change this value
	%			by providing an input for the optional parameter 'min_channel_length', input should be in map units and greater than the default.
	%		'channel_heads' - will use a supplied list of coordinates of channel heads to select and iterate through streams of interest. If this
	%			option is used, the user must provide an input for the optional 'channel_head_list' parameter.
	%	pick_method ['chi'] - choice of how you want to pick stream segments. The diagram within which to pick based on your selection will be 
	%			outline in red. Valid inputs are:
	%		'chi' - select segments on a chi - z plot (recommended and default)
	%		'stream' - select segments on a longitudinal profile
	%		'slope_area' - select segments on a slope area plot
	%	junction_method ['check'] - choice of how to deal with stream junctions:
	%		'check' - after each choice, will check whether downstream portions of the selected stream have already been fit, and if it has,
	%			the already fit portion of the stream will not be displayed or refit
	%		'ignore' - each stream will be displayed from its head to mouth independent of whether portions of the same stream network have 
	%			been fit
	%	concavity_method ['ref']- options for concavity:
	%		'ref' - uses a reference concavity, the user can specify this value with the reference concavity option (see below)
	%		'auto' - function finds a best fit concavity for each selected stream, if used in conjunction with 'junction_method','check'
	%			this means that short sections of streams picked will auto fit concavity that may differ from downstream portions of the same
	%			streams
	%
	%%% Input Method Options 
	%	min_channel_length [] - minimum channel length for consideration when using the 'all_streams' method of input, provide in map units.
	%	channel_head_list [] - m x 2 array of x and y coordinates of channel heads OR the name / location of a point shapefile of channel heads, 
	%			one of these is required when using 'channel_heads' method of input, must be in the same coordinate system as the input DEM etc. 
	%			The code will attempt to find the nearest channel head to the coordinates you provided, so the closer the provided user coordinates
	%			are to channel heads, the more accurate this selection method will be.
	%	min_length_to_extract [] - minimum stream length (in map units) to extract streams if 'input_method' is set to 'stream_length'.
	%
	%%% Redefine Threshold Area Options
	%	redefine_threshold [false] - logical flag to initiate an extra step for each stream where you manually define the hillslope-fluvial 
	%			transition (this will result in overriding the threshold area you used to generate the supplied STREAMobj, and it will also produce
	%			a STREAMobj with a variable threshold area for channel definition). See additional optional input 'rd_pick_method'.
	%	rd_pick_method ['slope_area'] - plot to use to choose new threshold area if 'redefine_threshold' is set to true. Valid inputs are 
	%			'slopearea' and 'chi'.
	%
	%%% Stream Network Modification Options
	%	complete_networks_only [false] - if true, the code will filter out portions of the stream network that are incomplete prior to choosing
	%			streams
	%	min_elev [] - minimum elevation below which the code stops extracting channel information (no action if left empty)
	%	max_area [] - maximum drainage area above which the code stops extracting channel information (in square map units, no action if left empty)
	%
	%%% Hydrological Conditioning Options
	%	conditioned_DEM [] - option to provide a hydrologically conditioned DEM for use in this function (do not provide a conditoned DEM
	%			for the main required DEM input!) which will be used for extracting elevations. See 'ConditionDEM' function for options 
	%			for making a hydrological conditioned DEM. If no input is provided the code defaults to using the mincosthydrocon function.
	%	interp_value [0.1] - value (between 0 and 1) used for interpolation parameter in mincosthydrocon (not used if user provides a 
	%			conditioned DEM). Values closer to 0 tend to 'carve' more, whereas values closer to 1 tend to fill. See info for 
	%			'mincosthydrocon'
	%
	%%% Display Options
	%	display_slope_area [false] - logical flag to display slope area plots. Some people love slope area plots (like one of the authors of
	%			the supporting paper), some people hate slope area plots (like the other author of the supporting paper), so you can either 
	%			not draw them at all (false - default) or include them (true). This will automatically be set to true if you select 'slope_area'
	%			as the 'pick_method'.
	%	plot_type ['vector'] - expects either 'vector' or 'grid', default is 'vector'. Controls whether all streams are drawn as individual 
	%			lines ('vector') or if the stream network is plotted as a grid and downsampled ('grid'). The 'grid' option is much faster on 
	%			large datasets, but can result in inaccurate channel head selection. The 'vector' option is easier to see, but can be very 
	%			slow to load and interact with on large datasets.	
	%
	%%% Constants
	%	ref_concavity [0.50] - refrence concavity used if 'theta_method' is set to 'ref'
	%	smooth_distance [1000] - distance in map units over which to smooth ksn measures when converting to shapefile
	%	max_ksn [250] - maximum  ksn used for the color scale, will not effect actual results, for display purposes only
	%	threshold_area [1e6] - used to redraw downsampled stream network if 'plot_type' is set to 'grid' 
	%
	%%% Output Options
	%	stack_method ['stack'] - if 'junction_method' is set to 'ignore', this parameter will control how the function deals with overlapping sections
	%			of stream networks when generating the shapefile. Valid inputs are 'stack' (default) and 'average'. If set to 'stack', the output shapefile 
	%			will have multiple stacked polylines in overlapping portions of networks. This is similar to how Profiler51 worked. If set to 'average', the
	%			function will average overlapping portions of networks on a node by node basis. Note that if 'junction_method' is set to 'check', then this 
	%			parameter is ignored.
	%	shape_name ['ksn'] - name for the shapefile to be export, must have no spaces to be a valid name for ArcGIS and should NOT include the '.shp'
	%	save_figures [false] - logical flag to either save figures showing ksn fits (true) or to not (false - default)	
	%
	%%%%%%%%%%	
	% Outputs:
	%	knl - n x 12 matrix of node list for selected stream segments, columns are x coordinate, y coordinate, drainage area, ksn, negative ksn error,
	%		positive ksn error, reference concavity, best fit concavity, mininum threshold area, gradient, fit residual, and an identifying number. Note
	%		that if using the code in 'concavity_method','auto' mode then the reference concavity and best fit concavity columns will be the same.
	%	ksn_master - identical to knl but as a cell array where individual cells are individual selected channels
	%	bnd_list - n x 4 matrix of selected bounds for fitting ksn, columns are x coordinate, y coordinate, elevation, and the stream identifying number, 
	%		 also output as a seperate shapefile ('_bounds.shp'). If x y and z values appear as NaN, this indicates that bounds for this stream were not selected. 
	%		
	%	kn_list - n x 4 matrix of extra identified knickpoints (not bounds),columns are x coordinate, y coordinate, elevation, and the stream identifying number,
	%		 also output as a seperate shapefile ('_knicks.shp'. If x y and z values appear as NaN, this indicates that knicks for this stream were not selected.
	%	Sc - STREAMobj of selected streams
	%
	% Examples:
	%	[knl,ksn_master,bnd_list,kn_list,Sc]=KSN_Profiler(DEM,FD,A,S);
	%	[knl,ksn_master,bnd_list,kn_list,Sc]=KSN_Profiler(DEM,FD,A,S,'junction_method','ignore','ref_concavity',0.65,'max_ksn',500);
	%	[knl,ksn_master,bnd_list,kn_list,Sc]=KSN_Profiler(DEM,FD,A,S,'input_method','channel_heads','channel_head_list',channel_head_array);
	%	[knl,ksn_master,bnd_list,kn_list,Sc]=KSN_Profiler(DEM,FD,A,S,'input_method','channel_heads','channel_head_list','channel_heads.shp');
	%	Restart Examples:
	%	[knl,ksn_master,bnd_list,kn_list,Sc]=KSN_Profiler(DEM,FD,A,S,'restart','continue'); % Continues where you left off
	%	[knl,ksn_master,bnd_list,kn_list,Sc]=KSN_Profiler(DEM,FD,A,S,'restart','skip'); % Skips next stream in non interactive sequence
	%	
	%
	% Note:
	%	-If no boundaries/knickpoints are selected for any of the streams selected, then a '_bounds.shp'/'_knicks.shp' shapefile will not be produced.
	%	-The '*_profiler.mat' that is saved out contains additional files besides the formal outputs of the code. These additional variables
	%		are necessary to be able to restart a run using the 'restart' option.
	%	-If you have set 'save_figures' to true, DO NOT close figures manually as this will cause the code to error.
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 01/09/20 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'KsnProfiler';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'FD',@(x) isa(x,'FLOWobj'));
	addRequired(p,'A',@(x) isa(x,'GRIDobj'));
	addRequired(p,'S',@(x) isa(x,'STREAMobj'));

	addParameter(p,'shape_name','ksn',@(x) ischar(x));
	addParameter(p,'smooth_distance',1000,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'concavity_method','ref',@(x) ischar(validatestring(x,{'ref','auto'})));
	addParameter(p,'complete_networks_only',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'pick_method','chi',@(x) ischar(validatestring(x,{'chi','stream','slope_area'})));
	addParameter(p,'junction_method','check',@(x) ischar(validatestring(x,{'check','ignore'})));	
	addParameter(p,'ref_concavity',0.50,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'redefine_threshold',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'rd_pick_method','slope_area',@(x) ischar(validatestring(x,{'chi','slope_area'})));
	addParameter(p,'display_slope_area',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'max_ksn',250,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'min_elev',[],@(x) isscalar(x) && isnumeric(x) || isempty(x));
	addParameter(p,'max_area',[],@(x) isscalar(x) && isnumeric(x) || isempty(x));
	addParameter(p,'plot_type','vector',@(x) ischar(validatestring(x,{'vector','grid'})));
	addParameter(p,'threshold_area',1e6,@(x) isnumeric(x));
	addParameter(p,'input_method','interactive',@(x) ischar(validatestring(x,{'interactive','channel_heads','all_streams','stream_length'})));
	addParameter(p,'channel_head_list',[],@(x) isnumeric(x) && size(x,2)==2 || regexp(x,regexptranslate('wildcard','*.shp')) || isempty(x));
	addParameter(p,'min_length_to_extract',[],@(x) isnumeric(x) && isscalar(x) || isempty(x));
	addParameter(p,'min_channel_length',[],@(x) isnumeric(x) && isscalar(x) || isempty(x));
	addParameter(p,'conditioned_DEM',[],@(x) isa(x,'GRIDobj') || isempty(x));
	addParameter(p,'interp_value',0.1,@(x) isnumeric(x) && x>=0 && x<=1);
	addParameter(p,'save_figures',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'restart',[],@(x) ischar(validatestring(x,{'continue','skip'})) || isempty(x));
	addParameter(p,'stack_method','stack',@(x) ischar(validatestring(x,{'average','stack'})));
	addParameter(p,'restart_loc',[],@(x) ischar(x) || isempty(x));

	parse(p,DEM,FD,A,S,varargin{:});
	DEM=p.Results.DEM;
	FD=p.Results.FD;
	S=p.Results.S;
	A=p.Results.A;

	shape_name=p.Results.shape_name;
	smooth_distance=p.Results.smooth_distance;
	theta_method=p.Results.concavity_method;
	cno=p.Results.complete_networks_only;
	pick_method=p.Results.pick_method;
	junction_method=p.Results.junction_method;
	ref_theta=p.Results.ref_concavity;
	display_slope_area=p.Results.display_slope_area;
	plot_type=p.Results.plot_type;
	threshold_area=p.Results.threshold_area;
	input_method=p.Results.input_method;
	chl=p.Results.channel_head_list;
	mlte=p.Results.min_length_to_extract;
	min_channel_length=p.Results.min_channel_length;
	iv=p.Results.interp_value;
	DEMc=p.Results.conditioned_DEM;
	min_elev=p.Results.min_elev;
	max_area=p.Results.max_area;
	save_figures=p.Results.save_figures;
	redefine_thresh=p.Results.redefine_threshold;
	rd_pick_method=p.Results.rd_pick_method;
	mksn=p.Results.max_ksn;
	restart=p.Results.restart;
	stack_method=p.Results.stack_method;
	restart_loc=p.Results.restart_loc; % Hidden parameter for GUI deployed versions

	wtb=waitbar(0,'Preparing inputs...');

	% Set restart flag
	if ~isempty(restart)
		rf=true;
	else
		rf=false;
	end

	% Check to see if restart if invoked and replace parameters if necessary
	if rf
		if ~isempty(restart_loc)
			restart_file=dir('*_restart.mat');
			out_file=dir('*_profiler.mat');
		else
			restart_file=dir(fullfile(restart_loc,'*_restart.mat'));
			out_file=dir(fullfile(restart_loc,'*_profiler.mat'));			
		end

		if numel(restart_file)>1
			if isdeployed 
				errordlg('Multiple restart files were found, please remove non-target restart files from active directory or search path')
			end
			error('Multiple restart files were found, please remove non-target restart files from active directory or search path');
		elseif isempty(restart_file) & ~isempty(out_file)		
			if numel(out_file)>1
				if isdeployed
					errordlg('Multiple profiler output mat files were found, please remove non-target restart files from active directory or search path')
				end
				error('Multiple profiler output mat files were found, please remove non-target restart files from active directory or search path');
			end
			load(out_file(1,1).name,'input_params');
			r_type='c';
		elseif isempty(restart_file) & isempty(out_file)
			if isdeployed
				errordlg('No previous run files were found, do not provide an entry to "restart" if this is your first time running KsnProfiler for these data')
			end
			error('No previous run files were found, do not provide an entry to "restart" if this is your first time running KsnProfiler for these data')
		else
			load(restart_file(1,1).name,'input_params');
			r_type='r';
		end

		% Load in parameters from previous run
		shape_name=input_params.shape_name;
		smooth_distance=input_params.smooth_distance;
		theta_method=input_params.concavity_method;
		cno=input_params.complete_networks_only;
		pick_method=input_params.pick_method;
		junction_method=input_params.junction_method;
		ref_theta=input_params.ref_concavity;
		display_slope_area=input_params.display_slope_area;
		plot_type=input_params.plot_type;
		threshold_area=input_params.threshold_area;
		input_method=input_params.input_method;
		chl=input_params.channel_head_list;
		mlte=input_params.min_length_to_extract;
		min_channel_length=input_params.min_channel_length;
		iv=input_params.interp_value;
		DEMc=input_params.conditioned_DEM;
		min_elev=input_params.min_elev;
		max_area=input_params.max_area;
		save_figures=input_params.save_figures;
		redefine_thresh=input_params.redefine_threshold;
		rd_pick_method=input_params.rd_pick_method;
		mksn=input_params.max_ksn;	
	end	

	waitbar(1/4,wtb);

	% Store out parameters in both final file and restart file
	out_mat_name=[shape_name '_profiler.mat'];
	out_restart_name=[shape_name '_restart.mat'];
	if rf
		save(out_mat_name,'input_params','-append');
		if exist(out_restart_name)==2
			save(out_restart_name,'input_params','-append');
		else
			save(out_restart_name,'input_params','-v7.3');
		end
	else
		input_params=p.Results;
		save(out_mat_name,'input_params','-v7.3');
		save(out_restart_name,'input_params','-v7.3');
	end

    % Remove edges if flag is thrown
    if cno
    	S=removeedgeeffects(S,FD,DEM);
    end

	% Find channel heads
	[ch]=streampoi(S,'channelheads','xy');

	% Create master KSN colormap
	KSN_col=ksncolor(100);

	waitbar(2/4,wtb);

	% Perform some checks and reassign values as needed
	if strcmp(input_method,'channel_heads')
		if isempty(chl)
			if isdeployed
				errordlg('Selection method is "channel_heads", must provide an input for the "channel_head_list" parameter')
			end
			error('Selection method is "channel_heads", must provide an input for the "channel_head_list" parameter');
		end

		if ischar(chl)
			% Load in shapefile if provided
			if logical(regexp(chl,regexptranslate('wildcard','*.shp')))
				ch_ms=shaperead(chl);
				ch_t=struct2table(ch_ms);
				if ~strcmp(ch_t.Geometry(1),'Point')
					if isdeployed
						errordlg('Shapefile provided as "channel_heads" does not appear to be a point shapefile')
					end
					error('Shapefile provided as "channel_heads" does not appear to be a point shapefile');
				end
				xi=ch_t.X;
				yi=ch_t.Y;
				chl=[xi yi];
			end
		end

		% Snap to nearest channel heads
		dists=pdist2(chl,ch);
		[~,s_ch_ix]=min(dists,[],2);
		s_ch=ch(s_ch_ix,:);
		num_ch=size(s_ch,1);

		plot_type='none';
		input_method='preselected';
	elseif strcmp(input_method,'all_streams')
		s_ch=ch;
		num_ch=size(s_ch,1);

		if isempty(min_channel_length);
			min_channel_length=sqrt((4*DEM.cellsize)^2 + (4*DEM.cellsize)^2);
		end

		plot_type='none';
		input_method='preselected';
	elseif strcmp(input_method,'stream_length')
		if isempty(mlte)
			if isdeployed
				errordlg('Selection method is "stream_length", must provide an input for "mean_length_to_extract" parameter')
			end
			error('Selection method is "stream_length", must provide an input for "mean_length_to_extract" parameter');
		end

		FlowD=flowdistance(FD);
		ix=coord2ind(DEM,ch(:,1),ch(:,2));
		idx=FlowD.Z(ix)>=mlte;

		s_ch=ch(idx,:);

		if isempty(s_ch)
			if isdeployed
				errordlg('Input to "mean_length_to_extract" resulted in no streams being selected, reduce the length and retry')
			end
			error('Input to "mean_length_to_extract" resulted in no streams being selected, reduce the length and retry')
		end

		num_ch=size(s_ch,1);

		if isempty(min_channel_length);
			min_channel_length=sqrt((4*DEM.cellsize)^2 + (4*DEM.cellsize)^2);
		end
		plot_type='none';
		input_method='preselected';
	end

	if strcmp(pick_method,'slope_area')
		display_slope_area=true;
	end

	if ~isempty(min_elev) && ~isempty(max_area)
		if isdeployed
			errordlg('Providing values to both "min_elev" and "max_area" is not permitted, please only provide values for one of these parameters')
		end
		error('Providing values to both "min_elev" and "max_area" is not permitted, please only provide values for one of these parameters')
	end

	% Hydrologically condition dem
	if isempty(DEMc)
		zc=mincosthydrocon(S,DEM,'interp',iv);
		DEMc=GRIDobj(DEM);
		DEMc.Z(DEMc.Z==0)=NaN;
		DEMc.Z(S.IXgrid)=zc;
	elseif ~isempty(DEMc) & redefine_thresh
		if isdeployed
			warndlg(['Supplying a Conditioned DEM and redefining the drainage area threshold may produce unexpected results as it may be necessary ' ...
				'to recondition portions of the DEM using a different method. To avoid this, make sure that a very low threshold area was used to ' ...
				'produce the supplied Conditioned DEM.'])
		else
			warning(['Supplying a Conditioned DEM and redefining the drainage area threshold may produce unexpected results as it may be necessary ' ...
				'to recondition portions of the DEM using a different method. To avoid this, make sure that a very low threshold area was used to ' ...
				'produce the supplied Conditioned DEM.'])
		end
	end

	% Make gradient
	G=gradient8(DEMc);
	% Make drainage area
	DA=A.*A.cellsize^2;

	waitbar(3/4,wtb);

	% Modify provided stream network if minimum elevation or maximum drainage area options are included
	if ~isempty(min_elev)
		zel=getnal(S,DEMc);
		idx=zel>min_elev;
		new_ix=S.IXgrid(idx);
		W=GRIDobj(DEMc,'logical');
		W.Z(new_ix)=true;
		S=STREAMobj(FD,W);
	elseif ~isempty(max_area)
		DA=A.*(A.cellsize^2);
		zda=getnal(S,DA);
		idx=zda<max_area;
		new_ix=S.IXgrid(idx);
		W=GRIDobj(DEMc,'logical');
		W.Z(new_ix)=true;
		S=STREAMobj(FD,W);		
	end

	% Generate flow distance grid if redefining threshold
	if redefine_thresh
		FLUS=flowdistance(FD);
	end

	% Generate Map Figures for interactive picking
	switch plot_type
	case 'grid'	

		disp('Downsampling datasets for display purposes')
		% Redo flow direction	
		DEMr=resample(DEM,DEM.cellsize*4);
		FDr=FLOWobj(DEMr,'preprocess','carve');
		% True outlets
		out_T_xy=streampoi(S,'outlets','xy');
		% Downsampled total stream network
		Sr_temp=STREAMobj(FDr,'minarea',threshold_area,'unit','mapunits');
		out_D_xy=streampoi(Sr_temp,'outlets','xy');
		out_D_ix=streampoi(Sr_temp,'outlets','ix');
		% Find if outlet list is different
		dists=pdist2(out_T_xy,out_D_xy);
		[~,s_out_ix]=min(dists,[],2);
		out_D_ix=out_D_ix(s_out_ix);
		% Rebuild downsampled network
		Sr=STREAMobj(FDr,'minarea',threshold_area,'unit','mapunits','outlets',out_D_ix);
		% Turn it into a grid
		SG=STREAMobj2GRIDobj(Sr);

		% Initiate Map Figure

		f1=figure(1);
		set(f1,'Visible','off');
		hold on
		[RGB]=imageschs(DEMr,SG,'colormap','gray');
		[~,R]=GRIDobj2im(DEMr);
		imshow(flipud(RGB),R);
		axis xy
		colormap(KSN_col);
		caxis([0 mksn])
		c1=colorbar;
		ylabel(c1,'Channel Steepness')
		if ~verLessThan('matlab','9.5')
	        disableDefaultInteractivity(gca);
	    end	
		hold off
		set(f1,'Visible','on','Units','normalized','Position',[0.05 0.1 0.5 0.5],'renderer','painters');

	case 'vector'

		% Initiate Map Figure;

		f1=figure(1);
		set(f1,'Visible','off');

		[RGB]=imageschs(DEM,DEM,'colormap','gray');
		[~,R]=GRIDobj2im(DEM);

		imshow(flipud(RGB),R);
		axis xy
		hold on
		colormap(KSN_col);
		plot(S,'-w');
		caxis([0 mksn])
		c1=colorbar;
		ylabel(c1,'Channel Steepness')
		if ~verLessThan('matlab','9.5')
	        disableDefaultInteractivity(gca);
	    end			
		hold off
		set(f1,'Visible','on','Units','normalized','Position',[0.05 0.1 0.5 0.5],'renderer','painters');	
	end

	waitbar(1,wtb);
	close(wtb);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Main switch for graphical selection vs list of channel heads %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	switch input_method
	case 'interactive'
		% Initiate counters and while loop values
		str1='N';
		str2='Y';

		if rf & strcmp(r_type,'c')
			% Load in data from previous run
			load(out_mat_name,'count','ksn_master','bnd_master','kn_master','res_master','Sc');
			ii=count+1;
			% Regenerate plotted streams
			km=vertcat(ksn_master{:});
			kix=coord2ind(DEM,km(:,1),km(:,2));
			K=GRIDobj(DEM);
			K.Z(kix)=km(:,4);
			figure(1)
			hold on
			plotc(Sc,K);
			hold off
			clear km kix K;
		elseif rf & strcmp(r_type,'r')
			% Load in data from failed or aborted run
			load(out_restart_name,'count','ksn_master','bnd_master','kn_master','res_master','Sc');
			ii=count+1;
			% Regenerate plotted streams
			km=vertcat(ksn_master{:});
			kix=coord2ind(DEM,km(:,1),km(:,2));
			K=GRIDobj(DEM);
			K.Z(kix)=km(:,4);
			figure(1)
			hold on
			plotc(Sc,K);
			hold off
			clear km kix K;
		else
			ii=1;
		end

		if strcmp(theta_method,'ref')
			% Autocalculate ksn for comparison purposes
			[auto_ksn]=KSN_Quick(DEM,A,S,ref_theta);
		end

		% Begin picking streams
		while strcmpi(str2,'Y');
			% Select channel to fit
			while strcmpi(str1,'N');	
				str3='R'; %Reset redo flag

	            figure(1)
	            hold on
	            title('Zoom or pan to area of interest and then press enter');
	            hold off
				pause()

	            figure(1)
	            hold on
	            title('Choose point near channel head of interest');
	            hold off

				[x,y]=ginput(1);
				pOI=[x y];

				% Find nearest channel head
				distance=sqrt(sum(bsxfun(@minus, ch, pOI).^2,2));
				chOI=ch(distance==min(distance),:);

				% Build logical raster
				ix=coord2ind(DEM,chOI(:,1),chOI(:,2));
				IX=GRIDobj(DEM,'logical');
				IX.Z(ix)=true;

				if redefine_thresh

					% Extract stream of interest
					Sn=modify(S,'downstreamto',IX);

					figure(f1)
					hold on
					p1=plot(Sn,'-b','LineWidth',2);
					hold off

					qa1=questdlg('Is this the stream segment you wanted?','Stream Selection','Yes','No','Yes');
					switch qa1
					case 'Yes'

						delete(p1);

						[Sn]=RedefineThreshold(DEM,FD,A,Sn,FLUS,ref_theta,rd_pick_method,smooth_distance,ii,save_figures,shape_name);
						% Update DEMc
						if any(isnan(getnal(Sn,DEMc)));
							zc=mincosthydrocon(Sn,DEM,'interp',iv);
							DEMc.Z(Sn.IXgrid)=zc;
						end
						% Recalculate auto_ksn
						if strcmp(theta_method,'ref')
							[auto_ksn]=KSN_Quick(DEM,A,Sn,ref_theta);
						end

						str1='Y';

						if strcmp(junction_method,'check')
							if ii>1
								[IIXX,~,~,Si]=intersectlocs(Sc,Sn);
								if isempty(IIXX)
									Sn=Sn;
									Sct=union(Sn,Sc,FD);
								else
									Sn=Si;
									Sct=union(Sn,Sc,FD);
								end
							else
								Sct=Sn;
							end
						elseif strcmp(junction_method,'ignore')					
							if ii>1
								Sct=union(Sn,Sc,FD);
							else
								Sct=Sn;
							end
						end

						% Plot updated stream
						figure(f1)
						hold on
						p1=plot(Sn,'-b','LineWidth',2);
						hold off

						Sc=Sct;
					case 'No'
						str1='N';
						delete(p1);
					end
				else
					% Extract stream of interest
					Sn=modify(S,'downstreamto',IX);

					% Build composite stream network of picked streams
					if strcmp(junction_method,'check')
						if ii>1
							[IIXX,~,~,Si]=intersectlocs(Sc,Sn);
							if isempty(IIXX)
								Sn=Sn;
								Sct=union(Sn,Sc,FD);
							else
								Sn=Si;
								Sct=union(Sn,Sc,FD);
							end
						else
							Sct=Sn;
						end
					elseif strcmp(junction_method,'ignore')					
						if ii>1
							Sct=union(Sn,Sc,FD);
						else
							Sct=Sn;
						end
					end

					figure(f1)
					hold on
					p1=plot(Sn,'-b','LineWidth',2);
					hold off

					qa1=questdlg('Is this the stream segment you wanted?','Stream Selection','Yes','No','Yes');
					switch qa1
					case 'Yes'
						str1='Y';
						Sc=Sct;
					case 'No'
						str1='N';
						delete(p1);
					end					
				end
			end % End single channel select

			%% Extract threshold drainage area
			snchix=streampoi(Sn,'channelheads','ix');
			snda=DA.Z(snchix);

			%% Calculate chi and extract ksn data
			if strcmp(theta_method,'ref')
				C=ChiCalc(Sn,DEMc,A,1,ref_theta);
				ak=getnal(Sn,auto_ksn);
			elseif strcmp(theta_method,'auto')
				C=ChiCalc(Sn,DEMc,A,1);
				if redefine_thresh
					[auto_ksn]=KSN_Quick(DEM,A,Sn,C.mn);
				else
					[auto_ksn]=KSN_Quick(DEM,A,S,C.mn);
				end
				ak=getnal(Sn,auto_ksn);
			end

			%% Bin data
			[DAvg,KsnAvg]=BinAverage(Sn.distance,ak,smooth_distance);
			[~,CAvg]=BinAverage(C.distance,C.chi,smooth_distance);

			%% Begin fitting loop
			while strcmpi(str3,'R')
				%% Initiate figure to pick bounds
				f2=figure(2);
				set(f2,'Units','normalized','Position',[0.5 0.1 0.45 0.8],'renderer','painters');
				clf	

				%% Main swith for different pick methods
				if strcmp(pick_method,'chi')

					if display_slope_area
						[bs,ba,bc,bd,aa,ag,ad,ac]=sa(DEMc,Sn,A,C.chi,smooth_distance);

						ax4=subplot(4,1,4);
						hold on
						scatter(aa,ag,5,ac,'+');
						scatter(ba,bs,20,bc,'filled','MarkerEdgeColor','k');
						xlabel('Log Area');
						ylabel('Log Gradient');
						title('Slope-Area');
						set(ax4,'YScale','log','XScale','log','XDir','reverse');
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax4);
					    end							
						hold off

						ax3=subplot(4,1,3);
						hold on
						pl1=plotdz(Sn,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
						pl2=plotdz(Sn,DEMc,'dunit','km','Color','k');
						pl3=scatter((C.distance)./1000,C.elev,5,C.chi,'filled');
						xlabel('Distance from Mouth (km)')
						ylabel('Elevation (m)')
						legend([pl1 pl2 pl3],'Unconditioned DEM','Conditioned DEM','\chi','location','best');
						title('Long Profile')
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax3);
					    end							
						hold off

						ax2=subplot(4,1,2);
						hold on
						scatter(CAvg,KsnAvg,20,CAvg,'filled','MarkerEdgeColor','k');
						xlabel('\chi')
						ylabel('Auto k_{sn}');
						title('\chi - Auto k_{sn}');
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax2);
					    end							
						hold off

						ax1=subplot(4,1,1);
						hold on
						plot(C.chi,C.elev,'-k');
						scatter(C.chi,C.elev,10,C.chi,'filled');
						xlabel('\chi')
						ylabel('Elevation (m)')
						title(['\chi - Z : \theta = ' num2str(C.mn) ' : Pick Segments - Press Enter When Done'],'Color','r')
						ax1.XColor='Red';
						ax1.YColor='Red';
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax1);
					    end							
						hold off

						linkaxes([ax1,ax2],'x');
						colormap(ax1,'jet'); colormap(ax2,'jet'); colormap(ax3,'jet'); colormap(ax4,'jet');

					else
						ax3=subplot(3,1,3);
						hold on
						pl1=plotdz(Sn,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
						pl2=plotdz(Sn,DEMc,'dunit','km','Color','k');
						pl3=scatter((C.distance)./1000,C.elev,5,C.chi,'filled');
						xlabel('Distance from Mouth (km)')
						ylabel('Elevation (m)')
						legend([pl1 pl2 pl3],'Unconditioned DEM','Conditioned DEM','\chi','location','best');
						title('Long Profile')
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax3);
					    end							
						hold off

						ax2=subplot(3,1,2);
						hold on
						scatter(CAvg,KsnAvg,20,CAvg,'filled','MarkerEdgeColor','k');
						xlabel('\chi')
						ylabel('Auto k_{sn}');
						title('\chi - Auto k_{sn}');
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax2);
					    end							
						hold off

						ax1=subplot(3,1,1);
						hold on
						plot(C.chi,C.elev,'-k');
						scatter(C.chi,C.elev,10,C.chi,'filled');
						xlabel('\chi')
						ylabel('Elevation (m)')
						title(['\chi - Z : \theta = ' num2str(C.mn) ' : Pick Segments - Press Enter When Done'],'Color','r')
						ax1.XColor='Red';
						ax1.YColor='Red';
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax1);
					    end							
						hold off

						linkaxes([ax1,ax2],'x');
						colormap(ax1,'jet'); colormap(ax2,'jet'); colormap(ax3,'jet');
					end	

					[cv,~,bttn]=ginput;
					% Determine if there any non bound knicks
					bttn_idx=bttn~=1;
					if any(bttn_idx)
						% Parse out non bounds
						cv_kn=cv(bttn_idx);
						cv(bttn_idx)=[];
						% Convert to indices
						rc=C.chi; rx=C.x; ry=C.y;
						kn_ix=zeros(numel(cv_kn),1); 
						for jj=1:numel(cv_kn);
							chidist=sqrt(sum(bsxfun(@minus, rc, cv_kn(jj)).^2,2));
							[~,knbix]=min(chidist);
							knbx=rx(knbix);
							knby=ry(knbix);
							kn_ix(jj)=coord2ind(DEM,knbx,knby);
						end
					else
						kn_ix=NaN;
					end

					if isempty(cv)
						if strcmp(theta_method,'ref')
							Cbf=ChiCalc(Sn,DEMc,A,1);
							ksn_list=[C.x C.y C.area ones(numel(C.x),1)*C.ks ones(numel(C.x),1)*C.ks_neg ones(numel(C.x),1)*C.ks_pos ...
								ones(numel(C.x),1)*C.mn ones(numel(C.x),1)*Cbf.mn ones(numel(C.x),1)*snda];
						elseif strcmp(theta_method,'auto')
							ksn_list=[C.x C.y C.area ones(numel(C.x),1)*C.ks ones(numel(C.x),1)*C.ks_neg ones(numel(C.x),1)*C.ks_pos ...
								ones(numel(C.x),1)*C.mn ones(numel(C.x),1)*C.mn ones(numel(C.x),1)*snda];
						end

						% Determine where ksn value fits into color scale and plot
						ksn_val=C.ks;

						if ksn_val > mksn;
							figure(f1)
							hold on
							plot(Sn,'Color',KSN_col(end,:),'LineWidth',2);
							hold off
						else
							edges=linspace(0,mksn,10);
							n=histc(ksn_val,edges);
							figure(f1)
							hold on
							plot(Sn,'Color',KSN_col(logical(n),:),'LineWidth',2);
							hold off	
						end

						[~,lbix]=min(C.chi);
						elbl=C.elev(lbix);

						if display_slope_area
							figure(f2)
							subplot(4,1,1);
							hold on
							plot(C.chi,C.pred+elbl,'-k','LineWidth',2);
							hold off

							subplot(4,1,3);
							hold on
							pl4=plot((C.distance)/1000,C.pred+elbl,'-k','LineWidth',2);
							legend([pl1 pl2 pl3 pl4],'Unconditioned DEM','Conditioned DEM','\chi','Segment Fit','location','best');
							hold off

							subplot(4,1,4);
							hold on
							plot(C.area,ksn_val.*C.area.^(-C.mn),'-k','LineWidth',2);
							hold off

						else
							figure(f2)
							subplot(3,1,1);
							hold on
							plot(C.chi,C.pred+elbl,'-k','LineWidth',2);
							hold off

							subplot(3,1,3);
							hold on
							pl4=plot((C.distance)/1000,C.pred+elbl,'-k','LineWidth',2);
							legend([pl1 pl2 pl3 pl4],'Unconditioned DEM','Conditioned DEM','\chi','Segment Fit','location','best');
							hold off
						end

						res_list=[C.chi C.res];
						bnd_ix=NaN;;

					else
						% Sort knickpoint list and construct bounds list
						cvs=sortrows(cv);
						bnds=vertcat(0,cvs,C.chi(1));

						num_bnds=numel(bnds);
						rc=C.chi;
						rx=C.x;
						ry=C.y;
						rd=C.distance;
						for jj=1:num_bnds-1
							% Extract bounds
							lb=bnds(jj);
							rb=bnds(jj+1);

							% Clip out stream segment
							lb_chidist=sqrt(sum(bsxfun(@minus, rc, lb).^2,2));
							rb_chidist=sqrt(sum(bsxfun(@minus, rc, rb).^2,2));

							[~,lbix]=min(lb_chidist);
							[~,rbix]=min(rb_chidist);

							lbx=rx(lbix);
							lby=ry(lbix);

							rbx=rx(rbix);
							rby=ry(rbix);

							lix=coord2ind(DEM,lbx,lby);
							rix=coord2ind(DEM,rbx,rby);

							Seg=modify(Sn,'downstreamto',rix);
							Seg=modify(Seg,'upstreamto',lix);

							%Remake stream with downstream bound node added back in
							WSEG=GRIDobj(DEM,'logical');
							WSEG.Z(Seg.IXgrid)=true;
							WSEG.Z(lix)=true;
							% Add back in upstream node if it's the end of the stream
							if jj==num_bnds-1
								WSEG.Z(rix)=true;
							end
							Seg=STREAMobj(FD,WSEG);

							% Check length of stream, if it's less than two nodes,
							% move down stream until it's greater than two nodes
							lbix_new=lbix+1;
							first_time=true;
							while numel(Seg.IXgrid)<=2
								if first_time
									wrn_mssg=['Segment ' num2str(jj) ' of chosen segments was too short, segment bound was expanded downstream'];
									wd=warndlg(wrn_mssg);
									uiwait(wd);
								end
								lix_new=coord2ind(DEM,rx(lbix_new),ry(lbix_new));
								WSEG.Z(lix_new)=true;
								Seg=STREAMobj(FD,WSEG);
								lbix_new=lbix_new+1;
								first_time=false;
							end

							% Construct bound list
							if jj<num_bnds-1
								bnd_ix(jj,1)=rix;
							end

							% Calculate chi to find ksn and bestfit concavity 
							if strcmp(theta_method,'ref')
								Cseg=ChiCalc(Seg,DEMc,A,1,ref_theta);
							elseif strcmp(theta_method,'auto')
								Cseg=ChiCalc(Seg,DEMc,A,1,C.mn);
							end
							Cbfseg=ChiCalc(Seg,DEMc,A,1);								

							% Determine where ksn value fits into color scale and plot
							ksn_val=Cseg.ks;

							if ksn_val > mksn;
								figure(f1)
								hold on
								plot(Seg,'Color',KSN_col(end,:),'LineWidth',2);
								hold off
							else
								edges=linspace(0,mksn,10);
								n=histc(ksn_val,edges);
								figure(f1)
								hold on
								plot(Seg,'Color',KSN_col(logical(n),:),'LineWidth',2);
								hold off	
							end

							ksn_nodes{jj,1}=[Cseg.x Cseg.y Cseg.area ones(numel(Cseg.x),1)*Cseg.ks ones(numel(Cseg.x),1)*Cseg.ks_neg ones(numel(Cseg.x),1)*Cseg.ks_pos ...
								ones(numel(Cseg.x),1)*Cseg.mn ones(numel(Cseg.x),1)*Cbfseg.mn ones(numel(Cseg.x),1)*snda];

							% Plot linear fits	
							rchi=rc(rb_chidist==min(rb_chidist));
							lchi=rc(lb_chidist==min(lb_chidist));
							segChi=linspace(lchi,rchi,numel(Cseg.chi));
							seg0Chi=linspace(0,max(Cseg.chi),numel(Cseg.chi));
							[~,lbix]=min(Cseg.chi);
							elbl=Cseg.elev(lbix);

							if display_slope_area
								figure(f2)
								subplot(4,1,1);
								hold on
								plot(segChi,(seg0Chi.*ksn_val)+elbl,'-k','LineWidth',2);
								hold off

								seg_st=rd(lb_chidist==min(lb_chidist));
								subplot(4,1,3);
								hold on
								pl4=plot((Cseg.distance+seg_st)/1000,(Cseg.pred)+elbl,'-k','LineWidth',2);
								legend([pl1 pl2 pl3 pl4],'Unconditioned DEM','Conditioned DEM','\chi','Segment Fit','location','best');
								hold off

								subplot(4,1,4);
								hold on
								plot(Cseg.area,ksn_val.*Cseg.area.^(-1*Cseg.mn),'-k','LineWidth',2);
								hold off
							else								
								figure(f2)
								subplot(3,1,1);
								hold on
								plot(segChi,(seg0Chi.*ksn_val)+elbl,'-k','LineWidth',2);
								hold off

								seg_st=rd(lb_chidist==min(lb_chidist));
								subplot(3,1,3);
								hold on
								pl4=plot((Cseg.distance+seg_st)/1000,(Cseg.pred)+elbl,'-k','LineWidth',2);
								legend([pl1 pl2 pl3 pl4],'Unconditioned DEM','Conditioned DEM','\chi','Segment Fit','location','best');
								hold off
							end

							res_list{jj,1}=[Cseg.chi+lchi Cseg.res];	

						end

						ksn_list=vertcat(ksn_nodes{:});
						res_list=vertcat(res_list{:});
					end

					%% Plot result figure
					if display_slope_area
						figure(f2)
						subplot(4,1,2)
						hold on
						plot(res_list(:,1),ksn_list(:,4),'-k','LineWidth',2);
						hold off
					else
						figure(f2)
						subplot(3,1,2)
						hold on
						plot(res_list(:,1),ksn_list(:,4),'-k','LineWidth',2);
						hold off
					end

					f3=figure(3);
					set(f3,'Units','normalized','Position',[0.05 0.1 0.45 0.8],'renderer','painters');

					clf
					sbplt2=subplot(2,1,2);
					hold on
					s1=scatter(CAvg,KsnAvg,20,'k','filled');
					p1=plot(res_list(:,1),ksn_list(:,4)-ksn_list(:,5),':k');
					plot(res_list(:,1),ksn_list(:,4)+ksn_list(:,6),':k');
					p2=plot(res_list(:,1),ksn_list(:,4),'-k','LineWidth',2);
					[ksn_vals,~,ksn_ix]=unique(ksn_list(:,4));
					chi_means=accumarray(ksn_ix,res_list(:,1),[],@(x) mean(x,'omitnan'));
					for kk=1:numel(ksn_vals)
						text(chi_means(kk),ksn_vals(kk),['k_{sn} = ' num2str(ksn_vals(kk))],...
							'VerticalAlignment','bottom','HorizontalAlignment','center');
					end
					xlabel('\chi')
					ylabel('k_{sn}')
					legend([s1 p1 p2],{'Auto k_{sn}','k_{sn} uncertainty','k_{sn} of fit segments'},'location','best');
					title('k_{sn} - \chi')
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(sbplt2);
				    end						
					hold off

					sbplt1=subplot(2,1,1);
					hold on
					plot([min(res_list(:,1)) max(res_list(:,1))],[0 0],'-k');
					scatter(res_list(:,1),res_list(:,2),10,'k','filled');
					xlabel('\chi')
					ylabel('Residual (m)')
					title('Residual on k_{sn} fit')
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(sbplt1);
				    end						
					hold off				


				elseif strcmp(pick_method,'stream')	

					if display_slope_area
						[bs,ba,bc,bd,bk,aa,ag,ad,ac]=sa_ksn(DEMc,Sn,A,C.chi,ak,smooth_distance);	

						ax4=subplot(4,1,4);
						hold on
						scatter(aa,ag,5,ad./1000,'+');
						scatter(ba,bs,20,bd./1000,'filled','MarkerEdgeColor','k');
						xlabel('Log Area');
						ylabel('Log Gradient');
						title('Slope-Area');
						set(ax4,'XScale','log','YScale','log','XDir','reverse');
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax4);
					    end							
						hold off

						ax1=subplot(4,1,1);
						hold on
						plot(C.chi,C.elev,'Color',[0.5 0.5 0.5]);
						scatter(C.chi,C.elev,10,(C.distance)./1000,'filled');
						xlabel('\chi')
						ylabel('Elevation (m)')
						title(['\chi - Z : \theta = ' num2str(C.mn)])
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax1);
					    end							
						hold off

						ax2=subplot(4,1,2);
						hold on
						scatter(DAvg./1000,KsnAvg,20,DAvg./1000,'filled','MarkerEdgeColor','k');
						xlabel('Distance (km)')
						ylabel('Auto k_{sn}');
						title('Distance - Auto k_{sn}');
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax2);
					    end							
						hold off

						ax3=subplot(4,1,3);
						hold on
						pl1=plotdz(Sn,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
						pl2=plotdz(Sn,DEMc,'dunit','km','Color','k');
						pl3=scatter((C.distance)./1000,C.elev,5,(C.distance)./1000,'filled');
						xlabel('Distance from Mouth (km)')
						ylabel('Elevation (m)')
						legend([pl1 pl2 pl3],'Unconditioned DEM','Conditioned DEM','Stream Distance','location','best');
						title('Long Profile : Pick Segments - Press Enter When Done','Color','r')
						ax3.XColor='Red';
						ax3.YColor='Red';
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax3);
					    end							
						hold off

						linkaxes([ax2,ax3],'x');
						colormap(ax1,'jet'); colormap(ax2,'jet'); colormap(ax3,'jet'); colormap(ax4,'jet')
					else
						ax1=subplot(3,1,1);
						hold on
						plot(C.chi,C.elev,'Color',[0.5 0.5 0.5]);
						scatter(C.chi,C.elev,10,(C.distance)./1000,'filled');
						xlabel('\chi')
						ylabel('Elevation (m)')
						title(['\chi - Z : \theta = ' num2str(C.mn)])
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax1);
					    end							
						hold off

						ax2=subplot(3,1,2);
						hold on
						scatter(DAvg./1000,KsnAvg,20,DAvg./1000,'filled','MarkerEdgeColor','k');
						xlabel('Distance (km)')
						ylabel('Auto k_{sn}');
						title('Distance - Auto k_{sn}');
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax2);
					    end							
						hold off

						ax3=subplot(3,1,3);
						hold on
						pl1=plotdz(Sn,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
						pl2=plotdz(Sn,DEMc,'dunit','km','Color','k');
						pl3=scatter((C.distance)./1000,C.elev,5,(C.distance)./1000,'filled');
						xlabel('Distance from Mouth (km)')
						ylabel('Elevation (m)')
						legend([pl1 pl2 pl3],'Unconditioned DEM','Conditioned DEM','Stream Distance','location','best');
						title('Long Profile : Pick Segments - Press Enter When Done','Color','r')
						ax3.XColor='Red';
						ax3.YColor='Red';
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax3);
					    end							
						hold off

						linkaxes([ax2,ax3],'x');
						colormap(ax1,'jet'); colormap(ax2,'jet'); colormap(ax3,'jet');	
					end	

					[d,~,bttn]=ginput;
					d=d*1000; % convert back to meters;

					% Determine if there any non bound knicks
					bttn_idx=bttn~=1;
					if any(bttn_idx)
						% Parse out non bounds
						d_kn=d(bttn_idx);
						d(bttn_idx)=[];
						% Convert to indices
						rd=C.distance; rx=C.x; ry=C.y;
						kn_ix=zeros(numel(d_kn),1); 
						for jj=1:numel(d_kn);
							distdist=sqrt(sum(bsxfun(@minus, rd, d_kn(jj)).^2,2));
							[~,knbix]=min(distdist);
							knbx=rx(knbix);
							knby=ry(knbix);
							kn_ix(jj)=coord2ind(DEM,knbx,knby);
						end
					else
						kn_ix=NaN;
					end

					if isempty(d)
						if strcmp(theta_method,'ref')
							Cbf=ChiCalc(Sn,DEMc,A,1);
							ksn_list=[C.x C.y C.area ones(numel(C.x),1)*C.ks ones(numel(C.x),1)*C.ks_neg ones(numel(C.x),1)*C.ks_pos ...
								ones(numel(C.x),1)*C.mn ones(numel(C.x),1)*Cbf.mn ones(numel(C.x),1)*snda];
						elseif strcmp(theta_method,'auto')
							ksn_list=[C.x C.y C.area ones(numel(C.x),1)*C.ks ones(numel(C.x),1)*C.ks_neg ones(numel(C.x),1)*C.ks_pos ...
								ones(numel(C.x),1)*C.mn ones(numel(C.x),1)*C.mn ones(numel(C.x),1)*snda];
						end

						% Determine where ksn value fits into color scale and plot
						ksn_val=C.ks;

						if ksn_val > mksn;
							figure(f1)
							hold on
							plot(Sn,'Color',KSN_col(end,:),'LineWidth',2);
							hold off
						else
							edges=linspace(0,mksn,10);
							n=histc(ksn_val,edges);
							figure(f1)
							hold on
							plot(Sn,'Color',KSN_col(logical(n),:),'LineWidth',2);
							hold off	
						end

						[~,lbix]=min(C.chi);
						elbl=C.elev(lbix);
						if display_slope_area
							figure(f2)
							subplot(4,1,1);
							hold on
							plot(C.chi,C.pred+elbl,'-k','LineWidth',2);
							hold off

							subplot(4,1,3);
							hold on
							pl4=plot((C.distance)./1000,C.pred+elbl,'-k','LineWidth',2);
							legend([pl1 pl2 pl3 pl4],'Unconditioned DEM','Conditioned DEM','Stream Distance','Segment Fit','location','best');
							hold off

							subplot(4,1,4);
							hold on
							plot(C.area,ksn_val.*C.area.^(-C.mn),'-k','LineWidth',2);
							hold off							
						else
							figure(f2)
							subplot(3,1,1);
							hold on
							plot(C.chi,C.pred+elbl,'-k','LineWidth',2);
							hold off

							subplot(3,1,3);
							hold on
							pl4=plot((C.distance)./1000,C.pred+elbl,'-k','LineWidth',2);
							legend([pl1 pl2 pl3 pl4],'Unconditioned DEM','Conditioned DEM','Stream Distance','Segment Fit','location','best');
							hold off
						end

						res_list=[C.chi C.res];
						bnd_ix=NaN;;

					else
						% Sort knickpoint list and construct bounds list
						ds=sortrows(d);
						bnds=vertcat(0,ds,max(Sn.distance));

						num_bnds=numel(bnds);
						rd=C.distance; 
						rx=C.x;
						ry=C.y; 
						rc=C.chi;
						for jj=1:num_bnds-1
							% Extract bounds
							lb=bnds(jj);
							rb=bnds(jj+1);

							% Clip out stream segment
							lb_dist=sqrt(sum(bsxfun(@minus, rd, lb).^2,2));
							rb_dist=sqrt(sum(bsxfun(@minus, rd, rb).^2,2));

							[~,lbix]=min(lb_dist);
							[~,rbix]=min(rb_dist);

							lbx=rx(lbix);
							lby=ry(lbix);

							rbx=rx(rbix);
							rby=ry(rbix);	

							lix=coord2ind(DEM,lbx,lby);	
							rix=coord2ind(DEM,rbx,rby);

							Seg=modify(Sn,'downstreamto',rix);
							Seg=modify(Seg,'upstreamto',lix);

							%Remake stream with downstream bound node added back in
							WSEG=GRIDobj(DEM,'logical');
							WSEG.Z(Seg.IXgrid)=true;
							WSEG.Z(lix)=true;
							% Add back in upstream node if it's the end of the stream
							if jj==num_bnds-1
								WSEG.Z(rix)=true;
							end
							Seg=STREAMobj(FD,WSEG);

							% Check length of stream, if it's less than two nodes,
							% move down stream until it's greater than two nodes
							lbix_new=lbix+1;
							first_time=true;
							while numel(Seg.IXgrid)<=2
								if first_time
									wrn_mssg=['Segment ' num2str(jj) ' of chosen segments was too short, segment bound was expanded downstream'];
									wd=warndlg(wrn_mssg);
									uiwait(wd);
								end
								lix_new=coord2ind(DEM,rx(lbix_new),ry(lbix_new));
								WSEG.Z(lix_new)=true;
								Seg=STREAMobj(FD,WSEG);
								lbix_new=lbix_new+1;
								first_time=false;
							end

							% Construct bound list
							if jj<num_bnds-1
								bnd_ix(jj,1)=rix;
							end	

							% Calculate chi to find ksn and bestfit concavity 
							if strcmp(theta_method,'ref')
								Cseg=ChiCalc(Seg,DEMc,A,1,ref_theta);
							elseif strcmp(theta_method,'auto')
								Cseg=ChiCalc(Seg,DEMc,A,1,C.mn);
							end
							Cbfseg=ChiCalc(Seg,DEMc,A,1);

							% Determine where ksn value fits into color scale and plot
							ksn_val=Cseg.ks;

							if ksn_val > mksn;
								figure(f1)
								hold on
								plot(Seg,'Color',KSN_col(end,:),'LineWidth',2);
								hold off
							else
								edges=linspace(0,mksn,10);
								n=histc(ksn_val,edges);
								figure(f1)
								hold on
								plot(Seg,'Color',KSN_col(logical(n),:),'LineWidth',2);
								hold off	
							end

							ksn_nodes{jj,1}=[Cseg.x Cseg.y Cseg.area ones(numel(Cseg.x),1)*Cseg.ks ones(numel(Cseg.x),1)*Cseg.ks_neg ones(numel(Cseg.x),1)*Cseg.ks_pos ...
								ones(numel(Cseg.x),1)*Cseg.mn ones(numel(Cseg.x),1)*Cbfseg.mn ones(numel(Cseg.x),1)*snda];

							% Plot linear fits	
							rchi=rc(rb_dist==min(rb_dist));
							lchi=rc(lb_dist==min(lb_dist));
							ld=rd(lb_dist==min(lb_dist));
							segChi=linspace(lchi,rchi,numel(Cseg.chi));
							seg0Chi=linspace(0,max(Cseg.chi),numel(Cseg.chi));
							[~,lbix]=min(Cseg.chi);
							elbl=Cseg.elev(lbix);
							if display_slope_area
								figure(f2)
								subplot(4,1,1);
								hold on
								plot(segChi,(seg0Chi.*ksn_val)+elbl,'-k','LineWidth',2);
								hold off

								seg_st=rd(lb_dist==min(lb_dist));
								subplot(4,1,3);
								hold on
								pl4=plot((Cseg.distance+seg_st)/1000,(Cseg.pred)+elbl,'-k','LineWidth',2);
								legend([pl1 pl2 pl3 pl4],'Unconditioned DEM','Conditioned DEM','Stream Distance','Segment Fit','location','best');
								hold off

								subplot(4,1,4);
								hold on
								plot(Cseg.area,ksn_val.*Cseg.area.^(-Cseg.mn),'-k','LineWidth',2);
								hold off									
							else
								figure(f2)
								subplot(3,1,1);
								hold on
								plot(segChi,(seg0Chi.*ksn_val)+elbl,'-k','LineWidth',2);
								hold off

								seg_st=rd(lb_dist==min(lb_dist));
								subplot(3,1,3);
								hold on
								pl4=plot((Cseg.distance+seg_st)/1000,(Cseg.pred)+elbl,'-k','LineWidth',2);
								legend([pl1 pl2 pl3 pl4],'Unconditioned DEM','Conditioned DEM','Stream Distance','Segment Fit','location','best');
								hold off	
							end						

							res_list{jj,1}=[Cseg.distance+ld Cseg.res];
						end

						ksn_list=vertcat(ksn_nodes{:});
						res_list=vertcat(res_list{:});
					end

					%% Plot fit result
					if display_slope_area
						figure(f2)
						subplot(4,1,2)
						hold on
						plot(res_list(:,1)./1000,ksn_list(:,4),'-k','LineWidth',2);
						hold off
					else
						figure(f2)
						subplot(3,1,2)
						hold on
						plot(res_list(:,1)./1000,ksn_list(:,4),'-k','LineWidth',2);
						hold off
					end

					f3=figure(3);
					set(f3,'Units','normalized','Position',[0.05 0.1 0.45 0.8],'renderer','painters');

					clf
					sbplt1=subplot(2,1,2);
					hold on
					s1=scatter(DAvg./1000,KsnAvg,20,'k','filled');
					p1=plot(res_list(:,1)./1000,ksn_list(:,4)-ksn_list(:,5),':k');
					plot(res_list(:,1)./1000,ksn_list(:,4)+ksn_list(:,6),':k');
					p2=plot(res_list(:,1)./1000,ksn_list(:,4),'-k','LineWidth',2);
					[ksn_vals,~,ksn_ix]=unique(ksn_list(:,4));
					d_means=accumarray(ksn_ix,res_list(:,1),[],@(x) mean(x,'omitnan'));
					for kk=1:numel(ksn_vals)
						text(d_means(kk)./1000,ksn_vals(kk),['k_{sn} = ' num2str(ksn_vals(kk))],...
							'VerticalAlignment','bottom','HorizontalAlignment','center');
					end
					xlabel('Distance (km)')
					ylabel('k_{sn}')
					title('k_{sn} - Distance')
					legend([s1 p1 p2],{'Auto k_{sn}','k_{sn} uncertainty','k_{sn} of fit segments'},'location','best');
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(sbplt1);
				    end						
					hold off

					sbplt2=subplot(2,1,1);
					hold on
					plot([min(res_list(:,1)./1000) max(res_list(:,1)./1000)],[0 0],'-k');
					scatter(res_list(:,1)./1000,res_list(:,2),10,'k','filled');
					xlabel('Distance (km)')
					ylabel('Residual (m)')
					title('Residual on k_{sn} fit')
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(sbplt2);
				    end						
					hold off

				elseif strcmp(pick_method,'slope_area')

					[bs,ba,bc,bd,bk,aa,ag,ad,ac]=sa_ksn(DEMc,Sn,A,C.chi,ak,smooth_distance);

					ax3=subplot(4,1,3);
					hold on
					pl1=plotdz(Sn,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
					pl2=plotdz(Sn,DEMc,'dunit','km','Color','k');
					pl3=scatter((C.distance)./1000,C.elev,5,log10(C.area),'filled');
					xlabel('Distance from Mouth (km)')
					ylabel('Elevation (m)')
					legend([pl1 pl2 pl3],'Unconditioned DEM','Conditioned DEM','Log Drainage Area','location','best');
					title('Long Profile')
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(ax3);
				    end						
					hold off

					ax2=subplot(4,1,2);
					hold on
					scatter(ba,bk,20,log10(ba),'filled','MarkerEdgeColor','k');
					xlabel('Log Area')
					ylabel('Auto k_{sn}');
					title('Log Area - Auto k_{sn}');
					set(ax2,'XScale','log','XDir','reverse');
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(ax2);
				    end						
					hold off

					ax1=subplot(4,1,1);
					hold on
					plot(C.chi,C.elev,'-k');
					scatter(C.chi,C.elev,10,log10(C.area),'filled');
					xlabel('\chi')
					ylabel('Elevation (m)')
					title('\chi - Z')
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(ax1);
				    end						
					hold off

					ax4=subplot(4,1,4);
					hold on
					scatter(aa,ag,5,log10(aa),'+');
					scatter(ba,bs,20,log10(ba),'filled','MarkerEdgeColor','k');
					xlabel('Log Area');
					ylabel('Log Gradient');
					title(['Slope-Area: \theta = ' num2str(C.mn) ' : Pick Segments - Press Enter When Done'],'Color','r');
					set(ax4,'YScale','log','XScale','log','XDir','reverse');
					ax4.XColor='Red';
					ax4.YColor='Red';
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(ax4);
				    end						
					hold off

					linkaxes([ax4,ax2],'x');
					colormap(ax1,'jet'); colormap(ax2,'jet'); colormap(ax3,'jet'); colormap(ax4,'jet');

					[av,~,bttn]=ginput;
					% Determine if there any non bound knicks
					bttn_idx=bttn~=1;
					if any(bttn_idx)
						% Parse out non bounds
						av_kn=av(bttn_idx);
						av(bttn_idx)=[];
						% Convert to indices
						ra=C.area; rx=C.x; ry=C.y;
						kn_ix=zeros(numel(av_kn),1); 
						for jj=1:numel(av_kn);
							areadist=sqrt(sum(bsxfun(@minus, ra, av_kn(jj)).^2,2));
							[~,knbix]=min(areadist);
							knbx=rx(knbix);
							knby=ry(knbix);
							kn_ix(jj)=coord2ind(DEM,knbx,knby);
						end
					else
						kn_ix=NaN;
					end					

					if isempty(av)
						if strcmp(theta_method,'ref')
							Cbf=ChiCalc(Sn,DEMc,A,1);
							ksn_list=[C.x C.y C.area ones(numel(C.x),1)*C.ks ones(numel(C.x),1)*C.ks_neg ones(numel(C.x),1)*C.ks_pos ...
								ones(numel(C.x),1)*C.mn ones(numel(C.x),1)*Cbf.mn ones(numel(C.x),1)*snda];
						elseif strcmp(theta_method,'auto')
							ksn_list=[C.x C.y C.area ones(numel(C.x),1)*C.ks ones(numel(C.x),1)*C.ks_neg ones(numel(C.x),1)*C.ks_pos ...
								ones(numel(C.x),1)*C.mn ones(numel(C.x),1)*C.mn ones(numel(C.x),1)*snda];
						end

						% Determine where ksn value fits into color scale and plot
						ksn_val=C.ks;

						if ksn_val > mksn;
							figure(f1)
							hold on
							plot(Sn,'Color',KSN_col(end,:),'LineWidth',2);
							hold off
						else
							edges=linspace(0,mksn,10);
							n=histc(ksn_val,edges);
							figure(f1)
							hold on
							plot(Sn,'Color',KSN_col(logical(n),:),'LineWidth',2);
							hold off	
						end

						[~,lbix]=min(C.chi);
						elbl=C.elev(lbix);

						figure(f2)
						subplot(4,1,1);
						hold on
						plot(C.chi,C.pred+elbl,'-k','LineWidth',2);
						hold off

						subplot(4,1,3);
						hold on
						pl4=plot((C.distance)/1000,C.pred+elbl,'-k','LineWidth',2);
						legend([pl1 pl2 pl3 pl4],'Unconditioned DEM','Conditioned DEM','Log Drainage Area','Segment Fit','location','best');
						hold off

						subplot(4,1,4);
						hold on
						plot(C.area,ksn_val.*C.area.^(-C.mn),'-k','LineWidth',2);
						hold off

						res_list=[C.area C.res];
						bnd_ix=NaN;;

					else
						% Sort knickpoint list and construct bounds list
						avs=sortrows(av,'descend');
						bnds=vertcat(max(C.area,[],'omitnan'),avs,min(C.area,[],'omitnan'));

						num_bnds=numel(bnds);
						rc=C.chi;
						rx=C.x;
						ry=C.y;
						rd=C.distance;
						ra=C.area;
						for jj=1:num_bnds-1
							% Extract bounds
							lb=bnds(jj);
							rb=bnds(jj+1);

							% Clip out stream segment
							lb_dadist=sqrt(sum(bsxfun(@minus, ra, lb).^2,2));
							rb_dadist=sqrt(sum(bsxfun(@minus, ra, rb).^2,2));

							[~,lbix]=min(lb_dadist);
							[~,rbix]=min(rb_dadist);

							lbx=rx(lbix);
							lby=ry(lbix);

							rbx=rx(rbix);
							rby=ry(rbix);	

							lix=coord2ind(DEM,lbx,lby);
							rix=coord2ind(DEM,rbx,rby);

							Seg=modify(Sn,'downstreamto',rix);
							Seg=modify(Seg,'upstreamto',lix);

							%Remake stream with downstream bound node added back in
							WSEG=GRIDobj(DEM,'logical');
							WSEG.Z(Seg.IXgrid)=true;
							WSEG.Z(lix)=true;
							% Add back in upstream node if it's the end of the stream
							if jj==num_bnds-1
								WSEG.Z(rix)=true;
							end
							Seg=STREAMobj(FD,WSEG);

							% Check length of stream, if it's less than two nodes,
							% move down stream until it's greater than two nodes
							lbix_new=lbix+1;
							first_time=true;
							while numel(Seg.IXgrid)<=2
								if first_time
									wrn_mssg=['Segment ' num2str(jj) ' of chosen segments was too short, segment bound was expanded downstream'];
									wd=warndlg(wrn_mssg);
									uiwait(wd);
								end
								lix_new=coord2ind(DEM,rx(lbix_new),ry(lbix_new));
								WSEG.Z(lix_new)=true;
								Seg=STREAMobj(FD,WSEG);
								lbix_new=lbix_new+1;
								first_time=false;
							end							

							% Construct bound list
							if jj<num_bnds-1
								bnd_ix(jj,1)=rix;
							end

							% Calculate chi to find ksn and bestfit concavity 
							if strcmp(theta_method,'ref')
								Cseg=ChiCalc(Seg,DEMc,A,1,ref_theta);
							elseif strcmp(theta_method,'auto')
								Cseg=ChiCalc(Seg,DEMc,A,1,C.mn);
							end
							Cbfseg=ChiCalc(Seg,DEMc,A,1);								

							% Determine where ksn value fits into color scale and plot
							ksn_val=Cseg.ks;

							if ksn_val > mksn;
								figure(f1)
								hold on
								plot(Seg,'Color',KSN_col(end,:),'LineWidth',2);
								hold off
							else
								edges=linspace(0,mksn,10);
								n=histc(ksn_val,edges);
								figure(f1)
								hold on
								plot(Seg,'Color',KSN_col(logical(n),:),'LineWidth',2);
								hold off	
							end

							ksn_nodes{jj,1}=[Cseg.x Cseg.y Cseg.area ones(numel(Cseg.x),1)*Cseg.ks ones(numel(Cseg.x),1)*Cseg.ks_neg ones(numel(Cseg.x),1)*Cseg.ks_pos ...
								ones(numel(Cseg.x),1)*Cseg.mn ones(numel(Cseg.x),1)*Cbfseg.mn ones(numel(Cseg.x),1)*snda];

							% Plot linear fits	
							rchi=rc(rb_dadist==min(rb_dadist));
							lchi=rc(lb_dadist==min(lb_dadist));
							segChi=linspace(lchi,rchi,numel(Cseg.chi));
							seg0Chi=linspace(0,max(Cseg.chi),numel(Cseg.chi));
							[~,lbix]=min(Cseg.chi);
							elbl=Cseg.elev(lbix);

							figure(f2)
							subplot(4,1,1);
							hold on
							plot(segChi,(seg0Chi.*ksn_val)+elbl,'-k','LineWidth',2);
							hold off

							seg_st=rd(lb_dadist==min(lb_dadist));
							subplot(4,1,3);
							hold on
							pl4=plot((Cseg.distance+seg_st)/1000,(Cseg.pred)+elbl,'-k','LineWidth',2);
							legend([pl1 pl2 pl3 pl4],'Unconditioned DEM','Conditioned DEM','Log Drainage Area','Segment Fit','location','best');
							hold off

							subplot(4,1,4);
							hold on
							plot(Cseg.area,ksn_val.*Cseg.area.^(-1*Cseg.mn),'-k','LineWidth',2);
							hold off


							res_list{jj,1}=[Cseg.area Cseg.res];


						end

						ksn_list=vertcat(ksn_nodes{:});
						res_list=vertcat(res_list{:});
					end


					%% Plot result figure
					if display_slope_area
						figure(f2)
						subplot(4,1,2)
						hold on
						plot(res_list(:,1),ksn_list(:,4),'-k','LineWidth',2);
						hold off
					else
						figure(f2)
						subplot(3,1,2)
						hold on
						plot(res_list(:,1),ksn_list(:,4),'-k','LineWidth',2);
						hold off
					end

					f3=figure(3);
					set(f3,'Units','normalized','Position',[0.05 0.1 0.45 0.8],'renderer','painters');
					clf
					ax31=subplot(2,1,2);
					hold on
					s1=scatter(log10(ba),bk,20,'k','filled');
					p1=plot(log10(res_list(:,1)),ksn_list(:,4)-ksn_list(:,5),':k');
					plot(log10(res_list(:,1)),ksn_list(:,4)+ksn_list(:,6),':k');
					p2=plot(log10(res_list(:,1)),ksn_list(:,4),'-k','LineWidth',2);
					[ksn_vals,~,ksn_ix]=unique(ksn_list(:,4));
					a_means=accumarray(ksn_ix,res_list(:,1),[],@(x) mean(x,'omitnan'));
					for kk=1:numel(ksn_vals)
						text(log10(a_means(kk)),ksn_vals(kk),['k_{sn} = ' num2str(ksn_vals(kk))],...
							'VerticalAlignment','bottom','HorizontalAlignment','center');
					end
					xlabel('Log Area')
					ylabel('k_{sn}')
					title('k_{sn} - Log Area')
					legend([s1 p1 p2],{'Auto k_{sn}','k_{sn} uncertainty','k_{sn} of fit segments'},'location','best');
					set(ax31,'XScale','log','XDir','reverse');
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(ax31);
				    end						
					hold off

					ax32=subplot(2,1,1);
					hold on
					plot([log10(min(res_list(:,1))) log10(max(res_list(:,1)))],[0 0],'-k');
					scatter(log10(res_list(:,1)),res_list(:,2),10,'k','filled');
					xlabel('Log Area')
					ylabel('Residual (m)')
					title('Residual on k_{sn} fit')
					set(ax32,'XScale','log','XDir','Reverse');
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(ax32);
				    end						
					hold off

				end
			
				qa2=questdlg('What would you like to do?','Stream Fitting','Stop Picking','Redo Fit','Continue Picking','Continue Picking');

				switch qa2
				case 'Continue Picking'
					str2 = 'Y';
					str1 = 'N';
					str3 = 'C';
					% Add Gradient, Residual, and Riv Number to node list and clear NaNs
					kidx=isnan(ksn_list(:,1));
					ksn_list(kidx,:)=[];
					res_list(kidx,:)=[];
					gix=coord2ind(DEM,ksn_list(:,1),ksn_list(:,2));
					kmat=[ksn_list G.Z(gix) res_list(:,2) ones(size(ksn_list,1),1).*ii];
					ksn_master{ii,1}=kmat;
					bnd_master{ii,1}=bnd_ix;
					kn_master{ii,1}=kn_ix;
					res_master{ii,1}=res_list;
					count=ii;
					save(out_restart_name,'ksn_master','bnd_master','kn_master','res_master','Sc','count','-append');
					if save_figures
						f2_name=[shape_name '_stream_fits_' num2str(ii) '.pdf'];
						f3_name=[shape_name '_stream_rsds_' num2str(ii) '.pdf'];
						print(f2,f2_name,'-dpdf','-fillpage');
						print(f3,f3_name,'-dpdf','-fillpage');
					end
					clear ksn_list ksn_nodes res_list bnd_ix kn_ix kidx gix kmat;
					ii=ii+1;
				case 'Stop Picking'
					wtb=waitbar(0,'Cleaning up and generating outputs, do not close windows');
					str1 = 'N';
					str2 = 'N';
					str3 = 'C';
					% Add Gradient, Residual, and Riv Number to node list and clear NaNs
					kidx=isnan(ksn_list(:,1));
					ksn_list(kidx,:)=[];
					res_list(kidx,:)=[];
					gix=coord2ind(DEM,ksn_list(:,1),ksn_list(:,2));
					kmat=[ksn_list G.Z(gix) res_list(:,2) ones(size(ksn_list,1),1).*ii];
					ksn_master{ii,1}=kmat;
					bnd_master{ii,1}=bnd_ix;
					kn_master{ii,1}=kn_ix;
					res_master{ii,1}=res_list;
					count=ii;
					save(out_restart_name,'ksn_master','bnd_master','kn_master','res_master','Sc','count','-append');					
					if save_figures
						f2_name=[shape_name '_stream_fits_' num2str(ii) '.pdf'];
						f3_name=[shape_name '_stream_rsds_' num2str(ii) '.pdf'];
						print(f2,f2_name,'-dpdf','-fillpage');
						print(f3,f3_name,'-dpdf','-fillpage');
					end
					clear ksn_list ksn_nodes res_list bnd_ix kn_ix kidx gix kmat;
				case 'Redo Fit'
					str3 = 'R';
					str2 = 'Y';
					clear ksn_list ksn_nodes res_list bnd_ix kn_ix;
				end

				close figure 2
				close figure 3
			end
		end


	case 'preselected'
		% Initiate while loop value
		str1='R';
		break_flag=false;

		if strcmp(theta_method,'ref')
			% Autocalculate ksn for comparison purposes
			[auto_ksn]=KSN_Quick(DEM,A,S,ref_theta);
		end

		if rf & strcmp(r_type,'c') & strcmp(restart,'continue')
			load(out_mat_name,'count','ksn_master','bnd_master','kn_master','res_master','Sc');
			if count>=num_ch
				if isdeployed
					errordlg('Run appears to have already completed, cannot restart')
				end
				error('Run appears to have already completed, cannot restart');
			end
			rng=count+1:num_ch;
		elseif rf & strcmp(r_type,'c') & strcmp(restart,'skip')
			load(out_mat_name,'count','ksn_master','bnd_master','kn_master','res_master','Sc');
			if count>=num_ch
				if isdeployed
					errordlg('Run appears to have already completed, cannot restart')
				end
				error('Run appears to have already completed, cannot restart');
			end
			rng=count+2:num_ch;
		elseif rf & strcmp(r_type,'r') & strcmp(restart,'continue')
			load(out_restart_name,'count','ksn_master','bnd_master','kn_master','res_master','Sc');
			if count>=num_ch
				if isdeployed
					errordlg('Run appears to have already completed, cannot restart')
				end
				error('Run appears to have already completed, cannot restart');
			end
			rng=count+1:num_ch;
		elseif rf & strcmp(r_type,'r') & strcmp(restart,'skip')
			load(out_restart_name,'count','ksn_master','bnd_master','kn_master','res_master','Sc');
			if count>=num_ch
				if isdeployed
					errordlg('Run appears to have already completed, cannot restart')
				end
				error('Run appears to have already completed, cannot restart');
			end
			rng=count+2:num_ch;
		else
			rng=1:num_ch;
		end		

		for ii=rng

			while strcmpi(str1,'R');
	           
				chOI=s_ch(ii,:);

				% Build logical raster
				ix=coord2ind(DEM,chOI(:,1),chOI(:,2));
				IX=GRIDobj(DEM,'logical');
				IX.Z(ix)=true;

				% Extract stream of interest
				Sn=modify(S,'downstreamto',IX);

				if redefine_thresh
					[Sn]=RedefineThreshold(DEM,FD,A,Sn,FLUS,ref_theta,rd_pick_method,smooth_distance,ii,save_figures,shape_name);
					% Update DEMc
					if any(isnan(getnal(Sn,DEMc)));
						zc=mincosthydrocon(Sn,DEM,'interp',iv);
						DEMc.Z(Sn.IXgrid)=zc;
					end

					% Recalculate auto_ksn
					if strcmp(theta_method,'ref')
						[auto_ksn]=KSN_Quick(DEM,A,Sn,ref_theta);
					end
				end

				% Build composite stream network of picked streams
				if strcmp(junction_method,'check')
					if ii>1
						[IIXX,~,~,Si]=intersectlocs(Sc,Sn);
						if isempty(IIXX)
							Sn=Sn;
							Sct=union(Sn,Sc,FD);
						else
							Sn=Si;

							if ~isempty(min_channel_length) & max(Sn.distance)<min_channel_length
								disp(['Skipping channel head ' num2str(ii) ', stream segment is too short']);
								break;
							end

							Sct=union(Sn,Sc,FD);
						end
					else
						Sct=Sn;
					end
				elseif strcmp(junction_method,'ignore')
					if ii>1
						Sct=union(Sn,Sc,FD);
					else
						Sct=Sn;
					end
				end

				%% Extract threshold drainage area
				snchix=streampoi(Sn,'channelheads','ix');
				snda=DA.Z(snchix);

				%% Calculate chi and extract ksn data
				if strcmp(theta_method,'ref')
					C=ChiCalc(Sn,DEMc,A,1,ref_theta);
					ak=getnal(Sn,auto_ksn);
				elseif strcmp(theta_method,'auto')
					C=ChiCalc(Sn,DEMc,A,1);
					if redefine_thresh
						[auto_ksn]=KSN_Quick(DEM,A,Sn,C.mn);
					else
						[auto_ksn]=KSN_Quick(DEM,A,S,C.mn);
					end
					ak=getnal(Sn,auto_ksn);
				end

				%% Bin data
				[DAvg,KsnAvg]=BinAverage(Sn.distance,ak,smooth_distance);
				[~,CAvg]=BinAverage(C.distance,C.chi,smooth_distance);


				%% Initiate figure to pick bounds
				f2=figure(2);
				set(f2,'Units','normalized','Position',[0.5 0.1 0.45 0.8],'renderer','painters');
				clf	

				%% Main swith for different pick methods
				if strcmp(pick_method,'chi')

					if display_slope_area
						[bs,ba,bc,bd,aa,ag,ad,ac]=sa(DEMc,Sn,A,C.chi,smooth_distance);

						ax4=subplot(4,1,4);
						hold on
						scatter(aa,ag,5,ac,'+');
						scatter(ba,bs,20,bc,'filled','MarkerEdgeColor','k');
						xlabel('Log Area');
						ylabel('Log Gradient');
						title('Slope-Area');
						set(ax4,'YScale','log','XScale','log','XDir','reverse');
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax4);
					    end							
						hold off

						ax3=subplot(4,1,3);
						hold on
						pl1=plotdz(Sn,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
						pl2=plotdz(Sn,DEMc,'dunit','km','Color','k');
						pl3=scatter((C.distance)./1000,C.elev,5,C.chi,'filled');
						xlabel('Distance from Mouth (km)')
						ylabel('Elevation (m)')
						legend([pl1 pl2 pl3],'Unconditioned DEM','Conditioned DEM','\chi','location','best');
						title('Long Profile')
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax3);
					    end							
						hold off

						ax2=subplot(4,1,2);
						hold on
						scatter(CAvg,KsnAvg,20,CAvg,'filled','MarkerEdgeColor','k');
						xlabel('\chi')
						ylabel('Auto k_{sn}');
						title('\chi - Auto k_{sn}');
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax2);
					    end							
						hold off

						ax1=subplot(4,1,1);
						hold on
						plot(C.chi,C.elev,'-k');
						scatter(C.chi,C.elev,10,C.chi,'filled');
						xlabel('\chi')
						ylabel('Elevation (m)')
						title(['\chi - Z : \theta = ' num2str(C.mn) ' : Pick Segments - Press Enter When Done'],'Color','r')
						ax1.XColor='Red';
						ax1.YColor='Red';
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax1);
					    end							
						hold off

						linkaxes([ax1,ax2],'x');
						colormap(ax1,'jet'); colormap(ax2,'jet'); colormap(ax3,'jet'); colormap(ax4,'jet');

					else
						ax3=subplot(3,1,3);
						hold on
						pl1=plotdz(Sn,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
						pl2=plotdz(Sn,DEMc,'dunit','km','Color','k');
						pl3=scatter((C.distance)./1000,C.elev,5,C.chi,'filled');
						xlabel('Distance from Mouth (km)')
						ylabel('Elevation (m)')
						legend([pl1 pl2 pl3],'Unconditioned DEM','Conditioned DEM','\chi','location','best');
						title('Long Profile')
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax3);
					    end							
						hold off

						ax2=subplot(3,1,2);
						hold on
						scatter(CAvg,KsnAvg,20,CAvg,'filled','MarkerEdgeColor','k');
						xlabel('\chi')
						ylabel('Auto k_{sn}');
						title('\chi - Auto k_{sn}');
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax2);
					    end							
						hold off

						ax1=subplot(3,1,1);
						hold on
						plot(C.chi,C.elev,'-k');
						scatter(C.chi,C.elev,10,C.chi,'filled');
						xlabel('\chi')
						ylabel('Elevation (m)')
						title(['\chi - Z : \theta = ' num2str(C.mn) ' : Pick Segments - Press Enter When Done'],'Color','r')
						ax1.XColor='Red';
						ax1.YColor='Red';
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax1);
					    end							
						hold off

						linkaxes([ax1,ax2],'x');
						colormap(ax1,'jet'); colormap(ax2,'jet'); colormap(ax3,'jet');
					end	

					[cv,~,bttn]=ginput;
					% Determine if there any non bound knicks
					bttn_idx=bttn~=1;
					if any(bttn_idx)
						% Parse out non bounds
						cv_kn=cv(bttn_idx);
						cv(bttn_idx)=[];
						% Convert to indices
						rc=C.chi; rx=C.x; ry=C.y;
						kn_ix=zeros(numel(cv_kn),1); 
						for jj=1:numel(cv_kn);
							chidist=sqrt(sum(bsxfun(@minus, rc, cv_kn(jj)).^2,2));
							[~,knbix]=min(chidist);
							knbx=rx(knbix);
							knby=ry(knbix);
							kn_ix(jj)=coord2ind(DEM,knbx,knby);
						end
					else
						kn_ix=NaN;
					end

					if isempty(cv)
						if strcmp(theta_method,'ref')
							Cbf=ChiCalc(Sn,DEMc,A,1);
							ksn_list=[C.x C.y C.area ones(numel(C.x),1)*C.ks ones(numel(C.x),1)*C.ks_neg ones(numel(C.x),1)*C.ks_pos ...
								ones(numel(C.x),1)*C.mn ones(numel(C.x),1)*Cbf.mn ones(numel(C.x),1)*snda];
						elseif strcmp(theta_method,'auto')
							ksn_list=[C.x C.y C.area ones(numel(C.x),1)*C.ks ones(numel(C.x),1)*C.ks_neg ones(numel(C.x),1)*C.ks_pos ...
								ones(numel(C.x),1)*C.mn ones(numel(C.x),1)*C.mn ones(numel(C.x),1)*snda];
						end

						% Determine where ksn value fits into color scale and plot
						ksn_val=C.ks;

						[~,lbix]=min(C.chi);
						elbl=C.elev(lbix);
						if display_slope_area
							figure(f2)
							subplot(4,1,1);
							hold on
							plot(C.chi,C.pred+elbl,'-k','LineWidth',2);
							hold off

							subplot(4,1,3);
							hold on
							pl4=plot((C.distance)/1000,C.pred+elbl,'-k','LineWidth',2);
							legend([pl1 pl2 pl3 pl4],'Unconditioned DEM','Conditioned DEM','\chi','Segment Fit','location','best');
							hold off

							subplot(4,1,4);
							hold on
							plot(C.area,ksn_val.*C.area.^(-C.mn),'-k','LineWidth',2);
							hold off
						else
							figure(f2)
							subplot(3,1,1);
							hold on
							plot(C.chi,C.pred+elbl,'-k','LineWidth',2);
							hold off

							subplot(3,1,3);
							hold on
							pl4=plot((C.distance)/1000,C.pred+elbl,'-k','LineWidth',2);
							legend([pl1 pl2 pl3 pl4],'Unconditioned DEM','Conditioned DEM','\chi','Segment Fit','location','best');
							hold off
						end

						res_list=[C.chi C.res];
						bnd_ix=NaN;;

					else
						% Sort knickpoint list and construct bounds list
						cvs=sortrows(cv);
						bnds=vertcat(0,cvs,C.chi(1));

						num_bnds=numel(bnds);
						rc=C.chi;
						rx=C.x;
						ry=C.y;
						rd=C.distance;
						for jj=1:num_bnds-1
							% Extract bounds
							lb=bnds(jj);
							rb=bnds(jj+1);

							% Clip out stream segment
							lb_chidist=sqrt(sum(bsxfun(@minus, rc, lb).^2,2));
							rb_chidist=sqrt(sum(bsxfun(@minus, rc, rb).^2,2));

							[~,lbix]=min(lb_chidist);
							[~,rbix]=min(rb_chidist);

							lbx=rx(lbix);
							lby=ry(lbix);

							rbx=rx(rbix);
							rby=ry(rbix);

							lix=coord2ind(DEM,lbx,lby);
							rix=coord2ind(DEM,rbx,rby);

							Seg=modify(Sn,'downstreamto',rix);
							Seg=modify(Seg,'upstreamto',lix);

							%Remake stream with downstream bound node added back in
							WSEG=GRIDobj(DEM,'logical');
							WSEG.Z(Seg.IXgrid)=true;
							WSEG.Z(lix)=true;
							% Add back in upstream node if it's the end of the stream
							if jj==num_bnds-1
								WSEG.Z(rix)=true;
							end
							Seg=STREAMobj(FD,WSEG);

							% Check length of stream, if it's less than two nodes,
							% move down stream until it's greater than two nodes
							lbix_new=lbix+1;
							first_time=true;
							while numel(Seg.IXgrid)<=2
								if first_time
									wrn_mssg=['Segment ' num2str(jj) ' of chosen segments was too short, segment bound was expanded downstream'];
									wd=warndlg(wrn_mssg);
									uiwait(wd);
								end
								lix_new=coord2ind(DEM,rx(lbix_new),ry(lbix_new));
								WSEG.Z(lix_new)=true;
								Seg=STREAMobj(FD,WSEG);
								lbix_new=lbix_new+1;
								first_time=false;
							end

							% Construct bound list
							if jj<num_bnds-1
								bnd_ix(jj,1)=rix;
							end

							% Calculate chi to find ksn and bestfit concavity 
							if strcmp(theta_method,'ref')
								Cseg=ChiCalc(Seg,DEMc,A,1,ref_theta);
							elseif strcmp(theta_method,'auto')
								Cseg=ChiCalc(Seg,DEMc,A,1,C.mn);
							end
							Cbfseg=ChiCalc(Seg,DEMc,A,1);								

							% Determine where ksn value fits into color scale and plot
							ksn_val=Cseg.ks;

							ksn_nodes{jj,1}=[Cseg.x Cseg.y Cseg.area ones(numel(Cseg.x),1)*Cseg.ks ones(numel(Cseg.x),1)*Cseg.ks_neg ones(numel(Cseg.x),1)*Cseg.ks_pos ...
								ones(numel(Cseg.x),1)*Cseg.mn ones(numel(Cseg.x),1)*Cbfseg.mn ones(numel(Cseg.x),1)*snda];

							% Plot linear fits	
							rchi=rc(rb_chidist==min(rb_chidist));
							lchi=rc(lb_chidist==min(lb_chidist));
							segChi=linspace(lchi,rchi,numel(Cseg.chi));
							seg0Chi=linspace(0,max(Cseg.chi),numel(Cseg.chi));
							[~,lbix]=min(Cseg.chi);
							elbl=Cseg.elev(lbix);

							if display_slope_area
								figure(f2)
								subplot(4,1,1);
								hold on
								plot(segChi,(seg0Chi.*ksn_val)+elbl,'-k','LineWidth',2);
								hold off

								seg_st=rd(lb_chidist==min(lb_chidist));
								subplot(4,1,3);
								hold on
								pl4=plot((Cseg.distance+seg_st)/1000,(Cseg.pred)+elbl,'-k','LineWidth',2);
								legend([pl1 pl2 pl3 pl4],'Unconditioned DEM','Conditioned DEM','\chi','Segment Fit','location','best');
								hold off

								subplot(4,1,4);
								hold on
								plot(Cseg.area,ksn_val.*Cseg.area.^(-Cseg.mn),'-k','LineWidth',2);
								hold off
							else
								figure(f2)
								subplot(3,1,1);
								hold on
								plot(segChi,(seg0Chi.*ksn_val)+elbl,'-k','LineWidth',2);
								hold off

								seg_st=rd(lb_chidist==min(lb_chidist));
								subplot(3,1,3);
								hold on
								pl4=plot((Cseg.distance+seg_st)/1000,(Cseg.pred)+elbl,'-k','LineWidth',2);
								legend([pl1 pl2 pl3 pl4],'Unconditioned DEM','Conditioned DEM','\chi','Segment Fit','location','best');
								hold off
							end

							res_list{jj,1}=[Cseg.chi+lchi Cseg.res];	

						end

						ksn_list=vertcat(ksn_nodes{:});
						res_list=vertcat(res_list{:});
					end

					%% Plot result figure
					if display_slope_area
						figure(f2)
						subplot(4,1,2)
						hold on
						plot(res_list(:,1),ksn_list(:,4),'-k','LineWidth',2);
						hold off
					else
						figure(f2)
						subplot(3,1,2)
						hold on
						plot(res_list(:,1),ksn_list(:,4),'-k','LineWidth',2);
						hold off
					end

					f3=figure(3);
					set(f3,'Units','normalized','Position',[0.05 0.1 0.45 0.8],'renderer','painters');

					clf
					sbplt2=subplot(2,1,2);
					hold on
					s1=scatter(CAvg,KsnAvg,20,'k','filled');
					p1=plot(res_list(:,1),ksn_list(:,4)-ksn_list(:,5),':k');
					plot(res_list(:,1),ksn_list(:,4)+ksn_list(:,6),':k');
					p2=plot(res_list(:,1),ksn_list(:,4),'-k','LineWidth',2);
					[ksn_vals,~,ksn_ix]=unique(ksn_list(:,4));
					chi_means=accumarray(ksn_ix,res_list(:,1),[],@(x) mean(x,'omitnan'));
					for kk=1:numel(ksn_vals)
						text(chi_means(kk),ksn_vals(kk),['k_{sn} = ' num2str(ksn_vals(kk))],...
							'VerticalAlignment','bottom','HorizontalAlignment','center','Color','r');
					end
					xlabel('\chi')
					ylabel('k_{sn}')
					legend([s1 p1 p2],{'Auto k_{sn}','k_{sn} uncertainty','k_{sn} of fit segments'},'location','best');
					title('k_{sn} - \chi')
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(sbplt2);
				    end						
					hold off

					sbplt1=subplot(2,1,1);
					hold on
					plot([min(res_list(:,1)) max(res_list(:,1))],[0 0],'-k');
					scatter(res_list(:,1),res_list(:,2),10,'k','filled');
					xlabel('\chi')
					ylabel('Residual (m)')
					title('Residual on k_{sn} fit')
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(sbplt1);
				    end						
					hold off				

				elseif strcmp(pick_method,'stream')	

					if display_slope_area
						[bs,ba,bc,bd,aa,ag,ad,ac]=sa(DEMc,Sn,A,C.chi,smooth_distance);	

						ax4=subplot(4,1,4);
						hold on
						scatter(aa,ag,5,ad./1000,'+');
						scatter(ba,bs,20,bd./1000,'filled','MarkerEdgeColor','k');
						xlabel('Log Area');
						ylabel('Log Gradient');
						title('Slope-Area');
						set(ax4,'XScale','log','YScale','log','XDir','reverse');
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax4);
					    end							
						hold off

						ax1=subplot(4,1,1);
						hold on
						plot(C.chi,C.elev,'Color',[0.5 0.5 0.5]);
						scatter(C.chi,C.elev,10,(C.distance)./1000,'filled');
						xlabel('\chi')
						ylabel('Elevation (m)')
						title(['\chi - Z : \theta = ' num2str(C.mn)])
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax1);
					    end							
						hold off

						ax2=subplot(4,1,2);
						hold on
						scatter(DAvg./1000,KsnAvg,20,DAvg./1000,'filled','MarkerEdgeColor','k');
						xlabel('Distance (km)')
						ylabel('Auto k_{sn}');
						title('Distance - Auto k_{sn}');
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax2);
					    end							
						hold off

						ax3=subplot(4,1,3);
						hold on
						pl1=plotdz(Sn,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
						pl2=plotdz(Sn,DEMc,'dunit','km','Color','k');
						pl3=scatter((C.distance)./1000,C.elev,5,(C.distance)./1000,'filled');
						xlabel('Distance from Mouth (km)')
						ylabel('Elevation (m)')
						legend([pl1 pl2 pl3],'Unconditioned DEM','Conditioned DEM','Stream Distance','location','best');
						title('Long Profile : Pick Segments - Press Enter When Done','Color','r')
						ax3.XColor='Red';
						ax3.YColor='Red';
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax3);
					    end							
						hold off

						linkaxes([ax2,ax3],'x');
						colormap(ax1,'jet'); colormap(ax2,'jet'); colormap(ax3,'jet'); colormap(ax4,'jet')
					else
						ax1=subplot(3,1,1);
						hold on
						plot(C.chi,C.elev,'Color',[0.5 0.5 0.5]);
						scatter(C.chi,C.elev,10,(C.distance)./1000,'filled');
						xlabel('\chi')
						ylabel('Elevation (m)')
						title(['\chi - Z : \theta = ' num2str(C.mn)])
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax1);
					    end							
						hold off

						ax2=subplot(3,1,2);
						hold on
						scatter(DAvg./1000,KsnAvg,20,DAvg./1000,'filled','MarkerEdgeColor','k');
						xlabel('Distance (km)')
						ylabel('Auto k_{sn}');
						title('Distance - Auto k_{sn}');
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax2);
					    end							
						hold off

						ax3=subplot(3,1,3);
						hold on
						pl1=plotdz(Sn,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
						pl2=plotdz(Sn,DEMc,'dunit','km','Color','k');
						pl3=scatter((C.distance)./1000,C.elev,5,(C.distance)./1000,'filled');
						xlabel('Distance from Mouth (km)')
						ylabel('Elevation (m)')
						legend([pl1 pl2 pl3],'Unconditioned DEM','Conditioned DEM','Stream Distance','location','best');
						title('Long Profile : Pick Segments - Press Enter When Done','Color','r')
						ax3.XColor='Red';
						ax3.YColor='Red';
						if ~verLessThan('matlab','9.5')
					        disableDefaultInteractivity(ax3);
					    end							
						hold off

						linkaxes([ax2,ax3],'x');
						colormap(ax1,'jet'); colormap(ax2,'jet'); colormap(ax3,'jet');	
					end	

					linkaxes([ax2,ax3],'x');
					colormap(ax1,'jet'); colormap(ax2,'jet'); colormap(ax3,'jet');

					[d,~,bttn]=ginput;
					d=d*1000; % convert back to meters;

					% Determine if there any non bound knicks
					bttn_idx=bttn~=1;
					if any(bttn_idx)
						% Parse out non bounds
						d_kn=d(bttn_idx);
						d(bttn_idx)=[];
						% Convert to indices
						rd=C.distance; rx=C.x; ry=C.y;
						kn_ix=zeros(numel(d_kn),1); 
						for jj=1:numel(d_kn);
							distdist=sqrt(sum(bsxfun(@minus, rd, d_kn(jj)).^2,2));
							[~,knbix]=min(distdist);
							knbx=rx(knbix);
							knby=ry(knbix);
							kn_ix(jj)=coord2ind(DEM,knbx,knby);
						end
					else
						kn_ix=NaN;
					end

					if isempty(d)
						if strcmp(theta_method,'ref')
							Cbf=ChiCalc(Sn,DEMc,A,1);
							ksn_list=[C.x C.y C.area ones(numel(C.x),1)*C.ks ones(numel(C.x),1)*C.ks_neg ones(numel(C.x),1)*C.ks_pos ...
								ones(numel(C.x),1)*C.mn ones(numel(C.x),1)*Cbf.mn ones(numel(C.x),1)*snda];
						elseif strcmp(theta_method,'auto')
							ksn_list=[C.x C.y C.area ones(numel(C.x),1)*C.ks ones(numel(C.x),1)*C.ks_neg ones(numel(C.x),1)*C.ks_pos ...
								ones(numel(C.x),1)*C.mn ones(numel(C.x),1)*C.mn ones(numel(C.x),1)*snda];
						end

						% Determine where ksn value fits into color scale and plot
						ksn_val=C.ks;

						[~,lbix]=min(C.chi);
						elbl=C.elev(lbix);
						if display_slope_area
							figure(f2)
							subplot(4,1,1);
							hold on
							plot(C.chi,C.pred+elbl,'-k','LineWidth',2);
							hold off

							subplot(4,1,3);
							hold on
							pl4=plot((C.distance)./1000,C.pred+elbl,'-k','LineWidth',2);
							legend([pl1 pl2 pl3 pl4],'Unconditioned DEM','Conditioned DEM','Stream Distance','Segment Fit','location','best');
							hold off

							subplot(4,1,4);
							hold on
							plot(C.area,ksn_val.*C.area.^(-C.mn),'-k','LineWidth',2);
							hold off							
						else
							figure(f2)
							subplot(3,1,1);
							hold on
							plot(C.chi,C.pred+elbl,'-k','LineWidth',2);
							hold off

							subplot(3,1,3);
							hold on
							pl4=plot((C.distance)./1000,C.pred+elbl,'-k','LineWidth',2);
							legend([pl1 pl2 pl3 pl4],'Unconditioned DEM','Conditioned DEM','Stream Distance','Segment Fit','location','best');
							hold off
						end

						res_list=[C.chi C.res];
						bnd_ix=NaN;;

					else
						% Sort knickpoint list and construct bounds list
						ds=sortrows(d);
						bnds=vertcat(0,ds,max(Sn.distance));

						num_bnds=numel(bnds);
						rd=C.distance; 
						rx=C.x;
						ry=C.y; 
						rc=C.chi;
						for jj=1:num_bnds-1
							% Extract bounds
							lb=bnds(jj);
							rb=bnds(jj+1);

							% Clip out stream segment
							lb_dist=sqrt(sum(bsxfun(@minus, rd, lb).^2,2));
							rb_dist=sqrt(sum(bsxfun(@minus, rd, rb).^2,2));

							[~,lbix]=min(lb_dist);
							[~,rbix]=min(rb_dist);

							lbx=rx(lbix);
							lby=ry(lbix);

							rbx=rx(rbix);
							rby=ry(rbix);	

							lix=coord2ind(DEM,lbx,lby);	
							rix=coord2ind(DEM,rbx,rby);

							Seg=modify(Sn,'downstreamto',rix);
							Seg=modify(Seg,'upstreamto',lix);

							%Remake stream with downstream bound node added back in
							WSEG=GRIDobj(DEM,'logical');
							WSEG.Z(Seg.IXgrid)=true;
							WSEG.Z(lix)=true;
							% Add back in upstream node if it's the end of the stream
							if jj==num_bnds-1
								WSEG.Z(rix)=true;
							end
							Seg=STREAMobj(FD,WSEG);

							% Check length of stream, if it's less than two nodes,
							% move down stream until it's greater than two nodes
							lbix_new=lbix+1;
							first_time=true;
							while numel(Seg.IXgrid)<=2
								if first_time
									wrn_mssg=['Segment ' num2str(jj) ' of chosen segments was too short, segment bound was expanded downstream'];
									wd=warndlg(wrn_mssg);
									uiwait(wd);
								end
								lix_new=coord2ind(DEM,rx(lbix_new),ry(lbix_new));
								WSEG.Z(lix_new)=true;
								Seg=STREAMobj(FD,WSEG);
								lbix_new=lbix_new+1;
								first_time=false;
							end

							% Calculate chi to find ksn and bestfit concavity 
							if strcmp(theta_method,'ref')
								Cseg=ChiCalc(Seg,DEMc,A,1,ref_theta);
							elseif strcmp(theta_method,'auto')
								Cseg=ChiCalc(Seg,DEMc,A,1,C.mn);
							end
							Cbfseg=ChiCalc(Seg,DEMc,A,1);

							% Determine where ksn value fits into color scale and plot
							ksn_val=Cseg.ks;

							ksn_nodes{jj,1}=[Cseg.x Cseg.y Cseg.area ones(numel(Cseg.x),1)*Cseg.ks ones(numel(Cseg.x),1)*Cseg.ks_neg ones(numel(Cseg.x),1)*Cseg.ks_pos ...
								ones(numel(Cseg.x),1)*Cseg.mn ones(numel(Cseg.x),1)*Cbfseg.mn ones(numel(Cseg.x),1)*snda];

							% Plot linear fits	
							rchi=rc(rb_dist==min(rb_dist));
							lchi=rc(lb_dist==min(lb_dist));
							ld=rd(lb_dist==min(lb_dist));
							segChi=linspace(lchi,rchi,numel(Cseg.chi));
							seg0Chi=linspace(0,max(Cseg.chi),numel(Cseg.chi));
							[~,lbix]=min(Cseg.chi);
							elbl=Cseg.elev(lbix);

							if display_slope_area
								figure(f2)
								subplot(4,1,1);
								hold on
								plot(segChi,(seg0Chi.*ksn_val)+elbl,'-k','LineWidth',2);
								hold off

								seg_st=rd(lb_dist==min(lb_dist));
								subplot(4,1,3);
								hold on
								pl4=plot((Cseg.distance+seg_st)/1000,(Cseg.pred)+elbl,'-k','LineWidth',2);
								legend([pl1 pl2 pl3 pl4],'Unconditioned DEM','Conditioned DEM','Stream Distance','Segment Fit','location','best');
								hold off

								subplot(4,1,4);
								hold on
								plot(Cseg.area,ksn_val.*Cseg.area.^(-Cseg.mn),'-k','LineWidth',2);
								hold off
							else
								figure(f2)
								subplot(3,1,1);
								hold on
								plot(segChi,(seg0Chi.*ksn_val)+elbl,'-k','LineWidth',2);
								hold off

								seg_st=rd(lb_dist==min(lb_dist));
								subplot(3,1,3);
								hold on
								pl4=plot((Cseg.distance+seg_st)/1000,(Cseg.pred)+elbl,'-k','LineWidth',2);
								legend([pl1 pl2 pl3 pl4],'Unconditioned DEM','Conditioned DEM','Stream Distance','Segment Fit','location','best');
								hold off
							end

							res_list{jj,1}=[Cseg.distance+ld Cseg.res];
						end

						ksn_list=vertcat(ksn_nodes{:});
						res_list=vertcat(res_list{:});
					end

					%% Plot fit result
					if display_slope_area
						figure(f2)
						subplot(4,1,2)
						hold on
						plot(res_list(:,1)./1000,ksn_list(:,4),'-k','LineWidth',2);
						hold off
					else
						figure(f2)
						subplot(3,1,2)
						hold on
						plot(res_list(:,1)./1000,ksn_list(:,4),'-k','LineWidth',2);
						hold off
					end

					f3=figure(3);
					set(f3,'Units','normalized','Position',[0.05 0.1 0.45 0.8],'renderer','painters');

					clf
					sbplt2=subplot(2,1,2);
					hold on
					s1=scatter(DAvg./1000,KsnAvg,20,'k','filled');
					p1=plot(res_list(:,1)./1000,ksn_list(:,4)-ksn_list(:,5),':k');
					plot(res_list(:,1)./1000,ksn_list(:,4)+ksn_list(:,6),':k');
					p2=plot(res_list(:,1)./1000,ksn_list(:,4),'-k','LineWidth',2);
					[ksn_vals,~,ksn_ix]=unique(ksn_list(:,4));
					d_means=accumarray(ksn_ix,res_list(:,1),[],@(x) mean(x,'omitnan'));
					for kk=1:numel(ksn_vals)
						text(d_means(kk)./1000,ksn_vals(kk),['k_{sn} = ' num2str(ksn_vals(kk))],...
							'VerticalAlignment','bottom','HorizontalAlignment','center');
					end
					xlabel('Distance (km)')
					ylabel('k_{sn}')
					title('k_{sn} - Distance')
					legend([s1 p1 p2],{'Auto k_{sn}','k_{sn} uncertainty','k_{sn} of fit segments'},'location','best');
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(sbplt2);
				    end						
					hold off

					sbplt1=subplot(2,1,1);
					hold on
					plot([min(res_list(:,1)./1000) max(res_list(:,1)./1000)],[0 0],'-k');
					scatter(res_list(:,1)./1000,res_list(:,2),10,'k','filled');
					xlabel('Distance (km)')
					ylabel('Residual (m)')
					title('Residual on k_{sn} fit')
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(sbplt1);
				    end						
					hold off

				elseif strcmp(pick_method,'slope_area')

					[bs,ba,bc,bd,bk,aa,ag,ad,ac]=sa_ksn(DEMc,Sn,A,C.chi,ak,smooth_distance);

					ax3=subplot(4,1,3);
					hold on
					pl1=plotdz(Sn,DEM,'dunit','km','Color',[0.5 0.5 0.5]);
					pl2=plotdz(Sn,DEMc,'dunit','km','Color','k');
					pl3=scatter((C.distance)./1000,C.elev,5,log10(C.area),'filled');
					xlabel('Distance from Mouth (km)')
					ylabel('Elevation (m)')
					legend([pl1 pl2 pl3],'Unconditioned DEM','Conditioned DEM','Log Drainage Area','location','best');
					title('Long Profile')
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(ax3);
				    end						
					hold off

					ax2=subplot(4,1,2);
					hold on
					scatter(ba,bk,20,log10(ba),'filled','MarkerEdgeColor','k');
					xlabel('Log Area')
					ylabel('Auto k_{sn}');
					title('Log Area - Auto k_{sn}');
					set(ax2,'XScale','log','XDir','reverse');
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(ax2);
				    end						
					hold off

					ax1=subplot(4,1,1);
					hold on
					plot(C.chi,C.elev,'-k');
					scatter(C.chi,C.elev,10,log10(C.area),'filled');
					xlabel('\chi')
					ylabel('Elevation (m)')
					title('\chi - Z')
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(ax1);
				    end						
					hold off

					ax4=subplot(4,1,4);
					hold on
					scatter(aa,ag,5,log10(aa),'+');
					scatter(ba,bs,20,log10(ba),'filled','MarkerEdgeColor','k');
					xlabel('Log Area');
					ylabel('Log Gradient');
					title(['Slope-Area: \theta = ' num2str(C.mn) ' : Pick Segments - Press Enter When Done'],'Color','r');
					set(ax4,'YScale','log','XScale','log','XDir','reverse');
					ax4.XColor='Red';
					ax4.YColor='Red';
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(ax4);
				    end						
					hold off

					linkaxes([ax4,ax2],'x');
					colormap(ax1,'jet'); colormap(ax2,'jet'); colormap(ax3,'jet'); colormap(ax4,'jet');

					[av,~,bttn]=ginput;
					% Determine if there any non bound knicks
					bttn_idx=bttn~=1;
					if any(bttn_idx)
						% Parse out non bounds
						av_kn=av(bttn_idx);
						av(bttn_idx)=[];
						% Convert to indices
						ra=C.area; rx=C.x; ry=C.y;
						kn_ix=zeros(numel(av_kn),1); 
						for jj=1:numel(av_kn);
							areadist=sqrt(sum(bsxfun(@minus, ra, av_kn(jj)).^2,2));
							[~,knbix]=min(areadist);
							knbx=rx(knbix);
							knby=ry(knbix);
							kn_ix(jj)=coord2ind(DEM,knbx,knby);
						end
					else
						kn_ix=NaN;
					end	

					if isempty(av)
						if strcmp(theta_method,'ref')
							Cbf=ChiCalc(Sn,DEMc,A,1);
							ksn_list=[C.x C.y C.area ones(numel(C.x),1)*C.ks ones(numel(C.x),1)*C.ks_neg ones(numel(C.x),1)*C.ks_pos ...
								ones(numel(C.x),1)*C.mn ones(numel(C.x),1)*Cbf.mn ones(numel(C.x),1)*snda];
						elseif strcmp(theta_method,'auto')
							ksn_list=[C.x C.y C.area ones(numel(C.x),1)*C.ks ones(numel(C.x),1)*C.ks_neg ones(numel(C.x),1)*C.ks_pos ...
								ones(numel(C.x),1)*C.mn ones(numel(C.x),1)*C.mn ones(numel(C.x),1)*snda];
						end

						% Determine where ksn value fits into color scale and plot
						ksn_val=C.ks;

						[~,lbix]=min(C.chi);
						elbl=C.elev(lbix);

						figure(f2)
						subplot(4,1,1);
						hold on
						plot(C.chi,C.pred+elbl,'-k','LineWidth',2);
						hold off

						subplot(4,1,3);
						hold on
						pl4=plot((C.distance)/1000,C.pred+elbl,'-k','LineWidth',2);
						legend([pl1 pl2 pl3 pl4],'Unconditioned DEM','Conditioned DEM','Log Drainage Area','Segment Fit','location','best');
						hold off

						subplot(4,1,4);
						hold on
						plot(C.area,ksn_val.*C.area.^(-C.mn),'-k','LineWidth',2);
						hold off

						res_list=[C.area C.res];
						bnd_ix=NaN;;

					else
						% Sort knickpoint list and construct bounds list
						avs=sortrows(av,'descend');
						bnds=vertcat(max(C.area,[],'omitnan'),avs,min(C.area,[],'omitnan'));

						num_bnds=numel(bnds);
						rc=C.chi;
						rx=C.x;
						ry=C.y;
						rd=C.distance;
						ra=C.area;
						for jj=1:num_bnds-1
							% Extract bounds
							lb=bnds(jj);
							rb=bnds(jj+1);

							% Clip out stream segment
							lb_dadist=sqrt(sum(bsxfun(@minus, ra, lb).^2,2));
							rb_dadist=sqrt(sum(bsxfun(@minus, ra, rb).^2,2));

							[~,lbix]=min(lb_dadist);
							[~,rbix]=min(rb_dadist);

							lbx=rx(lbix);
							lby=ry(lbix);

							rbx=rx(rbix);
							rby=ry(rbix);	

							lix=coord2ind(DEM,lbx,lby);
							rix=coord2ind(DEM,rbx,rby);

							Seg=modify(Sn,'downstreamto',rix);
							Seg=modify(Seg,'upstreamto',lix);

							%Remake stream with downstream bound node added back in
							WSEG=GRIDobj(DEM,'logical');
							WSEG.Z(Seg.IXgrid)=true;
							WSEG.Z(lix)=true;
							% Add back in upstream node if it's the end of the stream
							if jj==num_bnds-1
								WSEG.Z(rix)=true;
							end
							Seg=STREAMobj(FD,WSEG);

							% Check length of stream, if it's less than two nodes,
							% move down stream until it's greater than two nodes
							lbix_new=lbix+1;
							first_time=true;
							while numel(Seg.IXgrid)<=2
								if first_time
									wrn_mssg=['Segment ' num2str(jj) ' of chosen segments was too short, segment bound was expanded downstream'];
									wd=warndlg(wrn_mssg);
									uiwait(wd);
								end
								lix_new=coord2ind(DEM,rx(lbix_new),ry(lbix_new));
								WSEG.Z(lix_new)=true;
								Seg=STREAMobj(FD,WSEG);
								lbix_new=lbix_new+1;
								first_time=false;
							end

							% Construct bound list
							if jj<num_bnds-1
								bnd_ix(jj,1)=rix;
							end

							% Calculate chi to find ksn and bestfit concavity 
							if strcmp(theta_method,'ref')
								Cseg=ChiCalc(Seg,DEMc,A,1,ref_theta);
							elseif strcmp(theta_method,'auto')
								Cseg=ChiCalc(Seg,DEMc,A,1,C.mn);
							end
							Cbfseg=ChiCalc(Seg,DEMc,A,1);								

							% Determine where ksn value fits into color scale and plot
							ksn_val=Cseg.ks;

							ksn_nodes{jj,1}=[Cseg.x Cseg.y Cseg.area ones(numel(Cseg.x),1)*Cseg.ks ones(numel(Cseg.x),1)*Cseg.ks_neg ones(numel(Cseg.x),1)*Cseg.ks_pos ...
								ones(numel(Cseg.x),1)*Cseg.mn ones(numel(Cseg.x),1)*Cbfseg.mn ones(numel(Cseg.x),1)*snda];

							% Plot linear fits	
							rchi=rc(rb_dadist==min(rb_dadist));
							lchi=rc(lb_dadist==min(lb_dadist));
							segChi=linspace(lchi,rchi,numel(Cseg.chi));
							seg0Chi=linspace(0,max(Cseg.chi),numel(Cseg.chi));
							[~,lbix]=min(Cseg.chi);
							elbl=Cseg.elev(lbix);

							figure(f2)
							subplot(4,1,1);
							hold on
							plot(segChi,(seg0Chi.*ksn_val)+elbl,'-k','LineWidth',2);
							hold off

							seg_st=rd(lb_dadist==min(lb_dadist));
							subplot(4,1,3);
							hold on
							pl4=plot((Cseg.distance+seg_st)/1000,(Cseg.pred)+elbl,'-k','LineWidth',2);
							legend([pl1 pl2 pl3 pl4],'Unconditioned DEM','Conditioned DEM','Log Drainage Area','Segment Fit','location','best');
							hold off

							subplot(4,1,4);
							hold on
							plot(Cseg.area,ksn_val.*Cseg.area.^(-1*Cseg.mn),'-k','LineWidth',2);
							hold off


							res_list{jj,1}=[Cseg.area Cseg.res];


						end

						ksn_list=vertcat(ksn_nodes{:});
						res_list=vertcat(res_list{:});
					end


					%% Plot result figure
					if display_slope_area
						figure(f2)
						subplot(4,1,2)
						hold on
						plot(res_list(:,1),ksn_list(:,4),'-k','LineWidth',2);
						hold off
					else
						figure(f2)
						subplot(3,1,2)
						hold on
						plot(res_list(:,1),ksn_list(:,4),'-k','LineWidth',2);
						hold off
					end

					f3=figure(3);
					set(f3,'Units','normalized','Position',[0.05 0.1 0.45 0.8],'renderer','painters');
					clf
					ax31=subplot(2,1,2);
					hold on
					s1=scatter(log10(ba),bk,20,'k','filled');
					p1=plot(log10(res_list(:,1)),ksn_list(:,4)-ksn_list(:,5),':k');
					plot(log10(res_list(:,1)),ksn_list(:,4)+ksn_list(:,6),':k');
					p2=plot(log10(res_list(:,1)),ksn_list(:,4),'-k','LineWidth',2);
					[ksn_vals,~,ksn_ix]=unique(ksn_list(:,4));
					a_means=accumarray(ksn_ix,res_list(:,1),[],@(x) mean(x,'omitnan'));
					for kk=1:numel(ksn_vals)
						text(log10(a_means(kk)),ksn_vals(kk),['k_{sn} = ' num2str(ksn_vals(kk))],...
							'VerticalAlignment','bottom','HorizontalAlignment','center');
					end
					xlabel('Log Area')
					ylabel('k_{sn}')
					title('k_{sn} - Log Area')
					legend([s1 p1 p2],{'Auto k_{sn}','k_{sn} uncertainty','k_{sn} of fit segments'},'location','best');
					set(ax31,'XScale','log','XDir','reverse');
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(ax31);
				    end						
					hold off

					ax32=subplot(2,1,1);
					hold on
					plot([log10(min(res_list(:,1))) log10(max(res_list(:,1)))],[0 0],'-k');
					scatter(log10(res_list(:,1)),res_list(:,2),10,'k','filled');
					xlabel('Log Area')
					ylabel('Residual (m)')
					title('Residual on k_{sn} fit')
					set(ax32,'XScale','log','XDir','Reverse');
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(ax32);
				    end						
					hold off
				% End pick method switch	
				end

				if ii<num_ch

					ignore_str=['Ignore ' num2str(num_ch-ii) ' Remaining Streams'];
					qa2=questdlg('What would you like to do?','Stream Fitting',ignore_str,'Redo Fit','Continue Picking','Continue Picking');

					switch qa2
					case 'Continue Picking'
						str2 = 'Y';
						str1 = [];
						% Add Gradient, Residual, and Riv Number to node list and clear NaNs
						kidx=isnan(ksn_list(:,1));
						ksn_list(kidx,:)=[];
						res_list(kidx,:)=[];
						gix=coord2ind(DEM,ksn_list(:,1),ksn_list(:,2));
						kmat=[ksn_list G.Z(gix) res_list(:,2) ones(size(ksn_list,1),1).*ii];
						ksn_master{ii,1}=kmat;
						bnd_master{ii,1}=bnd_ix;
						kn_master{ii,1}=kn_ix;
						res_master{ii,1}=res_list;		
						Sc=Sct;
						count=ii;
						save(out_restart_name,'ksn_master','bnd_master','kn_master','res_master','Sc','count','-append');						
						if save_figures
							f2_name=[shape_name '_stream_fits_' num2str(ii) '.pdf'];
							f3_name=[shape_name '_stream_rsds_' num2str(ii) '.pdf'];
							print(f2,f2_name,'-dpdf','-fillpage');
							print(f3,f3_name,'-dpdf','-fillpage');
						end
						clear ksn_list ksn_nodes res_list bnd_ix kn_ix kidx gix kmat;
					case 'Redo Fit'
						str2 = 'R';
						str1 = 'R';
						clear ksn_list ksn_nodes res_list bnd_ix kn_ix;
					case ignore_str
						wtb=waitbar(0,'Cleaning up and generating outputs, do not close windows');
						str1=[];
						str2=[];
						% Add Gradient, Residual, and Riv Number to node list and clear NaNs
						kidx=isnan(ksn_list(:,1));
						ksn_list(kidx,:)=[];
						res_list(kidx,:)=[];
						gix=coord2ind(DEM,ksn_list(:,1),ksn_list(:,2));
						kmat=[ksn_list G.Z(gix) res_list(:,2) ones(size(ksn_list,1),1).*ii];
						ksn_master{ii,1}=kmat;
						bnd_master{ii,1}=bnd_ix;
						kn_master{ii,1}=kn_ix;
						res_master{ii,1}=res_list;		
						Sc=Sct;
						count=ii;
						save(out_restart_name,'ksn_master','bnd_master','kn_master','res_master','Sc','count','-append');
						if save_figures
							f2_name=[shape_name '_stream_fits_' num2str(ii) '.pdf'];
							f3_name=[shape_name '_stream_rsds_' num2str(ii) '.pdf'];
							print(f2,f2_name,'-dpdf','-fillpage');
							print(f3,f3_name,'-dpdf','-fillpage');
						end
						clear ksn_list ksn_nodes res_list bnd_ix kn_ix kidx gix kmat;	
						break_flag=true;
					end

					close figure 2
					close figure 3
				else
					qa2=questdlg('What would you like to do?','Stream Fitting','Redo Fit','Complete Routine','Complete Routine');

					switch qa2
					case 'Complete Routine'
						wtb=waitbar(0,'Cleaning up and generating outputs, do not close windows');
						str2 = 'Y';
						str1 = [];
						% Add Gradient, Residual, and Riv Number to node list and clear NaNs
						kidx=isnan(ksn_list(:,1));
						ksn_list(kidx,:)=[];
						res_list(kidx,:)=[];
						gix=coord2ind(DEM,ksn_list(:,1),ksn_list(:,2));
						kmat=[ksn_list G.Z(gix) res_list(:,2) ones(size(ksn_list,1),1).*ii];
						ksn_master{ii,1}=kmat;
						bnd_master{ii,1}=bnd_ix;
						kn_master{ii,1}=kn_ix;
						res_master{ii,1}=res_list;		
						Sc=Sct;
						count=ii;
						save(out_restart_name,'ksn_master','bnd_master','res_master','Sc','count','-append');
						if save_figures
							f2_name=[shape_name '_stream_fits_' num2str(ii) '.pdf'];
							f3_name=[shape_name '_stream_rsds_' num2str(ii) '.pdf'];
							print(f2,f2_name,'-dpdf','-fillpage');
							print(f3,f3_name,'-dpdf','-fillpage');
						end
						clear ksn_list ksn_nodes res_list bnd_ix kn_ix kidx gix kmat;
					case 'Redo Fit'
						str2 = 'R';
						str1 = 'R';
						clear ksn_list ksn_nodes res_list bnd_ix kn_ix;
					end

					close figure 2
					close figure 3					
				end


			end
			% Reset for next loop
			str1='R';
			
			if break_flag
				break
			end

		end

	% End switch between interactive vs provided
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Bundle, plot, and export %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	switch input_method
	case 'interactive'
		close figure 1
	end

	% Add stream numbers to bound list and convert bound list
	bnd_list=cell(numel(bnd_master),1);
	for jj=1:numel(bnd_master)
		bnd=bnd_master{jj,1};
		if ~isnan(bnd)
			% Extract elevations
			bz=DEM.Z(bnd);
			% Find X Y
			[bx,by]=ind2coord(DEM,bnd);
			% Preserve linear index for later indexing
			bix=bnd;
			% Add River Number
			brm=ones(size(bnd)).*jj;
			% Compile and store
			bnd_list{jj,1}=[bx by bz brm bix];
		elseif isempty(bnd)
			% bnd_list{jj,1}=[0 0 0 jj nan];
			bnd_list{jj,1}=[];
		else
			bnd_list{jj,1}=[0 0 0 jj bnd];
		end
	end

	% Collapse bnd_list
	bnd_list=vertcat(bnd_list{:});

	% Add stream numbers to knick list and convert knick list
	kn_list=cell(numel(kn_master),1);
	for jj=1:numel(kn_master)
		kn=kn_master{jj,1};
		if ~isnan(kn)
			knz=DEM.Z(kn);
			[knx,kny]=ind2coord(DEM,kn);
			knix=kn;
			krm=ones(size(kn)).*jj;
			kn_list{jj,1}=[knx kny knz krm knix];
		elseif isempty(kn)
			kn_list{jj,1}=[];
		else
			kn_list{jj,1}=[0 0 0 jj kn];
		end
	end
	kn_list=vertcat(kn_list{:});

	waitbar(1/5,wtb);

	% Build KSN Mapstructure
	if strcmp(junction_method,'check')
		% Bundle all observations into single node list
		knl=vertcat(ksn_master{:});
		ix=coord2ind(DEM,knl(:,1),knl(:,2));

		% Generate empty rasters
		ksnR=GRIDobj(DEM);
		ksnRn=GRIDobj(DEM);
		ksnRp=GRIDobj(DEM);
		thetaR=GRIDobj(DEM);
		segthetaR=GRIDobj(DEM);
		threshaR=GRIDobj(DEM);
		resR=GRIDobj(DEM);
		rivnumR=GRIDobj(DEM);

		ksnR.Z(ix)=knl(:,4);
		ksnRn.Z(ix)=knl(:,5);
		ksnRp.Z(ix)=knl(:,6);
		thetaR.Z(ix)=knl(:,7);
		segthetaR.Z(ix)=knl(:,8);
		threshaR.Z(ix)=knl(:,9);
		resR.Z(ix)=knl(:,11);
		rivnumR.Z(ix)=knl(:,12);

		% Create KSN map structure and export shapefile
		KSN=STREAMobj2mapstruct(Sc,'seglength',smooth_distance,'attributes',...
			{'ksn' ksnR @mean 'ksn_neg' ksnRn @mean 'ksn_pos' ksnRp @mean 'uparea' (A.*(A.cellsize^2)) @mean...
			'gradient' G @mean 'theta' thetaR @mean 'seg_theta' segthetaR @mean 'thrsh_ar' threshaR @mean...
			'resid' resR @mean 'riv_num' rivnumR @median});
	elseif strcmp(junction_method,'ignore') & strcmp(stack_method,'stack')
		num_streams=numel(ksn_master);
		KSN=cell(num_streams,1);
		for ii=1:num_streams
			knl=ksn_master{ii};
			ix=coord2ind(DEM,knl(:,1),knl(:,2));
			WW=GRIDobj(DEM,'logical');
			WW.Z(ix)=true;
			ScT=STREAMobj(FD,WW);

			% Generate empty rasters
			ksnR=GRIDobj(DEM);
			ksnRn=GRIDobj(DEM);
			ksnRp=GRIDobj(DEM);
			thetaR=GRIDobj(DEM);
			segthetaR=GRIDobj(DEM);
			threshaR=GRIDobj(DEM);
			resR=GRIDobj(DEM);
			rivnumR=GRIDobj(DEM);

			ksnR.Z(ix)=knl(:,4);
			ksnRn.Z(ix)=knl(:,5);
			ksnRp.Z(ix)=knl(:,6);
			thetaR.Z(ix)=knl(:,7);
			segthetaR.Z(ix)=knl(:,8);
			threshaR.Z(ix)=knl(:,9);
			resR.Z(ix)=knl(:,11);
			rivnumR.Z(ix)=knl(:,12);

			KSN{ii}=STREAMobj2mapstruct(ScT,'seglength',smooth_distance,'attributes',...
				{'ksn' ksnR @mean 'ksn_neg' ksnRn @mean 'ksn_pos' ksnRp @mean 'uparea' (A.*(A.cellsize^2)) @mean...
				'gradient' G @mean 'theta' thetaR @mean 'seg_theta' segthetaR @mean 'thrsh_ar' threshaR @mean...
				'resid' resR @mean 'riv_num' rivnumR @median});	
		end
		KSN=vertcat(KSN{:});
	elseif strcmp(junction_method,'ignore') & strcmp(stack_method,'average')
		% Bundle all observations into single node list
		knl=vertcat(ksn_master{:});
		ix=coord2ind(DEM,knl(:,1),knl(:,2));

		% Generate empty rasters
		ksnR=GRIDobj(DEM);
		ksnRn=GRIDobj(DEM);
		ksnRp=GRIDobj(DEM);
		thetaR=GRIDobj(DEM);
		segthetaR=GRIDobj(DEM);
		threshaR=GRIDobj(DEM);
		resR=GRIDobj(DEM);
		rivnumR=GRIDobj(DEM);

		% Accumulate values
		ksnr=accumarray(ix,knl(:,4),[],@mean);
		ksnrn=accumarray(ix,knl(:,5),[],@mean);
		ksnrp=accumarray(ix,knl(:,6),[],@mean);
		thetar=accumarray(ix,knl(:,7),[],@mean);
		segthetar=accumarray(ix,knl(:,8),[],@mean);
		threshar=accumarray(ix,knl(:,9),[],@mean);
		resr=accumarray(ix,knl(:,11),[],@mean);
		rivnumr=accumarray(ix,knl(:,12),[],@mode);

		% Filter values
		ksnr=ksnr(ix);
		ksnrn=ksnrn(ix);
		ksnrp=ksnrp(ix);
		thetar=thetar(ix);
		segthetar=segthetar(ix);
		threshar=threshar(ix);
		resr=resr(ix);
		rivnumr=rivnumr(ix);

		% Populate grids
		ksnR.Z(ix)=ksnr;
		ksnRn.Z(ix)=ksnrn;
		ksnRp.Z(ix)=ksnrp;
		thetaR.Z(ix)=thetar;
		segthetaR.Z(ix)=segthetar;
		threshaR.Z(ix)=threshar;
		resR.Z(ix)=resr;
		rivnumR.Z(ix)=rivnumr;

		% Create KSN map structure and export shapefile
		KSN=STREAMobj2mapstruct(Sc,'seglength',smooth_distance,'attributes',...
			{'ksn' ksnR @mean 'ksn_neg' ksnRn @mean 'ksn_pos' ksnRp @mean 'uparea' (A.*(A.cellsize^2)) @mean...
			'gradient' G @mean 'theta' thetaR @mean 'seg_theta' segthetaR @mean 'thrsh_ar' threshaR @mean...
			'resid' resR @mean 'riv_num' rivnumR @mode});
	end	

	waitbar(2/5,wtb);	

	% Create knickpoint map structure and prepare bound output
	idx=~isnan(bnd_list(:,5));
	bnd_strc=bnd_list(idx,:);
	bnd_list(~idx,1)=NaN; bnd_list(~idx,2)=NaN; bnd_list(~idx,3)=NaN;
	bnd_list=bnd_list(:,[1:4]);

	if ~isempty(bnd_strc)
		KNK=struct;
		for jj=1:numel(bnd_strc(:,1));
			KNK(jj,1).Geometry='Point';
			KNK(jj,1).X=double(bnd_strc(jj,1));
			KNK(jj,1).Y=double(bnd_strc(jj,2));
			KNK(jj,1).Elev=double(bnd_strc(jj,3));
			KNK(jj,1).StrNum=double(bnd_strc(jj,4));
		end
		out_bound_name=[shape_name '_bounds.shp'];
		shapewrite(KNK,out_bound_name);
	end

	idx=~isnan(kn_list(:,5));
	kn_strc=kn_list(idx,:);
	kn_list(~idx,1)=NaN; kn_list(~idx,2)=NaN; kn_list(~idx,3)=NaN;
	kn_list=kn_list(:,[1:4]);

	if ~isempty(kn_strc)
		XKNK=struct;
		for jj=1:numel(kn_strc(:,1))
			XKNK(jj,1).Geometry='Point';
			XKNK(jj,1).X=double(kn_strc(jj,1));
			XKNK(jj,1).Y=double(kn_strc(jj,2));
			XKNK(jj,1).Elev=double(kn_strc(jj,3));
			XKNK(jj,1).StrNum=double(kn_strc(jj,4));
		end
		out_knick_name=[shape_name '_knicks.shp'];
		shapewrite(XKNK,out_knick_name);
	end

	waitbar(3/5,wtb);

	out_shape_name=[shape_name '.shp'];
	shapewrite(KSN,out_shape_name);

	waitbar(4/5,wtb);

	% Save out file
	save(out_mat_name,'knl','ksn_master','bnd_list','kn_list','Sc','bnd_master','kn_master','res_master','count','-append');
	% Delete restart file after successful completion
	delete(out_restart_name);

	waitbar(5/5,wtb);
	close(wtb);

%FUNCTION END
end

function [OUT]=ChiCalc(S,DEM,A,a0,varargin)
	% Modified version of chiplot function by Wolfgang Schwanghart to remove unused options/outputs and 
	% to interpolate chi-z relationship so that chi values are equally spaced to avoid biasing of ksn fit
	% by clustering at high drainage areas 

	% nr of nodes in the entire stream network
	nrc = numel(S.x);
	M   = sparse(double(S.ix),double(S.ixc),true,nrc,nrc);
	% find outlet
	outlet = sum(M,2) == 0 & sum(M,1)'~=0;
	if nnz(outlet)>1
	    % there must not be more than one outlet (constraint could be removed
	    % in the future).
	    if isdeployed
	    	errordlg('The stream network must not have more than one outlet')
	    end
	    error('The stream network must not have more than one outlet');
	end


	% elevation values at nodes
	zx   = double(DEM.Z(S.IXgrid));
	% elevation at outlet
	if nnz(outlet)==0
		zb=min(zx);
	else
		zb   = double(DEM.Z(S.IXgrid(outlet)));
	end
	% a is the term inside the brackets of equation 6b 
	a    = double(a0./(A.Z(S.IXgrid)*(A.cellsize.^2)));
	% x is the cumulative horizontal distance in upstream direction
	x    = S.distance;
	Lib = true(size(x));

	% Find bestfit concavity if no concavity provided or use provided concavity
	if isempty(varargin)
	    mn0  = 0.5; % initial value
	    mn   = fminsearch(@mnfit,mn0);
	else
		mn=varargin{1};
	end

	% calculate chi
	chi = netcumtrapz(x,a.^mn,S.ix,S.ixc);

	% Resample chi-elevation relationship using cubic spline interpolation
	chiF=chi(Lib);
	zabsF=zx(Lib)-zb;

	% The splining and fitting generates lots of warning for small basins so turning warnings off
	warning off

	chiS=linspace(0,max(chiF),numel(chiF)).';
	try
		zS=spline(chiF,zabsF,chiS);
	catch
		% cubic spline will fail if segments is nearly straight, skip spline fit
		% in this case to avoid erroring out
		zS=zabsF;
		chiS=chiF;
	end

	

	OUT=struct;
	try
		ft=fittype('a*x');
		fobj=fit(chiS,zS,ft,'StartPoint',chiS\zS);
		BETA=coeffvalues(fobj);
		BETA_UNC=confint(fobj);
		OUT.ks   = BETA*a0^mn;
		OUT.ks_neg = (BETA*a0^mn)-(min(BETA_UNC)*a0^mn);
		OUT.ks_pos = (max(BETA_UNC)*a0^mn)-(BETA*a0^mn);
		OUT.mn   = mn;
	catch
		BETA = chiS\(zS);
		OUT.ks   = BETA*a0^mn;
		OUT.ks_neg = 0;
		OUT.ks_pos = 0;
		OUT.mn   = mn;
	end

	warning on

	[OUT.x,...
	 OUT.y,...
	 OUT.chi,...
	 OUT.elev,...
	 OUT.elevbl,...
	 OUT.distance,...
	 OUT.pred,...
	 OUT.area] = STREAMobj2XY(S,chi,DEM,zx-zb,S.distance,BETA*chi,A.*(A.cellsize^2));
	 OUT.res = OUT.elevbl - OUT.pred;

	 	% Nested function for calculating mn ratio
	function sqres = mnfit(mn)
		% calculate chi with a given mn ratio
		% and integrate in upstream direction
		CHI = netcumtrapz(x(Lib),a(Lib).^mn,S.ix,S.ixc);%*ab.^mn
		% normalize both variables
		CHI = CHI ./ max(CHI);
		z   = zx(Lib)-zb;
		z   = z./max(z);
        sqres = sum((CHI - z).^2);
	end

end

function z = netcumtrapz(x,y,ix,ixc)
	% cumtrapz along upward direction in a directed tree network

	z = zeros(size(x));
	for lp = numel(ix):-1:1;
	    z(ix(lp)) = z(ixc(lp)) + (y(ixc(lp))+(y(ix(lp))-y(ixc(lp)))/2) *(abs(x(ixc(lp))-x(ix(lp))));
	end
end

function [Xavg,Yavg]=BinAverage(X,Y,bin_size);

	ix=~isnan(X);
	X=X(ix); Y=Y(ix);

	minX=min(X);
	maxX=max(X);

	b=[minX:bin_size:maxX+bin_size];

	try
		[idx]=discretize(X,b);
	catch
		[~,idx]=histc(X,b);
	end

	Xavg=accumarray(idx(:),X,[],@mean);
	Yavg=accumarray(idx(:),Y,[],@mean);
end

function [ksn]=KSN_Quick(DEM,A,S,theta_ref)
	zc=mincosthydrocon(S,DEM,'interp',0.1);

	g=gradient(S,zc);
	G=GRIDobj(DEM);
	G.Z(G.Z==0)=NaN;
	G.Z(S.IXgrid)=g;

	ksn=G./(A.*(A.cellsize^2)).^(-theta_ref);
end

function [bs,ba,bc,bd,a,g,d,C]=sa(DEM,S,A,C,bin_size)
	% Modified slope area function that uses the smooth length to
	%	to determine the number of bins and uses those same bins
	%	to find mean values of chi and distance for plotting
	%	purposes

	minX=min(S.distance);
	maxX=max(S.distance);
	b=[minX:bin_size:maxX+bin_size];

	numbins=round(max([numel(b) numel(S.IXgrid)/10]));

	an=getnal(S,A.*A.cellsize^2);
	z=getnal(S,DEM);
	gn=gradient(S,z,'unit','tangent','method','robust','drop',20);
	gn=smooth(gn,3);

	% Run through STREAMobj2XY so chi and everything else are same size
	[~,~,a,g,d]=STREAMobj2XY(S,an,gn,S.distance);
	% Remove NaNs
	% a(isnan(a))=[];
	% g(isnan(g))=[];
	% d(isnan(d))=[];
	% C(isnan(C))=[];

	% This avoids erroneous selection of confluence points at end of selected streams
	rlIDX= ~isnan(a) & ~isnan(g) & ~isnan(d) & ~isnan(C);
	a = a(rlIDX); 
	g = g(rlIDX);
	d = d(rlIDX);
	C = C(rlIDX):

	mina=min(a);
	maxa=max(a);

    edges = logspace(log10(mina-0.1),log10(maxa+1),numbins+1);
    try
    	% histc is deprecated
    	[ix]=discretize(a,edges);
    catch
	    [~,ix] = histc(a,edges);
	end

	ba=accumarray(ix,a,[numbins 1],@median,nan);
	bs=accumarray(ix,g,[numbins 1],@(x) mean(x(~isnan(x))),nan);
	bd=accumarray(ix,d,[numbins 1],@mean,nan);
	bc=accumarray(ix,C,[numbins 1],@mean,nan);

	% Filter negatives
	idx=bs>=0 & ba>=0 & bc>=0 & bd>=0;
	bs=bs(idx);
	ba=ba(idx);
	bc=bc(idx);
	bd=bd(idx);

	idx=a>=0 & g>=0 & d>=0 & C>=0;
	a=a(idx);
	g=g(idx);
	d=d(idx);
	C=C(idx);
end

function [bs,ba,bc,bd,bk,a,g,d,C]=sa_ksn(DEM,S,A,C,ak,bin_size);
	% Modified slope area function that uses the smooth length to
	%	to determine the number of bins and uses those same bins
	%	to find mean values of chi and distance for plotting
	%	purposes

	minX=min(S.distance);
	maxX=max(S.distance);
	b=[minX:bin_size:maxX+bin_size];

	numbins=round(max([numel(b) numel(S.IXgrid)/10]));

	an=getnal(S,A.*A.cellsize^2);
	z=getnal(S,DEM);
	gn=gradient(S,z,'unit','tangent','method','robust','drop',20);
	gn=smooth(gn,3);

	% Run through STREAMobj2XY so chi and everything else are same size
	[~,~,a,g,d,k]=STREAMobj2XY(S,an,gn,S.distance,ak);
	% Remove NaNs
	a(isnan(a))=[];
	g(isnan(g))=[];
	d(isnan(d))=[];
	C(isnan(C))=[];
	k(isnan(k))=[];

	mina=min(a);
	maxa=max(a);

    edges = logspace(log10(mina-0.1),log10(maxa+1),numbins+1);
    try
    	% histc is deprecated
    	[ix]=discretize(a,edges);
    catch
	    [~,ix] = histc(a,edges);
	end

	ba=accumarray(ix,a,[numbins 1],@median,nan);
	bs=accumarray(ix,g,[numbins 1],@(x) mean(x(~isnan(x))),nan);
	bd=accumarray(ix,d,[numbins 1],@mean,nan);
	bc=accumarray(ix,C,[numbins 1],@mean,nan);
	bk=accumarray(ix,k,[numbins 1],@mean,nan);

	% Filter negatives
	idx=bs>=0 & ba>=0 & bc>=0 & bd>=0 & bk>=0;
	bs=bs(idx);
	ba=ba(idx);
	bc=bc(idx);
	bd=bd(idx);
	bk=bk(idx);

	idx=a>=0 & g>=0 & d>=0 & C>=0 & k>=0;
	a=a(idx);
	g=g(idx);
	d=d(idx);
	C=C(idx);
	k=k(idx);
end

function [Sn]=RedefineThreshold(DEM,FD,A,S,FLUS,ref_theta,pick_method,bin_size,count,figure_flag,shape_name)

	% Find channel head and flow distances
	chix=streampoi(S,'channelheads','ix');
	DA=A.*(DEM.cellsize^2);


	UP=dependencemap(FD,chix);
	FLDSt=DEM.*UP;

	[~,ix]=max(FLDSt);

	IX=influencemap(FD,ix);

	ST=STREAMobj(FD,IX);
	z=mincosthydrocon(ST,DEM,'interp',0.1);

	C=chiplot(ST,z,A,'a0',1,'mn',ref_theta,'plot',false);
	[bs,ba,bc,bd,aa,ag,ad,ac]=sa(DEM,ST,A,C.chi,bin_size);


	% % Filter negatives
	% idx=bs>=0 & ba>=0 & bc>=0 & bd>=0;
	% bs=bs(idx);
	% ba=ba(idx);
	% bc=bc(idx);
	% bd=bd(idx);

	% idx=aa>=0 & ag>=0 & ad>=0 & ac>=0;
	% aa=aa(idx);
	% ag=ag(idx);
	% ad=ad(idx);
	% ac=ac(idx);

	str11='R';

	while strcmp(str11,'R');
		f4=figure(4);
		set(f4,'Units','normalized','Position',[0.5 0.1 0.45 0.8],'renderer','painters');
		clf

		colormap(jet);

		switch pick_method
		case 'chi'

			ax2=subplot(2,1,2);
			hold on 
			scatter(aa,ag,5,ac,'+');
			scatter(ba,bs,20,bc,'filled','MarkerEdgeColor','k');
			xlabel('Log Drainage Area');
			ylabel('Log Gradient');
			caxis([0 max(C.chi)]);
			set(ax2,'YScale','log','XScale','log','XDir','reverse');
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(ax2);
		    end	
			hold off

			ax1=subplot(2,1,1);
			hold on
			plot(C.chi,C.elev,'-k');
			scatter(C.chi,C.elev,10,C.chi,'filled');
			xlabel('\chi');
			ylabel('Elevation (m)');
			title('Choose hillslope to channel transition');
			caxis([0 max(C.chi)]);
			ax1.XColor='Red';
			ax1.YColor='Red';
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(ax1);
		    end				
			hold off

			% Find user selected threshold area
			[c,~]=ginput(1);
			if ~isempty(c)
				[~,cix]=min(abs(C.chi-c),[],'omitnan');
				a=C.area(cix);
			else
				a=1e6;
			end

			chi_idx=C.area<a;
			sl_idx=ba<a;
			aa_idx=aa<a;

			subplot(2,1,2)
			hold on
			scatter(aa(aa_idx),ag(aa_idx),10,'k','+');
			scatter(ba(sl_idx),bs(sl_idx),30,'k','filled');
			hold off

			subplot(2,1,1)
			hold on
			scatter(C.chi(chi_idx),C.elev(chi_idx),20,'k','filled');
			title('Black points will be excluded from stream definition');
			hold off

		case 'slope_area'

			ax2=subplot(2,1,2);
			hold on
			plot(C.chi,C.elev,'-k');
			scatter(C.chi,C.elev,10,C.chi,'filled');
			xlabel('\chi');
			ylabel('Elevation (m)');
			caxis([0 max(C.chi)]);
			xlim([0 max(C.chi)+0.5]);
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(ax2);
		    end				
			hold off

			ax1=subplot(2,1,1);
			hold on 
			scatter(aa,ag,5,ac,'+');
			scatter(ba,bs,20,bc,'filled','MarkerEdgeColor','k');
			xlabel('Log Drainage Area');
			ylabel('Log Gradient');
			title('Choose hillslope to channel transition');
			caxis([0 max(C.chi)]);
			set(ax1,'YScale','log','XScale','log','XDir','reverse');
			ax1.XColor='Red';
			ax1.YColor='Red';
			if ~verLessThan('matlab','9.5')
		        disableDefaultInteractivity(ax1);
		    end				
			hold off

			% Find user selected threshold area
			[a,~]=ginput(1);
			if isempty(a)
				a=1e6;
			end

			chi_idx=C.area<a;
			sl_idx=ba<a;
			aa_idx=aa<a;

			subplot(2,1,2)
			hold on
			scatter(C.chi(chi_idx),C.elev(chi_idx),20,'k','filled');
			hold off

			subplot(2,1,1)
			hold on
			scatter(aa(aa_idx),ag(aa_idx),10,'k','+');
			scatter(ba(sl_idx),bs(sl_idx),30,'k','filled');
			title('Black points will be excluded from stream definition');
			hold off

		end

		qa3=questdlg('Accept new threshold area?','Set Threshold','No, Redo','Yes','Yes');
		switch qa3
		case 'Yes'
			str11='C';
		case 'No, Redo'
			str11='R';
		end
	end

	if figure_flag
		f4_name=[shape_name '_stream_thresh_' num2str(count) '.pdf'];
		print(f4,f4_name,'-dpdf','-fillpage');
	end

	close(f4);

	da=getnal(ST,A.*A.cellsize^2);
	nix=ST.IXgrid(da>=a);
	IX=GRIDobj(DEM,'logical');
	IX.Z(nix)=true;

	Sn=STREAMobj(FD,IX);
end

