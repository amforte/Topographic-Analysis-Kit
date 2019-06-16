function TAK(func_name,wdir,varargin)
	% Bundles all compiled functions together and deals with the conversion of input 
	% variables to their proper format

	% Number of arguments supplied to function minus the function name call and working
	% directory arguments
	nf_args=nargin-2;

	switch func_name
	case 'MakeStreams'
		num_req=3;
		req_args=varargin(1:num_req);
		if nf_args==num_req
			cmpMakeStreams(wdir,req_args{1},str2num(req_args{2}),req_args{3});
		else
			opt_args=varargin(num_req+1:end);
			num_var=numel(opt_args);
			odds=1:2:num_var;
			evens=2:2:num_var;
			for ii=1:num_var/2
				pn=opt_args{odds(ii)};
				pv=opt_args{evens(ii)};
				switch pn
				case 'min_flat_area'
					opt_args{evens(ii)}=str2num(pv);
				case 'new_cellsize'
					opt_args{evens(ii)}=str2num(pv);
				case 'resample_grid'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end
				end
			end
			cmpMakeStreams(wdir,req_args{1},str2num(req_args{2}),req_args{3},opt_args);
		end
		disp('Function Successfully Completed')
	case 'ConditionDEM'
		num_req=2;
		req_args=varargin(1:num_req);
		if nf_args==num_req
			cmpConditionDEM(wdir,req_args{1},req_args{2});
		else
			opt_args=varargin(num_req+1:end);
			num_var=numel(opt_args);
			odds=1:2:num_var;
			evens=2:2:num_var;
			for ii=1:num_var/2
				pn=opt_args{odds(ii)};
				pv=opt_args{evens(ii)};
				switch pn
				case 'fillp'
					opt_args{evens(ii)}=str2num(pv);
				case 'ming'
					opt_args{evens(ii)}=str2num(pv);
				case 'tau'
					opt_args{evens(ii)}=str2num(pv);
				case 'split'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end
				case 'stiffness'
					opt_args{evens(ii)}=str2num(pv);
				case 'stiff_tribs'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end	
				case 'positive'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end
				case 'imposemin'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end
				case 'attachtomin'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end
				case 'attachheads'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end
				case 'discardflats'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end
				case 'maxcurvature'
					opt_args{evens(ii)}=str2num(pv);
				end
			end
			cmpConditionDEM(wdir,req_args{1},req_args{2},opt_args);
		end
		disp('Function Successfully Completed')
	case 'RemoveFlats'
		num_req=3;
		req_args=varargin;
		cmpRemoveFlats(wdir,req_args{1},str2num(req_args{2}),req_args{3});
		disp('Function Successfully Completed')
	case 'FindThreshold'
		num_req=2;
		req_args=varargin(1:num_req);
		if nf_args==num_req
			if strcmp(req_args{2},'all')
				cmpFindThreshold(wdir,req_args{1},req_args{2});
			elseif strcmp(req_args{2},'auto')
				cmpFindThreshold(wdir,req_args{1},req_args{2});
			else
				cmpFindThreshold(wdir,req_args{1},str2num(req_args{2}));
			end
		else
			opt_args=varargin(num_req+1:end);
			num_var=numel(opt_args);
			odds=1:2:num_var;
			evens=2:2:num_var;
			for ii=1:num_var/2
				pn=opt_args{odds(ii)};
				pv=opt_args{evens(ii)};
				switch pn
				case 'ref_concavity'
					opt_args{evens(ii)}=str2num(pv);
				case 'max_threshold'
					opt_args{evens(ii)}=str2num(pv);
				case 'remake_network'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end
				end
			end

			if strcmp(req_args{2},'all')
				cmpFindThreshold(wdir,req_args{1},req_args{2},opt_args);
			elseif strcmp(req_args{2},'auto')
				cmpFindThreshold(wdir,req_args{1},req_args{2},opt_args);
			else
				cmpFindThreshold(wdir,req_args{1},str2num(req_args{2}),opt_args);
			end
		end
		disp('Function Successfully Completed')
	case 'SegmentPicker'
		num_req=2;
		req_args=varargin(1:num_req);
		if nf_args==num_req
			cmpSegmentPicker(wdir,req_args{1},str2num(req_args{2}));
		else
			opt_args=varargin(num_req+1:end);
			num_var=numel(opt_args);
			odds=1:2:num_var;
			evens=2:2:num_var;
			for ii=1:num_var/2
				pn=opt_args{odds(ii)};
				pv=opt_args{evens(ii)};
				switch pn
				case 'ref_concavity'
					opt_args{evens(ii)}=str2num(pv);
				case 'min_elev'
					opt_args{evens(ii)}=str2num(pv);
				case 'max_area'
					opt_args{evens(ii)}=str2num(pv);
				case 'threshold_area'
					opt_args{evens(ii)}=str2num(pv);
				case 'interp_value'
					opt_args{evens(ii)}=str2num(pv);
				case 'bin_size'
					opt_args{evens(ii)}=str2num(pv);
				case 'complete_networks_only'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end
				case 'recalc'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end
				end
			end
			cmpSegmentPicker(wdir,req_args{1},str2num(req_args{2}),opt_args);
		end
		disp('Function Successfully Completed')
	case 'SegmentPlotter'
		num_req=1;
		req_args=varargin(1:num_req);
		if nf_args==num_req
			cmpSegmentPlotter(wdir,str2num(req_args{1}));
		else
			opt_args=varargin(num_req+1:end);
			num_var=numel(opt_args);
			odds=1:2:num_var;
			evens=2:2:num_var;
			for ii=1:num_var/2
				pn=opt_args{odds(ii)};
				pv=opt_args{evens(ii)};
				switch pn
				case 'subset'
					opt_args{evens(ii)}=str2num(pv);
				case 'separate'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end
				case 'label'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end
				end
			end
			cmpSegmentPlotter(wdir,str2num(req_args{1}),opt_args);
		end
		disp('Function Successfully Completed')
	case 'SegmentProjector'
		num_req=1;
		req_args=varargin(1:num_req);
		if nf_args==num_req
			cmpSegmentProjector(wdir,req_args{1});
		else
			opt_args=varargin(num_req+1:end);
			num_var=numel(opt_args);
			odds=1:2:num_var;
			evens=2:2:num_var;
			for ii=1:num_var/2
				pn=opt_args{odds(ii)};
				pv=opt_args{evens(ii)};
				switch pn
				case 'ref_concavity'
					opt_args{evens(ii)}=str2num(pv);
				case 'interp_value'
					opt_args{evens(ii)}=str2num(pv);
				case 'save_figures'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end
				case 'refit_streams'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end
				end
			end
			cmpSegmentProjector(wdir,req_args{1},opt_args);
		end
		disp('Function Successfully Completed')
	case 'BasinPicker'
		num_req=1;
		req_args=varargin(1:num_req);
		if nf_args==num_req
			cmpBasinPicker(wdir,req_args{1});
		else
			opt_args=varargin(num_req+1:end);
			num_var=numel(opt_args);
			odds=1:2:num_var;
			evens=2:2:num_var;
			for ii=1:num_var/2
				pn=opt_args{odds(ii)};
				pv=opt_args{evens(ii)};
				switch pn
				case 'ref_concavity'
					opt_args{evens(ii)}=str2num(pv);
				case 'rlf_radius'
					opt_args{evens(ii)}=str2num(pv);
				case 'interp_value'
					opt_args{evens(ii)}=str2num(pv);
				case 'threshold_area'
					opt_args{evens(ii)}=str2num(pv);
				end
			end
			cmpBasinPicker(wdir,req_args{1},opt_args);
		end
		disp('Function Successfully Completed')
	case 'KsnChiBatch'
		num_req=2;
		req_args=varargin(1:num_req);
		if nf_args==num_req
			cmpKsnChiBatch(wdir,req_args{1},req_args{2});
		else
			opt_args=varargin(num_req+1:end);
			num_var=numel(opt_args);
			odds=1:2:num_var;
			evens=2:2:num_var;
			for ii=1:num_var/2
				pn=opt_args{odds(ii)};
				pv=opt_args{evens(ii)};
				switch pn
				case 'segment_length'
					opt_args{evens(ii)}=str2num(pv);
				case 'smooth_distance'
					opt_args{evens(ii)}=str2num(pv);
				case 'ref_concavity'
					opt_args{evens(ii)}=str2num(pv);
				case 'min_order'
					opt_args{evens(ii)}=str2num(pv);
				case 'min_elevation'
					opt_args{evens(ii)}=str2num(pv);
				case 'interp_value'
					opt_args{evens(ii)}=str2num(pv);
				case 'radius'
					opt_args{evens(ii)}=str2num(pv);
				case 'complete_networks_only'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end
				end
			end
			cmpKsnChiBatch(wdir,req_args{1},req_args{2},opt_args);
		end
		disp('Function Successfully Completed')
	case 'KsnProfiler'
		num_req=1;
		req_args=varargin(1:num_req);
		if nf_args==num_req
			cmpKsnProfiler(wdir,req_args{1});
		else
			opt_args=varargin(num_req+1:end);
			num_var=numel(opt_args);
			odds=1:2:num_var;
			evens=2:2:num_var;
			for ii=1:num_var/2
				pn=opt_args{odds(ii)};
				pv=opt_args{evens(ii)};
				switch pn
				case 'min_channel_length'
					opt_args{evens(ii)}=str2num(pv);
				case 'min_length_to_extract'
					opt_args{evens(ii)}=str2num(pv);
				case 'min_elev'
					opt_args{evens(ii)}=str2num(pv);
				case 'max_area'
					opt_args{evens(ii)}=str2num(pv);
				case 'interp_value'
					opt_args{evens(ii)}=str2num(pv);
				case 'ref_concavity'
					opt_args{evens(ii)}=str2num(pv);
				case 'smooth_distance'
					opt_args{evens(ii)}=str2num(pv);
				case 'max_ksn'
					opt_args{evens(ii)}=str2num(pv);
				case 'threshold_area'
					opt_args{evens(ii)}=str2num(pv);
				case 'redefine_threshold'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end
				case 'complete_networks_only'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end
				case 'display_slope_area'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end
				case 'save_figures'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end
				end
			end
			cmpKsnProfiler(wdir,req_args{1},opt_args);
		end
		disp('Function Successfully Completed')
	case 'ClassifyKnicks'
		num_req=2;
		req_args=varargin(1:num_req);
		if nf_args==num_req
			cmpClassifyKnicks(wdir,req_args{1},req_args{2});
		else
			opt_args=varargin(num_req+1:end);
			num_var=numel(opt_args);
			cmpClassifyKnicks(wdir,req_args{1},req_args{2},opt_args);
		end
		disp('Function Successfully Completed')
	case 'ProcessRiverBasins'
		num_req=3;
		req_args=varargin(1:num_req);
		if nf_args==num_req
			if ~isempty(regexp(req_args{2},regexptranslate('wildcard','*.shp'))) 
				cmpProcessRiverBasins(wdir,req_args{1},req_args{2},req_args{3});
			elseif ~isempty(regexp(req_args{2},regexptranslate('wildcard','*.txt'))) 
				cmpProcessRiverBasins(wdir,req_args{1},req_args{2},req_args{3});
			else
				cmpProcessRiverBasins(wdir,req_args{1},str2num(req_args{2}),req_args{3});
			end
		else
			opt_args=varargin(num_req+1:end);
			num_var=numel(opt_args);
			odds=1:2:num_var;
			evens=2:2:num_var;
			for ii=1:num_var/2
				pn=opt_args{odds(ii)};
				pv=opt_args{evens(ii)};
				switch pn
				case 'ref_concavity'
					opt_args{evens(ii)}=str2num(pv);
				case 'relief_radii'
					opt_args{evens(ii)}=str2num(pv);
				case 'segment_length'
					opt_args{evens(ii)}=str2num(pv);
				case 'interp_value'
					opt_args{evens(ii)}=str2num(pv);
				case 'threshold_area'
					opt_args{evens(ii)}=str2num(pv);
				case 'min_order'
					opt_args{evens(ii)}=str2num(pv);
				case 'ksn_radius'
					opt_args{evens(ii)}=str2num(pv);
				case 'write_arc_files'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end
				case 'calc_relief'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end
				end
			end

			if ~isempty(regexp(req_args{2},regexptranslate('wildcard','*.shp'))) 
				cmpProcessRiverBasins(wdir,req_args{1},req_args{2},req_args{3},opt_args);
			elseif ~isempty(regexp(req_args{2},regexptranslate('wildcard','*.txt'))) 
				cmpProcessRiverBasins(wdir,req_args{1},req_args{2},req_args{3},opt_args);
			else
				cmpProcessRiverBasins(wdir,req_args{1},str2num(req_args{2}),req_args{3},opt_args);
			end
		end
		disp('Function Successfully Completed')
	case 'PrepareAddGrids'
		num_req=1;
		req_args=varargin(1:num_req);
		opt_args=varargin(num_req+1:end);
		cmpPrepareAddGrids(wdir,req_args{1},opt_args);
		disp('Function Successfully Completed')
	case 'PrepareAddCatGrids'
		num_req=2;
		req_args=varargin(1:num_req);
		opt_args=varargin(num_req+1:end);
		cmpPrepareCatAddGrids(wdir,req_args{1},req_args{2},opt_args);
		disp('Function Successfully Completed')
	case 'SubDivideBigBasins'
		num_req=3;
		req_args=varargin(1:num_req);
		if nf_args==num_req
			cmpSubDivideBigBasins(wdir,req_args{1},str2num(req_args{2}),req_args{3});
		else
			opt_args=varargin(num_req+1:end);
			num_var=numel(opt_args);
			odds=1:2:num_var;
			evens=2:2:num_var;
			for ii=1:num_var/2
				pn=opt_args{odds(ii)};
				pv=opt_args{evens(ii)};
				switch pn
				case 'ref_concavity'
					opt_args{evens(ii)}=str2num(pv);
				case 'threshold_area'
					opt_args{evens(ii)}=str2num(pv);
				case 'min_basin_size'
					opt_args{evens(ii)}=str2num(pv);
				case 's_order'
					opt_args{evens(ii)}=str2num(pv);
				case 'recursive'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end
				case 'write_arc_files'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end
				case 'no_nested'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end
				end
			end
			cmpSubDivideBigBasins(wdir,req_args{1},str2num(req_args{2}),req_args{3},opt_args);
		end
		disp('Function Successfully Completed')
	case 'FindBasinKnicks'
		num_req=3;
		req_args=varargin(1:num_req);
		if strcmpi(req_args{3},'false')
			req_args{3}=false;
		elseif strcmpi(req_args{3},'true')
			req_args{3}=true;
		else
			req_args{3}=logical(str2num(pv));
		end

		if nf_args==num_req
			cmpFindBasinKnicks(wdir,req_args{1},req_args{2},req_args{3});
		else
			opt_args=varargin(num_req+1:end);
			num_var=numel(opt_args);
			odds=1:2:num_var;
			evens=2:2:num_var;
			for ii=1:num_var/2
				pn=opt_args{odds(ii)};
				pv=opt_args{evens(ii)};
				switch pn
				case 'ref_concavity'
					opt_args{evens(ii)}=str2num(pv);
				case 'classify_knicks'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end
				end
			end
			cmpFindBasinKnicks(wdir,req_args{1},req_args{2},req_args{3},opt_args);
		end
		disp('Function Successfully Completed')
	case 'BasinStatsPlots'
		num_req=2;
		req_args=varargin(1:num_req);
		if nf_args==num_req
			cmpBasinStatsPlots(wdir,req_args{1},req_args{2});
		else
			opt_args=varargin(num_req+1:end);
			num_var=numel(opt_args);
			odds=1:2:num_var;
			evens=2:2:num_var;
			for ii=1:num_var/2
				pn=opt_args{odds(ii)};
				pv=opt_args{evens(ii)};
				switch pn
				case 'use_filtered'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end
				case 'save_figure'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end
				case 'define_region'
					opt_args{evens(ii)}=str2num(pv);
				case 'rlf_radius'
					opt_args{evens(ii)}=str2num(pv);
				case 'basin_num'
					opt_args{evens(ii)}=str2num(pv);
				case 'fit_grd_ksn'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end
				case 'start_diffusivity'
					opt_args{evens(ii)}=str2num(pv);
				case 'start_erodibility'
					opt_args{evens(ii)}=str2num(pv);
				case 'start_threshold_gradient'
					opt_args{evens(ii)}=str2num(pv);
				case 'n_val'
					opt_args{evens(ii)}=str2num(pv);
				case 'fit_rlf_ksn'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end
				case 'fit_filtered'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end
				end
			end
			cmpBasinStatsPlots(wdir,req_args{1},req_args{2},opt_args);
		end
		disp('Function Successfully Completed')
	case 'PlotIndividualBasins'
		num_req=1;
		req_args=varargin(1:num_req);
		if nf_args==num_req
			cmpPlotIndividualBasins(wdir,req_args{1});
		else
			opt_args=varargin(num_req+1:end);
			num_var=numel(opt_args);
			odds=1:2:num_var;
			evens=2:2:num_var;
			for ii=1:num_var/2
				pn=opt_args{odds(ii)};
				pv=opt_args{evens(ii)};
				switch pn
				case 'bin_size'
					opt_args{evens(ii)}=str2num(pv);
				end
			end
			cmplotIndividualBasins(wdir,req_args{1},opt_args);
		end
		disp('Function Successfully Completed')
	case 'Basin2Raster'
		num_req=3;
		req_args=varargin(1:num_req);
		if nf_args==num_req
			cmpBasin2Raster(wdir,req_args{1},req_args{2},req_args{3});
		else
			opt_args=varargin(num_req+1:end);
			num_var=numel(opt_args);
			odds=1:2:num_var;
			evens=2:2:num_var;
			for ii=1:num_var/2
				pn=opt_args{odds(ii)};
				pv=opt_args{evens(ii)};
				switch pn
				case 'relief_radius'
					opt_args{evens(ii)}=str2num(pv);
				end
			end
			cmpBasin2Raster(wdir,req_args{1},req_args{2},req_args{3},opt_args);
		end
		disp('Function Successfully Completed')
	case 'Basin2Shape'
		num_req=2;
		req_args=varargin(1:num_req);
		if nf_args==num_req
			cmpBasin2Shape(wdir,req_args{1},req_args{2});
		else
			opt_args=varargin(num_req+1:end);
			num_var=numel(opt_args);
			odds=1:2:num_var;
			evens=2:2:num_var;
			for ii=1:num_var/2
				pn=opt_args{odds(ii)};
				pv=opt_args{evens(ii)};
				switch pn
				case 'new_concavity'
					opt_args{evens(ii)}=str2num(pv);	
				case 'populate_categories'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end
				end
			end
			cmpBasin2Shape(wdir,req_args{1},req_args{2},opt_args);
		end
		disp('Function Successfully Completed')
	case 'CompileBasinStats'
		num_req=1;
		req_args=varargin(1:num_req);
		if nf_args==num_req
			cmpCompileBasinStats(wdir,req_args{1});
		else
			opt_args=varargin(num_req+1:end);
			num_var=numel(opt_args);
			odds=1:2:num_var;
			evens=2:2:num_var;
			for ii=1:num_var/2
				pn=opt_args{odds(ii)};
				pv=opt_args{evens(ii)};
				switch pn
				case 'dist_along_azimuth'
					opt_args{evens(ii)}=str2num(pv);
				case 'new_concavity'
					opt_args{evens(ii)}=str2num(pv);
				case 'filter_by_category'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end
				case 'populate_categories'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end
				end
			end
			cmpCompileBasinStats(wdir,req_args{1},opt_args);
		end
		disp('Function Successfully Completed')
	case 'MakeTopoSwath'
		num_req=3;
		req_args=varargin(1:num_req);
		if nf_args==num_req
			cmpMakeTopoSwath(wdir,req_args{1},req_args{2},str2num(req_args{3}));
		else
			opt_args=varargin(num_req+1:end);
			num_var=numel(opt_args);
			odds=1:2:num_var;
			evens=2:2:num_var;
			for ii=1:num_var/2
				pn=opt_args{odds(ii)};
				pv=opt_args{evens(ii)};
				switch pn
				case 'sample'
					opt_args{evens(ii)}=str2num(pv);
				case 'smooth'
					opt_args{evens(ii)}=str2num(pv);
				case 'vex'
					opt_args{evens(ii)}=str2num(pv);
				case 'save_figure'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end
				case 'plot_as_points'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end
				case 'plot_as_heatmap'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end
				end
			end
			cmpMakeTopoSwath(wdir,req_args{1},req_args{2},str2num(req_args{3}),opt_args);
		end
		disp('Function Successfully Completed')
	case 'MakeCombinedSwath'
		num_req=6;
		req_args=varargin(1:num_req);
		if nf_args==num_req
			cmpMakeCombinedSwath(wdir,req_args{1},req_args{2},str2num(req_args{3}),req_args{4},...
				req_args{5},str2num(req_args{6}));
		else
			opt_args=varargin(num_req+1:end);
			num_var=numel(opt_args);
			odds=1:2:num_var;
			evens=2:2:num_var;
			for ii=1:num_var/2
				pn=opt_args{odds(ii)};
				pv=opt_args{evens(ii)};
				switch pn
				case 'sample'
					opt_args{evens(ii)}=str2num(pv);
				case 'smooth'
					opt_args{evens(ii)}=str2num(pv);
				case 'vex'
					opt_args{evens(ii)}=str2num(pv);
				case 'plot_map'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end
				case 'save_figure'
					if strcmpi(pv,'false')
						opt_args{evens(ii)}=false;
					elseif strcmpi(pv,'true')
						opt_args{evens(ii)}=true;
					else
						opt_args{evens(ii)}=logical(str2num(pv));
					end
				end
			end
			cmpMakeCombinedSwath(wdir,req_args{1},req_args{2},str2num(req_args{3}),req_args{4},...
				req_args{5},str2num(req_args{6}),opt_args);
		end
		disp('Function Successfully Completed')
	case 'DippingBedFinder'
		num_req=6;
		req_args=varargin(1:num_req);
		if nf_args==num_req
			cmpDippingBedFinder(wdir,req_args{1},str2num(req_args{2}),str2num(req_args{3}),...
				str2num(req_args{4}),str2num(req_args{5}),str2num(req_args{6}));
		end
		disp('Function Successfully Completed')
	case 'Mat2Arc'
		num_req=2;
		req_args=varargin(1:num_req);
		if nf_args==num_req
			cmpMat2Arc(wdir,req_args{1},req_args{2});
		else
			opt_args=varargin(num_req+1:end);
			cmpMat2Arc(wdir,req_args{1},req_args{2},opt_args);
		end
		disp('Function Successfully Completed')
	case 'PlotKsn'
		num_req=2;
		req_args=varargin(1:num_req);
		if nf_args==num_req
			cmpPlotKsn(wdir,req_args{1},req_args{2});
		else
			opt_args=varargin(num_req+1:end);
			cmpPlotKsn(wdir,req_args{1},req_args{2},opt_args);
		end
		disp('Function Successfully Completed')
	case 'CatPoly2GRIDobj'
		disp('There is no compiled version of the CatPoly2GRIDobj function, use the compiled PrepareAddCatGrids function')
	case 'ksncolor'
		disp('There is no compiled version of the ksncolor function')
	case 'ProjectOntoSwath'
		disp('There is no compiled version of the ProjectOntoSwath function, use the compiled MakeCombinedSwath function')
	otherwise
		disp([func_name ' is not a recognized function name within TAK'])
	%Main switch end	
	end

	% Call colormaps to make sure they're included
	[col]=ksncolor(20);
	[col]=flowcolor(20);
	[col]=landcolor(20);
	[col]=magmacolor(20);


% Function end
end