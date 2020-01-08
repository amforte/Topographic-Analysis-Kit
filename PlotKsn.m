function PlotKsn(DEM,FD,ksn,varargin)
	%
	% Usage:
	%	PlotKsn(DEM,FD,ksn);
	%	PlotKsn(DEM,FD,ksn,bnd_list);
	%	PlotKsn(DEM,FD,ksn,'knicks.shp');
	%
	% Description:
	% 	Function to plot a map of normalized channel steepness on a hillshade colored
	% 	by elevation.
	%
	% Required Inputs:
	%	DEM -  Digital Elevation as a GRIDobj used to produce the provided ksn data
	%	FD - Flow routing as a FLOWobj used to proudce the provided ksn data
	%	ksn - ksn data either as a shapefile (as ouput from KsnProfiler, ProcessRiverBasins
	%		KsnChiBatch), an ascii file (for continuous Ksn values as output from KsnChiBatch),
	%		a GRIDobj of continous ksn (as output from KsnChiBatch), or a mapstructure 
	%		(as output from ProcessRiverBasins or KsnChiBatch)
	% 
	% Optional Inputs:
	%	knicks [] - location of knickpoints as an array (as output from FindBasinKnicks
	%		or KsnProfiler 'bnd_list') or shapefile (as output by FindBasinKnicks or KsnProfiler)
	%	ksn_lim [] - 1 x n vector setting the min and max for the color scaling for ksn. If left blank
	%		will default to 0 and the maximum in the dataset
	%
	% Examples:
	%	PlotKsn(DEMoc,FDc,MSNc); % Plot ksn map of a basin from ProcessRiverBasins
	%	PlotKsn(DEM,FD,'ksn.shp'); % Plot ksn map from a shapefile
	%	PlotKsn(DEM,FD,'ksn.shp','knicks','ksn_knicks.shp'); % Include knickpoint locations from KsnProfiler
	%	PlotKsn(DEMoc,FDc,MSNc,,'knicks',KnickPoints); %Include knickpoints output from FindBasinKnicks
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 06/18/18 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if ischar(ksn)
		[~,~,ext]=fileparts(ksn);
	else
		ext=' ';
	end

    % Parse Inputs
    p = inputParser;
    p.FunctionName = 'PlotKsn';
    addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
    addRequired(p,'FD',@(x) isa(x,'FLOWobj'));
    addRequired(p,'ksn',@(x) isstruct(x) || strcmpi(ext,'.shp') || strcmpi(ext,'.txt') || isa(x,'GRIDobj'));

    addParameter(p,'knicks',[],@(x) isnumeric(x) || regexp(x,regexptranslate('wildcard','*.shp')) || istable(x));
    addParameter(p,'ksn_lim',[],@(x) isnumeric(x) && numel(x)==2);

    parse(p,DEM,FD,ksn,varargin{:});
    DEM=p.Results.DEM;
    FD=p.Results.FD;
    ksn=p.Results.ksn;

    knks=p.Results.knicks;
    ksn_lim=p.Results.ksn_lim;


	if ischar(ksn) & logical(regexp(ksn,regexptranslate('wildcard','*.shp')))
		ksn=shaperead(ksn);
		grid_flag=false;
	elseif ischar(ksn) & logical(regexp(ksn,regexptranslate('wildcard','*.txt')))
		ksn=GRIDobj(ksn);
		if ~validatealignment(DEM,ksn)
			ksn=resample(ksn,DEM);
		end
		grid_flag=true;
	elseif isstruct(ksn)
		ksn=ksn;
		grid_flag=false;
	elseif isa(ksn,'GRIDobj');
		if ~validatealignment(DEM,ksn)
			ksn=resample(ksn,DEM);
		end
		grid_flag=true;
	else
		if isdeployed
			errordlg('Input to "ksn" not recognized as a shapefile or a mapstructure')
		end
		error('Input to "ksn" not recognized as a shapefile or a mapstructure');
	end
	
	if grid_flag
		f1=figure(1);
		set(f1,'Units','normalized','Position',[0.05 0.1 0.8 0.8],'renderer','painters');
		hold on
		if isempty(ksn_lim)
			imageschs(DEM,ksn,'colormap','ksncolor');
		else
			imageschs(DEM,ksn,'colormap','ksncolor','caxis',ksn_lim);
		end

        if ~verLessThan('matlab','9.5')
            disableDefaultInteractivity(gca);
        end 
		hold off
	else	
		num_seg=numel(ksn);

		sx=cell(num_seg,1);
		sy=cell(num_seg,1);
		sk=cell(num_seg,1);
		for ii=1:num_seg
			sx{ii,1}=ksn(ii,1).X(:);
			sy{ii,1}=ksn(ii,1).Y(:);
			if isfield(ksn,'ksn')
				sk{ii,1}=ones(numel(sx{ii,1}),1)*ksn(ii,1).ksn;
			elseif isfield(ksn,'fit_ksn')
				sk{ii,1}=ones(numel(sx{ii,1}),1)*ksn(ii,1).fit_ksn;
			else
				if isdeployed
					errordlg('There is no valid field in the provided shapefile or mapstructure named "ksn"')
				end
				error('There is no valid field in the provided shapefile or mapstructure named "ksn"')
			end
		end

		sx=vertcat(sx{:});
		sy=vertcat(sy{:});
		sk=vertcat(sk{:});

		ix=coord2ind(DEM,sx,sy);
		idx=isnan(ix);

		ix(idx)=[];
		sk(idx)=[];

		W=GRIDobj(DEM,'logical');
		W.Z(ix)=true;
		S=STREAMobj(FD,W);

		[~,loc,~]=unique(ix);
		sk=sk(loc);

		f1=figure(1);
		set(f1,'Visible','off');

		[RGB]=imageschs(DEM,DEM,'colormap','gray');
		[~,R]=GRIDobj2im(DEM);

		imshow(flipud(RGB),R);
		axis xy
		hold on
		colormap(ksncolor(20));
		plotc(S,sk);
		if isempty(ksn_lim)
			caxis([0 max(sk)]);
		else
			caxis([min(ksn_lim) max(ksn_lim)]);
		end
		c1=colorbar;
		ylabel(c1,'Normalized Channel Steepness')
		if ~isempty(knks)
			if ischar(knks) & logical(regexp(knks,regexptranslate('wildcard','*.shp')))
				knk=shaperead(knks);
				knkx=[knk.X];
				knky=[knk.Y];
				scatter(knkx,knky,100,'w','p','filled','MarkerEdgeColor','k');
			elseif istable(knks)
				knkx=knks.x_coord;
				knky=knks.y_coord;
				scatter(knkx,knky,100,'w','p','filled','MarkerEdgeColor','k');
			else
				scatter(knks(:,1),knks(:,2),100,'w','p','filled','MarkerEdgeColor','k');
			end
		end
        if ~verLessThan('matlab','9.5')
            disableDefaultInteractivity(gca);
        end 
		hold off
		set(f1,'Visible','on','Units','normalized','Position',[0.05 0.1 0.8 0.8],'renderer','painters');
	end
end