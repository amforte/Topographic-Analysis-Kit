function PlotChi(DEM,S,chi,chi_type,varargin)
	%
	% Usage:
	%	PlotChi(DEM,S,chi,'chi_map');
	%	PlotChi(DEM,S,chi,'chi_grid');
	%	PlotChi(DEM,S,chi,'chi_map','name',value,...);
	%
	% Description:
	% 	Function to plot a map of normalized channel steepness on a hillshade colored
	% 	by elevation.
	%
	% Required Inputs:
	%	DEM -  Digital Elevation as a GRIDobj used to produce the provided chi data
	%	S - STREAMobj used when producing chi values
	%	chi - chi data either as an ascii file or a GRIDobj (e.g. as output from KsnChiBatch) 
	%	chi_type - value indicating if the provided chi data is a 'chimap' (along streams only) or
	%		a continuous grid 'chigrid'
	% 
	% Optional Inputs:
	%	chi_lim [] - 1 x n vector setting the min and max for the color scaling for chi. If left blank
	%		will default to 0 and the maximum in the dataset
	%
	% Examples:
	%	PlotChi(DEM,S,chimap,'chimap'); % Plot chi from an output chimap (as a GRIDobj)
	%	PlotChi(DEM,S,'Topo_chigrid.txt','chigrid','chi_lim',[0 10]);
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 09/30/19 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   % Parse Inputs
    p = inputParser;
    p.FunctionName = 'PlotChi';
    addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
    addRequired(p,'S',@(x) isa(x,'STREAMobj'));
    addRequired(p,'chi',@(x) isa(x,'GRIDobj') || regexp(x,regexptranslate('wildcard','*.txt')));
	addRequired(p,'chi_type',@(x) ischar(validatestring(x,{'chimap','chigrid'})));

    addParameter(p,'chi_lim',[],@(x) isnumeric(x) && numel(x)==2);
    addParameter(p,'override_resample',false,@(x) isscalar(x) && islogical(x)); % Hidden option for GUIs

    parse(p,DEM,S,chi,chi_type,varargin{:});
    DEM=p.Results.DEM;
    S=p.Results.S;
    chi=p.Results.chi;
    chi_type=p.Results.chi_type;

   	chi_lim=p.Results.chi_lim;
   	os=p.Results.override_resample;


	if ischar(chi) & logical(regexp(chi,regexptranslate('wildcard','*.txt')))
		chi=GRIDobj(chi);
		if ~validatealignment(DEM,chi) && ~os
			chi=resample(chi,DEM);
		elseif ~validatealignment(DEM,chi) && os
			chi.refmat=DEM.refmat;
			chi.georef=DEM.georef;
		end
	elseif isa(chi,'GRIDobj');
		if ~validatealignment(DEM,chi)
			chi=resample(chi,DEM);
		end
	else
		error('Input to "chi" not recognized as a valid ascii or a GRIDobj');
	end
	
	switch chi_type
	case 'chigrid'
		f1=figure(1);
		set(f1,'Units','normalized','Position',[0.05 0.1 0.8 0.8],'renderer','painters');
		hold on
		if isempty(chi_lim)
			imageschs(DEM,chi,'colormap','jet');
		else
			imageschs(DEM,chi,'colormap','jet','caxis',chi_lim);
		end
        if ~verLessThan('matlab','9.5')
            disableDefaultInteractivity(gca);
        end 
		hold off
	case 'chimap'
		nal=getnal(S,chi);

		f1=figure(1);
		set(f1,'Visible','off');

		[RGB]=imageschs(DEM,DEM,'colormap','gray');
		[~,R]=GRIDobj2im(DEM);

		imshow(flipud(RGB),R);
		axis xy
		hold on
		colormap(jet);
		plotc(S,nal);
		if isempty(chi_lim)
			caxis([0 max(nal)]);
		else
			caxis([min(chi_lim) max(chi_lim)])
		end
		c1=colorbar;
		ylabel(c1,'\chi');

        if ~verLessThan('matlab','9.5')
            disableDefaultInteractivity(gca);
        end 
		hold off
		set(f1,'Visible','on','Units','normalized','Position',[0.05 0.1 0.8 0.8],'renderer','painters');
	end



end