function BasinStatsPlots(basin_table,plots,varargin)
	% Function to take the outputs from 'ProcessRiverBasins' and 'SubDivideBigBasins' and produce various plots
	%	of aggregated basin values.


	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'BasinStatsPlots';
	addRequired(p,'basin_table',@(x) isa(x,'table'));
	addRequired(p,'plots',@(x) ischar(validatestring(x,{'grd_ksn','xy'})));

	addParamValue(p,'uncertainty','se',@(x) ischar(validatestring(x,{'se','std'})));
	addParamValue(p,'color_by',[],@(x) ischar(x));
	addParamValue(p,'cmap','jet',@(x) ischar(x) || isnumeric(x) & size(x,2)==3);
	addParamValue(p,'xval',[],@(x) ischar(x));
	addParamValue(p,'yval',[],@(x) ischar(x));

	parse(p,basin_table,plots,varargin{:});
	T=p.Results.basin_table;
	plts=p.Results.plots;

	uncertainty=p.Results.uncertainty;
	color_by=p.Results.color_by;
	cmap=p.Results.cmap;
	xval=p.Results.xval;
	yval=p.Results.yval;


	% Generate Plots

	switch plts
	case 'grd_ksn'
		g=T.mean_gradient;
		k=T.mean_ksn;

		if ismember('std_ksn',T.Properties.VariableNames) & strcmp(uncertainty,'std')
			sk=T.std_ksn;
			sg=T.std_gradient;
		elseif ismember('se_ksn',T.Properties.VariableNames) & strcmp(uncertainty,'se')
			sk=T.se_ksn;
			sg=T.se_gradient;
		elseif ismember('std_ksn',T.Properties.VariableNames) & ~ismember('se_ksn',T.Properties.VariableNames)
			sk=T.std_ksn;
			sg=T.std_gradient;
		elseif ismember('se_ksn',T.Properties.VariableNames) & ~ismember('std_ksn',T.Properties.VariableNames)
			sk=T.se_ksn;
			sg=T.se_gradient;
		end

		f1=figure(1);
		hold on 
		errorbar(k,g,sg,sg,sk,sk,'.k','CapSize',0);

		if ~isempty(color_by);
			colormap(cmap);
			scatter(k,g,30,T.(color_by),'filled');
			cb=colorbar;
			% Remove any underscores
			color_by_label=strrep(color_by,'_',' ');
			ylabel(cb,color_by_label);
		else
			scatter(k,g,30,'k','filled');
		end

		xlabel('Mean Basin k_{sn}');
		ylabel('Mean Basin Gradient');
		hold off
	case 'xy'
		if isempty(xval) | isempty(yval)
			error('For "xy" plot, entries for "xval" and "yval" are required');
		end

		x=T.(xval);
		y=T.(yval);

		sxN=strrep(xval,'mean_',''); 
		syN=strrep(yval,'mean_','');

		if ismember(['std_' sxN],T.Properties.VariableNames) & strcmp(uncertainty,'std')
			sx=T.(['std_' sxN]);
		elseif ismember(['se_' sxN],T.Properties.VariableNames) & strcmp(uncertainty,'se')
			sx=T.(['se_' sxN]);
		elseif ismember(['std_' sxN],T.Properties.VariableNames) & ~ismember(['se_' sxN],T.Properties.VariableNames)
			sx=T.(['std_' sxN]);
		elseif ismember(['se_' sxN],T.Properties.VariableNames) & ~ismember(['std_' sxN],T.Properties.VariableNames)
			sx=T.(['se_' sxN]);
		else
			sx=[];
		end

		if ismember(['std_' syN],T.Properties.VariableNames) & strcmp(uncertainty,'std')
			sy=T.(['std_' syN]);
		elseif ismember(['se_' syN],T.Properties.VariableNames) & strcmp(uncertainty,'se')
			sy=T.(['se_' syN]);
		elseif ismember(['std_' syN],T.Properties.VariableNames) & ~ismember(['se_' syN],T.Properties.VariableNames)
			sy=T.(['std_' syN]);
		elseif ismember(['se_' syN],T.Properties.VariableNames) & ~ismember(['std_' syN],T.Properties.VariableNames)
			sy=T.(['se_' syN]);
		else
			sy=[];
		end

		f1=figure(1);
		hold on 

		if ~isempty(sx) & isempty(sy)
			errobar(x,y,sx,'horizontal','.k','CapSize',0);
		elseif isempty(sy) & ~isempty(sx)
			errobar(x,y,sy,'vertical','.k','CapSize',0);
		elseif ~isempty(sx) & ~isempty(sy)
			errorbar(x,y,sy,sy,sx,sx,'.k','CapSize',0);
		end

		if ~isempty(color_by);
			colormap(cmap);
			scatter(x,y,30,T.(color_by),'filled');
			cb=colorbar;
			% Remove any underscores
			color_by_label=strrep(color_by,'_',' ');
			ylabel(cb,color_by_label);
		else
			scatter(x,y,30,'k','filled');
		end

		xlabel(['Mean ' sxN]);
		ylabel(['Mean ' syN]);
		hold off


	end
