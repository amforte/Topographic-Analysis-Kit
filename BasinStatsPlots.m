function BasinStatsPlots(basin_table,plots,varargin)
	%
	% Usage:
	%	BasinStatsPlots(basin_table,plots);
	%	BasinStatsPlots(basin_table,plots,'name',value,...);
	%
	% Description:
	% 	Function to take the complied outputs from 'ProcessRiverBasins' and 'SubDivideBigBasins' and produce various plots
	%	of aggregated basin values. 
	%
	% Required inputs:
	%	basin_table - Table output from 'CompileBasinStats'
	%	plots - Type of plot you want to produce, valid inputs are:
	%		'grd_ksn' - plot of mean basin gradient vs mean basin channel steepness (e.g. see Forte et al, 2016, Earth and Planetary 
	%					Science Letters for discussion of use of these plots)
	%		'grd_rlf' - similar to 'grd_ksn' but uses local relief instead of ksn, requires that relief was calculated when running
	%					ProcessRiverBasins. Assumes relief radius is 2500 (can set alternative radii with 'rlf_radius' optional parameter)
	%		'rlf_ksn' - plot of mean basin relief vs mean basin channel steepness
	%		'compare_filtered' - plot comparing mean values vs filtered mean values if you ran 'CompileBasinStats' and filtered by a category
	%		'category_mean_hist' - if you calculated 'means_by_category' when running 'CompileBasinStats', you can plot distributions of the
	%					means by category as histograms using this option. Requires an input to 'cat_mean1'
	%		'category_mean_compare' -if you calculated 'means_by_category' for more than one value (e.g. both gradient and ksn), you can compare
	%					the mean values by category using this plot. Requires inputs to both 'cat_mean1' (value that will be plotted on x axis) 
	%					and 'cat_mean2' (value that will be plotted on y axis)
	%		'stacked_hypsometry' - plot hypsometries for the basins
	%		'compare_mean_and_dist' - plots a histogram of values within a selected basin or across all basins for a statistic of interest to 
	%					compare to the mean value, accepts an input for 'statistic_of_interest' and 'basin_num'.
	%		'scatterplot_matrix' - matrix of scatterplots and histograms, designed to be sort of similar to 'lattice' plots in R. Providing a table 
	%					for which you calculated filtered means and leaving 'use_filtered' set to false may produce a large matrix 
	%		'xy' - generic plot, requires entries to optional 'xval' and 'yval' inputs
	%
	%%%%%%%%%%%%%%%%%%
	% Optional Inputs:
	%
	%%% General Parameters
	%	uncertianty ['se'] - uncertainty to value use for plots, valid options are 'se' (standard error), 'std' (standard deviation), or 'none'. 
	%		Providing 'none' indicates you do not want to plot errorbars. Behavior of this option will depend on how you ran ProcessRiverBasins, 
	%		e.g. if you only calculated standard deviations when running ProcessRiverBasins but supply 'se'	here, the code will ignore your choice
	%		and use the standard deviation values.
	%	use_filtered [false] - logical flag to use filtered values for 'grd_ksn', 'grd_rlf', 'rlf_ksn', or 'scatterplot_matrix'. Will only work if 
	%		you calculated filtered values when running 'CompileBasinStats'.
	%	color_by [] - value to color points by, valid for 'grd_ksn','grd_rlf','rlf_ksn', and 'xy', either the name of a column in the provided table
	%		or a m x 1 array of numeric values the same length as the provided table
	%	cmap [] - colormap to use if an entry is provided to 'color_by', can be the name of a standard colormap or a nx3 array of rgb values
	%		to use as a colormap. 
	%	define_region [] - set of coordinates to define a rectangular region to draw data from, expects a four element matrix (row or column) that define
	%		the minimum x, maximum x, minimum y, and maximum y coordinate to include OR define as true to bring up a plot of all basin centers
	%		for you to select a region by drawing a rectangle. Works with all plots.
	%	rlf_radius [2500] - radius of relief used when plotting relief related values
	%	save_figure [false] - logical flag to save pdfs of all figures produced
	%
	%%% xy plot
	%	xval [] - value to plot on x axis for plot type 'xy' provided as name of column as it appears in the provided table or a a m x 1 array of numeric 
	%		values the same length as the provided table 
	%	yval [] - value to plot on y axis for plot type 'xy' provided as name of column as it appears in the provided table or a a m x 1 array of numeric 
	%		values the same length as the provided table 
	%
	%%% compare_mean_and_dist plot
	%	statistic_of_interest ['ksn'] - statistic of interest for plotting histogram to compare with mean value. Valid inputs are 'ksn', 'gradient',
	%		'elevation', 'relief' (if you provide relief, the code will look for relief calculated at the radius specified with the optional 'rlf_radius' parameter),
	%		or the name of an additional grid provided to 'ProcessRiverBasins', e.g. if you provided a precipitation grid and provided the name 'precip'
	%		and a column named 'mean_precip' exists in the table, then 'precip' would be a valid input to this parameter.
	%	basin_num [] - number of basin (as it appears in the ID column of the table) to use for 'compare_mean_and_dist', if empty, 'compare_mean_and_dist'
	%		will use all basins.
	%
	%%% category_mean_hist OR category_mean_compare 
	%	cat_mean1 [] - category to use for plotting, see 'category_mean_hist' or 'category_mean_compare', valid inputs are 'ksn', 'rlf', 'gradient', or 
	%		the name of an additional grid provided to ProcessRiverMeans.
	%	cat_mean2 [] - category to use for plotting, 'category_mean_compare' , valid inputs are 'ksn', 'rlf', 'gradient', or the name of an additional grid 
	%		provided to ProcessRiverMeans.
	%
	%%% Fit Gradient-Ksn Relationship
	%	fit_grd_ksn [false] - logical flag to initiate fitting of gradient - ksn relationship. Setting this flag to true only produces a result if the 
	%		plot type is set to 'grd_ksn'. The relationship is a fit using a power law relationship between erosion rate and channel steepness and erosion
	%		rate and mean hillslope gradient. See Forte et al, 2016, Earth and Planetary Science Letters for further discussion. The fit optimizes
	%		values of hillslope diffusivity (D), fluvial erodibility (K), and threshold gradient (Sc). The best fit values for these will be printed
	%		to the console.
	%	start_diffusivity [0.01] - starting value for optimization of hillslope diffusivity parameter.
	%	start_erodibility [1e-7] - starting value for optimization of fluvial erodibility parameter.
	%	start_threshold_gradient [0.8] - starting value for optimzation of threshold hillslope gradient parameter.
	%	n_val [2] - n value on slope parameter, this is not a free parameter in the fit.
	%
	%%% Fit Relief-Ksn Relationship
	%	fit_rlf_ksn [false] - logical flag to initiate simple linear fit of relief-ksn relationship (expectation is a linear relationship). Setting this
	%		flag to true only produces a result if the plot type is set to 'rlf_ksn'.
	%
	%%% Fit Filtered Data
	%	fit_filtered [false] - logical flat to initiate a simple linear fit to filtered data. Setting this to true only produces a result if plot type
	%		is set to 'compare_filtered'.
	%
	%%%%%%%%%%%%%%%%%%
	% Examples:
	%
	%	% Plot of mean basin gradient vs 2500 m^2 relief using default relief radius
	%	BasinStatsPlots(T,'grd_rlf');
	%
	%	% Plot of mean basin gradient vs mean basin channel steepness colored by mean elevation
	%	BasinStatsPlots(T,'grd_ksn','color_by','mean_el');
	%
	%	% Plot of mean basin gradient vs mean basin channel steepenss colored by mode geology, where colormap has been scaled to provide unique colors for units
	%	cmap=colorcube(numel(unique(T.mode_geology)));
	%	BasinStatsPlots(T,'grd_ksn','color_by','mode_geology','cmap',cmap);
	%
	%	% Plot of mean basin channel steepenss vs basin drainage area
	%	BasinStatsPlots(T,'xy','xval','drainage_area','yval','mean_ksn');
	%
	%	% Histograms of mean 2500 m^2 relief by individual categories (e.g. geology)
	%	BasinStatsPlots(T,'category_mean_hist','cat_mean1','rlf');
	%
	%	% Plots of mean basin gradient vs mean basin channel steepness within individual categories
	%	BasinStatsPlots(T,'category_mean_compare','cat_mean1','ksn','cat_mean2','gradient');
	% 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 06/18/18 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'BasinStatsPlots';
	addRequired(p,'basin_table',@(x) isa(x,'table'));
	addRequired(p,'plots',@(x) ischar(validatestring(x,{'grd_ksn','grd_rlf','rlf_ksn',...
		'compare_filtered','category_mean_hist','category_mean_compare','xy','stacked_hypsometry',...
		'compare_mean_and_dist','scatterplot_matrix'})));

	addParameter(p,'uncertainty','se',@(x) ischar(validatestring(x,{'se','std','none'})));
	addParameter(p,'use_filtered',false,@(x) islogical(x) && isscalar(x));
	addParameter(p,'color_by',[],@(x) ischar(x) || isnumeric(x) & size(x,2)==1 || isempty(x));
	addParameter(p,'cmap',[],@(x) ischar(x) || isnumeric(x) & size(x,2)==3);
	addParameter(p,'xval',[],@(x) ischar(x) || isnumeric(x) & size(x,2)==1);
	addParameter(p,'yval',[],@(x) ischar(x) || isnumeric(x) & size(x,2)==1);
	addParameter(p,'define_region',[],@(x) isnumeric(x) & numel(x)==4 || islogical(x));
	addParameter(p,'statistic_of_interest','ksn',@(x) ischar(x));
	addParameter(p,'basin_num',[],@(x) isnumeric(x) && isscalar(x));
	addParameter(p,'rlf_radius',2500,@(x) isnumeric(x) && isscalar(x));
	addParameter(p,'cat_mean1',[],@(x) ischar(x));
	addParameter(p,'cat_mean2',[],@(x) ischar(x));
	addParameter(p,'fit_grd_ksn',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'start_diffusivity',0.01,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'start_erodibility',1e-7,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'start_threshold_gradient',0.8,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'n_val',2,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'fit_rlf_ksn',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'fit_filtered',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'save_figure',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'out_dir',[],@(x) isdir(x));

	parse(p,basin_table,plots,varargin{:});
	T=p.Results.basin_table;
	plts=p.Results.plots;

	uncertainty=p.Results.uncertainty;
	use_filtered=p.Results.use_filtered;
	color_by=p.Results.color_by;
	cmap=p.Results.cmap;
	xval=p.Results.xval;
	yval=p.Results.yval;
	regionOI=p.Results.define_region;
	basin_num=p.Results.basin_num;
	stOI=p.Results.statistic_of_interest;
	rr=p.Results.rlf_radius;
	cm1=p.Results.cat_mean1;
	cm2=p.Results.cat_mean2;
	fit_grd_ksn=p.Results.fit_grd_ksn;
	fit_rlf_ksn=p.Results.fit_rlf_ksn;
	fit_filtered=p.Results.fit_filtered;
	D_in=p.Results.start_diffusivity;
	K_in=p.Results.start_erodibility;
	s_in=p.Results.start_threshold_gradient;
	n_val=p.Results.n_val;
	save_figure=p.Results.save_figure;
	out_dir=p.Results.out_dir;

	if isempty(out_dir)
		out_dir=pwd;
	end
 
	if isempty(cmap);
		cmap=parula(50);
	end

	% Deal with variable inputs
	if ~isempty(color_by) & isnumeric(color_by)
		cval=color_by;
		color_by='color_by';
		T.color_by=cval;
	end

	if ~isempty(xval) & isnumeric(xval);
		xv=xval;
		xval='xval';
		T.xval=xv;
	end

	if ~isempty(yval) & isnumeric(yval);
		yv=yval;
		yval='yval';
		T.yval=yv;
	end	

	% Deal with region if specified
	if ~isempty(regionOI) && ~islogical(regionOI)
		rIDX=T.center_x>=regionOI(1) & T.center_x<=regionOI(2) & T.center_y>=regionOI(3) & T.center_y<=regionOI(4);
		T=T(rIDX,:);
		if isempty(T)
			error('Provided region has eliminated all entries from the table, check that the coordinates are correct');
		end
	elseif ~isempty(regionOI) && islogical(regionOI) && regionOI

		f1=figure(1);
		clf
		hold on
		scatter(T.center_x,T.center_y,20,'k','filled');
		title('Draw a rectangle around the data you would like to select');
		hold off

		rgn=getrect;

		close(f1);

		regionOI=zeros(4,1);
		regionOI(1)=rgn(1);
		regionOI(2)=rgn(1)+rgn(3);
		regionOI(3)=rgn(2);
		regionOI(4)=rgn(2)+rgn(4);

		rIDX=T.center_x>=regionOI(1) & T.center_x<=regionOI(2) & T.center_y>=regionOI(3) & T.center_y<=regionOI(4);
		T=T(rIDX,:);
	end



	% Generate Plots
	switch plts
	case 'grd_ksn'

		if use_filtered
			g=T.mean_gradient_f;
			k=T.mean_ksn_f;
		else
			g=T.mean_gradient;
			k=T.mean_ksn;
		end

		if use_filtered
			if ismember('std_ksn_f',T.Properties.VariableNames) & strcmp(uncertainty,'std')
				sk=T.std_ksn_f;
				sg=T.std_gradient_f;
			elseif ismember('se_ksn_f',T.Properties.VariableNames) & strcmp(uncertainty,'se')
				sk=T.se_ksn_f;
				sg=T.se_gradient_f;
			elseif ismember('std_ksn_f',T.Properties.VariableNames) & ~ismember('se_ksn_f',T.Properties.VariableNames)
				sk=T.std_ksn_f;
				sg=T.std_gradient_f;
			elseif ismember('se_ksn_f',T.Properties.VariableNames) & ~ismember('std_ksn_f',T.Properties.VariableNames)
				sk=T.se_ksn_f;
				sg=T.se_gradient_f;
			end
		else
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
		end


		if fit_grd_ksn
			[pf]=KGF(k,g,D_in,K_in,s_in,n_val);
			[mg,mk,er,l]=KGG(pf.D,pf.K,pf.Sc,pf.n,k);
		end			


		f=figure(1);
		clf
		set(f,'Units','normalized','Position',[0.05 0.1 0.5 0.5],'renderer','painters');
		hold on 

		if ~strcmp(uncertainty,'none')
			errorbar(k,g,sg,sg,sk,sk,'.k','CapSize',0);
		end

		if ~isempty(color_by) & isnumeric(T.(color_by));
			colormap(cmap);
			scatter(k,g,30,T.(color_by),'filled');
			cb=colorbar;
			% Remove any underscores
			color_by_label=strrep(color_by,'_',' ');
			ylabel(cb,color_by_label);
		elseif ~isempty(color_by) & isa(T.(color_by),'cell');
			colormap(cmap);
			scatter(k,g,30,categorical(T.(color_by)),'filled');
			cb=colorbar('Ticks',[1:numel(unique(T.(color_by)))],'YTickLabel',unique(T.(color_by)));
			% Remove any underscores
			color_by_label=strrep(color_by,'_',' ');
			ylabel(cb,color_by_label);
		else
			scatter(k,g,30,'k','filled');
		end

		if fit_grd_ksn
			plot(mk,mg,'-r','LineWidth',2);
			disp(['Fluvial Erodibility (K) = ' num2str(pf.K)]);
			disp(['Diffusivity (D) = ' num2str(pf.D)]);
			disp(['Threshold Gradient (Sc) = ' num2str(pf.Sc)]);
		end

		xlabel('Mean Basin k_{sn}');
		ylabel('Mean Basin Gradient');
		hold off
	case 'grd_rlf'
		% Validate Relief Entry
		if use_filtered
			m_rlfN=['mean_rlf' num2str(rr) '_f'];
			se_rlfN=['se_rlf' num2str(rr) '_f'];
			std_rlfN=['std_rlf' num2str(rr) '_f'];
			if ismember(m_rlfN,T.Properties.VariableNames)
				r=T.(m_rlfN);
			else
				error('Relief radius is not recognized, confirm that you calculated local relief at this radius when running "ProcessRiverBasins"')
			end

			g=T.mean_gradient_f;
		else
			m_rlfN=['mean_rlf' num2str(rr)];
			se_rlfN=['se_rlf' num2str(rr)];
			std_rlfN=['std_rlf' num2str(rr)];
			if ismember(m_rlfN,T.Properties.VariableNames)
				r=T.(m_rlfN);
			else
				error('Relief radius is not recognized, confirm that you calculated local relief at this radius when running "ProcessRiverBasins"')
			end

			g=T.mean_gradient;
		end

		if use_filtered
			if ismember('std_gradient_f',T.Properties.VariableNames) & strcmp(uncertainty,'std')
				sr=T.(std_rlfN);
				sg=T.std_gradient_f;
			elseif ismember('se_gradient_f',T.Properties.VariableNames) & strcmp(uncertainty,'se')
				sr=T.(se_rlfN);
				sg=T.se_gradient_f;
			elseif ismember('std_gradient_f',T.Properties.VariableNames) & ~ismember('se_gradient_f',T.Properties.VariableNames)
				sr=T.(std_rlfN);
				sg=T.std_gradient_f;
			elseif ismember('se_gradient_f',T.Properties.VariableNames) & ~ismember('std_gradient_f',T.Properties.VariableNames)
				sr=T.(se_rlfN);
				sg=T.se_gradient_f;
			end
		else
			if ismember('std_gradient',T.Properties.VariableNames) & strcmp(uncertainty,'std')
				sr=T.(std_rlfN);
				sg=T.std_gradient;
			elseif ismember('se_gradient',T.Properties.VariableNames) & strcmp(uncertainty,'se')
				sr=T.(se_rlfN);
				sg=T.se_gradient;
			elseif ismember('std_gradient',T.Properties.VariableNames) & ~ismember('se_gradient',T.Properties.VariableNames)
				sr=T.(std_rlfN);
				sg=T.std_gradient;
			elseif ismember('se_gradient',T.Properties.VariableNames) & ~ismember('std_gradient',T.Properties.VariableNames)
				sr=T.(se_rlfN);
				sg=T.se_gradient;
			end
		end

		f=figure(1);
		set(f,'Units','normalized','Position',[0.05 0.1 0.5 0.5],'renderer','painters');
		clf
		hold on 

		if ~strcmp(uncertainty,'none')
			errorbar(r,g,sg,sg,sr,sr,'.k','CapSize',0);
		end

		if ~isempty(color_by) & isnumeric(T.(color_by));
			colormap(cmap);
			scatter(r,g,30,T.(color_by),'filled');
			cb=colorbar;
			% Remove any underscores
			color_by_label=strrep(color_by,'_',' ');
			ylabel(cb,color_by_label);
		elseif ~isempty(color_by) & isa(T.(color_by),'cell');
			colormap(cmap);
			scatter(r,g,30,categorical(T.(color_by)),'filled');
			cb=colorbar('Ticks',[1:numel(unique(T.(color_by)))],'YTickLabel',unique(T.(color_by)));
			% Remove any underscores
			color_by_label=strrep(color_by,'_',' ');
			ylabel(cb,color_by_label);
		else
			scatter(r,g,30,'k','filled');
		end

		xlabel(['Mean Basin ' num2str(rr) ' m^2 Relief']);
		ylabel('Mean Basin Gradient');
		hold off
		
	case 'rlf_ksn'
			% Validate Relief Entry
		if use_filtered
			m_rlfN=['mean_rlf' num2str(rr) '_f'];
			se_rlfN=['se_rlf' num2str(rr) '_f'];
			std_rlfN=['std_rlf' num2str(rr) '_f'];
			if ismember(m_rlfN,T.Properties.VariableNames)
				r=T.(m_rlfN);
			else
				error('Relief radius is not recognized, confirm that you calculated local relief at this radius when running "ProcessRiverBasins"')
			end

			k=T.mean_ksn_f;
		else
			m_rlfN=['mean_rlf' num2str(rr)];
			se_rlfN=['se_rlf' num2str(rr)];
			std_rlfN=['std_rlf' num2str(rr)];
			if ismember(m_rlfN,T.Properties.VariableNames)
				r=T.(m_rlfN);
			else
				error('Relief radius is not recognized, confirm that you calculated local relief at this radius when running "ProcessRiverBasins"')
			end

			k=T.mean_ksn;
		end

		if use_filtered
			if ismember('std_ksn_f',T.Properties.VariableNames) & strcmp(uncertainty,'std')
				sk=T.std_ksn_f;
				sr=T.(std_rlfN);
			elseif ismember('se_ksn_f',T.Properties.VariableNames) & strcmp(uncertainty,'se')
				sk=T.se_ksn_f;
				sr=T.(se_rlfN);
			elseif ismember('std_ksn_f',T.Properties.VariableNames) & ~ismember('se_ksn_f',T.Properties.VariableNames)
				sk=T.std_ksn_f;
				sr=T.(std_rlfN);
			elseif ismember('se_ksn_f',T.Properties.VariableNames) & ~ismember('std_ksn_f',T.Properties.VariableNames)
				sk=T.se_ksn_f;
				sr=T.(se_rlfN);
			end		
		else
			if ismember('std_ksn',T.Properties.VariableNames) & strcmp(uncertainty,'std')
				sk=T.std_ksn;
				sr=T.(std_rlfN);
			elseif ismember('se_ksn',T.Properties.VariableNames) & strcmp(uncertainty,'se')
				sk=T.se_ksn;
				sr=T.(se_rlfN);
			elseif ismember('std_ksn',T.Properties.VariableNames) & ~ismember('se_ksn',T.Properties.VariableNames)
				sk=T.std_ksn;
				sr=T.(std_rlfN);
			elseif ismember('se_ksn',T.Properties.VariableNames) & ~ismember('std_ksn',T.Properties.VariableNames)
				sk=T.se_ksn;
				sr=T.(se_rlfN);
			end
		end

		if fit_rlf_ksn
			ft=fittype('a*x');
			[fobj,gof]=fit(k,r,ft,'StartPoint',k\r);
			krs=coeffvalues(fobj);
			krs_unc=confint(fobj);
		end	

		f=figure(1);
		set(f,'Units','normalized','Position',[0.05 0.1 0.5 0.5],'renderer','painters');
		clf
		hold on 

		if ~strcmp(uncertainty,'none')
			errorbar(k,r,sr,sr,sk,sk,'.k','CapSize',0);
		end

		if ~isempty(color_by) & isnumeric(T.(color_by));
			colormap(cmap);
			scatter(k,r,30,T.(color_by),'filled');
			cb=colorbar;
			% Remove any underscores
			color_by_label=strrep(color_by,'_',' ');
			ylabel(cb,color_by_label);
		elseif ~isempty(color_by) & isa(T.(color_by),'cell');
			colormap(cmap);
			scatter(k,r,30,categorical(T.(color_by)),'filled');
			cb=colorbar('Ticks',[1:numel(unique(T.(color_by)))],'YTickLabel',unique(T.(color_by)));
			% Remove any underscores
			color_by_label=strrep(color_by,'_',' ');
			ylabel(cb,color_by_label);
		else
			scatter(k,r,30,'k','filled');
		end

		if fit_rlf_ksn
			ksn_vec=linspace(0,max(k)+50,100);
			rlf_vec=krs*ksn_vec;
			rlf_vec_pos=(krs_unc(2))*ksn_vec;
			rlf_vec_neg=(krs_unc(1))*ksn_vec;	
			plot(ksn_vec,rlf_vec,'-r','LineWidth',2);
			plot(ksn_vec,rlf_vec_pos,'--r');
			plot(ksn_vec,rlf_vec_neg,'--r');
			disp(['Relationship slope = ' num2str(krs)]);
			disp(['Uncertainty on slope = ' num2str(krs_unc(1)) ' : ' num2str(krs_unc(2))]);
			disp(['r-square = ' num2str(gof.rsquare)]);
		end		

		ylabel(['Mean Basin ' num2str(rr) ' m^2 Relief']);
		xlabel('Mean Basin k_{sn}');
		hold off

	case 'stacked_hypsometry'
		hypsCell=T.hypsometry;
		HI=T.hyp_integral;

		numHyps=numel(hypsCell);
		normEl=zeros(100,numHyps);
		normF=zeros(100,numHyps);

		El=normEl;
		F=normF;
		for ii=1:numHyps
			normEl(:,ii)=(hypsCell{ii}(:,2)-min(hypsCell{ii}(:,2)))/(max(hypsCell{ii}(:,2))-min(hypsCell{ii}(:,2)));
			normF(:,ii)=(hypsCell{ii}(:,1))/100;

			El(:,ii)=hypsCell{ii}(:,2);
			F(:,ii)=hypsCell{ii}(:,1);
		end

		histIm=zeros(100,100);
		jj=100:-1:1;
		for ii=1:100
			[N,~]=histcounts(normF(ii,:),linspace(0,1,101));
			histIm(jj(ii),:)=log10(N);
		end

		if ~isempty(color_by) & isnumeric(T.(color_by));
			col_val=T.(color_by);
			f(1)=figure(1);
			set(f(1),'Units','normalized','Position',[0.05 0.5 0.4 0.4],'renderer','painters');
			clf
			hold on
			colormap(cmap);
			cm=colormap;
			num_col=size(cm,1);
			[cix,ed]=discretize(col_val,num_col);
			for ii=1:num_col
				idx=cix==ii;
				plot(normF(:,idx),normEl(:,idx),'Color',cm(ii,:));
			end
			caxis([min(ed) max(ed)]);
			cb=colorbar;
			% Remove any underscores
			color_by_label=strrep(color_by,'_',' ');
			ylabel(cb,color_by_label);
			axis equal
			xlabel('Normalized Area');
			ylabel('Normalized Elevation');
			xlim([0 1]);
			ylim([0 1]);
			hold off			
		else
			f(1)=figure(1);
			set(f(1),'Units','normalized','Position',[0.05 0.5 0.4 0.4],'renderer','painters');
			clf
			hold on
			plot(normF,normEl,'-k');
			axis equal
			xlabel('Normalized Area');
			ylabel('Normalized Elevation');
			xlim([0 1]);
			ylim([0 1]);
			hold off
		end


		if ~isempty(color_by) & isnumeric(T.(color_by));
			col_val=T.(color_by);
			f(2)=figure(2);
			set(f(2),'Units','normalized','Position',[0.05 0.05 0.4 0.4],'renderer','painters');
			clf
			hold on
			colormap(cmap);
			cm=colormap;
			num_col=size(cm,1);
			[cix,ed]=discretize(col_val,num_col);
			for ii=1:num_col
				idx=cix==ii;
				plot(F(:,idx),El(:,idx),'Color',cm(ii,:));
			end
			caxis([min(ed) max(ed)]);			
			cb=colorbar;
			% Remove any underscores
			color_by_label=strrep(color_by,'_',' ');
			ylabel(cb,color_by_label);
			axis square
			xlabel('Percentage Area');
			ylabel('Elevation');
			hold off
		else
			f(2)=figure(2);
			set(f(2),'Units','normalized','Position',[0.05 0.05 0.4 0.4],'renderer','painters');
			clf
			hold on
			plot(F,El,'-k');
			axis square
			xlabel('Percentage Area');
			ylabel('Elevation');
			hold off
		end

		num_bins=20;
		[hix,hed]=discretize(HI,linspace(0,1,num_bins+1));

		f(3)=figure(3);
		set(f(3),'Units','normalized','Position',[0.5 0.5 0.4 0.4],'renderer','painters');
		clf 

		for ii=1:num_bins
			subplot(4,5,ii)
			hold on
			idx=hix==ii;
			perc(ii,1)=round((nnz(idx)/numel(idx))*100,1);
			nF=normF(:,idx);
			nE=normEl(:,idx);
			mF{ii,1}=mean(nF,2);
			mE{ii,1}=mean(nE,2);
			plot(nF,nE,'LineWidth',0.5,'Color',[0.4 0.4 0.4]);
			plot(mF{ii,1},mE{ii,1},'-r','LineWidth',2);
			if ii<=num_bins/2
				text(0.75,0.9,[num2str(perc(ii,1)) '%']);
			else
				text(0.1,0.1,[num2str(perc(ii,1)) '%']);
			end
			axis equal
			title(['HI ' num2str(hed(ii)) ' to ' num2str(hed(ii+1))])
			xlabel('Normalized Area');
			ylabel('Normalized Elevation');
			xlim([0 1]);
			ylim([0 1]);
			hold off
		end				

		f(4)=figure(4);
		set(f(4),'Units','normalized','Position',[0.5 0.05 0.4 0.4]);
		clf 
		X=linspace(0,1,100);
		Y=linspace(0,1,100);
		hold on
		colormap(jet);
		imagesc(X,Y,histIm);
		axis equal
		xlabel('Normalized Area');
		ylabel('Normalized Elevation');
		xlim([0 1]);
		ylim([0 1]);
		hold off

		f(5)=figure(5);
		set(f(5),'Units','normalized','Position',[0.25 0.25 0.4 0.4],'renderer','painters');
		clf 
		hold on
		colormap(jet);
		idx=perc>0;
		pf=perc(idx);
		jc=jet(100);
		for ii=1:num_bins
			ix=round((perc(ii,1)/max(pf))*100);
			if ix==0
				cl=[1 1 1];
			else
				cl=jc(ix,:);
			end
			plot(mF{ii,1},mE{ii,1},'Color',cl,'LineWidth',2);
		end
		caxis([min(pf) max(pf)]);
		cb=colorbar;
		ylabel(cb,'Percentage of Basins')
		axis equal
		xlabel('Normalized Area');
		ylabel('Normalized Elevation');
		xlim([0 1]);
		ylim([0 1]);
		hold off

	case 'scatterplot_matrix'
		% Find values with means
		if use_filtered
			VN=T.Properties.VariableNames;
			ix=regexp(VN,regexptranslate('wildcard','mean_*_f'));
			ix=cellfun(@any,ix);
			VNoi=VN(ix);
		else
			VN=T.Properties.VariableNames;
			ix=regexp(VN,regexptranslate('wildcard','mean_*'));
			ix=cellfun(@any,ix);
			VNoi=VN(ix);
		end

		% Add azimuth if it was calculated
		aix=regexp(VN,regexptranslate('wildcard','dist_along*'));
		aix=cellfun(@any,aix);
		if ~isempty(aix);
			VNoi=horzcat(VNoi,VN(aix));
		end

		num_sc=numel(VNoi);

		f=figure(1);
		set(f,'Units','normalized','Position',[0.05 0.5 0.9 0.9],'renderer','painters');		
		clf
		pos=1;

		for ii=1:num_sc
			yval=T.(VNoi{ii});
			yname=VNoi{ii};
			yname=strrep(yname,'_',' ');
			for jj=1:num_sc
				subplot(num_sc,num_sc,pos)
				hold on
				if ii==jj
					histogram(yval,25,'FaceColor','k');
					if jj==num_sc & ii==num_sc
						xlabel(yname);
					end
					axis square
				else
					xval=T.(VNoi{jj});
					xname=VNoi{jj};
					xname=strrep(xname,'_',' ');
					scatter(xval,yval,5,'k','filled');
					warning off
					f=fit(xval,yval,'poly2');
					warning on
					xx=linspace(min(xval),max(xval),50);
					yy=f(xx);
					plot(xx,yy,'-r');
					if ii==num_sc
						xlabel(xname);
					end

					if jj==1
						ylabel(yname);
					end

					axis square
				end
				hold off
				pos=pos+1;
			end
		end

	case 'compare_filtered'

		VN=T.Properties.VariableNames;
		ix=regexp(VN,regexptranslate('wildcard','mean_*_f'));
		ix=cellfun(@any,ix);
		VNoi=VN(ix);
		if isempty(VNoi)
			error('No filtered values were found in provided table');
		end

		if fit_filtered
			ft=fittype('a*x');
		end		

		num_filt=numel(VNoi);
		for ii=1:num_filt
			fN=VNoi{ii};
			N=strrep(fN,'_f','');

			% Parse value
			if strcmp(N,'mean_ksn')
				t='Mean Channel Steepness';
			elseif strcmp(N,'mean_gradient');
				t='Mean Gradient';
			elseif regexp(N,regexptranslate('wildcard','mean_rlf*'));
				t=['Mean ' strrep(N,'mean_rlf','') ' m^2 Relief'];
			else
				t=['Mean ' strrep(N,'mean_','')];
			end

			max_val=max([max(T.(fN)) max(T.(N))]);
			max_vec=[0 max_val];

			slp=round(T.(N)./T.(fN),1);

			idx1=slp==1;
			idx2=slp<1;
			idx3=slp>1;

			if fit_filtered
				idx=isnan(T.(fN)) | isnan(T.(N));
				x=double(T.(fN)(~idx));
				y=double(T.(N)(~idx));
				[fobj,gof]=fit(x,y,ft,'StartPoint',x\y);
				cf=coeffvalues(fobj);
				cf_unc=confint(fobj);
			end	

			f(ii)=figure(ii);
			set(gcf,'Units','normalized','Position',[0.05 0.1 0.5 0.5],'renderer','painters');
			clf
			hold on
			plot(max_vec,max_vec,'-k');
			sp(1)=scatter(T.(fN)(idx1),T.(N)(idx1),30,'k','filled');
			sp(2)=scatter(T.(fN)(idx2),T.(N)(idx2),30,'r','filled');
			sp(3)=scatter(T.(fN)(idx3),T.(N)(idx3),30,'b','filled');

			if fit_filtered
				fvec=linspace(0,max(T.(fN)),100);
				vec=cf*fvec;
				vec_pos=(cf_unc(2))*fvec;
				vec_neg=(cf_unc(1))*fvec;	
				pl=plot(fvec,vec,'-g','LineWidth',2);
				plot(fvec,vec_pos,'--g');
				plot(fvec,vec_neg,'--g');
				disp(['Fits on ' t ':']);
				disp(['	Relationship slope = ' num2str(cf)]);
				disp(['	Uncertainty on slope = ' num2str(cf_unc(1)) ' : ' num2str(cf_unc(2))]);
				disp(['	r-square = ' num2str(gof.rsquare)]);
				legend([sp pl],{'Filtered = Unfiltered','Filtered > Unfiltered','Filtered < Unfiltered','Best Fit'},'location','northwest')
			else
				legend(sp,{'Filtered = Unfiltered','Filtered > Unfiltered','Filtered < Unfiltered'},'location','northwest')			
			end



			xlabel('Filtered Means');
			ylabel('Unfiltered Means');
			title(t);
			hold off
		end

	case 'category_mean_hist'

		if isempty(cm1)
			error('For plot option "category_mean_hist" you must provide an input for "cat_mean1"')
		elseif strcmp(cm1,'ksn')
			VN=T.Properties.VariableNames;
			ix=regexp(VN,regexptranslate('wildcard','mksn*'));
			ix=cellfun(@any,ix);
			VNoi=VN(ix);
			Cat_Names=cellfun(@(x) strrep(x,'mksn_',''),VNoi,'UniformOutput',false);
			Main_Title='Mean k_{sn} within ';
		elseif strcmp(cm1,'gradient')
			VN=T.Properties.VariableNames;
			ix=regexp(VN,regexptranslate('wildcard','mgrad*'));
			ix=cellfun(@any,ix);
			VNoi=VN(ix);
			Cat_Names=cellfun(@(x) strrep(x,'mgrad_',''),VNoi,'UniformOutput',false);
			Main_Title='Mean gradient within ';
		elseif strcmp(cm1,'rlf')
			VN=T.Properties.VariableNames;
			srch_strng=['mr' num2str(rr) '*'];
			ix=regexp(VN,regexptranslate('wildcard',srch_strng));
			ix=cellfun(@any,ix);
			VNoi=VN(ix);
			if isempty(VNoi)
				error('Entry for "cat_mean1" is not recognized, check that relief radius correct');
			end
			srch_strng=['mr' num2str(rr) '_'];
			Cat_Names=cellfun(@(x) strrep(x,srch_strng,''),VNoi,'UniformOutput',false);	
			Main_Title=['Mean ' num2str(rr) '_Relief within '];	
		else
			VN=T.Properties.VariableNames;
			srch_strng=['m' cm1 '*'];
			ix=regexp(VN,regexptranslate('wildcard',srch_strng));
			ix=cellfun(@any,ix);
			VNoi=VN(ix);
			if isempty(VNoi)
				error('Entry for "cat_mean1" is not recognized');
			end	
			srch_strng=['m' cm1 '_'];			
			Cat_Names=cellfun(@(x) strrep(x,srch_strng,''),VNoi,'UniformOutput',false);	
			Main_Title=['Mean ' cm1 ' within '];
		end	

		for ii=1:numel(VNoi)
			vals=T.(VNoi{ii});
			vals(isnan(vals))=[];

			if ~isempty(vals)
				val_list{ii,1}=vals;
			end
		end
		val_list=vertcat(val_list{:});
		[~,edges]=discretize(val_list,100);

		for ii=1:numel(VNoi)
			vals=T.(VNoi{ii});
			vals(isnan(vals))=[];

			if ~isempty(vals)
				f(ii)=figure(ii);
				set(gcf,'Units','normalized','Position',[0.05 0.1 0.5 0.5],'renderer','painters');
				clf
				hold on
				histogram(vals,edges);
				title([Main_Title Cat_Names{ii}]);
				hold off
			end
		end

	case 'category_mean_compare'

		if isempty(cm1)
			error('For plot option "category_mean_compare" you must provide an input for "cat_mean1"')
		elseif strcmp(cm1,'ksn')
			VN=T.Properties.VariableNames;
			ix=regexp(VN,regexptranslate('wildcard','mksn*'));
			ix=cellfun(@any,ix);
			VNoi1=VN(ix);
			Cat_Names=cellfun(@(x) strrep(x,'mksn_',''),VNoi1,'UniformOutput',false);
			axis1='Mean k_{sn} within ';
		elseif strcmp(cm1,'gradient')
			VN=T.Properties.VariableNames;
			ix=regexp(VN,regexptranslate('wildcard','mgrad*'));
			ix=cellfun(@any,ix);
			VNoi1=VN(ix);
			Cat_Names=cellfun(@(x) strrep(x,'mgrad_',''),VNoi1,'UniformOutput',false);
			axis1='Mean gradient within ';
		elseif strcmp(cm1,'rlf')
			VN=T.Properties.VariableNames;
			srch_strng=['mr' num2str(rr) '*'];
			ix=regexp(VN,regexptranslate('wildcard',srch_strng));
			ix=cellfun(@any,ix);
			VNoi1=VN(ix);
			if isempty(VNoi1)
				error('Entry for "cat_mean1" is not recognized, check that relief radius correct');
			end
			srch_strng=['mr' num2str(rr) '_'];
			Cat_Names=cellfun(@(x) strrep(x,srch_strng,''),VNoi1,'UniformOutput',false);	
			axis1=['Mean ' num2str(rr) '_Relief within '];	
		else
			VN=T.Properties.VariableNames;
			srch_strng=['m' cm1 '*'];
			ix=regexp(VN,regexptranslate('wildcard',srch_strng));
			ix=cellfun(@any,ix);
			VNoi1=VN(ix);
			if isempty(VNoi1)
				error('Entry for "cat_mean1" is not recognized');
			end	
			srch_strng=['m' cm1 '_'];			
			Cat_Names=cellfun(@(x) strrep(x,srch_strng,''),VNoi1,'UniformOutput',false);	
			axis1=['Mean ' cm1 ' within '];
		end			

		if isempty(cm2)
			error('For plot option "category_mean_compare" you must provide an input for "cat_mean2"')
		elseif strcmp(cm2,'ksn')
			VN=T.Properties.VariableNames;
			ix=regexp(VN,regexptranslate('wildcard','mksn*'));
			ix=cellfun(@any,ix);
			VNoi2=VN(ix);
			Cat_Names=cellfun(@(x) strrep(x,'mksn_',''),VNoi2,'UniformOutput',false);
			axis2='Mean k_{sn} within ';
		elseif strcmp(cm2,'gradient')
			VN=T.Properties.VariableNames;
			ix=regexp(VN,regexptranslate('wildcard','mgrad*'));
			ix=cellfun(@any,ix);
			VNoi2=VN(ix);
			Cat_Names=cellfun(@(x) strrep(x,'mgrad_',''),VNoi2,'UniformOutput',false);
			axis2='Mean gradient within ';
		elseif strcmp(cm2,'rlf')
			VN=T.Properties.VariableNames;
			srch_strng=['mr' num2str(rr) '*'];
			ix=regexp(VN,regexptranslate('wildcard',srch_strng));
			ix=cellfun(@any,ix);
			VNoi2=VN(ix);
			if isempty(VNoi2)
				error('Entry for "cat_mean2" is not recognized, check that relief radius correct');
			end
			srch_strng=['mr' num2str(rr) '_'];
			Cat_Names=cellfun(@(x) strrep(x,srch_strng,''),VNoi2,'UniformOutput',false);	
			axis2=['Mean ' num2str(rr) '_Relief within '];	
		else
			VN=T.Properties.VariableNames;
			srch_strng=['m' cm2 '*'];
			ix=regexp(VN,regexptranslate('wildcard',srch_strng));
			ix=cellfun(@any,ix);
			VNoi2=VN(ix);
			if isempty(VNoi2)
				error('Entry for "cat_mean2" is not recognized');
			end	
			srch_strng=['m' cm2 '_'];			
			Cat_Names=cellfun(@(x) strrep(x,srch_strng,''),VNoi2,'UniformOutput',false);	
			axis2=['Mean ' cm2 ' within '];
		end	

		rng_v1=zeros(numel(VNoi1),2);
		rng_v2=zeros(numel(VNoi1),2);		

		for ii=1:numel(VNoi1)
			vals1=T.(VNoi1{ii});
			vals2=T.(VNoi2{ii});

			idx=~isnan(vals1) & ~isnan(vals2);
			vals1=vals1(idx);
			vals2=vals2(idx);

			if ~isempty(vals1) & ~isempty(vals2)
				rng_v1(ii,1)=min(vals1,[],'omitnan');
				rng_v2(ii,1)=min(vals2,[],'omitnan');
				rng_v1(ii,2)=max(vals1,[],'omitnan');
				rng_v2(ii,2)=max(vals2,[],'omitnan');				
			end
		end

		rng_v1=[min(rng_v1(:,1),[],'omitnan') max(rng_v1(:,2),[],'omitnan')];
		rng_v2=[min(rng_v2(:,1),[],'omitnan') max(rng_v2(:,2),[],'omitnan')];

		for ii=1:numel(VNoi1)
			vals1=T.(VNoi1{ii});
			vals2=T.(VNoi2{ii});

			idx=~isnan(vals1) & ~isnan(vals2);
			vals1=vals1(idx);
			vals2=vals2(idx);


			if ~isempty(vals1) & ~isempty(vals2)
				f(ii)=figure(ii);
				set(gcf,'Units','normalized','Position',[0.05 0.1 0.5 0.5],'renderer','painters');
				clf
				hold on
				scatter(vals1,vals2,30,'k','filled');
				xlim(rng_v1);
				ylim(rng_v2);
				xlabel([axis1 Cat_Names{ii}])
				ylabel([axis2 Cat_Names{ii}])
				hold off
			end
		end

	case 'compare_mean_and_dist'

		% Validate Input
		if strcmp(stOI,'ksn')
			val=T.mean_ksn;
			m_valN='mean_ksn';
			load_flag=1;
			title_str='Channel Steepness';
		elseif strcmp(stOI,'gradient')
			val=T.mean_gradient;
			m_valN='mean_gradient';
			load_flag=2;
			title_str='Gradient';
		elseif strcmp(stOI,'elevation')
			val=T.mean_el;
			m_valN='mean_el';
			load_flag=3;
			title_str='Elevation';
		elseif strcmp(stOI,'relief')
			m_valN=['mean_rlf' num2str(rr)];
			if ismember(m_valN,T.Properties.VariableNames)
				val=T.(m_valN);	
				load_flag=4;
			else
				error('Relief radius is not recognized, confirm that you calculated local relief at this radius when running "ProcessRiverBasins"')
			end
			title_str=[num2str(rr) ' m^2 Relief'];
		else
			m_valN=['mean_' stOI];
			if ismember(m_valN,T.Properties.VariableNames)
				val=T.(m_valN);
				load_flag=5;
			else
				error(['There is not a column named "mean_' stOI '" in the supplied table, confirm name of additional grid provided to "ProcessRiverBasins"'])
			end
			title_str=[upper(stOI(1)) stOI(2:end)];
		end			

		% Determine type of plot
		if isempty(basin_num)
			out=cell(numel(val),1);

			w1=waitbar(0,'Compiling Statistics');
			for ii=1:numel(val)
				if load_flag==1
					load(T.file_path{ii,1},'MSNc');
					out{ii,1}=[MSNc.ksn]';
				elseif load_flag==2
					load(T.file_path{ii,1},'Goc');
					g=Goc.Z(:);
					g(isnan(g))=[];
					out{ii,1}=g;
				elseif load_flag==3
					load(T.file_path{ii,1},'DEMoc');
					d=DEMoc.Z(:);
					d(isnan(d))=[];
					out{ii,1}=d;
				elseif load_flag==4
					load(T.file_path{ii,1},'rlf');
					ix=find(cell2mat(rlf(:,2))==rr);
					r=rlf{ix,1}.Z(:);		
					r(isnan(r))=[];
					out{ii,1}=r;	
				elseif load_flag==5
					load(T.file_path{ii,1},'AGc');
					ix=find(strcmp(AGc(:,2),stOI));
					a=AGc{ix,1}.Z(:);
					a(isnan(a))=[];
					out{ii,1}=a;
				end				
				waitbar(ii/numel(val));
			end
			close(w1);

			out=vertcat(out{:});

			f=figure(1);
			set(f,'Units','normalized','Position',[0.05 0.1 0.5 0.5],'renderer','painters');
			clf

			% Filter
			[N,ed]=histcounts(out,100);
			N=N/max(N);
			bin_list=[1:100]';
			[IX]=discretize(out,ed);
			idx=N>0.01;
			keep_bins=bin_list(idx);
			idx=ismember(IX,keep_bins);
			out_f=out(idx);

			subplot(1,2,1);
			hold on 
			[nF,edgesF]=histcounts(out_f,100);
			[nV,~]=histcounts(val,edgesF);
			nF=nF/max(nF);
			nV=nV/max(nV);
			histogram('BinEdges',edgesF,'BinCounts',nF);
			histogram('BinEdges',edgesF,'BinCounts',nV);
			title('<1% Bins Removed');
			legend('All Basins','Basin Means','location','best');
			xlabel(title_str);
			ylabel('Normalized Counts');
			hold off

			subplot(1,2,2)
			hold on 
			[n,edges]=histcounts(out,100);
			[nV,~]=histcounts(val,edges);
			n=n/max(n);
			nV=nV/max(nV);
			histogram('BinEdges',edges,'BinCounts',n);
			histogram('BinEdges',edges,'BinCounts',nV);
			title('All Data')
			legend('All Basins','Basin Means','location','best');			
			xlabel(title_str);			
			hold off

		else
			
			ii=find(T.river_mouth==basin_num);
			if isempty(ii)
				error('Basin number does not appear in provided table')
			end

			if load_flag==1
				load(T.file_path{ii,1},'MSNc');
				out=[MSNc.ksn]';
			elseif load_flag==2
				load(T.file_path{ii,1},'Goc');
				g=Goc.Z(:);
				g(isnan(g))=[];
				out=g;
			elseif load_flag==3
				load(T.file_path{ii,1},'DEMoc');
				d=DEMoc.Z(:);
				d(isnan(d))=[];
				out=d;
			elseif load_flag==4
				load(T.file_path{ii,1},'rlf');
				ix=find(cell2mat(rlf(:,2))==rr);
				r=rlf{ix,1}.Z(:);		
				r(isnan(r))=[];
				out=r;	
			elseif load_flag==5
				load(T.file_path{ii,1},'AGc');
				ix=find(strcmp(AGc(:,2),stOI));
				a=AGc{ix,1}.Z(:);
				a(isnan(a))=[];
				out=a;
			end	

			[N,~]=histcounts(out,100);

			f=figure(1);
			set(f,'Units','normalized','Position',[0.05 0.1 0.5 0.5],'renderer','painters');
			clf

			hold on 
			histogram(out,100);
			plot([T.(m_valN)(ii,1),T.(m_valN)(ii,1)],[0,max(N)],'-k','LineWidth',2);
			title(title_str);
			xlabel(['Values from Basin ' num2str(basin_num)]);
			hold off
		end		

	case 'xy'

		if isempty(xval) | isempty(yval)
			error('For "xy" plot, entries for "xval" and "yval" are required');
		end

		VN=T.Properties.VariableNames;

		if ~any(strcmp(VN,xval))
			error('Entry to provided to "xval" is not regonized as the name of a column in the provided table');
		elseif ~any(strcmp(VN,yval))
			error('Entry to provided to "yval" is not regonized as the name of a column in the provided table');
		end

		x=T.(xval);
		y=T.(yval);

		if ~isnumeric(x) | ~isnumeric(y)
			error('Values in both x and y must be numeric')
		end

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

		f=figure(1);
		set(f,'Units','normalized','Position',[0.05 0.1 0.5 0.5],'renderer','painters');
		clf
		hold on 

		if ~strcmp(uncertainty,'none') & ~isempty(sx) & isempty(sy)
			errorbar(x,y,sx,'horizontal','.k','CapSize',0);
		elseif ~strcmp(uncertainty,'none') & isempty(sy) & ~isempty(sx)
			errorbar(x,y,sy,'vertical','.k','CapSize',0);
		elseif ~strcmp(uncertainty,'none') & ~isempty(sx) & ~isempty(sy)
			errorbar(x,y,sy,sy,sx,sx,'.k','CapSize',0);
		end

		if ~isempty(color_by) & isnumeric(T.(color_by));
			colormap(cmap);
			scatter(x,y,30,T.(color_by),'filled');
			cb=colorbar;
			% Remove any underscores
			color_by_label=strrep(color_by,'_',' ');
			ylabel(cb,color_by_label);
		elseif ~isempty(color_by) & isa(T.(color_by),'cell');
			colormap(cmap);
			scatter(x,y,30,categorical(T.(color_by)),'filled');
			cb=colorbar('Ticks',[1:numel(unique(T.(color_by)))],'YTickLabel',unique(T.(color_by)));
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


	if save_figure
		num_figs=numel(f);

		if num_figs>1
			for ii=1:num_figs
				orient(f(ii),'landscape');
				print(f(ii),'-dpdf','-fillpage',fullfile(out_dir,['Figure_' num2str(ii) '.pdf']));
			end
		else % WEIRD SIMULINK ERROR FROM HELL
			current=gcf;
			orient 'landscape';
			print(current,'-dpdf','-fillpage',fullfile(out_dir,['Figure_1.pdf']));
		end

	end
end

function [param_fit]=KGF(x0,y0,D_in,K_in,s_in,n)
	% n is not a free parameter in this version

    model=@RC;
    stp=[D_in,K_in,s_in];
    lb=[0,0,0];
    ub=[1,1,2];
	opt=optimset('Display','iter'); %,'MaxFunEvals',1000,'MaxIter',1000);

	disp('Fitting...')
    est=fmincon(model,stp,[],[],[],[],lb,ub);

    param_fit.D=est(1);
    param_fit.K=est(2);
    param_fit.Sc=est(3);
    param_fit.n=n;

    % Define least absolute deviation solver
    function [lad]=RC(params)
    	D=params(1);
    	K=params(2);
    	s=params(3);

    	er=ErKsn(x0,K,n);
    	er=er.';

		% Calculate hillslope lengths
		for ii=1:numel(er);
			l(1,ii)=CharLh(K,D,er(ii),n,s);
		end

		% Calculate mean slopes
		[mean_gradient]=NonLinearHillslope(er,s,D,l);
		res=(y0.')-mean_gradient;
		lad=sum(abs(res));
	end
end

function [me_grd,me_ksn,er,l]=KGG(D,K,s,n,ksn)

	% Define Erosion Rate Range
	% % er - erosion rate (m/yr)
	% er=logspace(-5,-1,100);
	ksn_vec=linspace(min(ksn),max(ksn)+50,100);
	er=ErKsn(ksn_vec,K,n);

	% Calculate hillslope lengths
	% l - hillslope length (m)
	for ii=1:numel(er);
		l(1,ii)=CharLh(K,D,er(ii),n,s);
	end

	% Calculate mean gradients
	[me_grd]=NonLinearHillslope(er,s,D,l);

	% Calcualte mean ksn
	% me_ksn - basin averaged channel steepness
	[me_ksn]=KsnEr(er,K,n);
end

function [lh]=CharLh(K,D,er,n,s)
	% er - erosion rate
	% s - threshold gradient
	% b - sediment density ratio
	% D - diffusivity
	% K - fluvial erodiblity
	b=2;
	m=n/2;

	fun=@(lh) (((D*(s^2))/(b*er*lh))*(sqrt(1+((b*er*lh)/(D*s))^2)-1))-((er/(2*K*(lh^(2*m))))^(1/n));

	lh=fzero(fun,100);
end

function [ksn]=KsnEr(er,K,n)
	% K = erosivity
	% n = slope exponent
		ksn=(er./K).^(1/n);
end

function [er]=ErKsn(ksn,K,n)
	% K = erosivity
	% n = slope exponent
	er=K.*(ksn.^n);
end

function [mean_gradient_r]=NonLinearHillslope(er,s,D,l)
	% er - erosion rate
	% s - threshold gradient
	% b - sediment density ratio
	% D - diffusivity
	% l - hillslope length

	b=2;

	[z_0]=el_func_x(er,s,b,D,0);
	[z_l]=el_func_x(er,s,b,D,l);

	mean_gradient_r=(z_0-z_l)./l;
end

function [z]=el_func_x(er,s,b,D,xx)
	% Helper function for NonLinearHillslope, elevation as a function of x

	out_term=(-(s^2)./(2*b.*er));

	rep_term=(2*b.*er)./s;

	first_inner=sqrt((D^2)+(rep_term.*xx).^2);

	second_inner=D * log((first_inner+D)./rep_term);

	z=out_term.*(first_inner-second_inner);

end
