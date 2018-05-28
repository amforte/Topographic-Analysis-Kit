function BasinStatsPlots(basin_table,plots,varargin)
	% Function to take the complied outputs from 'ProcessRiverBasins' and 'SubDivideBigBasins' and produce various plots
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
	%		'xy' - generic plot, requires entries to optional 'xval' and 'yval' inputs
	%
	% Optional Inputs:
	%	uncertianty ['se'] - uncertainty to value use for plots, valid options are 'se' (standard error), 'std' (standard deviation), or 'none'. 
	%		Providing 'none' indicates you do not want to plot errorbars. Behavior of this option will depend on how you ran ProcessRiverBasins, 
	%		e.g. if you only calculated standard deviations when running ProcessRiverBasins but supply 'se'	here, the code will ignore your choice
	%		and use the standard deviation values.
	%	color_by [] - value to color points by, valid for 'grd_ksn','grd_rlf','rlf_ksn', and 'xy'
	%	cmap [jet] - colormap to use if an entry is provided to 'color_by', can be the name of a standard colormap or a nx3 array of rgb values
	%		to use as a colormap. 
	%	xval [] - value to plot on x axis (name of column as it appears in the provided table) for plot type 'xy'
	%	yval [] - value to plot on y axis (name of column as it appears in the provided table) for plot type 'xy'
	%	rlf_radius [2500] - radius of relief used when plotting relief related values
	%	cat_mean1 [] - category to use for plotting, see 'category_mean_hist' or 'category_mean_compare', valid inputs are 'ksn', 'rlf', 'gradient', or 
	%		the name of an additional grid provided to ProcessRiverMeans.
	%	cat_mean2 [] - category to use for plotting, 'category_mean_compare' , valid inputs are 'ksn', 'rlf', 'gradient', or the name of an additional grid 
	%		provided to ProcessRiverMeans.
	%	only_positive [false] - filter out negative values when using either 'category_mean_hist' or 'category_mean_compare'	
	%	save_figure [false] - logical flag to save pdfs of all figures produced
	%
	% Examples:
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
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Last Revised Spring 2018 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'BasinStatsPlots';
	addRequired(p,'basin_table',@(x) isa(x,'table'));
	addRequired(p,'plots',@(x) ischar(validatestring(x,{'grd_ksn','grd_rlf','rlf_ksn','compare_filtered','category_mean_hist','category_mean_compare','xy'})));

	addParamValue(p,'uncertainty','se',@(x) ischar(validatestring(x,{'se','std','none'})));
	addParamValue(p,'color_by',[],@(x) ischar(x));
	addParamValue(p,'cmap','jet',@(x) ischar(x) || isnumeric(x) & size(x,2)==3);
	addParamValue(p,'xval',[],@(x) ischar(x));
	addParamValue(p,'yval',[],@(x) ischar(x));
	addParamValue(p,'rlf_radius',2500,@(x) isnumeric(x) && isscalar(x));
	addParamValue(p,'cat_mean1',[],@(x) ischar(x));
	addParamValue(p,'cat_mean2',[],@(x) ischar(x));
	addParamValue(p,'only_positive',false,@(x) isscalar(x) && islogical(x));
	addParamValue(p,'save_figure',false,@(x) isscalar(x) && islogical(x));

	parse(p,basin_table,plots,varargin{:});
	T=p.Results.basin_table;
	plts=p.Results.plots;

	uncertainty=p.Results.uncertainty;
	color_by=p.Results.color_by;
	cmap=p.Results.cmap;
	xval=p.Results.xval;
	yval=p.Results.yval;
	rr=p.Results.rlf_radius;
	cm1=p.Results.cat_mean1;
	cm2=p.Results.cat_mean2;
	op=p.Results.only_positive;
	save_figure=p.Results.save_figure;


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

		f=figure(1);
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

		xlabel('Mean Basin k_{sn}');
		ylabel('Mean Basin Gradient');
		hold off
	case 'grd_rlf'
		% Validate Relief Entry
		m_rlfN=['mean_rlf' num2str(rr)];
		se_rlfN=['se_rlf' num2str(rr)];
		std_rlfN=['std_rlf' num2str(rr)];
		if ismember(m_rlfN,T.Properties.VariableNames)
			r=T.(m_rlfN);
		else
			error('Relief radius is not recognized, confirm that you calculated local relief at this radius when running "ProcessRiverBasins"')
		end

		g=T.mean_gradient;

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
		m_rlfN=['mean_rlf' num2str(rr)];
		se_rlfN=['se_rlf' num2str(rr)];
		std_rlfN=['std_rlf' num2str(rr)];
		if ismember(m_rlfN,T.Properties.VariableNames)
			r=T.(m_rlfN);
		else
			error('Relief radius is not recognized, confirm that you calculated local relief at this radius when running "ProcessRiverBasins"')
		end

		k=T.mean_ksn;

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
			scatter(k,g,30,'k','filled');
		end

		ylabel(['Mean Basin ' num2str(rr) ' m^2 Relief']);
		xlabel('Mean Basin k_{sn}');
		hold off

	case 'compare_filtered'
		VN=T.Properties.VariableNames;
		ix=regexp(VN,regexptranslate('wildcard','mean_*_f'));
		ix=cellfun(@any,ix);
		VNoi=VN(ix);
		if isempty(VNoi)
			error('No filtered values were found in provided table');
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

			f(ii)=figure(ii);
			set(gcf,'Units','normalized','Position',[0.05 0.1 0.5 0.5],'renderer','painters');
			clf
			hold on
			plot(max_vec,max_vec,'-k');
			sp(1)=scatter(T.(fN)(idx1),T.(N)(idx1),30,'k','filled');
			sp(2)=scatter(T.(fN)(idx2),T.(N)(idx2),30,'r','filled');
			sp(3)=scatter(T.(fN)(idx3),T.(N)(idx3),30,'b','filled');
			legend(sp,{'Filtered = Unfiltered','Filtered > Unfiltered','Filtered < Unfiltered'},'location','northwest')
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
			if op
				vals(vals<0)=[];
			end

			if ~isempty(vals)
				val_list{ii,1}=vals;
			end
		end
		val_list=vertcat(val_list{:});
		[~,edges]=discretize(val_list,100);

		for ii=1:numel(VNoi)
			vals=T.(VNoi{ii});
			vals(isnan(vals))=[];
			if op
				vals(vals<0)=[];
			end

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

			if op
				idx= vals1>=0 & vals2>=0;
				vals1=vals1(idx);
				vals2=vals2(idx);
			end

			if ~isempty(vals1) & ~isempty(vals2)
				rng_v1(ii,1)=nanmin(vals1);
				rng_v2(ii,1)=nanmin(vals2);
				rng_v1(ii,2)=nanmax(vals1);
				rng_v2(ii,2)=nanmax(vals2);				
			end
		end

		rng_v1=[nanmin(rng_v1(:,1)) nanmax(rng_v1(:,2))];
		rng_v2=[nanmin(rng_v2(:,1)) nanmax(rng_v2(:,2))];

		for ii=1:numel(VNoi1)
			vals1=T.(VNoi1{ii});
			vals2=T.(VNoi2{ii});

			idx=~isnan(vals1) & ~isnan(vals2);
			vals1=vals1(idx);
			vals2=vals2(idx);

			if op
				idx= vals1>=0 & vals2>=0;
				vals1=vals1(idx);
				vals2=vals2(idx);
			end

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

	case 'xy'
		if isempty(xval) | isempty(yval)
			error('For "xy" plot, entries for "xval" and "yval" are required');
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
			errobar(x,y,sx,'horizontal','.k','CapSize',0);
		elseif ~strcmp(uncertainty,'none') & isempty(sy) & ~isempty(sx)
			errobar(x,y,sy,'vertical','.k','CapSize',0);
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
			scatter(k,g,30,'k','filled');
		end

		xlabel(['Mean ' sxN]);
		ylabel(['Mean ' syN]);
		hold off
	end


	if save_figure
		num_figs=numel(f);

		for ii=1:num_figs
			orient(f(ii),'landscape');
			print(f(ii),'-dpdf','-fillpage',['Figure_' num2str(ii) '.pdf']);
		end
	end
