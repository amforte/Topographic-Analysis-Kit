function cmpSegmentPlotter(wdir,basin_nums,varargin)
	% Description:
	% 	Function to plot all of the chi-Z relationships, longitudinal profiles, and slope area plots from a series of picked segments of river networks 
	%	that result from the 'SegmentPicker' function. 
	%
	% Required Input:
	%	wdir - full path of working directory
	% 	basin_nums - row or column vector of basin numbers used for the SegmentPicker you wish to plot together. Code expects that
	%		the mat files saved from cmpSegmentPicker are in the present working directory.
	% 
	% Optionl Input:
	%	separate [false] - logical flag to plot all segments as separate figures 
	%	subset [] - list of specific river numbers (i.e. the third column of either the 'Heads' or the 'Outlets' variable) that you wish to include
	%				in the plot. Only valid if you have only provided a single basin number for 'basin_nums'.
	%	label [false]  - logical flag to either label individual streams with the river number (true) or not label them (false, default). If 'separate' flag
	%					is true then the input for label is ignored as the stream number will be in the title of the plots
	%	names [] - option to add an identifying name for streams when 'label' is set to true.
	%
	% Ouptupt:
	%	saves pdfs of all figures produced
	%
   	% Examples if running for the command line, minus OS specific way of calling main TAK function:
    %   SegmentPlotter /path/to/wdir 1
    %   SegmentPlotter /path/to/wdir [1 2 3 6] separate true
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
	p.FunctionName = 'cmpSegmentPlotter';
	addRequired(p,'wdir',@(x) ischar(x));
	addRequired(p,'basin_nums',@(x) isnumeric(x));

	addParameter(p,'separate',false,@(x) islogical(x));
	addParameter(p,'subset',[],@(x) isnumeric(x));
	addParameter(p,'label',false,@(x) islogical(x));
	addParameter(p,'trunks_only',false,@(x) islogical(x));
	addParameter(p,'names','',@(x) ischar(x) | iscell(x));

	parse(p,wdir,basin_nums,varargin{:});
	wdir=p.Results.wdir;
	basin_nums=p.Results.basin_nums;
	sep=p.Results.separate;
	subset=p.Results.subset;
	lab=p.Results.label;
	nm=p.Results.names;

	num_basins=numel(basin_nums);

	% Check for errors
	if ~isempty(subset) & num_basins>1
		warning('Providing a list of specific rivers via "subset" is only recognized if only one basin number is provided, ignoring "subset" input')
		su=0;
	elseif ~isempty(subset) & num_basins==1
		su=1;
	elseif isempty(subset)
		su=0;
	end


	switch su
	case 0
		if ~sep
			cMap=jet(num_basins);

			f1=figure(1);
			set(f1,'Units','normalized','Position',[0.05 0.1 0.8 0.8],'renderer','painters');

			for ii=1:num_basins
				basin_num=basin_nums(ii);
				fileName=fullfile(wdir,['PickedSegments_' num2str(basin_num) '.mat']);
				load(fileName,'ChiSgmnts','SlpAreaSgmnts');

				vInfo=who('-file',fileName);
				if ismember('Heads',vInfo)
					load(fileName,'Heads');
					seg_num=Heads(:,3);
				elseif ismember('Outlets',vInfo)
					load(fileName,'Outlets');
					seg_num=Outlets(:,3);
				end

				num_seg=numel(ChiSgmnts);
					
				for jj=1:num_seg
					C=ChiSgmnts{jj};
					bSA=SlpAreaSgmnts{jj,1};
					aSA=SlpAreaSgmnts{jj,2};

					subplot(3,1,1);
					hold on
					p1(ii)=plot(C.chi,C.elev,'Color',cMap(ii,:));
					if lab
						[mc,ix]=max(C.chi);
						me=C.elev(ix);
						text(mc+mc*0.01,me,[nm ' ' num2str(seg_num(jj))],'Color',cMap(ii,:));
					end
					hold off

					subplot(3,1,2);
					hold on
					p2(ii)=plot(C.distance,C.elev,'Color',cMap(ii,:));
					if lab
						[md,ix]=max(C.distance);
						me=C.elev(ix);
						text(md+md*0.01,me,[nm ' ' num2str(seg_num(jj))],'Color',cMap(ii,:));
					end
					hold off

					ax3=subplot(3,1,3);
					hold on
					scatter(aSA(:,2),aSA(:,1),5,cMap(ii,:),'+');
					p3(ii)=scatter(bSA(:,2),bSA(:,1),20,'MarkerFaceColor',cMap(ii,:),'MarkerEdgeColor','k');
					set(ax3,'Xscale','log','Yscale','log','XDir','reverse');
					hold off
	
				end
			LegendEnt{ii}=['Basin ' num2str(basin_num)];
			end

			subplot(3,1,1);
			hold on
				xlabel('\chi')
				ylabel('Elevation (m)')
				title('\chi - Z')
				legend(p1,LegendEnt,'location','best');
			hold off

			subplot(3,1,2);
			hold on
				xlabel('Distance from Mouth (m)')
				ylabel('Elevation (m)')
				title('Long Profile')
				legend(p2,LegendEnt,'location','best');
			hold off

			subplot(3,1,3);
			hold on
				xlabel('Log Drainage Area');
				ylabel('Log Gradient');
				legend(p3,LegendEnt,'location','best');
				title('Slope-Area');
			hold off

			figFile=fullfile(wdir,['StreamSegments_' num2str(basin_num) '.pdf']);
			print(f1,'-dpdf',figFile,'-fillpage');

		elseif sep

			fig_num=1;

			for ii=1:num_basins
				basin_num=basin_nums(ii);
				fileName=fullfile(wdir,['PickedSegments_' num2str(basin_num) '.mat']);
				load(fileName,'ChiSgmnts','SlpAreaSgmnts');
				vInfo=who('-file',fileName);
				if ismember('Heads',vInfo)
					load(fileName,'Heads');
					seg_num=Heads(:,3);
				elseif ismember('Outlets',vInfo)
					load(fileName,'Outlets');
					seg_num=Outlets(:,3);
				end

				num_seg=numel(ChiSgmnts);

				for jj=1:num_seg
					f=figure(fig_num);
					set(f,'Units','normalized','Position',[0.05 0.1 0.8 0.8],'renderer','painters');
					C=ChiSgmnts{jj};
					bSA=SlpAreaSgmnts{jj,1};
					aSA=SlpAreaSgmnts{jj,2};

					subplot(3,1,1);
					hold on
					plot(C.chi,C.elev,'-k');
					xlabel('Chi')
					ylabel('Elevation (m)')
					title(['Basin ' num2str(basin_num) '- Stream ' num2str(seg_num(jj))]);
					hold off

					subplot(3,1,2);
					hold on
					plot(C.distance,C.elev,'-k');
					xlabel('Distance from Mouth (m)')
					ylabel('Elevation (m)')
					hold off

					ax3=subplot(3,1,3);
					hold on
					scatter(aSA(:,2),aSA(:,1),5,[0.5 0.5 0.5],'+');
					scatter(bSA(:,2),bSA(:,1),20,'k','filled');
					set(ax3,'Xscale','log','Yscale','log','XDir','reverse');
					xlabel('Log Drainage Area');
					ylabel('Log Gradient');
					hold off					

					drawnow

					figFile=fullfile(wdir,['StreamSegments_' num2str(basin_num) '_' num2str(fig_num) '_.pdf']);
					print(f,'-dpdf',figFile,'-fillpage');

					fig_num=fig_num+1;
				end
			end
		end

	case 1

		if ~sep

			f1=figure(1);
			set(f1,'Units','normalized','Position',[0.05 0.1 0.8 0.8],'renderer','painters');


			basin_num=basin_nums;
			fileName=fullfile(wdir,['PickedSegments_' num2str(basin_num) '.mat']);
			load(fileName,'ChiSgmnts','SlpAreaSgmnts');

			vInfo=who('-file',fileName);
			if ismember('Heads',vInfo)
				load(fileName,'Heads');
				seg_num=Heads(:,3);
			elseif ismember('Outlets',vInfo)
				load(fileName,'Outlets');
				seg_num=Outlets(:,3);
			end

			idx=ismember(seg_num,subset);
			ChiSgmnts=ChiSgmnts(idx);
			SlpAreaSgmnts=SlpAreaSgmnts(idx,:);
			seg_num=seg_num(idx);

			num_seg=numel(ChiSgmnts);

			for jj=1:num_seg
				C=ChiSgmnts{jj};
				bSA=SlpAreaSgmnts{jj,1};
				aSA=SlpAreaSgmnts{jj,2};

				subplot(3,1,1);
				hold on
				plot(C.chi,C.elev,'-k');
				if lab
					[mc,ix]=max(C.chi);
					me=C.elev(ix);
					text(mc+mc*0.01,me,[nm ' ' num2str(seg_num(jj))]);
				end
				hold off

				subplot(3,1,2);
				hold on
				plot(C.distance,C.elev,'-k');
				if lab
					[md,ix]=max(C.distance);
					me=C.elev(ix);
					text(md+md*0.01,me,[nm ' ' num2str(seg_num(jj))]);
				end
				hold off

				ax3=subplot(3,1,3);
				hold on
				scatter(aSA(:,2),aSA(:,1),5,[0.5 0.5 0.5],'+');
				scatter(bSA(:,2),bSA(:,1),20,'k','filled');
				set(ax3,'Xscale','log','Yscale','log','XDir','reverse');
				hold off

			end

			subplot(3,1,1);
			hold on
				xlabel('Chi')
				ylabel('Elevation (m)')
				title('Chi - Z')
			hold off

			subplot(3,1,2);
			hold on
				xlabel('Distance from Mouth (m)')
				ylabel('Elevation (m)')
				title('Long Profile')
			hold off

			subplot(3,1,3);
			hold on
				xlabel('Log Drainage Area');
				ylabel('Log Gradient');
				title('Slope-Area');
			hold off	

			figFile=fullfile(wdir,['StreamSegments_' num2str(basin_num) '.pdf']);
			print(f1,'-dpdf',figFile,'-fillpage');		

		elseif sep

			fig_num=1;

			for ii=1:num_basins
				basin_num=basin_nums(ii);
				fileName=fullfile(wdir,['PickedSegments_' num2str(basin_num) '.mat']);
				load(fileName,'ChiSgmnts','SlpAreaSgmnts');
				vInfo=who('-file',fileName);
				if ismember('Heads',vInfo)
					load(fileName,'Heads');
					seg_num=Heads(:,3);
				elseif ismember('Outlets',vInfo)
					load(fileName,'Outlets');
					seg_num=Outlets(:,3);
				end

				idx=ismember(seg_num,subset);
				ChiSgmnts=ChiSgmnts(idx);
				SlpAreaSgmnts=SlpAreaSgmnts(idx,:);
				seg_num=seg_num(idx);

				num_seg=numel(ChiSgmnts);

				for jj=1:num_seg
					f=figure(fig_num);
					set(f,'Units','normalized','Position',[0.05 0.1 0.8 0.8],'renderer','painters');
					C=ChiSgmnts{jj};
					bSA=SlpAreaSgmnts{jj,1};
					aSA=SlpAreaSgmnts{jj,2};

					subplot(3,1,1);
					hold on
					plot(C.chi,C.elev,'-k');
					xlabel('Chi')
					ylabel('Elevation (m)')
					title(['Basin ' num2str(basin_num) '- Stream ' num2str(seg_num(jj))]);
					hold off

					subplot(3,1,2);
					hold on
					plot(C.distance,C.elev,'-k');
					xlabel('Distance from Mouth (m)')
					ylabel('Elevation (m)')
					hold off

					ax3=subplot(3,1,3);
					hold on
					scatter(aSA(:,2),aSA(:,1),5,[0.5 0.5 0.5],'+');
					scatter(bSA(:,2),bSA(:,1),20,'k','filled');
					set(ax3,'Xscale','log','Yscale','log','XDir','reverse');
					xlabel('Log Drainage Area');
					ylabel('Log Gradient');
					hold off

					figFile=fullfile(wdir,['StreamSegments_' num2str(basin_num) '_' num2str(fig_num) '_.pdf']);
					print(f,'-dpdf',figFile,'-fillpage');

					drawnow
					fig_num=fig_num+1;

				end
			end
		end
	end

end