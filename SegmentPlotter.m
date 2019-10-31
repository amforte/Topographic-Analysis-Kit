function SegmentPlotter(basin_nums,varargin)
	%
	% Usage:
	%	SegmentPloter(basin_nums);
	%	SegmentPloter(basin_nums,'name',value,...):
	%
	% Description:
	% 	Function to plot all of the chi-Z relationships, longitudinal profiles, and slope area plots from a series of picked segments of river networks 
	%	that result from the 'SegmentPicker' function. 
	%
	% Required Input:
	% 	basin_nums - row or column vector of basin numbers used for the SegmentPicker you wish to plot together. Code expects that
	%		the mat files saved from SegmentPicker are in the active directory on your matlab path.
	% 
	% Optionl Input:
	%	separate [false] - logical flag to plot all segments as separate figures 
	%	subset [] - list of specific river numbers (i.e. the third column of either the 'Heads' or the 'Outlets' variable) that you wish to include
	%				in the plot. Only valid if you have only provided a single basin number for 'basin_nums'.
	%	label [false]  - logical flag to either label individual streams with the river number (true) or not label them (false, default). If 'separate' flag
	%					is true then the input for label is ignored as the stream number will be in the title of the plots
	%	names [] - option to add an identifying name for streams in a basin when 'label' is set to true. If you have input one basins data, the argument provided
	%				to 'names' should be a string, if you've provided several basins, it should be a cell array of strings in the order of the basins
	%
	% Examples:
	%	SegmentPlotter([4,5],'label',true,'names',{'big','little'});
	%	SegmentPlotter(4,'label',true,'names','big');
	%	SegmentPlotter(1,'subset',[3,5,7,8]);
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 06/18/18 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'SegmentPlotter';
	addRequired(p,'basin_nums',@(x) isnumeric(x));

	addParameter(p,'separate',false,@(x) islogical(x));
	addParameter(p,'subset',[],@(x) isnumeric(x) || isempty(x));
	addParameter(p,'label',false,@(x) islogical(x));
	addParameter(p,'trunks_only',false,@(x) islogical(x));
	addParameter(p,'names',[],@(x) ischar(x) || iscell(x) || isempty(x));
	addParameter(p,'in_dir',[],@(x) isdir(x));
	addParameter(p,'out_dir',[],@(x) isdir(x));
	addParameter(p,'save_fig',false,@(x) islogical(x));

	parse(p,basin_nums,varargin{:});
	basin_nums=p.Results.basin_nums;
	sep=p.Results.separate;
	subset=p.Results.subset;
	lab=p.Results.label;
	nm=p.Results.names;
	in_dir=p.Results.in_dir;
	out_dir=p.Results.out_dir;
	save_fig=p.Results.save_fig;

	num_basins=numel(basin_nums);

	% Check for errors
	if ~isempty(subset) & num_basins>1
		if isdeployed
			warndlg('Providing a list of specific rivers via "subset" is only recognized if only one basin number is provided, ignoring "subset" input')
		else
			warning('Providing a list of specific rivers via "subset" is only recognized if only one basin number is provided, ignoring "subset" input')
		end
		su=0;
	elseif ~isempty(subset) & num_basins==1
		su=1;
	elseif isempty(subset)
		su=0;
	end

	if isempty(in_dir)
		in_dir=pwd;
	end

	if isempty(out_dir)
		out_dir=pwd;
	end

	switch su
	case 0
		if ~sep
			cMap=jet(num_basins);

			f1=figure(1);
			set(f1,'Units','inches','Position',[1.0 1.5 10 10],'renderer','painters','PaperSize',[10 10],'PaperPositionMode','auto');

			for ii=1:num_basins
				basin_num=basin_nums(ii);
				fileName=fullfile(in_dir,['PickedSegments_' num2str(basin_num) '.mat']);
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

				if ischar(nm) || isempty(nm)
					pre='';
				elseif iscell(nm);
					pre=nm{ii};
				end
					

				for jj=1:num_seg
					C=ChiSgmnts{jj};
					bSA=SlpAreaSgmnts{jj,1};
					aSA=SlpAreaSgmnts{jj,2};

					sbplt1=subplot(3,1,1);
					hold on
					p1(ii)=plot(C.chi,C.elev,'Color',cMap(ii,:));
					if lab
						[mc,ix]=max(C.chi);
						me=C.elev(ix);
						text(mc+mc*0.01,me,[pre ' ' num2str(seg_num(jj))],'Color',cMap(ii,:));
					end
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(sbplt1);
				    end						
					hold off

					sbplt2=subplot(3,1,2);
					hold on
					p2(ii)=plot(C.distance,C.elev,'Color',cMap(ii,:));
					if lab
						[md,ix]=max(C.distance);
						me=C.elev(ix);
						text(md+md*0.01,me,[pre ' ' num2str(seg_num(jj))],'Color',cMap(ii,:));
					end
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(sbplt1);
				    end						
					hold off

					ax3=subplot(3,1,3);
					hold on
					scatter(aSA(:,2),aSA(:,1),5,cMap(ii,:),'+');
					p3(ii)=scatter(bSA(:,2),bSA(:,1),20,'MarkerFaceColor',cMap(ii,:),'MarkerEdgeColor','k');
					set(ax3,'Xscale','log','Yscale','log','XDir','reverse');
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(ax3);
				    end						
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

			if save_fig
				figFile=fullfile(out_dir,['StreamSegments_' num2str(basin_num) '.pdf']);
				print(f1,'-dpdf',figFile,'-fillpage');
			end

		elseif sep

			fig_num=1;

			for ii=1:num_basins
				basin_num=basin_nums(ii);
				fileName=fullfile(in_dir,['PickedSegments_' num2str(basin_num) '.mat']);
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
					set(f,'Units','inches','Position',[1.0 1.5 10 10],'renderer','painters','PaperSize',[10 10],'PaperPositionMode','auto');
					C=ChiSgmnts{jj};
					bSA=SlpAreaSgmnts{jj,1};
					aSA=SlpAreaSgmnts{jj,2};

					sbplt1=subplot(3,1,1);
					hold on
					plot(C.chi,C.elev,'-k');
					xlabel('Chi')
					ylabel('Elevation (m)')
					title(['Basin ' num2str(basin_num) '- Stream ' num2str(seg_num(jj))]);
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(sbplt1);
				    end						
					hold off

					sbplt2=subplot(3,1,2);
					hold on
					plot(C.distance,C.elev,'-k');
					xlabel('Distance from Mouth (m)')
					ylabel('Elevation (m)')
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(sbplt2);
				    end						
					hold off

					ax3=subplot(3,1,3);
					hold on
					scatter(aSA(:,2),aSA(:,1),5,[0.5 0.5 0.5],'+');
					scatter(bSA(:,2),bSA(:,1),20,'k','filled');
					set(ax3,'Xscale','log','Yscale','log','XDir','reverse');
					xlabel('Log Drainage Area');
					ylabel('Log Gradient');
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(ax3);
				    end						
					hold off		

					if save_fig
						figFile=fullfile(out_dir,['StreamSegments_' num2str(basin_num) '_' num2str(fig_num) '_.pdf']);
						print(f1,'-dpdf',figFile,'-fillpage');	
					end		

					fig_num=fig_num+1;
					drawnow
				end
			end
		end

	case 1

		if ~sep

			f1=figure(1);
			set(f1,'Units','inches','Position',[1.0 1.5 10 10],'renderer','painters','PaperSize',[10 10],'PaperPositionMode','auto');


			basin_num=basin_nums;
			fileName=fullfile(in_dir,['PickedSegments_' num2str(basin_num) '.mat']);
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

				sbplt1=subplot(3,1,1);
				hold on
				plot(C.chi,C.elev,'-k');
				if lab
					[mc,ix]=max(C.chi);
					me=C.elev(ix);
					text(mc+mc*0.01,me,[nm ' ' num2str(seg_num(jj))]);
				end
				if ~verLessThan('matlab','9.5')
			        disableDefaultInteractivity(sbplt1);
			    end					
				hold off

				sbplt2=subplot(3,1,2);
				hold on
				plot(C.distance,C.elev,'-k');
				if lab
					[md,ix]=max(C.distance);
					me=C.elev(ix);
					text(md+md*0.01,me,[nm ' ' num2str(seg_num(jj))]);
				end
				if ~verLessThan('matlab','9.5')
			        disableDefaultInteractivity(sbplt1);
			    end					
				hold off

				ax3=subplot(3,1,3);
				hold on
				scatter(aSA(:,2),aSA(:,1),5,[0.5 0.5 0.5],'+');
				scatter(bSA(:,2),bSA(:,1),20,'k','filled');
				set(ax3,'Xscale','log','Yscale','log','XDir','reverse');
				if ~verLessThan('matlab','9.5')
			        disableDefaultInteractivity(ax3);
			    end					
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

			if save_fig
				figFile=fullfile(out_dir,['StreamSegments_' num2str(basin_num) '.pdf']);
				print(f1,'-dpdf',figFile,'-fillpage');
			end			

		elseif sep

			fig_num=1;

			for ii=1:num_basins
				basin_num=basin_nums(ii);
				fileName=fullfile(in_dir,['PickedSegments_' num2str(basin_num) '.mat']);
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
					set(f,'Units','inches','Position',[1.0 1.5 10 10],'renderer','painters','PaperSize',[10 10],'PaperPositionMode','auto');
					C=ChiSgmnts{jj};
					bSA=SlpAreaSgmnts{jj,1};
					aSA=SlpAreaSgmnts{jj,2};

					sbplt1=subplot(3,1,1);
					hold on
					plot(C.chi,C.elev,'-k');
					xlabel('Chi')
					ylabel('Elevation (m)')
					title(['Basin ' num2str(basin_num) '- Stream ' num2str(seg_num(jj))]);
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(sbplt1);
				    end						
					hold off

					sbplt2=subplot(3,1,2);
					hold on
					plot(C.distance,C.elev,'-k');
					xlabel('Distance from Mouth (m)')
					ylabel('Elevation (m)')
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(sbplt2);
				    end						
					hold off

					ax3=subplot(3,1,3);
					hold on
					scatter(aSA(:,2),aSA(:,1),5,[0.5 0.5 0.5],'+');
					scatter(bSA(:,2),bSA(:,1),20,'k','filled');
					set(ax3,'Xscale','log','Yscale','log','XDir','reverse');
					xlabel('Log Drainage Area');
					ylabel('Log Gradient');
					if ~verLessThan('matlab','9.5')
				        disableDefaultInteractivity(ax3);
				    end						
					hold off

					if save_fig
						figFile=fullfile(out_dir,['StreamSegments_' num2str(basin_num) '_' num2str(fig_num) '_.pdf']);
						print(f1,'-dpdf',figFile,'-fillpage');	
					end

					fig_num=fig_num+1;
					drawnow
				end
			end
		end
	end

end