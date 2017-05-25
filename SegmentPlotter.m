function SegmentPlotter(basin_nums,varargin)
	% Function to plot all of the chi-Z relationships from a series of picked segments of river networks that result
	% from the 'SegmentPicker' function. 
	%
	% Required Input:
	% 	basin_nums - row or column vector of basin numbers used for the SegmentPicker you wish to plot together.
	% 
	% Optionl Input:
	%	separate [false] - logical flag to plot all segments as separate figures 
	%	subset [] - list of specific river numbers (i.e. the first column of either the 'Heads' or the 'Outlets' variable) that you wish to include
	%				in the plot. Only valid if you have only provided a single basin number for 'basin_nums'.
	%	label [false]  - logical flag to either label individual streams with the river number (true) or not label them (false, default). If 'separate' flag
	%					is true then the input for label is ignored as the stream number will be in the title of the plots
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Last Revised Summer 2017 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'SegmentPlotter';
	addRequired(p,'basin_nums',@(x) isnumeric(x));

	addParamValue(p,'separate',false,@(x) islogical(x));
	addParamValue(p,'subset',[],@(x) isnumeric(x));
	addParamValue(p,'label',false,@(x) islogical(x));

	parse(p,basin_nums,varargin{:});
	basin_nums=p.Results.basin_nums;
	sep=p.Results.separate;
	sub=p.Results.subset;
	lab=p.Results.label;

	num_basins=numel(basin_nums);

	% Check for errors
	if ~isempty(sub) & num_basins>1
		warning('Providing a list of specific rivers via "subset" is only recognized if only one basin number is provided, ignoring "subset" input')
		su=0;
	elseif ~isempty(sub) & num_basins==1
		su=1;
	elseif isempty(sub)
		su=0;
	end

	switch su
	case 0
		if ~sep
			cMap=jet(num_basins);

			f1=figure(1);
			set(f1,'Units','inches','Position',[1.0 1.5 10 10],'renderer','painters','PaperSize',[10 10],'PaperPositionMode','auto');

			for ii=1:num_basins
				basin_num=basin_nums(ii);
				fileName=['PickedSegments_' num2str(basin_num) '.mat'];
				load(fileName,'ChiSgmnts');

				vInfo=who('-file',fileName);
				if ismember('Heads',vInfo)
					load(fileName,'Heads');
					seg_num=Heads(:,1);
				elseif ismember('Outlets',vInfo)
					load(filename,'Outlets');
					seg_num=Outlets(:,1);
				end

				num_seg=numel(ChiSgmnts);

				for jj=1:num_seg
					C=ChiSgmnts{jj};

					subplot(2,1,1);
					hold on
					p1(ii)=plot(C.chi,C.elev,'Color',cMap(ii,:));
					if lab
						[mc,ix]=max(C.chi);
						me=C.elev(ix);
						text(mc+mc*0.01,me,num2str(seg_num(jj)));
					end
					hold off

					subplot(2,1,2);
					hold on
					p2(ii)=plot(C.distance,C.elev,'Color',cMap(ii,:));
					if lab
						[md,ix]=max(C.distance);
						me=C.elev(ix);
						text(md+md*0.01,me,num2str(seg_num(jj)));
					end
					hold off
				end
			LegendEnt{ii}=['Basin ' num2str(basin_num)];
			end

			subplot(2,1,1);
			hold on
				xlabel('Chi')
				ylabel('Elevation (m)')
				title('Chi - Z')
				legend(p1,LegendEnt);
			hold off

			subplot(2,1,2);
			hold on
				xlabel('Distance from Mouth (m)')
				ylabel('Elevation (m)')
				title('Long Profile')
				legend(p2,LegendEnt);
			hold off
		elseif sep

			fig_num=1;

			for ii=1:num_basins
				basin_num=basin_nums(ii);
				fileName=['PickedSegments_' num2str(basin_num) '.mat'];
				load(fileName,'ChiSgmnts');
				vInfo=who('-file',fileName);
				if ismember('Heads',vInfo)
					load(fileName,'Heads');
					seg_num=Heads(:,1);
				elseif ismember('Outlets',vInfo)
					load(filename,'Outlets');
					seg_num=Outlets(:,1);
				end

				num_seg=numel(ChiSgmnts);

				for jj=1:num_seg
					f=figure(fig_num);
					set(f,'Units','inches','Position',[1.0 1.5 10 10],'renderer','painters','PaperSize',[10 10],'PaperPositionMode','auto');
					C=ChiSgmnts{jj};

					subplot(2,1,1);
					hold on
					plot(C.chi,C.elev,'-k');
					xlabel('Chi')
					ylabel('Elevation (m)')
					title(['Basin ' num2str(basin_num) '- Stream ' num2str(seg_num(jj))]);
					hold off

					subplot(2,1,2);
					hold on
					plot(C.distance,C.elev,'-k');
					xlabel('Distance from Mouth (m)')
					ylabel('Elevation (m)')
					hold off


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
			fileName=['PickedSegments_' num2str(basin_num) '.mat'];
			load(fileName,'ChiSgmnts');

			vInfo=who('-file',fileName);
			if ismember('Heads',vInfo)
				load(fileName,'Heads');
				seg_num=Heads(:,1);
			elseif ismember('Outlets',vInfo)
				load(filename,'Outlets');
				seg_num=Outlets(:,1);
			end

			idx=ismember(seg_num,sub);
			ChiSgmnts=ChiSgmnts(idx);
			seg_num=seg_num(idx);

			num_seg=numel(ChiSgmnts);

			for jj=1:num_seg
				C=ChiSgmnts{jj};

				subplot(2,1,1);
				hold on
				plot(C.chi,C.elev,'-k');
				if lab
					[mc,ix]=max(C.chi);
					me=C.elev(ix);
					text(mc+mc*0.01,me,num2str(seg_num(jj)));
				end
				hold off

				subplot(2,1,2);
				hold on
				plot(C.distance,C.elev,'-k');
				if lab
					[md,ix]=max(C.distance);
					me=C.elev(ix);
					text(md+md*0.01,me,num2str(seg_num(jj)));
				end
				hold off
			end

			subplot(2,1,1);
			hold on
				xlabel('Chi')
				ylabel('Elevation (m)')
				title('Chi - Z')
			hold off

			subplot(2,1,2);
			hold on
				xlabel('Distance from Mouth (m)')
				ylabel('Elevation (m)')
				title('Long Profile')
			hold off

		elseif sep

			fig_num=1;

			for ii=1:num_basins
				basin_num=basin_nums(ii);
				fileName=['PickedSegments_' num2str(basin_num) '.mat'];
				load(fileName,'ChiSgmnts');
				vInfo=who('-file',fileName);
				if ismember('Heads',vInfo)
					load(fileName,'Heads');
					seg_num=Heads(:,1);
				elseif ismember('Outlets',vInfo)
					load(filename,'Outlets');
					seg_num=Outlets(:,1);
				end

				idx=ismember(seg_num,sub);
				ChiSgmnts=ChiSgmnts(idx);
				seg_num=seg_num(idx);

				num_seg=numel(ChiSgmnts);

				for jj=1:num_seg
					f=figure(fig_num);
					set(f,'Units','inches','Position',[1.0 1.5 10 10],'renderer','painters','PaperSize',[10 10],'PaperPositionMode','auto');
					C=ChiSgmnts{jj};

					subplot(2,1,1);
					hold on
					plot(C.chi,C.elev,'-k');
					xlabel('Chi')
					ylabel('Elevation (m)')
					title(['Basin ' num2str(basin_num) '- Stream ' num2str(seg_num(jj))]);
					hold off

					subplot(2,1,2);
					hold on
					plot(C.distance,C.elev,'-k');
					xlabel('Distance from Mouth (m)')
					ylabel('Elevation (m)')
					hold off


					fig_num=fig_num+1;
					drawnow
				end
			end
		end
	end

end