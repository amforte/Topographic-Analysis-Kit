function SegmentPlotter(basin_nums)
	% Function to plot all of the chi-Z relationships from a series of picked segments of river networks that result
	% from the 'SegmentPicker' function. 
	%
	% Inputs:
	% basin_nums - row or column vector of basin numbers used for the SegmentPicker you wish to plot together.
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Last Revised Spring 2016 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	num_basins=numel(basin_nums);

	cMap=jet(num_basins);

	f1=figure(1);
	set(f1,'Units','inches','Position',[1.0 1.5 10 10],'renderer','painters','PaperSize',[10 10],'PaperPositionMode','auto');

	for ii=1:num_basins
		basin_num=basin_nums(ii);
		fileName=['PickedSegments_' num2str(basin_num) '.mat'];
		load(fileName,'ChiSgmnts');

		num_seg=numel(ChiSgmnts);

		for jj=1:num_seg
			C=ChiSgmnts{jj};

			subplot(2,1,1);
			hold on
			p1(ii)=plot(C.chi,C.elev,'Color',cMap(ii,:));
			hold off

			subplot(2,1,2);
			hold on
			p2(ii)=plot(C.distance,C.elev,'Color',cMap(ii,:));
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

end