function [KnickTable]=FindBasinKnicks(Basin_Data_File,plot_result,varargin)
	%
	% Usage:
	%	[KnickTable]=FindBasinKnicks(Basin_Data_File,plot_result);
	%	[KnickTable]=FindBasinKnicks(Basin_Data_File,plot_result,'name',value,...);
	%
	% Description:
	% 	Function for manually selecting knickpoints within a Basin_Data_File (i.e. result of ProcessRiverBasins). 
	% 	Choose knickpoints on Chi-Elevation plot with mouse clicks and press return when you have selected
	% 	all the knickpoints for a given stream segment. As you progress through, knickpoints you have already picked 
	% 	(i.e. on shared portions of river profiles) will be displayed as red dots. If you're interested in trying out
	%	an automated method of finding knickpoints, try 'knickpointfinder' included with TopoToolbox. If you choose to 
	%	classify knickpoints ('classify_knicks' = true) The code expects you to input a number or character 
	%	to categorize the knickpoint higlighted in red. You must be consistent in your choice (i.e. you must either use 
	%	numbers for all of the classifications or characters for all the classifications within a given run), mixing numbers 
	%	and characters will result in an error at the end of the run. For entering characters, it's recommended you keep these 
	%	short strings without spaces (i.e. entries supported into a shapefile a attribute table), e.g. knick or bound 
	%
	% Required Inputs:
	% 	Basin_Data_File - full file path to a saved result from the ProcessRiverBasins script
	% 	plot_result - logical flag to either plot the results (true) or not (false) 
	%
	% Optional Inputs
	%	classify_knicks [false] - logical flag to provide a classification for each chosen knickpoint
	% 	ref_concavity [0.5] - reference concavity for chi calculation
	%	save_mat [true] - logical flag to save output mat file containing the KnickPoints array. The name of the file will 
	%		be 'Knicks_NUM.mat' where NUM is the river number. Do not change the file name if you want to plot knickpoints
	%		using 'MakeCombinedSwath'.
	%	shape_name [] - character string to name output shapefile (without .shp), if no input is provided then
	%		no shapefile is output
	%
	% Outputs:
	%	KnickTable - table with one row for each selected knickpoints. If classify_knicks is false, will have with columns x_coord, 
	%		y_coord, elevation, distance, and chi. If classify_knicks is true, will have a sixth column containing the classification
	%		of the knickpoints.
	%	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 07/02/18 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'FindBasinKnicks';
	addRequired(p,'Basin_Data_File',@(x) ischar(x));
	addRequired(p,'plot_result',@(x) islogical(x));

	addParameter(p,'classify_knicks',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'ref_concavity',0.5,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'shape_name',[],@(x) ischar(x));
	addParameter(p,'save_mat',true,@(x) isscalar(x) && islogical(x));

	parse(p,Basin_Data_File,plot_result,varargin{:});
	Basin_Data_File=p.Results.Basin_Data_File;
	plot_result=p.Results.plot_result;

	classify_knicks=p.Results.classify_knicks;
	theta_ref=p.Results.ref_concavity;
	shape_name=p.Results.shape_name;
	save_mat=p.Results.save_mat;

	% Load in File Contents
	load(Basin_Data_File);
    
	% De-Densify Network
	if drainage_area>20
		S=modify(Sc,'streamorder','>1');
		if isempty(S.x)
			S=Sc;
		end
	else
		S=Sc;
	end

	% Find Channel Heads of Channel Network
	ChXY=streampoi(S,'channelheads','xy');
	ChIX=coord2ind(DEMoc,ChXY(:,1),ChXY(:,2));
	NumHeads=numel(ChIX);

	% Create Logical Seed Raster
	SEED=GRIDobj(DEMoc,'logical');

	f1=figure(1);
	set(f1,'Units','normalized','Position',[0.05 0.1 0.45 0.8],'renderer','painters');
	
	f2=figure(2);
	set(f2,'Units','normalized','Position',[0.5 0.1 0.45 0.8],'renderer','painters');


	% Loop Through and Create Individual Channel Paths
	Channels=cell(NumHeads,1);
	Chi=cell(NumHeads,1);
	Knicks=cell(NumHeads,1);

	% disp('Choose knickpoints with mouse clicks and press return when done')
	uiwait(msgbox('Choose knickpoints with mouse clicks on  chi - elevation plot and press return when done picking for that stream'));

	for ii=1:NumHeads
		ChOI=ChIX(ii);
		ChR=SEED;

		ChR.Z(ChOI)=true;

		SS=modify(S,'downstreamto',ChR);

		Channels{ii}=SS;

		C=chiplot(SS,DEMcc,Ac,'a0',1,'mn',theta_ref,'plot',false);
		Chi{ii}=C;

		chi=C.chi;
		elev=C.elev;
		x=C.x;
		y=C.y;
		d=C.distance;

		% Plot Working River Within Network
		figure(f2);
		clf
		hold on
		imageschs(DEMoc,DEMoc,'colormap','gray');
		plot(S,'-w');
		plot(SS,'-r');
		hold off

		% Chi-Z Plot
		figure(f1);
		clf
		hold on
		scatter(C.chi,C.elev,10,'k','filled');

        KC=vertcat(Knicks{1:ii-1});
        
		if ii>1 && isempty(KC)~=1
			KC=vertcat(Knicks{1:ii-1});
			xx=KC(:,1); yy=KC(:,2); cc=KC(:,5); ee=KC(:,3);
			currentRiv=[x y];
			pastKnicks=[xx yy];
			cidx=ismember(pastKnicks,currentRiv,'rows');
			cc=cc(cidx); ee=ee(cidx);
			scatter(cc,ee,30,'r','filled');
		else
			;
		end

		xlabel('\chi');
		ylabel('Elevation (m)');
		title([num2str(NumHeads-ii) ' channels remaining']);
		hold off

		[c,e]=ginput;

		if numel(c)>=1;
			% Find point closest to pick
			for jj=1:numel(c);
				coi=c(jj);
				[~,idx]=min(abs(chi-coi));
				knp(jj,:)=[x(idx) y(idx) elev(idx) d(idx) chi(idx)];
			end

			Knicks{ii}=knp;
		end

		if classify_knicks
			if numel(c)>=1
				figure(f1);
				for jj=1:numel(c);
					hold on
					s1=scatter(knp(jj,5),knp(jj,3),40,'g');
					cl=inputdlg('Enter the classification for the selected knickpoint:','Knickpoint Classification');
					hold off
					delete(s1);

					cl_num=str2num(cl{1});
					if isempty(cl_num)
						knpClass{jj,1}=cl{1};
					else 
						knpClass{jj,1}=double(cl_num);
					end

				end	
			end	
			knpClasses{ii}=knpClass;
		end

	end

	KnickPoints=vertcat(Knicks{:});
	[KnickPoints,idx,~]=unique(KnickPoints,'rows');

	KnickTable=array2table(KnickPoints,'VariableNames',{'x_coord','y_coord','elevation','distance','chi'});

	if classify_knicks
		if exist('knpClasses')
			knpClasses=vertcat(knpClasses{:});
			knpClasses=knpClasses(idx);
			KnickTable.classification=knpClasses;
		end
	end

	close all

	if plot_result
		f1=figure(1);
		set(f1,'Units','normalized','Position',[0.05 0.1 0.8 0.8],'renderer','painters');
       
		subplot(1,2,1)
		hold on
		[RGB]=imageschs(DEMoc,DEMoc,'colormap','gray');
		[~,R]=GRIDobj2im(DEMoc);
		imshow(flipud(RGB),R);
		axis xy
		colormap(jet);
        caxis([0 max(KnickPoints(:,3))]);
		plot(S,'-w');
		scatter(KnickPoints(:,1),KnickPoints(:,2),20,KnickPoints(:,3),'filled');
		hold off

		subplot(1,2,2)
		hold on
        plotdz(S,DEMoc,'Color',[0.5 0.5 0.5]);
		plotdz(S,DEMcc,'color','k');
        caxis([0 max(KnickPoints(:,3))]);
		scatter(KnickPoints(:,4),KnickPoints(:,3),20,KnickPoints(:,3),'filled');
		c1=colorbar;
		ylabel(c1,'Knickpoint Elevation (m)');
		hold off
	end

	if save_mat
		save(['Knicks_' num2str(RiverMouth(:,3)) '.mat'],'KnickTable','-v7.3');
	end

	if ~isempty(shape_name)
		MS=struct;
		for ii=1:size(KnickPoints,1)
			MS(ii,1).ID=ii;
			MS(ii,1).Geometry='Point';
			MS(ii,1).X=KnickPoints(ii,1);
			MS(ii,1).Y=KnickPoints(ii,2);
			MS(ii,1).elev=KnickPoints(ii,3);
			MS(ii,1).dist=KnickPoints(ii,4);
			MS(ii,1).chi=KnickPoints(ii,5);
			if classify_knicks
				MS(ii,1).class=knpClasses{ii};
			end
		end
		shp_out=[shape_name '.shp'];
		shapewrite(MS,shp_out);
	end

% Function End
end


