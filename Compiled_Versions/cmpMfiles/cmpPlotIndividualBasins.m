function cmpPlotIndividualBasins(wdir,location_of_data_files,varargin)
	% Description:
	% 	Function takes outputs from 'ProcessRiverBasins' function and makes and saves plots for each basin with stream profiles, chi-z, and slope area
	%
	% Required Inputs:
	%	wdir - full path of working directory
	% 	location_of_data_files - name folder within working directory that contains the mat files from 'ProcessRiverBasins'
	%
	% Optional Inputs:
	%	location_of_subbasins ['SubBasins'] - name of folder that contains subbasins of interest (if you created subbasins using
	%		"SubDivideBigBasins"), expected to be within the main Basin folder provided with "location_of_data_files"
	%	bin_size [500] - bin size (in map units) for binning slope area data.
	%
   	% Examples if running for the command line, minus OS specific way of calling main TAK function:
    %   PlotIndividualBasins /path/to/wdir Basins
    %   PlotIndividualBasins /path/to/wdir Basins locations_of_subbasins MySubBasins
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 06/18/18 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if isdeployed
		if ~isempty(varargin)
			varargin=varargin{1};
		end
	end

	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'cmpPlotIndividualBasins';
	addRequired(p,'wdir',@(x) ischar(x));
	addRequired(p,'location_of_data_files',@(x) ischar(x));

	addParameter(p,'location_of_subbasins','SubBasins',@(x) ischar(x));
	addParameter(p,'bin_size',500,@(x) isscalar(x) && isnumeric(x));

	parse(p,wdir,location_of_data_files,varargin{:});
	wdir=p.Results.wdir;
	location_of_data_files=p.Results.location_of_data_files;

	location_of_subbasins=p.Results.location_of_subbasins;
	bin_size=p.Results.bin_size;

	
	%% Build File List
	% Get Basin Numbers
	AllFullFiles=dir(fullfile(wdir,location_of_data_files,'*_Data.mat'));
	num_basins=numel(AllFullFiles);
	basin_nums=zeros(num_basins,1);
	for jj=1:num_basins
		fileName=AllFullFiles(jj,1).name;
		basin_nums(jj)=sscanf(fileName,'%*6s %i'); %%%
	end

	FileCell=cell(num_basins,1);
	for kk=1:num_basins
		basin_num=basin_nums(kk);
		SearchAllString=fullfile(wdir,location_of_data_files,['*_' num2str(basin_num) '_Data.mat']);
		SearchSubString=fullfile(wdir,location_of_data_files,location_of_subbasins,['*_' num2str(basin_num) '_DataSubset*.mat']);

		if numel(dir(SearchSubString))>0
			Files=dir(SearchSubString);
		else
			Files=dir(SearchAllString);
		end

		FileCell{kk}=Files;
	end
	fileList=vertcat(FileCell{:});
	num_files=numel(fileList);

	for ii=1:num_files

		FileName=fullfile(fileList(ii,1).folder,fileList(ii,1).name);

		load(FileName,'DEMcc','ChiOBJc','Sc','Ac','RiverMouth','drainage_area');
		[bs,ba,~,~,aa,ag,~]=sa(DEMcc,Sc,Ac,ChiOBJc,bin_size);

		f1=figure(1);
		set(f1,'Units','inches','Position',[1.0 1.5 10 10],'renderer','painters','PaperSize',[10 10],'PaperPositionMode','auto');

		clf
		subplot(3,1,1);
		hold on
		title(['Basin Number: ' num2str(RiverMouth(:,3)) ' - Drainage Area: ' num2str(drainage_area)]);
		plotdz(Sc,DEMcc,'dunit','km','color','k');
		xlabel('Distance (km)');
		ylabel('Elevation (m)');
		hold off

		subplot(3,1,2);
		hold on
		plotdz(Sc,DEMcc,'distance',getnal(Sc,ChiOBJc),'color','k');
		xlabel('Chi');
		ylabel('Elevation (m)');
		hold off

		a1=subplot(3,1,3);
		hold on
		scatter(aa,ag,5,'k','+');
		scatter(ba,bs,'o','MarkerFaceColor','b','MarkerEdgeColor','k');
		set(a1,'XScale','log','YScale','log','XDir','reverse');	
		xlabel('Log Area');
		ylabel('Log Slope');
		hold off

		fileName=fullfile(wdir,location_of_data_files,['BasinPlot_' num2str(RiverMouth(:,3)) '.pdf']);
		print(f1,'-dpdf',fileName);
		close all
	end
end

function [bs,ba,bc,bd,a,g,C]=sa(DEM,S,A,Cg,bin_size)
	% Modified slope area function that uses the smooth length to
	%	to determine the number of bins and uses those same bins
	%	to find mean values of chi and distance for plotting
	%	purposes

	minX=min(S.distance);
	maxX=max(S.distance);
	b=[minX:bin_size:maxX+bin_size];

	numbins=round(max([numel(b) numel(S.IXgrid)/10]));

	an=getnal(S,A.*A.cellsize^2);
	z=getnal(S,DEM);
	cn=getnal(S,Cg);
	gn=gradient(S,z,'unit','tangent','method','robust','drop',20);
	gn=smooth(gn,3);

	% Run through STREAMobj2XY so chi and everything else are same size
	[~,~,a,g,C,d]=STREAMobj2XY(S,an,gn,cn,S.distance);
	% Remove NaNs
	a(isnan(a))=[];
	g(isnan(g))=[];
	d(isnan(d))=[];
	C(isnan(C))=[];

	mina=min(a);
	maxa=max(a);

    edges = logspace(log10(mina-0.1),log10(maxa+1),numbins+1);
    try
    	% histc is deprecated
    	[ix]=discretize(a,edges);
    catch
	    [~,ix] = histc(a,edges);
	end

	ba=accumarray(ix,a,[numbins 1],@median,nan);
	bs=accumarray(ix,g,[numbins 1],@(x) mean(x(~isnan(x))),nan);
	bd=accumarray(ix,d,[numbins 1],@mean,nan);
	bc=accumarray(ix,C,[numbins 1],@mean,nan);

	% Filter negatives
	idx=bs>=0 & ba>=0 & bc>=0 & bd>=0;
	bs=bs(idx);
	ba=ba(idx);
	bc=bc(idx);
	bd=bd(idx);

	idx=a>=0 & g>=0 & d>=0 & C>=0;
	a=a(idx);
	g=g(idx);
	d=d(idx);
	C=C(idx);
end