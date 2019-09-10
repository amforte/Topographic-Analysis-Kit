function cmpPlotKsn(wdir,MatFile,ksn,varargin)
	% Description:
	% 	Function to plot a map of normalized channel steepness on a hillshade colored
	% 	by elevation.
	%
	% Required Inputs:
	%	DEM -  Digital Elevation as a GRIDobj used to produce the provided ksn data
	%	FD - Flow routing as a FLOWobj used to proudce the provided ksn data
	%	ksn - ksn data as a shapefile (as ouput from KsnProfiler, ProcessRiverBasins
	%		KsnChiBatch)
	% 
	% Optional Inputs:
	%	knicks [] - location of knickpoints as shapefile (as output by FindBasinKnicks or KsnProfiler)
	%	ksn_lim [] - 1 x n vector setting the min and max for the color scaling for ksn. If left blank
	%		will default to 0 and the maximum in the dataset
	%
   	% Examples if running for the command line, minus OS specific way of calling main TAK function:
    %   PlotKsn /path/to/wdir Topo.mat ksn.shp
    %   PlotKsn /path/to/wdir Topo.mat ksn.shp knicks knicks.shp
   	%   PlotKsn /path/to/wdir Topo.mat ksn.shp knicks knicks.shp ksn_lim [100 200]
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
    p.FunctionName = 'cmpPlotKsn';
	addRequired(p,'wdir',@(x) ischar(x));
	addRequired(p,'MatFile',@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.mat'))));
    addRequired(p,'ksn',@(x) regexp(x,regexptranslate('wildcard','*.shp')));

    addParameter(p,'knicks',[],@(x) regexp(x,regexptranslate('wildcard','*.shp')));
    addParameter(p,'ksn_lim',[],@(x) isnumeric(x) && numel(x)==2);

    parse(p,wdir,MatFile,ksn,varargin{:});
    wdir=p.Results.wdir;
    MatFile=p.Results.MatFile;
    ksn=p.Results.ksn;

    knks=p.Results.knicks;
    ksn_lim=p.Results.ksn_lim;	

    % Determine the type of input
    MatFile=fullfile(wdir,MatFile);
    D=whos('-file',MatFile);
    VL=cell(numel(D),1);
    for ii=1:numel(D);
        VL{ii}=D(ii,1).name;
    end

    if any(strcmp(VL,'DEM')) & any(strcmp(VL,'FD')) & any(strcmp(VL,'A')) & any(strcmp(VL,'S'))
        load(MatFile,'DEM','FD');
    elseif any(strcmp(VL,'DEMoc')) & any(strcmp(VL,'FDc')) & any(strcmp(VL,'Ac')) & any(strcmp(VL,'Sc'))
        load(MatFile,'DEMoc','FDc');
        DEM=DEMoc;
        FD=FDc;
    end 

	if ischar(ksn) & ~isempty(regexp(ksn,regexptranslate('wildcard','*.shp')))
		ksn=shaperead(fullfile(wdir,ksn));
	else
		error('Input to "ksn" not recognized as the name of a shapefile');
	end
		
	num_seg=numel(ksn);

	sx=cell(num_seg,1);
	sy=cell(num_seg,1);
	sk=cell(num_seg,1);
	for ii=1:num_seg
		sx{ii,1}=ksn(ii,1).X(:);
		sy{ii,1}=ksn(ii,1).Y(:);
		if isfield(ksn,'ksn')
			sk{ii,1}=ones(numel(sx{ii,1}),1)*ksn(ii,1).ksn;
		elseif isfield(ksn,'fit_ksn')
			sk{ii,1}=ones(numel(sx{ii,1}),1)*ksn(ii,1).fit_ksn;
		else
			error('There is no valid field in the provided shapefile named "ksn"')
		end
	end

	sx=vertcat(sx{:});
	sy=vertcat(sy{:});
	sk=vertcat(sk{:});

	ix=coord2ind(DEM,sx,sy);
	idx=isnan(ix);

	ix(idx)=[];
	sk(idx)=[];

	W=GRIDobj(DEM,'logical');
	W.Z(ix)=true;
	S=STREAMobj(FD,W);

	[~,loc,~]=unique(ix);
	sk=sk(loc);

	f1=figure(1);
	set(f1,'Visible','off');

	[RGB]=imageschs(DEM,DEM,'colormap','gray');
	[~,R]=GRIDobj2im(DEM);

	imshow(flipud(RGB),R);
	axis xy
	hold on
	colormap(ksncolor(20));
	plotc(S,sk);
	if isempty(ksn_lim)
		caxis([0 max(sk)]);
	else
		caxis([min(ksn_lim) max(ksn_lim)])
	end
	c1=colorbar;
	ylabel(c1,'Normalized Channel Steepness')
	if ~isempty(knks)
		if ischar(knks) & logical(regexp(knks,regexptranslate('wildcard','*.shp')))
			knk=shaperead(fullfile(wdir,knks));
			knkx=[knk.X];
			knky=[knk.Y];
			scatter(knkx,knky,100,'w','p','filled','MarkerEdgeColor','k');
		end
	end

	hold off
	set(f1,'Visible','on','Units','normalized','Position',[0.05 0.1 0.8 0.8]);
	print(f1,'-dpdf',fullfile(wdir,'KsnMap.pdf'),'-bestfit');
	close(f1);

end