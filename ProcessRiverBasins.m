function ProcessRiverBasins(DEM,FD,S,river_mouths,varargin)
	% Function takes grid object outputs from MakeStreams script (DEM,FD,A,S), a series of x,y coordinates of river mouths,
	% and outputs clipped dem, stream network, variout topographic metrics, and river values (ks, ksn, chi)
	%
	% Required Inputs:
	% 		DEM - GRIDobj of the digital elevation model of your area loaded into the workspace
	% 		FD - FLOWobj of the flow direction of your area loaded into the workspace
	% 		S - STREAMobj of the stream network of your area loaded into the workspace	
	% 		river_mouths - nx3 matrix of river mouths with x, y, and a number identifying the stream/basin of interest.
	%					Needs to be in the same projection as the input DEM
	% Optional Inputs:
	% 		threshold_area [1e6] - minimum accumulation area to define streams in meters squared
	% 		smooth_distance [1000] - smoothing distance in meters for averaging along ksn, suggested value is 1000 meters
	% 		ref_concavity [0.45] - reference concavity for calculating ksn, suggested value is 0.45
	% 		write_arc_files [false] - set value to true to output a geotiff of the DEM and a shapefile of the ksn, false to not output arc files
	%
	% Examples:
	%		ProcessRiverBasins(DEM,FD,S,RiverMouths);
	%		ProcessRiverBasins(DEM,FD,S,RiverMouths,'ref_concavity',0.5,'write_arc_files',true);
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Last Revised Spring 2016 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'ProcessRiverBasins';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'FD',@(x) isa(x,'FLOWobj'));
	addRequired(p,'S',@(x) isa(x,'STREAMobj'));
	addRequired(p,'river_mouths',@(x) isnumeric(x) && size(x,2)==3);

	addParamValue(p,'ref_concavity',0.45,@(x) isscalar(x) && isnumeric(x));
	addParamValue(p,'threshold_area',1e6,@(x) isscalar(x) && isnumeric(x));
	addParamValue(p,'smooth_distance',1000,@(x) isscalar(x) && isnumeric(x));
	addParamValue(p,'write_arc_files',false,@(x) isscalar(x));

	parse(p,DEM,FD,S,river_mouths,varargin{:});
	DEM=p.Results.DEM;
	FD=p.Results.FD;
	S=p.Results.S;
	river_mouths=p.Results.river_mouths;

	ref_concavity=p.Results.ref_concavity;
	threshold_area=p.Results.threshold_area;
	smooth_distance=p.Results.smooth_distance;
	write_arc_files=p.Results.write_arc_files;

	current=pwd;

	% Clippin gout NaNs, max elev set arbitrarily large but below the 32,768 internal NaN value.
	max_elev=10000;
	min_elev=-200;
	IDX=DEM<max_elev & DEM>min_elev;
	DEM=crop(DEM,IDX,nan);
	cd(current);

	disp('Snapping river mouths to stream network')
	xi=river_mouths(:,1);
	yi=river_mouths(:,2);
	riv_nums=river_mouths(:,3);

	num_basins=numel(xi);

	[xn,yn]=snap2stream(S,xi,yi);

	RM=[xn yn riv_nums];

	RM=sortrows(RM,1);

	num_basins=numel(xn);

	for ii=1:num_basins
		disp(['Working on Basin Number ' num2str(ii) ' of ' num2str(num_basins) ' total basins'])
		xx=RM(ii,1);
		yy=RM(ii,2);
		basin_num=RM(ii,3);

		RiverMouth=[xx yy basin_num];

		% Build dependenc map and clip out drainage basins
		I=dependencemap(FD,xx,yy);
		DEMoc=crop(DEM,I,nan);

		% Calculate drainage area
		dep_map=GRIDobj2mat(I);
		num_pix=sum(sum(dep_map));
		drainage_area=(num_pix*DEMoc.cellsize*DEMoc.cellsize)/(1e6);

		% Find weighted centroid of drainage basin
		[Cx,Cy]=FindCentroid(DEMoc);
		Centroid=[Cx Cy];

		% Generate new stream map
		FDc=FLOWobj(DEMoc,'preprocess','carve');
		Ac=flowacc(FDc);
		DEMc=imposemin(FDc,DEMoc,.001); % Hyrdologically condition DEM
		DEMc_res=DEMc.cellsize;
		min_area=floor(threshold_area/(DEMc_res*DEMc_res));
		isstream=Ac>min_area;
		Sc=STREAMobj(FDc,isstream);

		% Extract largest drainage
		SLc=klargestconncomps(Sc,1);	

		% Calculate chi
		Chic=chiplot(SLc,DEMc,Ac,'a0',1,'mn',ref_concavity,'plot',false);
		% Make chi obj
		[ChiOBJc]=MakeChiMapInd(Chic,DEMc);

		% Calculate ksn
		SAc=slopearea(Sc,DEMc,Ac,'plot',false);
		Goc=gradient8(DEMoc);
		[GcREF,~]=CompGrad(DEMoc,FDc,Ac,1,ref_concavity);
		[GcSAc,~]=CompGrad(DEMoc,FDc,Ac,1,-1*SAc.theta);

		% Calculate channel steepness using concavity value from slope area fit
		KSc=GcSAc./(Ac.*(Ac.cellsize^2)).^SAc.theta;
		MSc=STREAMobj2mapstruct(Sc,'seglength',smooth_distance,'attributes',...
			{'ks' KSc @mean 'uparea' (Ac.*(Ac.cellsize^2)) @mean 'gradient' GcSAc @mean});

		% Calculate normalized channel steepness using reference concavity supplied by user
		KSNc=GcREF./(Ac.*(Ac.cellsize^2)).^(-ref_concavity);
		MSNc=STREAMobj2mapstruct(Sc,'seglength',smooth_distance,'attributes',...
			{'ksn' KSNc @mean 'uparea' (Ac.*(Ac.cellsize^2)) @mean 'gradient' GcREF @mean});

		% Calculate basin wide ksn statistics
		min_ksn=min([MSNc.ksn]);
		mean_ksn=mean([MSNc.ksn]);
		max_ksn=max([MSNc.ksn]);
		std_ksn=std([MSNc.ksn]);
		se_ksn=std_ksn/sqrt(numel(MSNc)); % Standard error

		% Calculate basin wide gradient statistics
		min_grad=min(nanmin(Goc.Z));
		mean_grad=mean(nanmean(Goc.Z));
		max_grad=max(nanmax(Goc.Z));
		std_grad=nanstd(vertcat(Goc.Z(:)));
		se_grad=std_grad/sqrt(sum(sum(~isnan(Goc.Z)))); % Standard error

		KSNc_stats=[mean_ksn se_ksn std_ksn min_ksn max_ksn];
		Gc_stats=[mean_grad se_grad std_grad min_grad max_grad];

		FileName=['Basin_' num2str(basin_num) '_Data.mat'];
		disp('	Writing MAT files')
		save(FileName,'RiverMouth','DEMc','DEMoc','drainage_area','FDc','Ac','Sc','SLc','Chic','SAc','Goc','GcREF','KSc','MSc','KSNc','MSNc','KSNc_stats','Gc_stats','Centroid','ChiOBJc');

		switch write_arc_files
		case true
			disp('	Writing GeoTIFF and ShapeFiles')

			% Replace NaNs in DEM with -32768
			Didx=isnan(DEMoc.Z);
			DEMoc_temp=DEMoc;
			DEMoc_temp.Z(Didx)=-32768;

			DEMFileName=['Basin_' num2str(basin_num) '_DEM.tif'];
			GRIDobj2geotiff(DEMoc_temp,DEMFileName);
			CHIFileName=['Basin_' num2str(basin_num) '_CHI.txt'];
			GRIDobj2ascii(ChiOBJc,CHIFileName);
			KSNFileName=['Basin_' num2str(basin_num) '_KSN.shp'];
			shapewrite(MSNc,KSNFileName);
		case false
			;
		end


	end
end

function [ChiOBJ]=MakeChiMapInd(ChiOI,DEM);
% Function to create a chi map from a precalculated chi structure (result of chiplot) and
% the DEM GRIDobj associated with that chi structure

	% Grab x, y, and chi
	cx=ChiOI.x;
	cy=ChiOI.y;
	cc=ChiOI.chi;

	% Strip out NaNs
	idx=~isnan(cx);
	cx=cx(idx); cy=cy(idx); cc=cc(idx);

	cix=coord2ind(DEM,cx,cy);

	% Make Empty GRIDobj
	ChiOBJ=GRIDobj(DEM);
	idx=ChiOBJ.Z==0;
	ChiOBJ.Z(idx)=-32768;

	% Populate Grid
	ChiOBJ.Z(cix)=cc;

end

function [G,GminIX] = CompGrad(DEM,FD,A,min_ksn,ref_concavity)

	% Convert Area in Pixels to Drainage Area in Sq Meters
	DA=A.*(A.cellsize^2);

	% Find Max Drainage Area Exponent
	maxDA=max(max(DA));
	ex=ceil(log10(maxDA));

	% Build Range
	bins=logspace(4,ex,ex-4+1); % Min area is 1e4
	num_bins=numel(bins);

	% Build Masks and Gradients
	for ii=1:num_bins
		da=bins(ii);
		% Build Indices
		if ii==1
			IX=DA<=da;
		else
			IX=DA>bins(ii-1) & DA<=da;
		end

		% Calculate and Impose Min Gradient
		minG=min_ksn*(da^-ref_concavity);
		Gtemp=gradient8(imposemin(FD,DEM,minG));
		GIX=Gtemp.*IX;

		% Layer Composite
		if ii==1
			G=GIX;
		else
			G=G+GIX;
		end
	end

	% Build map of where minimum gradient has been imposed
	Gn=gradient8(DEM);
	mIX=Gn.Z~=G.Z;
	[X,Y]=getcoordinates(DEM);
	GminIX=GRIDobj(X,Y,double(mIX));
end
