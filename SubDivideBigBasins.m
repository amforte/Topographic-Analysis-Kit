function SubDivideBigBasins(location_of_data_files,max_basin_size,varargin)
	% Function takes outputs from 'ProcessRiverBasins' function and subdvides any basin with a drainage area above a specified size and
	% outputs clipped dem, stream network, variout topographic metrics, and river values (ks, ksn, chi)
	%
	% Required Inputs:
	% 		location_of_data_files - full path of folder which contains the mat files from 'ProcessRiverBasins'
	% 		max_basin_size - size above which drainage basins will be subdivided in square kilometers
	%
	% Optional Inputs:
	% 		threshold_area [1e6] - minimum accumulation area to define streams in meters squared
	% 		smooth_distance [1000] - smoothing distance in meters for averaging along ksn, suggested value is 1000 meters
	% 		ref_concavity [0.45] - reference concavity for calculating ksn, suggested value is 0.45
	% 		write_arc_files [false] - set value to true to output a geotiff of the DEM and a shapefile of the ksn, false to not output arc files
	%
	% Examples:
	%		SubdivideBigBasins('/Users/JoeBlow/Project',100);
	%		SubdivideBigBasins('/Users/JoeBlow/Project',100,'threshold_area',1e5,'write_arc_files',true);
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Last Revised Fall 2015 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'SubDivideBigBasins';
	addRequired(p,'location_of_data_files',@(x) ischar(x));
	addRequired(p,'max_basin_size',@(x) isnumeric(x));

	addParamValue(p,'ref_concavity',0.45,@(x) isscalar(x) && isnumeric(x));
	addParamValue(p,'threshold_area',1e6,@(x) isscalar(x) && isnumeric(x));
	addParamValue(p,'smooth_distance',1000,@(x) isscalar(x) && isnumeric(x));
	addParamValue(p,'write_arc_files',false,@(x) isscalar(x));

	parse(p,DEM,FD,S,river_mouths,varargin{:});
	location_of_data_files=p.Results.location_of_data_files;
	max_basin_size=p.Results.max_basin_size;

	ref_concavity=p.Results.ref_concavity;
	threshold_area=p.Results.threshold_area;
	smooth_distance=p.Results.smooth_distance;
	write_arc_files=p.Results.write_arc_files;

	current=pwd;
	cd(location_of_data_files);

	FileList=dir('*Data.mat');
	num_files=numel(FileList);

	for ii=1:num_files;
		FileName=FileList(ii,1).name;
		load(FileName,'RiverMouth','drainage_area','DEMoc','Sc','FDc');

		DEM=DEMoc;
		S=Sc;
		FD=FDc;
		RM=RiverMouth;
		DA=drainage_area;
		basin_num=RM(:,3);

		% Check drainage area
		if DA>=max_basin_size
			disp(['Subdividing basin number ' num2str(basin_num)])
			% Extract 'outlets' of 3rd order streams
			Se=modify(S,'streamorder',3);
			outs=streampoi(Se,'outlets','xy');
			x=outs(:,1);
			y=outs(:,2);
			num_new_basins=numel(x);

			for jj=1:num_new_basins
				xx=x(jj);
				yy=y(jj);
				basin_string=sprintf([num2str(basin_num) '%03d'],jj);
				RiverMouth=[xx yy str2num(basin_string)];

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
				DEMc=imposemin(FDc,DEMoc,.001);
				DEMc_res=DEMc.cellsize;
				min_area=floor(threshold_area/(DEMc_res*DEMc_res));
				isstream=Ac>min_area;
				Sc=STREAMobj(FDc,isstream);

				%Extract largest drainage
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

				FileName=['Basin_' num2str(basin_num) '_DataSubset_' num2str(jj) '.mat'];
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
	end

	cd(current);
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
