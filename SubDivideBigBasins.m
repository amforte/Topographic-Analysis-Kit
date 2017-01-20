function SubDivideBigBasins(location_of_data_files,max_basin_size,divide_method,varargin)
	% Function takes outputs from 'ProcessRiverBasins' function and subdvides any basin with a drainage area above a specified size and
	% outputs clipped dem, stream network, variout topographic metrics, and river values (ks, ksn, chi)
	%
	% Required Inputs:
	% 		location_of_data_files - full path of folder which contains the mat files from 'ProcessRiverBasins'
	% 		max_basin_size - size above which drainage basins will be subdivided in square kilometers
	%		divide_method - method for subdividing basins, either 'order' to use the outlets of streams of a given order (that the user can specify with the optional
	%		's_order' parameter), 'confluences' to use the locations of confluences, or 'up_confluences' to use locations just upstream of confluences
	%
	% Optional Inputs:
	% 		threshold_area [1e6] - minimum accumulation area to define streams in meters squared
	% 		segment_length [1000] - smoothing distance in meters for averaging along ksn, suggested value is 1000 meters
	% 		theta_ref [0.45] - reference concavity for calculating ksn, suggested value is 0.45
	% 		write_arc_files [false] - set value to true to output a geotiff of the DEM and a shapefile of the ksn, false to not output arc files
	%		s_order [3] - stream order for defining stream outlets for subdividing if 'divide_method' is 'order'
	%		DEM_original [] - original DEM input to 'ProcessRiverBasins' (required if clip_method was set to 'segment' in 'ProcessRiverBasins')
	%		FD_original [] - original FD input to 'ProcessRiverBasins' (required if clip_method was set to 'segment' in 'ProcessRiverBasins')
	%		S_original [] - original S input to 'ProcessRiverBasins' (required if clip_method was set to 'segment' in 'ProcessRiverBasins')
	%
	% Examples:
	%		SubdivideBigBasins('/Users/JoeBlow/Project',100,'confluences');
	%		SubdivideBigBasins('/Users/JoeBlow/Project',100,'order','s_order',2,'threshold_area',1e5,'write_arc_files',true);
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Last Revised Winter 2017 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'SubDivideBigBasins';
	addRequired(p,'location_of_data_files',@(x) ischar(x));
	addRequired(p,'max_basin_size',@(x) isnumeric(x));
	addRequired(p,'divide_method',@(x) ischar(validatestring(x,{'order','confluences','up_confluences'})));

	addParamValue(p,'theta_ref',0.45,@(x) isscalar(x) && isnumeric(x));
	addParamValue(p,'threshold_area',1e6,@(x) isscalar(x) && isnumeric(x));
	addParamValue(p,'segment_length',1000,@(x) isscalar(x) && isnumeric(x));
	addParamValue(p,'write_arc_files',false,@(x) isscalar(x));
	addParamValue(p,'s_order',[3],@(x) isscalar(x));
	addParamValue(p,'DEM_Original',[],@(x) isa(x,'GRIDobj'));
	addParamValue(p,'FD_Original',[],@(x) isa(x,'FLOWobj'));
	addParamValue(p,'S_Original',[],@(x) isa(x,'STREAMobj'));

	parse(p,location_of_data_files,max_basin_size,divide_method,varargin{:});
	location_of_data_files=p.Results.location_of_data_files;
	max_basin_size=p.Results.max_basin_size;
	divide_method=p.Results.divide_method;

	theta_ref=p.Results.theta_ref;
	threshold_area=p.Results.threshold_area;
	segment_length=p.Results.segment_length;
	write_arc_files=p.Results.write_arc_files;
	s_order=p.Results.s_order;
	DEM0=p.Results.DEM_Original;
	FD0=p.Results.FD_Original;
	S0=p.Results.S_Original;

	current=pwd;
	cd(location_of_data_files);

	FileList=dir('*Data.mat');
	num_files=numel(FileList);


	% Load first file to grab clip method
	load(FileList(1,1).name,'clip_method');

	% If clip method is 'segment' perform checks and some pre processing
	if strcmp(clip_method,'segment')

		if isempty(DEM0) | isempty(FD0) | isempty(S0)
			error('Clip method used in ProcessRiverBasins was "segment" and as a result this code requires inputs for optional paramters "DEM0", "FD0", and "S0" ');
		end

		disp('Generating hydrologically conditioned DEM for the entire region')
		A0=flowacc(FD0);
		zcon=mincosthydrocon(S0,DEM0,'interp',0.1);
		DEMcon=GRIDobj(DEM0);
		DEMcon.Z(DEMcon.Z==0)=NaN;
		DEMcon.Z(S0.IXgrid)=zcon;
	end

	% Begin main file loop
	w1=waitbar(0,'Subdividing basins');
	for ii=1:num_files;
		FileName=FileList(ii,1).name;

		% Load in required basin files depending on clip method and rename
		if strcmp(clip_method,'clip')
			load(FileName,'RiverMouth','drainage_area','DEMoc','FDc','Sc','ksn_method');
			DEM=DEMoc;
			S=Sc;
			FD=FDc;
			RM=RiverMouth;
			DA=drainage_area;
			basin_num=RM(:,3);
		elseif strcmp(clip_method,'segment')
			load(FileName,'RiverMouth','drainage_area','DEMoc','Sc','ksn_method');
			DEM=DEMoc;
			S=Sc;
			RM=RiverMouth;
			DA=drainage_area;
			basin_num=RM(:,3);
		end

		% Check drainage area to determine if this basin will be processed
		if DA>=max_basin_size
			waitbar(ii/num_files,w1,['Subdividing basin number ' num2str(basin_num)]);

			switch divide_method
			case 'order'
				% Extract 'outlets' of 3rd order streams
				Se=modify(S,'streamorder',s_order);
				outs=streampoi(Se,'outlets','xy');
				x=outs(:,1);
				y=outs(:,2);
				num_new_basins=numel(x);
			case 'confluences'
				cons=streampoi(S,'confluences','xy');
				x=cons(:,1);
				y=cons(:,2);
				num_new_basins=numel(x);
			case 'up_confluences'
				cons=streampoi(S,'bconfluences','xy');
				x=cons(:,1);
				y=cons(:,2);
				num_new_basins=numel(x);
			end

			switch clip_method
			case 'clip'
				for jj=1:num_new_basins
					xx=x(jj);
					yy=y(jj);
					basin_string=sprintf([num2str(basin_num) '%03d'],jj);
					RiverMouth=[xx yy str2num(basin_string)];

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
					Sc=STREAMobj(FDc,'minarea',threshold_area,'unit','mapunits');

					% Calculate chi and create chi map
					Cc=chitransform(Sc,Ac.*(Ac.cellsize^2),'a0',1,'mn',theta_ref);
					ChiOBJc=GRIDobj(DEMoc);
					ChiOBJc.Z(Sc.IXgrid)=Cc;

					% Calculate gradient
					Goc=gradient8(DEMoc);

					% Find best fit concavity
					SLc=klargestconncomps(Sc,1);
					zcon=mincosthydrocon(SLc,DEMoc,'interp',0.1);
					DEMcc=GRIDobj(DEMoc);
					DEMcc.Z(DEMcc.Z==0)=NaN;
					DEMcc.Z(SLc.IXgrid)=zcon;
					Chic=chiplot(SLc,DEMcc,Ac,'a0',1,'plot',false);

					% Calculate ksn
					switch ksn_method
					case 'quick'
						[MSc]=KSN_Quick(DEMoc,Ac,Sc,Chic.mn,segment_length);
						[MSNc]=KSN_Quick(DEMoc,Ac,Sc,theta_ref,segment_length);
					case 'trib'
						% Overide choice if very small basin as KSN_Trib will fail for small basins
						if drainage_area>2.5
							[MSc]=KSN_Trib(DEMoc,FDc,Ac,Sc,Chic.mn,segment_length);
							[MSNc]=KSN_Trib(DEMoc,FDc,Ac,Sc,theta_ref,segment_length);
						else
							[MSc]=KSN_Quick(DEMoc,Ac,Sc,Chic.mn,segment_length);
							[MSNc]=KSN_Quick(DEMoc,Ac,Sc,theta_ref,segment_length);
						end
					end

					% Calculate basin wide ksn statistics
					min_ksn=min([MSNc.ksn]);
					mean_ksn=mean([MSNc.ksn]);
					max_ksn=max([MSNc.ksn]);
					std_ksn=std([MSNc.ksn]);
					se_ksn=std_ksn/sqrt(numel(MSNc)); % Standard error

					% Calculate basin wide gradient statistics
					min_grad=nanmin(Goc.Z(:));
					mean_grad=nanmean(Goc.Z(:));
					max_grad=nanmax(Goc.Z(:));
					std_grad=nanstd(Goc.Z(:));
					se_grad=std_grad/sqrt(sum(~isnan(Goc.Z(:)))); % Standard error

					KSNc_stats=[mean_ksn se_ksn std_ksn min_ksn max_ksn];
					Gc_stats=double([mean_grad se_grad std_grad min_grad max_grad]);

					SubFileName=['Basin_' num2str(basin_num) '_DataSubset_' num2str(jj) '.mat'];

					VarList=whos('-file',FileName);
					VarInd=find(strcmp(cellstr(char(VarList.name)),'AGc'));

					if ~isempty(VarInd)
						load(FileName,'AGc');
						AG=AGc;
						num_grids=size(AG,1);
						AGc=cell(size(AG));
						for kk=1:num_grids
							AGcOI=crop(AG{kk,1},I,nan);
							AGc{kk,1}=AGcOI;
							AGc{kk,2}=AG{kk,2};
							mean_AGc=nanmean(AGcOI.Z(:));
							min_AGc=nanmin(AGcOI.Z(:));
							max_AGc=nanmax(AGcOI.Z(:));
							std_AGc=nanstd(AGcOI.Z(:));
							se_AGc=std_AGc/sqrt(sum(~isnan(AGcOI.Z(:))));
							AGc_stats(kk,:)=[mean_AGc se_AGc std_AGc min_AGc max_AGc];
						end
						save(SubFileName,'RiverMouth','DEMcc','DEMoc','drainage_area','FDc','Ac','Sc','SLc','Chic','Goc','MSc','MSNc','KSNc_stats','Gc_stats','Centroid','ChiOBJc','ksn_method','clip_method','AGc','AGc_stats');
					elseif isempty(VarInd)
						save(SubFileName,'RiverMouth','DEMcc','DEMoc','drainage_area','FDc','Ac','Sc','SLc','Chic','Goc','MSc','MSNc','KSNc_stats','Gc_stats','Centroid','ChiOBJc','ksn_method','clip_method');
					end

					switch write_arc_files
					case true
						% Replace NaNs in DEM with -32768
						Didx=isnan(DEMoc.Z);
						DEMoc_temp=DEMoc;
						DEMoc_temp.Z(Didx)=-32768;

						DEMFileName=['Basin_' num2str(basin_num) '_DataSubset_' num2str(jj) '_DEM.tif'];
						GRIDobj2geotiff(DEMoc_temp,DEMFileName);
						CHIFileName=['Basin_' num2str(basin_num) '_DataSubset_' num2str(jj) '_CHI.txt'];
						GRIDobj2ascii(ChiOBJc,CHIFileName);
						KSNFileName=['Basin_' num2str(basin_num) '_DataSubset_' num2str(jj) '_KSN.shp'];
						shapewrite(MSNc,KSNFileName);
					end

				end % New basin loop end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%				
			case 'segment'
				for jj=1:num_new_basins
					xx=x(jj);
					yy=y(jj);
					basin_string=sprintf([num2str(basin_num) '%03d'],jj);
					RiverMouth=[xx yy str2num(basin_string)];

					% Build dependenc map and clip out drainage basins
					I=dependencemap(FD0,xx,yy);
					DEMoc=crop(DEM0,I,nan);

					% Calculate drainage area
					dep_map=GRIDobj2mat(I);
					num_pix=sum(sum(dep_map));
					drainage_area=(num_pix*DEMoc.cellsize*DEMoc.cellsize)/(1e6);

					% Find weighted centroid of drainage basin
					[Cx,Cy]=FindCentroid(DEMoc);
					Centroid=[Cx Cy];

					% Generate new stream map by segmenting original full stream network
					six=coord2ind(DEM0,xx,yy);
					Sc=STREAMobj(FD0,'minarea',threshold_area,'unit','mapunits','outlets',six);

					% Calculate chi and create chi map
					Cc=chitransform(Sc,A0.*(A0.cellsize^2),'a0',1,'mn',theta_ref);
					ChiOBJc=GRIDobj(DEMoc);
					ScIX=coord2ind(DEMoc,Sc.x,Sc.y);
					ChiOBJc.Z(ScIX)=Cc;

					% Calculate slope area
					SAc=slopearea(Sc,DEM0,A0,'plot',false);
					Goc=gradient8(DEMoc);

					% Calculate ksn
					switch ksn_method
					case 'quick'
						[MSc]=KSN_Quick(DEM0,A0,Sc,-1*(SAc.theta),segment_length);
						[MSNc]=KSN_Quick(DEM0,A0,Sc,theta_ref,segment_length);
					case 'trib'
						% Overide choice if very small basin as KSN_Trib will fail for small basins
						if drainage_area>2.5
							[MSc]=KSN_Trib(DEM0,FD0,A,Sc,-1*(SAc.theta),segment_length);
							[MSNc]=KSN_Trib(DEM0,FD0,A,Sc,theta_ref,segment_length);
						else
							[MSc]=KSN_Quick(DEM0,A0,Sc,-1*(SAc.theta),segment_length);
							[MSNc]=KSN_Quick(DEM0,A0,Sc,theta_ref,segment_length);				
						end
					end

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
					Gc_stats=double([mean_grad se_grad std_grad min_grad max_grad]);

					SubFileName=['Basin_' num2str(basin_num) '_DataSubset_' num2str(jj) '.mat'];

					VarList=whos('-file',FileName);
					VarInd=find(strcmp(cellstr(char(VarList.name)),'AGc'));

					if ~isempty(VarInd)
						load(FileName,'AGc');
						AG=AGc;
						num_grids=size(AG,1);
						AGc=cell(size(AG));
						for kk=1:num_grids
							AGcOI=crop(AG{kk,1},I,nan);
							AGc{kk,1}=AGcOI;
							AGc{kk,2}=AG{kk,2};
							mean_AGc=nanmean(AGcOI.Z(:));
							min_AGc=nanmin(AGcOI.Z(:));
							max_AGc=nanmax(AGcOI.Z(:));
							std_AGc=nanstd(AGcOI.Z(:));
							se_AGc=std_AGc/sqrt(sum(~isnan(AGcOI.Z(:))));
							AGc_stats(kk,:)=[mean_AGc se_AGc std_AGc min_AGc max_AGc];
						end
						save(SubFileName,'RiverMouth','DEMoc','drainage_area','Sc','SAc','Goc','MSc','MSNc','KSNc_stats','Gc_stats','Centroid','ChiOBJc','ksn_method','clip_method','AGc','AGc_stats');
					elseif isempty(VarInd)
						save(SubFileName,'RiverMouth','DEMoc','drainage_area','Sc','SAc','Goc','MSc','MSNc','KSNc_stats','Gc_stats','Centroid','ChiOBJc','ksn_method','clip_method');
					end

					switch write_arc_files
					case true
						% Replace NaNs in DEM with -32768
						Didx=isnan(DEMoc.Z);
						DEMoc_temp=DEMoc;
						DEMoc_temp.Z(Didx)=-32768;

						DEMFileName=['Basin_' num2str(basin_num) '_DataSubset_' num2str(jj) '_DEM.tif'];
						GRIDobj2geotiff(DEMoc_temp,DEMFileName);
						CHIFileName=['Basin_' num2str(basin_num) '_DataSubset_' num2str(jj) '_CHI.txt'];
						GRIDobj2ascii(ChiOBJc,CHIFileName);
						KSNFileName=['Basin_' num2str(basin_num) '_DataSubset_' num2str(jj) '_KSN.shp'];
						shapewrite(MSNc,KSNFileName);
					end
				end % New basin loop end

			end % Clip Method end
		end % Drainage Area check end
	end % Main Loop end
	close(w1);
	cd(current);
end % Main Function End



function [ksn_ms]=KSN_Quick(DEM,A,S,theta_ref,segment_length)

	zc=mincosthydrocon(S,DEM,'interp',0.1);
	DEMc=GRIDobj(DEM);
	DEMc.Z(DEMc.Z==0)=NaN;
	DEMc.Z(S.IXgrid)=zc;
	G=gradient8(DEMc);
	Z_RES=DEMc-DEM;

	ksn=G./(A.*(A.cellsize^2)).^(-theta_ref);
	
	ksn_ms=STREAMobj2mapstruct(S,'seglength',segment_length,'attributes',...
		{'ksn' ksn @mean 'uparea' (A.*(A.cellsize^2)) @mean 'gradient' G @mean 'cut_fill' Z_RES @mean});
end

function [ksn_ms]=KSN_Trib(DEM,FD,A,S,theta_ref,segment_length)

	% Define non-intersecting segments
	[as]=networksegment_slim(DEM,FD,S);
	seg_bnd_ix=as.ix;
	% Precompute values or extract values needed for later
	z=mincosthydrocon(S,DEM,'interp',0.1);
	zu=getnal(S,DEM);
	z_res=z-zu;
	c=chitransform(S,A.*(A.cellsize^2),'a0',1,'mn',theta_ref);
	d=S.distance;
	da=getnal(S,A.*(A.cellsize^2));
	ixgrid=S.IXgrid;
	% Extract ordered list of stream indices and find breaks between streams
	s_node_list=S.orderednanlist;
	streams_ix=find(isnan(s_node_list));
	streams_ix=vertcat(1,streams_ix);
	% Generate empty node attribute list for ksn values
	ksn_nal=zeros(size(d));
	% Begin main loop through channels
	num_streams=numel(streams_ix)-1;
	seg_count=1;
	for ii=1:num_streams
		% Extract node list for stream of interest
		if ii==1
			snlOI=s_node_list(streams_ix(ii):streams_ix(ii+1)-1);
		else
			snlOI=s_node_list(streams_ix(ii)+1:streams_ix(ii+1)-1);
		end

		% Determine which segments are within this stream
		[~,~,dn]=intersect(snlOI,seg_bnd_ix(:,1));
		[~,~,up]=intersect(snlOI,seg_bnd_ix(:,2));
		seg_ix=intersect(up,dn);

		num_segs=numel(seg_ix);
		dn_up=seg_bnd_ix(seg_ix,:);
		for jj=1:num_segs
			% Find positions within node list
			dnix=find(snlOI==dn_up(jj,1));
			upix=find(snlOI==dn_up(jj,2));
			% Extract segment indices of desired segment
			seg_ix_oi=snlOI(upix:dnix);
			% Extract flow distances and normalize
			dOI=d(seg_ix_oi);
			dnOI=dOI-min(dOI);
			num_bins=ceil(max(dnOI)/segment_length);
			bin_edges=[0:segment_length:num_bins*segment_length];
			% Loop through bins
			for kk=1:num_bins
				idx=dnOI>bin_edges(kk) & dnOI<=bin_edges(kk+1);
				bin_ix=seg_ix_oi(idx);
				cOI=c(bin_ix);
				zOI=z(bin_ix);
					if numel(cOI)>2
						[ksn_val]=Chi_Z_Spline(cOI,zOI);
						ksn_nal(bin_ix)=ksn_val;

						% Build mapstructure
						ksn_ms(seg_count).Geometry='Line';
						ksn_ms(seg_count).X=S.x(bin_ix);
						ksn_ms(seg_count).Y=S.y(bin_ix);
						ksn_ms(seg_count).ksn=ksn_val;
						ksn_ms(seg_count).cut_fill=mean(z_res(bin_ix));
						ksn_ms(seg_count).area=mean(da(bin_ix));
						seg_count=seg_count+1;
					end
			end
		end	
	end

end

function seg = networksegment_slim(DEM,FD,S)
	% Slimmed down version of 'networksegment' from main TopoToolbox library that also removes zero and single node length segments

	%% Identify channel heads, confluences, b-confluences and outlets
	Vhead = streampoi(S,'channelheads','logical');  ihead=find(Vhead==1);  IXhead=S.IXgrid(ihead);
	Vconf = streampoi(S,'confluences','logical');   iconf=find(Vconf==1);  IXconf=S.IXgrid(iconf);
	Vout = streampoi(S,'outlets','logical');        iout=find(Vout==1);    IXout=S.IXgrid(iout);
	Vbconf = streampoi(S,'bconfluences','logical'); ibconf=find(Vbconf==1);IXbconf=S.IXgrid(ibconf);

	%% Identify basins associated to b-confluences and outlets
	DB   = drainagebasins(FD,vertcat(IXbconf,IXout));DBhead=DB.Z(IXhead); DBbconf=DB.Z(IXbconf); DBconf=DB.Z(IXconf); DBout=DB.Z(IXout);

	%% Compute flowdistance
	D = flowdistance(FD);

	%% Identify river segments
	% links between channel heads and b-confluences
	[~,ind11,ind12]=intersect(DBbconf,DBhead);
	% links between confluences and b-confluences
	[~,ind21,ind22]=intersect(DBbconf,DBconf);
	% links between channel heads and outlets
	[~,ind31,ind32]=intersect(DBout,DBhead);
	% links between channel heads and outlets
	[~,ind41,ind42]=intersect(DBout,DBconf);
	% Connecting links into segments
	IX(:,1) = [ IXbconf(ind11)' IXbconf(ind21)' IXout(ind31)'  IXout(ind41)'  ];   ix(:,1)= [ ibconf(ind11)' ibconf(ind21)' iout(ind31)'  iout(ind41)'  ];
	IX(:,2) = [ IXhead(ind12)'  IXconf(ind22)'  IXhead(ind32)' IXconf(ind42)' ];   ix(:,2)= [ ihead(ind12)'  iconf(ind22)'  ihead(ind32)' iconf(ind42)' ];

	% Compute segment flow length
	flength=double(abs(D.Z(IX(:,1))-D.Z(IX(:,2))));

	% Remove zero and one node length elements
	idx=flength>=2*DEM.cellsize;
	seg.IX=IX(idx,:);
	seg.ix=ix(idx,:);
	seg.flength=flength(idx);

	% Number of segments
	seg.n=numel(IX(:,1));
end

function [KSN] = Chi_Z_Spline(c,z)

	% Resample chi-elevation relationship using cubic spline interpolation
	[~,minIX]=min(c);
	zb=z(minIX);
	chiF=c-min(c);
	zabsF=z-min(z);
	chiS=linspace(0,max(chiF),numel(chiF)).';
	zS=spline(chiF,zabsF,chiS);

	%Calculate beta
    BETA = chiS\(zS);

	KSN= BETA; %Beta.*a0^mn - if a0 set to 1, not needed
end


