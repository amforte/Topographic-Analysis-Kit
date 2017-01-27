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
	% 		segment_length [1000] - smoothing distance in meters for averaging along ksn, suggested value is 1000 meters
	% 		theta_ref [0.45] - reference concavity for calculating ksn, suggested value is 0.45
	%		ksn_method [quick] - switch between method to calculate ksn values, options are 'quick' and 'trib', the 'trib' method takes 3-4 times longer 
	%			than the 'quick' method. In most cases, the 'quick' method works well, but if values near tributary junctions are important, then 'trib'
	%			may be better as this calculates ksn values for individual channel segments individually
	% 		write_arc_files [false] - set value to true to output a geotiff of the DEM and a shapefile of the ksn, false to not output arc files
	%		clip_method ['clip'] - flag to determine how the code clips out stream networks, expects either 'clip' or 'segment'. The 'clip' option 
	%			(default) will clip out DEMs and rerun the flow algorithm on this clipped DEM and proceed from there to get a new stream network. 
	%			The 'segment' option clips out the DEM the same but then segments the original, full stream network based on the input river mouth,
	%			i.e. a new flow routing algorithm is not run on the clipped DEM. This 'segment' option should be tried if you are encountering errors
	%			with the 'clip' method, specifically warnings that STREAMobjs contain 0 connected components. This error often results if selected 
	%			watersheds are small and the DEMs are noisy, flow routing on the clipped DEMs in this case sometimes will route streams differently, 
	%			e.g. routing streams out the side of the basin. It is strongly recommended that you use the default 'clip' method as opposed to the 'segment'
	%			method, only using the 'segment' method if the clip method fails.
	%		add_grids [] - option to provide a cell array of additional grids to clip by selected river basins. The expected input is a m x 2 cell array,
	%			where the first column is a GRIDobj and the second column is a string identifying what this grid is (so you can remember what these grids
	%			are when looking at outputs later, but also used as the name of field values if you use 'Basin2Shape' on the output basins so these should be short 
	%			strings with no spaces). The code will perform a check on any input grid to determine if it is the same dimensions and cellsize as the input DEM, if
	%			it is not it will use the function 'resample' to transform the input grid. You can control the resampling method used with the 'resample_method' optional
	%			parameter (see below), but this method will be applied to all grids you provide, so if you want to use different resampling methods for different grids
	%			it is recommnended that you use the 'resample' function on the additional grids before you supply them to this function.
	%		resample_method ['nearest'] - method to use in the resample function on additional grids (if required). Acceptable inputs are 'nearest', 'bilinear', 
	%			or 'bicubic'. Method 'nearest' is appropriate if you do not want the resampling to interpolate between values (e.g. if an additinal grid has specific values
	%			that correlate to a property like rock type) and either 'bilinear' or 'bicubic' is appropriate if you want smooth variations between nodes. 
	%			
	%
	% Examples:
	%		ProcessRiverBasins(DEM,FD,S,RiverMouths);
	%		ProcessRiverBasins(DEM,FD,S,RiverMouths,'theta_ref',0.5,'write_arc_files',true);
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Last Revised Winter 2017 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'ProcessRiverBasins';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'FD',@(x) isa(x,'FLOWobj'));
	addRequired(p,'S',@(x) isa(x,'STREAMobj'));
	addRequired(p,'river_mouths',@(x) isnumeric(x) && size(x,2)==3);

	addParamValue(p,'theta_ref',0.45,@(x) isscalar(x) && isnumeric(x));
	addParamValue(p,'threshold_area',1e6,@(x) isscalar(x) && isnumeric(x));
	addParamValue(p,'segment_length',1000,@(x) isscalar(x) && isnumeric(x));
	addParamValue(p,'write_arc_files',false,@(x) isscalar(x));
	addParamValue(p,'clip_method','clip',@(x) ischar(validatestring(x,{'clip','segment'})));
	addParamValue(p,'ksn_method','quick',@(x) ischar(validatestring(x,{'quick','trib'})));
	addParamValue(p,'add_grids',[],@(x) isa(x,'cell') && size(x,2)==2);
	addParamValue(p,'resample_method','nearest',ischar(validatestring(x,{'nearest','bilinear','bicubic'})));

	parse(p,DEM,FD,S,river_mouths,varargin{:});
	DEM=p.Results.DEM;
	FD=p.Results.FD;
	S=p.Results.S;
	river_mouths=p.Results.river_mouths;

	theta_ref=p.Results.theta_ref;
	threshold_area=p.Results.threshold_area;
	segment_length=p.Results.segment_length;
	write_arc_files=p.Results.write_arc_files;
	clip_method=p.Results.clip_method;
	ksn_method=p.Results.ksn_method;
	AG=p.Results.add_grids;
	resample_method=p.Results.resample_method;

	% Clippin gout NaNs, max elev set arbitrarily large but below the 32,768 internal NaN value.
	max_elev=10000;
	min_elev=-200;
	IDX=DEM<max_elev & DEM>min_elev;
	DEM=crop(DEM,IDX,nan);

	% Perform check on dimensions and cellsize of additional grids and resample if necessary
	if ~isempty(AG)
		num_grids=size(AG,1);
		for jj=1:num_grids
			AGoi=AG{jj,1};
			if ~validatealignment(AGoi,DEM);
				disp(['Resampling ' AG{jj,2} ' GRIDobj to be the same resolution and dimensions as the input DEM by the ' resample_method ' method']);
				AG{jj,1}=resample(Agoi,DEM,resample_method);
			end
		end
	end

	disp('Snapping river mouths to stream network')
	xi=river_mouths(:,1);
	yi=river_mouths(:,2);
	riv_nums=river_mouths(:,3);

	num_basins=numel(xi);

	[xn,yn]=snap2stream(S,xi,yi);

	RM=[xn yn riv_nums];

	num_basins=numel(xn);

	switch clip_method
	case 'clip'
		w1=waitbar(0,['Working on Basin Number 0 of ' num2str(num_basins) ' total basins']);
		for ii=1:num_basins
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
			Sc=STREAMobj(FDc,'minarea',threshold_area,'unit','mapunits');

			% Check to make sure the stream object isn't empty
			if isempty(Sc.x)
				warning(['Input threshold drainage area is too large for basin ' num2str(basin_num) ' decreasing threshold area for this basin']);
				new_thresh=threshold_area;
				while isempty(Sc.x)
					new_thresh=new_thresh/2;
					Sc=STREAMobj(FDc,'minarea',new_thresh,'unit','mapunits');
				end
			end

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

			% Calculate basin wide elevation statistics
			min_z=nanmin(DEMoc.Z(:));
			mean_z=nanmean(DEMoc.Z(:));
			max_z=nanmax(DEMoc.Z(:));
			std_z=nanstd(DEMoc.Z(:));
			se_z=std_z/sqrt(sum(~isnan(DEMoc.Z(:)))); % Standard error

			KSNc_stats=[mean_ksn se_ksn std_ksn min_ksn max_ksn];
			Gc_stats=double([mean_grad se_grad std_grad min_grad max_grad]);
			Zc_stats=double([mean_z se_z std_z min_z max_z]);

			% Find outlet elevation
			out_ix=coord2ind(DEMoc,xx,yy);
			out_el=double(DEMoc.Z(out_ix));

			FileName=['Basin_' num2str(basin_num) '_Data.mat'];

			if ~isempty(AG)
				num_grids=size(AG,1);
				AGc=cell(size(AG));
				for jj=1:num_grids
					AGcOI=crop(AG{jj,1},I,nan);
					AGc{jj,1}=AGcOI;
					AGc{jj,2}=AG{jj,2};
					mean_AGc=nanmean(AGcOI.Z(:));
					min_AGc=nanmin(AGcOI.Z(:));
					max_AGc=nanmax(AGcOI.Z(:));
					std_AGc=nanstd(AGcOI.Z(:));
					se_AGc=std_AGc/sqrt(sum(~isnan(AGcOI.Z(:))));
					AGc_stats(jj,:)=[mean_AGc se_AGc std_AGc min_AGc max_AGc];
				end
				save(FileName,'RiverMouth','DEMcc','DEMoc','out_el','drainage_area','FDc','Ac','Sc','SLc','Chic','Goc','MSc','MSNc','KSNc_stats','Gc_stats','Zc_stats','Centroid','ChiOBJc','ksn_method','clip_method','AGc','AGc_stats');
			elseif isempty(AG)
				save(FileName,'RiverMouth','DEMcc','DEMoc','out_el','drainage_area','FDc','Ac','Sc','SLc','Chic','Goc','MSc','MSNc','KSNc_stats','Gc_stats','Zc_stats','Centroid','ChiOBJc','ksn_method','clip_method');
			end

			switch write_arc_files
			case true

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
			waitbar(ii/num_basins,w1,['Working on Basin Number ' num2str(ii) ' of ' num2str(num_basins) ' total basins'])
		end
		close(w1)

	case 'segment'

		% Generate flow accumulation and hydrologically conditioned DEM
		disp('Generating hydrologically conditioned DEM for the entire region')
		A=flowacc(FD);
		zcon=mincosthydrocon(S,DEM,'interp',0.1);
		DEMcon=GRIDobj(DEM);
		DEMcon.Z(DEMcon.Z==0)=NaN;
		DEMcon.Z(S.IXgrid)=zcon;
		w1=waitbar(0,['Working on Basin Number 0 of ' num2str(num_basins) ' total basins']);
		for ii=1:num_basins

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

			% Generate new stream map by segmenting original full stream network
			DEM_res=DEM.cellsize;
			six=coord2ind(DEM,xx,yy);
			Sc=STREAMobj(FD,'minarea',threshold_area,'unit','mapunits','outlets',six);

			% Check to make sure the stream object isn't empty
			if isempty(Sc.x)
				warning(['Input threshold drainage area is too large for basin ' num2str(basin_num) ' decreasing threshold area for this basin']);
				new_thresh=threshold_area;
				while isempty(Sc.x)
					new_thresh=new_thresh/2;
					Sc=STREAMobj(FD,'minarea',new_thresh,'unit','mapunits','outlets',six);
				end
			end

			% Calculate chi and create chi map
			Cc=chitransform(Sc,A.*(A.cellsize^2),'a0',1,'mn',theta_ref);
			ChiOBJc=GRIDobj(DEMoc);
			ScIX=coord2ind(DEMoc,Sc.x,Sc.y);
			ChiOBJc.Z(ScIX)=Cc;

			% Calculate slope area
			SAc=slopearea(Sc,DEM,A,'plot',false);
			Goc=gradient8(DEMoc);

			% Calculate ksn
			switch ksn_method
			case 'quick'
				[MSc]=KSN_Quick(DEM,A,Sc,-1*(SAc.theta),segment_length);
				[MSNc]=KSN_Quick(DEM,A,Sc,theta_ref,segment_length);
			case 'trib'
				% Overide choice if very small basin as KSN_Trib will fail for small basins
				if drainage_area>2.5
					[MSc]=KSN_Trib(DEM,FD,A,Sc,-1*(SAc.theta),segment_length);
					[MSNc]=KSN_Trib(DEM,FD,A,Sc,theta_ref,segment_length);
				else
					[MSc]=KSN_Quick(DEM,A,Sc,-1*(SAc.theta),segment_length);
					[MSNc]=KSN_Quick(DEM,A,Sc,theta_ref,segment_length);				
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

			% Calculate basin wide elevation statistics
			min_z=nanmin(DEMoc.Z(:));
			mean_z=nanmean(DEMoc.Z(:));
			max_z=nanmax(DEMoc.Z(:));
			std_z=nanstd(DEMoc.Z(:));
			se_z=std_z/sqrt(sum(~isnan(DEMoc.Z(:)))); % Standard error

			KSNc_stats=[mean_ksn se_ksn std_ksn min_ksn max_ksn];
			Gc_stats=double([mean_grad se_grad std_grad min_grad max_grad]);
			Zc_stats=double([mean_z se_z std_z min_z max_z]);

			% Find outlet elevation
			out_ix=coord2ind(DEMoc,xx,yy);
			out_el=double(DEMoc.Z(out_ix));

			FileName=['Basin_' num2str(basin_num) '_Data.mat'];

			if ~isempty(AG)
				num_grids=size(AG,1);
				AGc=cell(size(AG));
				for jj=1:num_grids
					AGcOI=crop(AG{jj,1},I,nan);
					AGc{jj,1}=AGcOI;
					AGc{jj,2}=AG{jj,2};
					mean_AGc=nanmean(AGcOI.Z(:));
					min_AGc=nanmin(AGcOI.Z(:));
					max_AGc=nanmax(AGcOI.Z(:));
					std_AGc=nanstd(AGcOI.Z(:));
					se_AGc=std_AGc/sqrt(sum(~isnan(AGcOI.Z(:))));
					AGc_stats(jj,:)=[mean_AGc se_AGc std_AGc min_AGc max_AGc];
				end
				save(FileName,'RiverMouth','DEMoc','out_el','drainage_area','Sc','SAc','Goc','MSc','MSNc','KSNc_stats','Gc_stats','Zc_stats','Centroid','ChiOBJc','ksn_method','clip_method','AGc','AGc_stats');
			elseif isempty(AG)
				save(FileName,'RiverMouth','DEMoc','out_el','drainage_area','Sc','SAc','Goc','MSc','MSNc','KSNc_stats','Gc_stats','Zc_stats','Centroid','ChiOBJc','ksn_method','clip_method');
			end


			switch write_arc_files
			case true

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

			waitbar(ii/num_basins,w1,['Working on Basin Number ' num2str(ii) ' of ' num2str(num_basins) ' total basins'])
		end
		close(w1)
	end
end


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
