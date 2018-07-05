function cmpSubDivideBigBasins(wdir,basin_dir,max_basin_size,divide_method,varargin)
	% Description:
	% 	Function takes outputs from 'ProcessRiverBasins' function and subdvides any basin with a drainage area above a specified size and
	% 	outputs clipped dem, stream network, variout topographic metrics, and river values (ks, ksn, chi)
	%
	% Required Inputs:
	%		wdir - full path of working directory
	% 		basin_dir - full path of folder which contains the mat files from 'ProcessRiverBasins'
	% 		max_basin_size - size above which drainage basins will be subdivided in square kilometers
	%		divide_method - method for subdividing basins, options are ('confluences' and 'up_confluences' is NOT recommended large datasets):
	%			'order' - use the outlets of streams of a given order that the user can specify with the optional 's_order' parameter 
	%			'confluences' - use the locations of confluences (WILL PRODUCE A LOT OF SUB BASINS!). There is an internal parameter to remove
	%				extremely short streams that would otherwise result in the code erroring out.
	%			'up_confluences' - use locations just upstream of confluences (WILL PRODUCE A LOT OF SUB BASINS!). There is an internal parameter
	%				to remove extremely short streams that otherwise result in the code erroring out.
	%			'filtered_confluences' - use locations of confluences if drainage basin above confluence is of a specified size that the user
	%				can specify with the optional 'min_basin_size'  
	%			'p_filtered_confluences' - similar to filtered confluences, but the user defines a percentage of the main basin area
	%				with the optional 'min_basin_size'
	%			'trunk' - uses the tributary junctions with the trunk stream within the main basin as pour points for subdivided basins. There is
	%				an internal parameter to remove extremely short streams that would otherwise result in the code erroring out.
	%			'filtered_trunk' - same as 'trunk' but will only include basins that are greater than the min_basin_size
	%			'p_filtered_trunk' - same as 'filtered_trunk' but 'min_basin_size' is interpreted as a percentage of the main basin area
	%
	% Optional Inputs:
	%		SBFiles_Dir ['SubBasins'] - name of folder (within the main Basins folder) to store the subbasin files. Subbasin files are now stored in
	%			a separate folder to aid in the creation of different sets of subbasins based on different requirements. 
	%		recursive [true] - logical flag to ensure no that no subbasins in the outputs exceed the 'max_basin_size' provided. If 'divide_method' is 
	%			one of the trunk varieties the code will continue redefining trunks and further split subbasins until no extracted basins are greater
	%			than the 'max_basin_size'. If the 'divide_method' is one of the confluence varities, subbasins greater than 'max_basin_size' will simply
	%			no be included in the output. The 'recursive' check is not implemented for the 'order' method.
	% 		threshold_area [1e6] - minimum accumulation area to define streams in meters squared
	% 		segment_length [1000] - smoothing distance in meters for averaging along ksn, suggested value is 1000 meters
	% 		ref_concavity [0.5] - reference concavity for calculating ksn
	% 		write_arc_files [false] - set value to true to output a ascii's of various grids and a shapefile of the ksn, false to not output arc files
	%		s_order [3] - stream order for defining stream outlets for subdividing if 'divide_method' is 'order' (lower number will result in more sub-basins)
	%		min_basin_size [10] - minimum basin size for auto-selecting sub basins. If 'divide_method' is 'filtered_confluences' this value is
	%			interpreted as a minimum drainage area in km^2. If 'divide_method' is 'p_filtered_confluences', this value is interpreted as
	%			the percentage of the input basin drainage area to use as a minimum drainage area, enter a value between 0 and 100 in this case.
	%		no_nested [false] - logical flag that when used in conjunction with either 'filtered_confluences' or 'p_filtered_confluences' will only extract
	%			subbasins if they are the lowest order basin that meets the drainage area requirements (this is to avoid producing nested basins)
	%
	% Notes:
	%	-Only the 'order', 'trunk', 'filtered_trunk', and 'p_filtered_trunk' divide methods will not produce nested subbasins.
	% 	-The interpolation necessary to produce a continous ksn grid will fail on extremely small basins. This will not cause the code to fail, but will result in
	%		no 'KsnOBJc' being saved for these basins.
	%	-Methods 'confluences', 'up_confluences', and 'trunk' can result in attempts to extract very small basins. There is an internal check on this that attempts to remove 
	%		these very small basins but it is not always effective and can occassionally result in errors. If you are encountering errors try running the drainage area
	%		filtered versions
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 07/02/18 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if isdeployed
		if ~isempty(varargin)
			varargin=varargin{1};
		end
	end
	
	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'cmpSubDivideBigBasins';
	addRequired(p,'wdir',@(x) ischar(x));
	addRequired(p,'basin_dir',@(x) ischar(x));
	addRequired(p,'max_basin_size',@(x) isnumeric(x));
	addRequired(p,'divide_method',@(x) ischar(validatestring(x,{'order','confluences','up_confluences','filtered_confluences','p_filtered_confluences','trunk','filtered_trunk','p_filtered_trunk'})));

	addParameter(p,'SBFiles_Dir','SubBasins',@(x) ischar(x));
	addParameter(p,'recursive',true,@(x) isscalar(x) && islogical(x));
	addParameter(p,'ref_concavity',0.5,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'threshold_area',1e6,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'segment_length',1000,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'write_arc_files',false,@(x) isscalar(x));
	addParameter(p,'s_order',[3],@(x) isscalar(x));
	addParameter(p,'min_basin_size',[10],@(x) isnumeric(x) & isscalar(x));
	addParameter(p,'no_nested',false,@(x) isscalar(x) && islogical(x));

	parse(p,wdir,basin_dir,max_basin_size,divide_method,varargin{:});
	wdir=p.Results.wdir;
	location_of_data_files=p.Results.basin_dir;
	max_basin_size=p.Results.max_basin_size;
	divide_method=p.Results.divide_method;

	SBFiles_Dir=p.Results.SBFiles_Dir;
	recursive=p.Results.recursive;
	theta_ref=p.Results.ref_concavity;
	threshold_area=p.Results.threshold_area;
	segment_length=p.Results.segment_length;
	write_arc_files=p.Results.write_arc_files;
	s_order=p.Results.s_order;
	min_basin_size=p.Results.min_basin_size;
	no_nested=p.Results.no_nested;


	FileList=dir(fullfile(wdir,location_of_data_files,'*Data.mat'));
	num_files=numel(FileList);

	sb_dir=fullfile(wdir,location_of_data_files,SBFiles_Dir);

	% Make Subbasin Directory if it doesn't exist
	if ~isdir(sb_dir)
		mkdir(sb_dir);
	end


	if strcmp(divide_method,'p_filtered_confluences') | strcmp(divide_method,'p_filtered_trunk') & min_basin_size>100 | min_basin_size<=0
		error('For divide_method "p_filtered_confluences" the entry to "min_basin_size" must be between 0 and 100')
	end

	% Begin main file loop
	w1=waitbar(0,'Subdividing basins');
	for ii=1:num_files;
		FileName=fullfile(FileList(ii,1).folder,FileList(ii,1).name);

		% Load in drainage area to check against input
		load(FileName,'drainage_area');
		DA=drainage_area;

		% Check drainage area to determine if this basin will be processed
		if DA>=max_basin_size
			
			% Load in required basin files and rename
			load(FileName,'RiverMouth','DEMoc','DEMcc','FDc','Ac','Sc','ksn_method','gradient_method');
			DEM=DEMoc;
			DEMhc=DEMcc;
			S=Sc;
			FD=FDc;
			A=Ac;
			RM=RiverMouth;	
			basin_num=RM(:,3);

			waitbar(ii/num_files,w1,['Subdividing basin number ' num2str(basin_num) ' - Determining number of subdivisions']);

			DAG=(A.*(A.cellsize^2))/1e6;

			switch divide_method
			case 'order'
				so=streamorder(S);
				if s_order<max(so) 
					Se=modify(S,'streamorder',s_order);
					outs=streampoi(Se,'outlets','xy');
					x=outs(:,1);
					y=outs(:,2);
					num_new_basins=numel(x);
				elseif s_order>=max(so) & max(so)>1
					s_order=s_order-1;
					Se=modify(S,'streamorder',s_order);
					outs=streampoi(Se,'outlets','xy');
					x=outs(:,1);
					y=outs(:,2);
					num_new_basins=numel(x);	
				else
					s_order=max(so);
					Se=modify(S,'streamorder',s_order);
					outs=streampoi(Se,'outlets','xy');
					x=outs(:,1);
					y=outs(:,2);
					num_new_basins=numel(x);
				end			
			case 'confluences'
				S=removeshortstreams(S,DEM.cellsize*10);	
				cons=streampoi(S,'confluences','xy');
				if recursive
					cons_ix=streampoi(S,'confluences','ix');
					idx=DAG.Z(cons_ix)<max_basin_size;
					x=cons(idx,1);
					y=cons(idx,2);
				else
					x=cons(:,1);
					y=cons(:,2);
				end
				num_new_basins=numel(x);
			case 'up_confluences'
				S=removeshortstreams(S,DEM.cellsize*10);
				cons=streampoi(S,'bconfluences','xy');
				if recursive
					cons_ix=streampoi(S,'bconfluences','ix');
					idx=DAG.Z(cons_ix)<max_basin_size;
					x=cons(idx,1);
					y=cons(idx,2);
				else
					x=cons(:,1);
					y=cons(:,2);
				end
				num_new_basins=numel(x);
			case 'filtered_confluences'
				if no_nested
					cons_ix=streampoi(S,'bconfluences','ix');
					if recursive
						da_idx=DAG.Z(cons_ix)>=min_basin_size & DAG.Z(cons_ix)<max_basin_size;
						cons_ix=cons_ix(da_idx);
					else
						da_idx=DAG.Z(cons_ix)>=min_basin_size;
						cons_ix=cons_ix(da_idx);
					end
					[x,y]=CheckUpstream(DEM,FD,cons_ix);
					num_new_basins=numel(x);
				else
					cons_ix=streampoi(S,'confluences','ix');
					cons=streampoi(S,'confluences','xy');
					if recursive
						da_idx=DAG.Z(cons_ix)>=min_basin_size & DAG.Z(cons_ix)<max_basin_size;
					else
						da_idx=DAG.Z(cons_ix)>=min_basin_size;
					end
					cons=cons(da_idx,:);
					x=cons(:,1);
					y=cons(:,2);
					num_new_basins=numel(x);
				end
			case 'p_filtered_confluences'
				if no_nested
					cons_ix=streampoi(S,'bconfluences','ix');
					da_cons=DAG.Z(cons_ix);
					mbz=DA*(min_basin_size/100);
					if recursive
						da_idx=da_cons>=mbz & da_cons<max_basin_size;
					else
						da_idx=da_cons>=mbz;
					end
					[x,y]=CheckUpstream(DEM,FD,cons_ix(da_idx));
					num_new_basins=numel(x);
				else
					cons_ix=streampoi(S,'confluences','ix');
					cons=streampoi(S,'confluences','xy');
					da_cons=DAG.Z(cons_ix);
					mbz=DA*(min_basin_size/100);
					if recursive
						da_idx=da_cons>=mbz & da_cons<max_basin_size;
					else
						da_idx=da_cons>=mbz;
					end
					cons=cons(da_idx,:);
					x=cons(:,1);
					y=cons(:,2);
					num_new_basins=numel(x);
				end
			case 'trunk'
				ST=trunk(klargestconncomps(S,1));
				S=removeshortstreams(S,DEM.cellsize*10);
				tix=streampoi(S,'bconfluences','ix');
				tix=ismember(ST.IXgrid,tix);
				ds=ST.distance;
				ds(~tix)=NaN;
				[~,tix]=max(ds);
				SupT=modify(S,'tributaryto',ST);
				cons=streampoi(SupT,'outlets','xy');
				cons_ix=streampoi(SupT,'outlets','ix');
				cons_ix=vertcat(cons_ix,ST.IXgrid(tix));
				x=cons(:,1); x=vertcat(x,ST.x(tix));
				y=cons(:,2); y=vertcat(y,ST.y(tix));
				num_new_basins=numel(x);

				if recursive
					rec_count=1;
					while any(DAG.Z(cons_ix)>=max_basin_size) & rec_count<=10;
						nidx=DAG.Z(cons_ix)>=max_basin_size;
						if any(nidx)
							x(nidx)=[];
							y(nidx)=[];
							ixs=cons_ix(nidx);
							for jj=1:numel(ixs)
								TIX=GRIDobj(DEM,'logical');
								TIX.Z(ixs(jj))=true;
								S_sub=modify(S,'upstreamto',TIX);
								S_sub=removeshortstreams(S_sub,DEM.cellsize*10);
								ST_sub=trunk(S_sub);
								tix=streampoi(S_sub,'bconfluences','ix');
								tix=ismember(ST_sub.IXgrid,tix);
								ds=ST_sub.distance;
								ds(~tix)=NaN;
								[~,tix]=max(ds);
								SupT_sub=modify(S_sub,'tributaryto',ST_sub);
								cons=streampoi(SupT_sub,'outlets','xy');
								cons_ix=streampoi(SupT_sub,'outlets','ix');
								cons_ix=vertcat(cons_ix,ST_sub.IXgrid(tix));
								xx=cons(:,1); xx=vertcat(xx,ST_sub.x(tix));
								yy=cons(:,2); yy=vertcat(yy,ST_sub.y(tix));
								x=vertcat(x,xx);
								y=vertcat(y,yy);
							end
						end
						rec_count=rec_count+1;
						if rec_count>10
							warning(['Subdivision of basin number ' num2str(basin_num) ' ended prematurely to avoid an infinite loop']);
						end
					end
					num_new_basins=numel(x);
				end

			case 'filtered_trunk'
				ST=trunk(klargestconncomps(S,1));
				S=removeshortstreams(S,DEM.cellsize*10);
				tix=streampoi(S,'bconfluences','ix');
				tix=ismember(ST.IXgrid,tix);
				ds=ST.distance;
				ds(~tix)=NaN;
				[~,tix]=max(ds);
				SupT=modify(S,'tributaryto',ST);				
				cons_ix=streampoi(SupT,'outlets','ix');
				cons_ix=vertcat(cons_ix,ST.IXgrid(tix));
				cons=streampoi(SupT,'outlets','xy');
				cons=vertcat(cons,[ST.x(tix) ST.y(tix)]);
				da_cons=DAG.Z(cons_ix);
				da_idx=da_cons>=min_basin_size;
				cons=cons(da_idx,:);
				cons_ix=cons_ix(da_idx);
				x=cons(:,1);
				y=cons(:,2);
				num_new_basins=numel(x);

				if recursive
					rec_count=1;
					while any(DAG.Z(cons_ix)>=max_basin_size) & rec_count<=10;
						nidx=DAG.Z(cons_ix)>=max_basin_size;
						if any(nidx)
							x(nidx)=[];
							y(nidx)=[];
							ixs=cons_ix(nidx);
							for jj=1:numel(ixs)
								TIX=GRIDobj(DEM,'logical');
								TIX.Z(ixs(jj))=true;
								S_sub=modify(S,'upstreamto',TIX);
								S_sub=removeshortstreams(S_sub,DEM.cellsize*10);
								ST_sub=trunk(S_sub);
								tix=streampoi(S_sub,'bconfluences','ix');
								tix=ismember(ST_sub.IXgrid,tix);
								ds=ST_sub.distance;
								ds(~tix)=NaN;
								[~,tix]=max(ds);
								SupT_sub=modify(S_sub,'tributaryto',ST_sub);
								cons=streampoi(SupT_sub,'outlets','xy');
								cons_ix=streampoi(SupT_sub,'outlets','ix');
								cons_ix=vertcat(cons_ix,ST_sub.IXgrid(tix));
								xx=cons(:,1); xx=vertcat(xx,ST_sub.x(tix));
								yy=cons(:,2); yy=vertcat(yy,ST_sub.y(tix));
								da_cons=DAG.Z(cons_ix);
								da_idx=da_cons>=min_basin_size;
								x=vertcat(x,xx(da_idx));
								y=vertcat(y,yy(da_idx));
							end
						end
						rec_count=rec_count+1;
						if rec_count>10
							warning(['Subdivision of basin number ' num2str(basin_num) ' ended prematurely to avoid an infinite loop']);
						end						
					end
					num_new_basins=numel(x);
				end

			case 'p_filtered_trunk'
				ST=trunk(klargestconncomps(S,1));
				S=removeshortstreams(S,DEM.cellsize*10);
				tix=streampoi(S,'bconfluences','ix');
				tix=ismember(ST.IXgrid,tix);
				ds=ST.distance;
				ds(~tix)=NaN;
				[~,tix]=max(ds);
				SupT=modify(S,'tributaryto',ST);
				cons_ix=streampoi(SupT,'confluences','ix');
				cons_ix=vertcat(cons_ix,ST.IXgrid(tix));
				cons=streampoi(SupT,'confluences','xy');
				cons=vertcat(cons,[ST.x(tix) ST.y(tix)]);
				da_cons=DAG.Z(cons_ix);
				mbz=DA*(min_basin_size/100);
				da_idx=da_cons>=mbz;
				cons=cons(da_idx,:);
				cons_ix=cons_ix(da_idx);
				x=cons(:,1);
				y=cons(:,2);
				num_new_basins=numel(x);

				if recursive
					rec_count=1;
					while any(DAG.Z(cons_ix)>=max_basin_size) & rec_count<=10;
						nidx=DAG.Z(cons_ix)>=max_basin_size;
						if any(nidx)
							x(nidx)=[];
							y(nidx)=[];
							ixs=cons_ix(nidx);
							for jj=1:numel(ixs)
								TIX=GRIDobj(DEM,'logical');
								TIX.Z(ixs(jj))=true;
								S_sub=modify(S,'upstreamto',TIX);
								S_sub=removeshortstreams(S_sub,DEM.cellsize*10);
								ST_sub=trunk(S_sub);
								tix=streampoi(S_sub,'bconfluences','ix');
								tix=ismember(ST_sub.IXgrid,tix);
								ds=ST_sub.distance;
								ds(~tix)=NaN;
								[~,tix]=max(ds);
								SupT_sub=modify(S_sub,'tributaryto',ST_sub);
								cons=streampoi(SupT_sub,'outlets','xy');
								cons_ix=streampoi(SupT_sub,'outlets','ix');
								cons_ix=vertcat(cons_ix,ST_sub.IXgrid(tix));
								xx=cons(:,1); xx=vertcat(xx,ST_sub.x(tix));
								yy=cons(:,2); yy=vertcat(yy,ST_sub.y(tix));
								da_cons=DAG.Z(cons_ix);
								da_idx=da_cons>=mbz;
								x=vertcat(x,xx(da_idx));
								y=vertcat(y,yy(da_idx));
							end
						end
						rec_count=rec_count+1;
						if rec_count>10
							warning(['Subdivision of basin number ' num2str(basin_num) ' ended prematurely to avoid an infinite loop']);
						end
					end
					num_new_basins=numel(x);
				end

			end

			% Nested Waitbar
			w2=waitbar(0,['Processing ' num2str(num_new_basins) ' new basins']);
			pos_w1=get(w1,'position');
			pos_w2=[pos_w1(1) pos_w1(2)-pos_w1(4) pos_w1(3) pos_w1(4)];
			set(w2,'position',pos_w2,'doublebuffer','on');

			waitbar(ii/num_files,w1,['Subdividing basin number ' num2str(basin_num)]);

			for jj=1:num_new_basins
				% waitbar(ii/num_files,w1,['Subdividing basin number ' num2str(basin_num) ' - Processing ' num2str(jj) ' of ' num2str(num_new_basins) ' new basins']);
				
				waitbar(jj/num_new_basins,w2,['Processing ' num2str(jj) ' of ' num2str(num_new_basins) ' new basins']);

				xx=x(jj);
				yy=y(jj);
				basin_string=sprintf([num2str(basin_num) '%03d'],jj);
				RiverMouth=[xx yy str2num(basin_string)];

				% Build dependenc map and clip out drainage basins
				I=dependencemap(FD,xx,yy);
				DEMoc=crop(DEM,I,nan);
				DEMcc=crop(DEMhc,I,nan);
				FDc=crop(FD,I);
				Ac=crop(A,I,nan);

				% Calculate drainage area
				dep_map=GRIDobj2mat(I);
				num_pix=sum(sum(dep_map));
				drainage_area=(num_pix*DEMoc.cellsize*DEMoc.cellsize)/(1e6);

				% Calculate hypsometry
				[rb,eb]=hypscurve(DEMoc,100);
				hyps=[rb eb];

				% Find weighted centroid of drainage basin
				[Cx,Cy]=FindCentroid(DEMoc);
				Centroid=[Cx Cy];

				% Generate new stream map
				Sc=STREAMobj(FDc,'minarea',threshold_area,'unit','mapunits');

				% Check to make sure the stream object isn't empty
				if isempty(Sc.x)
					warning(['Input threshold drainage area is too large for basin ' num2str(RiverMouth(:,3)) ' decreasing threshold area for this basin']);
					new_thresh=threshold_area;
					while isempty(Sc.x)
						new_thresh=new_thresh/2;
						Sc=STREAMobj(FDc,'minarea',new_thresh,'unit','mapunits');
					end
				end

				% Calculate chi and create chi map
				Cc=chitransform(Sc,Ac,'a0',1,'mn',theta_ref);
				ChiOBJc=GRIDobj(DEMoc);
				ChiOBJc.Z(Sc.IXgrid)=Cc;

				% Calculate gradient
				switch gradient_method
				case 'gradient8'
					Goc=gradient8(DEMoc);
				case 'arcslope'
					Goc=arcslope(DEMoc);
				end

				% Find best fit concavity
				SLc=klargestconncomps(Sc,1);
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
				min_ksn=nanmin([MSNc.ksn]);
				mean_ksn=nanmean([MSNc.ksn]);
				max_ksn=nanmax([MSNc.ksn]);
				std_ksn=nanstd([MSNc.ksn]);
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

				SubFileName=fullfile(sb_dir,['Basin_' num2str(basin_num) '_DataSubset_' num2str(jj) '.mat']);

				save(SubFileName,'RiverMouth','DEMcc','DEMoc','out_el','drainage_area','hyps','FDc','Ac','Sc','SLc','Chic','Goc','MSc','MSNc','KSNc_stats','Gc_stats','Zc_stats','Centroid','ChiOBJc','ksn_method','gradient_method','theta_ref');
				
				% Make interpolated ksn grid
				try 
					[KsnOBJc] = KsnInt(DEMoc,MSNc);
					save(FileName,'KsnOBJc','-append');
				catch
					warning(['Interpolation of KSN grid failed for basin ' num2str(RiverMouth(:,3))]);
				end

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
					save(SubFileName,'AGc','AGc_stats','-append');
				end

				VarInd=find(strcmp(cellstr(char(VarList.name)),'ACGc'));
				if ~isempty(VarInd)
					load(FileName,'ACGc');
					ACG=ACGc;
					num_grids=size(ACG,1);
					ACGc=cell(size(ACG));
					for kk=1:num_grids
						ACGcOI=crop(ACG{kk,1},I,nan);
						ACGc{kk,1}=ACGcOI;
						ACGc{kk,3}=ACG{kk,3};
						edg=ACG{kk,2}.Numbers;
						edg=edg+0.5;
						edg=vertcat(0.5,edg);
						[N,~]=histcounts(ACGcOI.Z(:),edg);
						T=ACG{kk,2};
						T.Counts=N';
						ACGc{kk,2}=T;
						ACGc_stats(kk,1)=[mode(ACGcOI.Z(:))];
					end
					save(SubFileName,'ACGc','ACGc_stats','-append');	
				end	

				VarInd=find(strcmp(cellstr(char(VarList.name)),'rlf'));
				if ~isempty(VarInd)
					load(FileName,'rlf');
					rlf_full=rlf; 
					num_rlf=size(rlf_full,1);
					rlf=cell(size(rlf_full));
					rlf_stats=zeros(num_rlf,6);
					for kk=1:num_rlf
						% Calculate relief
						radOI=rlf_full{kk,2};
						rlf{kk,2}=radOI;
						rlfOI=localtopography(DEMoc,radOI);
						rlf{kk,1}=rlfOI;
						% Calculate stats
						mean_rlf=nanmean(rlfOI.Z(:));
						min_rlf=nanmin(rlfOI.Z(:));
						max_rlf=nanmax(rlfOI.Z(:));
						std_rlf=nanstd(rlfOI.Z(:));
						se_rlf=std_rlf/sqrt(sum(~isnan(rlfOI.Z(:))));
						rlf_stats(kk,:)=[mean_rlf se_rlf std_rlf min_rlf max_rlf radOI];
					end
					save(SubFileName,'rlf','rlf_stats','-append');
				end					

				if write_arc_files
					% Replace NaNs in DEM with -32768
					Didx=isnan(DEMoc.Z);
					DEMoc_temp=DEMoc;
					DEMoc_temp.Z(Didx)=-32768;

					DEMFileName=fullfile(sb_dir,['Basin_' num2str(basin_num) '_DataSubset_' num2str(jj) '_DEM.txt']);
					GRIDobj2ascii(DEMoc_temp,DEMFileName);
					CHIFileName=fullfile(sb_dir,['Basin_' num2str(basin_num) '_DataSubset_' num2str(jj) '_CHI.txt']);
					GRIDobj2ascii(ChiOBJc,CHIFileName);
					KSNFileName=fullfile(sb_dir,['Basin_' num2str(basin_num) '_DataSubset_' num2str(jj) '_KSN.shp']);
					shapewrite(MSNc,KSNFileName);

					if calc_relief
						for kk=1:num_rlf
							RLFFileName=fullfile(sb_dir,['Basin_' num2str(basin_num) '_DataSubset_' num2str(jj) '_RLF_' num2str(rlf{kk,2}) '.txt']);
							GRIDobj2ascii(rlf{kk,1},RLFFileName);
						end
					end

					if ~isempty(AG);
						for kk=1:num_grids
							AGcFileName=fullfile(sb_dir,['Basin_' num2str(basin_num) '_DataSubset_' num2str(jj) '_' AGc{kk,2} '.txt']);
							GRIDobj2ascii(AGc{kk,1},AGcFileName);
						end
					end

					if ~isempty(ACG);
						for jj=1:num_grids
							ACGcFileName=fullfile(sb_dir,['Basin_' num2str(basin_num) '_' ACGc{jj,3} '.txt']);
							GRIDobj2ascii(ACGc{jj,1},ACGcFileName);
						end
					end
				end
			end % New basin loop end
			close(w2);
		end % Drainage Area check end
	end % Main Loop end
	close(w1);
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
	c=chitransform(S,A,'a0',1,'mn',theta_ref);
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

function [KSNGrid] = KsnInt(DEM,ksn_ms)
    [xx,yy]=getcoordinates(DEM);
    [X,Y]=meshgrid(xx,yy);

    ksn_cell=cell(numel(ksn_ms),1);
    for ii=1:numel(ksn_ms)
        ksn_cell{ii}=ones(numel(ksn_ms(ii).X),1)*ksn_ms(ii).ksn;
    end
    ksn_x=vertcat(ksn_ms.X); ksn_y=vertcat(ksn_ms.Y); ksn_ksn=vertcat(ksn_cell{:});
    idx=isnan(ksn_ksn);
    ksn_x(idx)=[];
    ksn_y(idx)=[];
    ksn_ksn(idx)=[];
    
    warning off
    Fk=scatteredInterpolant(ksn_x,ksn_y,ksn_ksn,'natural');
    ksn_int=Fk(X,Y);
    KSNGrid=GRIDobj(xx,yy,ksn_int);
    IDX=isnan(DEM);
    KSNGrid.Z(IDX.Z)=NaN;
    warning on
end

function [x,y] = CheckUpstream(DEM,FD,ix)
	% Build cell of influence list
	inflcs=cell(numel(ix),1);
	for ii=1:numel(ix)
	    IX=influencemap(FD,ix(ii));
	    inflcs{ii}=find(IX.Z);
	end
	    
	% Build index
	idx=zeros(numel(ix),1);
	idx=logical(idx);
	for ii=1:numel(ix)
	    inflcs_temp=inflcs;
	    inflcs_temp{ii}=[0];
	    up_member=cellfun(@(x) ismember(ix(ii),x),inflcs_temp);
	    if any(up_member)
	        idx(ii)=false;
	    else
	        idx(ii)=true;
	    end
	end

	[x,y]=ind2coord(DEM,ix(idx));
end