function ClassifyKnicks(DEM,FD,A,Sc,ksn_master,bnd_list,kn_list,varargin)
	%
	% Usage:
	%	ClassifyKnicks(DEM,FD,A,Sc,ksn_master,bnd_list,kn_list);
	%	ClassifyKnicks(DEM,FD,A,Sc,ksn_master,bnd_list,kn_list,'shape_name','name');
	%
	% Description:
	% 	Function to iterate through a set of bounds and knickpoints selected while running 'KsnProfiler'. The function 
	%	will display a long profile and chi - elevation plot for individual stream segments and will iterate through each
	%	bound point or knickpoint you selected in KsnProfiler. The code expects you to input a number or character (in the pop up dialog)
	%	to categorize the feature higlighted in red. You must be consistent in your choice (i.e. you must either use 
	%	numbers for all of the classifications or characters for all the classifications within a given run), mixing numbers 
	%	and characters will result in an error at the end of the run. For entering characters, it's recommended you keep these 
	%	short strings without spaces (i.e. entries supported into a shapefile a attribute table), e.g. knick or bound  
	%
	% Required Inputs:
	%	DEM - Digital Elevation used as input to KsnProfiler
	%	FD - Flow direction used as input to KsnProfiler
	%	A - Flow accumulation used as input to KsnProfiler
	%	Sc - Stream network output from KsnProfiler
	%	ksn_master - cell array of selected channels output from KsnProfiler
	%	bnd_list - matrix of bounds (i.e. knickpoints) output from KsnProfiler
	%
	% Optional Inputs:
	%	shape_name ['ksn'] - name for the shapefile to be export, must have no spaces to be a valid name for ArcGIS and should NOT include the '.shp'
	%
	% Outputs:
	%	saves a shapfile of knickpoints including the classification you assign using this tool
	%
	% Notes:
	%	You must provide both the 'bnd_list' and 'kn_list' outputs from KsnProfiler, even if one of them has no values that aren't NaN.
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 01/10/20 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'ClassifyKnicks';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'FD',@(x) isa(x,'FLOWobj'));
	addRequired(p,'A',@(x) isa(x,'GRIDobj'));
	addRequired(p,'Sc',@(x) isa(x,'STREAMobj'));
	addRequired(p,'ksn_master',@(x) iscell(x));
	addRequired(p,'bnd_list',@(x) ismatrix(x));
	addRequired(p,'kn_list',@(x) ismatrix(x));

	addParameter(p,'shape_name','ksn',@(x) ischar(x));

	parse(p,DEM,FD,A,Sc,ksn_master,bnd_list,kn_list,varargin{:});
	DEM=p.Results.DEM;
	FD=p.Results.FD;
	Sc=p.Results.Sc;
	A=p.Results.A;
	ksn_master=p.Results.ksn_master;
	bnd_list=p.Results.bnd_list;
	kn_list=p.Results.kn_list;

	shape_name=p.Results.shape_name;

	% Calculate flowdistances upstream of outlets of the streams
	OUTIX=GRIDobj(DEM,'logical');
	outix=streampoi(Sc,'outlets','ix');
	OUTIX.Z(outix)=true;
	FLDS=flowdistance(FD,OUTIX); 

	% Filter out streams with out bound knicks
	idx=isnan(bnd_list(:,1));
	bnd_list=bnd_list(~idx,:);

	bnd_x=bnd_list(:,1);
	bnd_y=bnd_list(:,2);
	bnd_ix=coord2ind(DEM,bnd_x,bnd_y);
	bnd_el=DEM.Z(bnd_ix);

	% Filter out streams with out knicks
	idx=isnan(kn_list(:,1));
	kn_list=kn_list(~idx,:);

	kn_x=kn_list(:,1);
	kn_y=kn_list(:,2);
	kn_ix=coord2ind(DEM,kn_x,kn_y);
	kn_el=DEM.Z(kn_ix);	

	% Find streams with knickpoints
	str_list=unique(vertcat(bnd_list(:,4),kn_list(:,4)));
	num_streams=numel(str_list);

	% Initiate arrays
	b_is_classified=logical(zeros(numel(bnd_ix),1));
	k_is_classified=logical(zeros(numel(kn_ix),1));

	for ii=1:num_streams

		% Generate entire stream to outlet
		ix=coord2ind(DEM,ksn_master{str_list(ii),1}(:,1),ksn_master{str_list(ii),1}(:,2));
		ref_theta=ksn_master{str_list(ii),1}(1,7);
		IX=influencemap(FD,ix);
		ST=STREAMobj(FD,IX);
		% Intersect with STREAMobj from KsnProfiler to correct for any modifications
		ST=intersect(ST,Sc);
		% Condition
		z=mincosthydrocon(ST,DEM,'interp',0.1);

		% Calculate chi
		c=chitransform(ST,A,'a0',1,'mn',ref_theta);
		C=GRIDobj(DEM);
		C.Z(ST.IXgrid)=c;

		% Determine the knicks on this stream
		bidx=ismember(bnd_ix,ST.IXgrid);
		kidx=ismember(kn_ix,ST.IXgrid);

		% Determine which have knicks have already been classified
		bidx2= bidx & ~b_is_classified;
		kidx2= kidx & ~k_is_classified;

		% Set is_classified for next loop
		b_is_classified(bidx2)=true;
		k_is_classified(kidx2)=true;

		% Select just the knickpoints of interest
		bnd_ixOI=bnd_ix(bidx2);
		bnd_xOI=bnd_x(bidx2);
		bnd_yOI=bnd_y(bidx2);
		bnd_class=zeros(numel(bnd_ixOI),1);

		kn_ixOI=kn_ix(kidx2);
		kn_xOI=kn_x(kidx2);
		kn_yOI=kn_y(kidx2);
		kn_class=ones(numel(kn_ixOI),1);

		% Combine bounds and knicks
		cm_ixOI=vertcat(bnd_ixOI,kn_ixOI);
		cm_xOI=vertcat(bnd_xOI,kn_xOI);
		cm_yOI=vertcat(bnd_yOI,kn_yOI);
		cm_class=vertcat(bnd_class,kn_class);

		if ~isempty(cm_ixOI)
			f1=figure(1);
			set(f1,'Units','normalized','Position',[0.5 0.1 0.5 0.5],'renderer','painters');
			clf
			sbplt1=subplot(2,1,1);
			hold on
			plotdz(ST,DEM,'Color',[0.5 0.5 0.5]);
			plotdz(ST,z,'Color','k');
			cnt=0;
			if ~isempty(bnd_ixOI)
				cnt=cnt+1;
				s1(cnt)=scatter(FLDS.Z(bnd_ixOI),DEM.Z(bnd_ixOI),30,'k','filled');
				leg{cnt,1}='Segment Bounds';
			end

			if ~isempty(kn_ixOI)
				cnt=cnt+1;
				s1(cnt)=scatter(FLDS.Z(kn_ixOI),DEM.Z(kn_ixOI),30,'k','s','LineWidth',2);
				leg{cnt,1}='Knickpoint';
			end

			if ~verLessThan('matlab','9.5')
				disableDefaultInteractivity(sbplt1);
			end	

			hold off

			sbplt2=subplot(2,1,2);
			hold on
			cvec=getnal(ST,C);
			[cvec,six]=sort(cvec);
			evec=z(six);
			plot(cvec,evec,'-k','LineWidth',2);
			if ~isempty(bnd_ixOI)
				scatter(C.Z(bnd_ixOI),DEM.Z(bnd_ixOI),30,'k','filled');
			end

			if ~isempty(kn_ixOI)
				scatter(C.Z(kn_ixOI),DEM.Z(kn_ixOI),30,'k','s','LineWidth',2);
			end
			if ~verLessThan('matlab','9.5')
				disableDefaultInteractivity(sbplt2);
			end				
			hold off

			hold off

			for jj=1:numel(cm_ixOI);
				subplot(2,1,1)
				hold on
				if cm_class(jj)==0
					p1=scatter(FLDS.Z(cm_ixOI(jj)),DEM.Z(cm_ixOI(jj)),50,'r','filled');
				elseif cm_class(jj)==1
					p1=scatter(FLDS.Z(cm_ixOI(jj)),DEM.Z(cm_ixOI(jj)),50,'r','filled','s');
				end

				% This is necessary because of 'leg' has identical inputs (i.e. {'Knickpoint','Knickpoint'})
				% it will error out
				if numel(unique(leg))==1
					legend(s1(1),leg{1},'location','best');
				else
					legend(s1,leg,'location','best');
				end
				hold off

				subplot(2,1,2)
				hold on
				if cm_class(jj)==0
					p2=scatter(C.Z(cm_ixOI(jj)),DEM.Z(cm_ixOI(jj)),50,'r','filled');
				elseif cm_class(jj)==1
					p2=scatter(C.Z(cm_ixOI(jj)),DEM.Z(cm_ixOI(jj)),50,'r','filled','s');
				end
				xlabel('\chi');
				ylabel('Elevation [m]');
				hold off

				c=inputdlg('Enter the classification for the selected feature:','Feature Classification');

				cn=str2num(c{1});
				if isempty(cn)
					char_flag=true;
				else
					char_flag=false;
				end

				new_bnds(jj,:)=[cm_xOI(jj) cm_yOI(jj) FLDS.Z(cm_ixOI(jj)) DEM.Z(cm_ixOI(jj)) str_list(ii)];
				if char_flag
					bnd_cat{jj,1}=c;
				else
					bnd_cat(jj,1)=cn;
				end
				delete(p1);
				delete(p2);
			end

			new_bnds_c{ii,1}=new_bnds;
			bnd_cats{ii,1}=bnd_cat;
		else
			new_bnds_c{ii,1}=NaN;
			bnd_cats{ii,1}=NaN;
		end
	end

	new_bnds=vertcat(new_bnds_c{:});
	[new_bnds,bix,~]=unique(new_bnds,'rows');
	bnd_cats=vertcat(bnd_cats{:});
	bnd_cats=bnd_cats(bix);

	if char_flag
		bnd_cats=vertcat(bnd_cats{:});
	end

	KNK=struct;
	for jj=1:numel(new_bnds(:,1));
		KNK(jj,1).Geometry='Point';
		KNK(jj,1).X=double(new_bnds(jj,1));
		KNK(jj,1).Y=double(new_bnds(jj,2));
		KNK(jj,1).Dist=double(new_bnds(jj,3));
		KNK(jj,1).Elev=double(new_bnds(jj,4));
		KNK(jj,1).StrNum=double(new_bnds(jj,5));
		if char_flag
			KNK(jj,1).Category=bnd_cats{jj,1};
		else 
			KNK(jj,1).Category=double(bnd_cats(jj,1));
		end
	end
	out_knick_name=[shape_name '_knicks_classified.shp'];
	shapewrite(KNK,out_knick_name);

	close(f1)
end



