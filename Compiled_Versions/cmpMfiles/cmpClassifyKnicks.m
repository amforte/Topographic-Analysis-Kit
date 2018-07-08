function cmpClassifyKnicks(wdir,MatFile,KsnProfilerMat,varargin)
	% Description:
	% 	Function to iterate through a set of bounds (i.e. knickpoints) selected while running 'KsnProfiler'. The function 
	%	will display a long profile and chi - elevation plot for individual stream segments and will iterate through each
	%	bound point you selected in KsnProfiler. The code expects you to input a number or character (at the command prompt)
	%	to categorize the knickpoint higlighted in red. You must be consistent in your choice (i.e. you must either use 
	%	numbers for all of the classifications or characters for all the classifications within a given run), mixing numbers 
	%	and characters will result in an error at the end of the run. For entering characters, it's recommended you keep these 
	%	short strings without spaces (i.e. entries supported into a shapefile a attribute table), e.g. knick or bound  
	%
	% Required Inputs:
	%	wdir - full path of working directory
	%	MatFile - Full path of matfile output from either 'cmpMakeStreams' or the name of a single basin mat file from 
	%		'cmpProcessRiverBasins' that was used to run 'cmpKsnProfiler'
	%	KsnProfilerMat - Full path of matfile output from 'cmpKsnProfiler'
	%
	% Optional Inputs:
	%	shape_name ['ksn'] - name for the shapefile to be export, must have no spaces to be a valid name for ArcGIS and should NOT include the '.shp'
	%
	% Outputs:
	%	saves a shapfile of knickpoints including the classification you assign using this tool
	%
	% Examples if running for the command line, minus OS specific way of calling main TAK function:
	%	ClassifyKnicks /path/to/wdir Topo.mat KsnProfiler.mat
	%	ClassifyKnicks /path/to/wdir Topo.mat KsnProfiler.mat shape_name my_knicks
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
	p.FunctionName = 'cmpClassifyKnicks';
	addRequired(p,'wdir',@(x) ischar(x));
	addRequired(p,'MatFile',@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.mat'))));
	addRequired(p,'KsnProfilerMat',@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.mat'))));

	addParameter(p,'shape_name','ksn',@(x) ischar(x));

	parse(p,wdir,MatFile,KsnProfilerMat,varargin{:});
	wdir=p.Results.wdir;
	MatFile=p.Results.MatFile;
	KsnProfilerMat=p.Results.KsnProfilerMat;

	shape_name=p.Results.shape_name;

	% Determine the type of input
	MatFile=fullfile(wdir,MatFile);
	D=whos('-file',MatFile);
	VL=cell(numel(D),1);
	for ii=1:numel(D);
		VL{ii}=D(ii,1).name;
	end

	if any(strcmp(VL,'DEM')) & any(strcmp(VL,'FD')) & any(strcmp(VL,'A')) & any(strcmp(VL,'S'))
		load(MatFile,'DEM','FD','A');
	elseif any(strcmp(VL,'DEMoc')) & any(strcmp(VL,'FDc')) & any(strcmp(VL,'Ac')) & any(strcmp(VL,'Sc'))
		load(MatFile,'FDc','Ac','Sc','DEMcc');
		DEM=DEMcc;
		FD=FDc;
		A=Ac;
	end	

	load(fullfile(wdir,KsnProfilerMat),'Sc','ksn_master','bnd_list');

	FLDS=flowdistance(FD);

	% Filter out streams with out knicks
	idx=isnan(bnd_list(:,1));
	bnd_list=bnd_list(~idx,:);

	bnd_x=bnd_list(:,1);
	bnd_y=bnd_list(:,2);
	bnd_ix=coord2ind(DEM,bnd_x,bnd_y);
	bnd_el=DEM.Z(bnd_ix);

	% Find streams with knickpoints
	str_list=unique(bnd_list(:,4));
	num_streams=numel(str_list);

	% Initiate arrays
	is_classified=logical(zeros(numel(bnd_ix),1));

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
		idx=ismember(bnd_ix,ST.IXgrid);
		% Determine which have knicks have already been classified
		idx2= idx & ~is_classified;

		% Set is_classified for next loop
		is_classified(idx2)=true;

		% Select just the knickpoints of interest
		bnd_ixOI=bnd_ix(idx2);
		bnd_xOI=bnd_x(idx2);
		bnd_yOI=bnd_y(idx);

		if ~isempty(bnd_ixOI)
			f1=figure(1);
			set(f1,'Units','normalized','Position',[0.5 0.1 0.5 0.5],'renderer','painters');
			clf
			subplot(2,1,1)
			hold on
			plotdz(ST,DEM,'Color',[0.5 0.5 0.5]);
			plotdz(ST,z,'Color','k');
			scatter(FLDS.Z(bnd_ixOI),DEM.Z(bnd_ixOI),20,'k','filled');
			hold off

			subplot(2,1,2)
			hold on
			cvec=getnal(ST,C);
			[cvec,six]=sort(cvec);
			evec=z(six);
			plot(cvec,evec,'-k','LineWidth',2);
			scatter(C.Z(bnd_ixOI),DEM.Z(bnd_ixOI),20,'k','filled');
			hold off

			hold off

			for jj=1:numel(bnd_ixOI);
				subplot(2,1,1)
				hold on
				p1=scatter(FLDS.Z(bnd_ixOI(jj)),DEM.Z(bnd_ixOI(jj)),40,'r','filled');
				hold off

				subplot(2,1,2)
				hold on
				p2=scatter(C.Z(bnd_ixOI(jj)),DEM.Z(bnd_ixOI(jj)),40,'r','filled');
				xlabel('\chi');
				ylabel('Elevation [m]');
				hold off

				c=inputdlg('Enter the classification for the selected knickpoint:','Knickpoint Classification');

				cn=str2num(c{1});
				if isempty(cn)
					char_flag=true;
				else
					char_flag=false;
				end

				new_bnds(jj,:)=[bnd_xOI(jj) bnd_yOI(jj) FLDS.Z(bnd_ixOI(jj)) DEM.Z(bnd_ixOI(jj)) str_list(ii)];
				if char_flag
					bnd_cat{jj,1}=c{1};
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
	out_knick_name=fullfile(wdir,[shape_name '_knicks_classified.shp']);
	shapewrite(KNK,out_knick_name);
	close(f1);
end