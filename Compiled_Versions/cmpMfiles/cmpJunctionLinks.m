function cmpJunctionLinks(wdir,MatFile,JunctionMatFile,varargin)
	% Description:
	%	Function for identifying the links within a stream network. Results will
	%	be similar (but not identical) to networksegments function. Classifies 
	%	links according to James & Krumbein, 1969, 'Frequency distributions of 
	%	stream link lengths', classifying links as either interior or exterior
	%	(exterior links have channel heads as their upstream ends) and either 
	%	cis or trans (cis links have the same handedness on both ends of the
	%	link, trans links have the oppostite handedness on the upstream and
	%	downstream links). Exterior links, interior links that connect to an
	%	outlet, links for which one or more of the junctions at the end
	%	have more than two upstream links, and interior links for which there is not
	% 	clear handedness (i.e. the incoming streams have the same Shreve order)
	%	will be undefined in terms of a cis or trans classification.
	%
	% Required Inputs:
	%	wdir - full path of working directory
	%	MatFile - Full path of matfile output from either 'cmpMakeStreams' or the name of a single basin mat 
	%		file from 'cmpProcessRiverBasins', in this case it must be the same matfile you provided to JunctionAngle
	%	JunctionMatFile - Full path of matfile output from JunctionAngle
	%
	% Optional Inputs:
	%	file_name_prefix ['topo'] - prefix for outputs shapefiles, textfiles, and mapfiles		
	%	make_shape [false] - logical flag to generate a shapefile containing the stream network broken
	%		into segments defined by junctions and categorized in terms of link_type and link_side
	%		(see outputs). Production of the mapstructure used to create shapefile can be time consuming
	%
	% Outputs:
	%	'*_junction_links.txt' - table containing the information for all the stream links:
	%		link_number - ID number of link
	%		downstream_x - x coordinate of downstream end of link
	%		donstream_y - y coordinate of downstream end of link
	%		downstream_IX - index into GRIDobj of downstream end of link
	%		upstream_x - x coordinate of upstream end of link
	%		upstream_y - y coordinate of upstream end of link
	%		upstream_IX - index into GRIDobj of upstream end of link
	%		link_type - classification of link is either Exterior (upstream end
	%			of link is channelhead) or Interior
	%		downstream_junc_num - junction number of downstream end of link referencing 
	%			the junction table provided to J, a NaN indicates that the downstream
	%			end of the link is not a junction
	%		upstream_junc_num - junction number of upstream end of link referencing
	%			the junction table provided to J, a NaN indicates that the upstream
	%			end of the link is not a junction
	%		link_side - classifcation of link as either Cis (both ends of link are 
	%			same handedness), Trans (ends of link are opposite handedness), or
	%			Undefined if one or both of the ends of the link do not have a handedness
	%			classification
	%
    % Examples if running for the command line, minus OS specific way of calling main TAK function:
	%	JunctionLinks /path/to/wdir Topo.mat Topo_junctions.mat
	%	JunctionLinks /path/to/wdir Topo.mat Topo_junctions.mat make_shape true	
	%
	% Related Functions:
	%	JunctionAngle InspectJunction networksegment
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 06/27/19 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if isdeployed
		if ~isempty(varargin)
			varargin=varargin{1};
		end
	end

	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'cmpJunctionLinks';
	addRequired(p,'wdir',@(x) ischar(x));
	addRequired(p,'MatFile',@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.mat'))));
	addRequired(p,'JunctionMatFile',@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.mat'))));

	addParameter(p,'file_name_prefix','topo',@(x) ischar(x));
	addParameter(p,'make_shape',false,@(x) isscalar(x) && islogical(x));

	parse(p,wdir,MatFile,JunctionMatFile,varargin{:});
	wdir=p.Results.wdir;
	MatFile=p.Results.MatFile;
	JMat=p.Results.JunctionMatFile;	

	file_name_prefix=p.Results.file_name_prefix;
	make_shape=p.Results.make_shape;

	% Determine the type of input
	MatFile=fullfile(wdir,MatFile);
	D=whos('-file',MatFile);
	VL=cell(numel(D),1);
	for ii=1:numel(D);
		VL{ii}=D(ii,1).name;
	end

	if any(strcmp(VL,'FD')) & any(strcmp(VL,'S'))
		load(MatFile,'FD','S');
	elseif any(strcmp(VL,'FDc')) & any(strcmp(VL,'Sc'))
		load(MatFile,'FDc','Sc');
		FD=FDc;
		S=Sc;
	end

	JMat=fullfile(wdir,JMat);
	load(JMat,'IX','junctions');
	J=junctions;

	% Grab handedness
	handedness=J.handedness;	

	% Find number of juctions (valid and invalid)
	num_con=size(IX,1);
	con=IX(:,1);
	% Determine junctions with invalid downstream nodes
	vld=~cellfun(@isempty,IX(:,2));

	% Generate original junction position list
	olist=1:num_con;
	olist=olist(vld)';

	% Populate lists of down- and up-stream indices
	dsix=zeros(num_con,1);
	us1ix=zeros(num_con,1);
	us2ix=zeros(num_con,1);

	for ii=1:nnz(vld)
		% Find nodes up and downstream of junctions
		JOI=IX(olist(ii),:);
		dsix(olist(ii),1)=JOI{1,2}(1);
		us1ix(olist(ii),1)=JOI{1,3}(1);
		us2ix(olist(ii),1)=JOI{1,4}(1);
	end

	%% Build links - Modified from networksegment.m
	% Build logical of head, up-stream confluences, and outlet
	% position in nal. 
	headIDX = streampoi(S,'channelheads','logical');      
	outIDX = streampoi(S,'outlets','logical');            
	bconfIDX = streampoi(S,'bconfluences','logical'); 

	% Find positions in nal of heads, up-stream confluences, down-stream confluences
	% and outlets
	ihead=find(headIDX==1);
	iout=find(outIDX==1);
	ibconf=find(bconfIDX==1);
	iconf=dsix(dsix>0);

	% Find GRID indices of heads, up-stream confluences, down-stream confluences and
	% outlets
	IXhead=S.IXgrid(ihead);
	IXout=S.IXgrid(iout);
	IXbconf=S.IXgrid(ibconf);
	IXconf=S.IXgrid(iconf);

	% Build drainage basin object establishing links
	DB=drainagebasins(FD,vertcat(IXbconf,IXout));

	% Extract link numbers for heads, up-stream confluences, down-stream confluences
	% and outlets
	DBhead=DB.Z(IXhead); 
	DBout=DB.Z(IXout);
	DBbconf=DB.Z(IXbconf); 
	DBconf=DB.Z(IXconf); 

	% Find links between channel heads and up-stream confluences
	[~,ind11,ind12]=intersect(DBbconf,DBhead);
	% Find links between down-stream confluences and up-stream confluencres
	[~,ind21,ind22]=intersect(DBbconf,DBconf);
	% Find inks between channel heads and outlets
	[~,ind31,ind32]=intersect(DBout,DBhead);
	% Find links between down-stream confluences and outlets
	[~,ind41,ind42]=intersect(DBout,DBconf);

	% Connect nodes into links
	ix(:,1)= [ ibconf(ind11)' ibconf(ind21)' iout(ind31)'  iout(ind41)'  ]; % Downstream end of link
	ix(:,2)= [ ihead(ind12)'  iconf(ind22)'  ihead(ind32)' iconf(ind42)' ]; % Upstream end of link

	% Build link table
	num_links=numel(ix(:,1));
	link_number=1:num_links;
	link_number=link_number(:);

	downstream_x=S.x(ix(:,1));
	downstream_y=S.y(ix(:,1));

	upstream_x=S.x(ix(:,2));
	upstream_y=S.y(ix(:,2));

	downstream_IX=S.IXgrid(ix(:,1));
	upstream_IX=S.IXgrid(ix(:,2));

	% Allocate arrays for junction numbers
	% downstream_junction=zeros(num_links,1);
	upstream_junc_num=zeros(num_links,1);
	downstream_junc_num=zeros(num_links,1);
	% Allocate link_side
	link_side=categorical(zeros(num_links,1));
	% Precompute data for determining interior vs exterior links
	inex=ismember(ix(:,2),ihead);
	link_type=categorical(zeros(num_links,1));
	for ii=1:num_links

		% Reset logical flags
		down1=false;
		down2=false;
		up=false;

		% Extract indices defining link
		down_ix=ix(ii,1);
		up_ix=ix(ii,2);

		% Label interior vs exterior links
		if inex(ii)
			link_type(ii,1)='Exterior';
			upstream_junc_num(ii,1)=NaN;
			link_side(ii,1)='Undefined';
		else
			link_type(ii,1)='Interior';
			% Find upstream junction number
			pos_in_ds=find(dsix==up_ix);
			if isempty(pos_in_ds)
				upstream_junc_num(ii,1)=NaN;
				link_side(ii,1)='Undefined';
			else
				upstream_junc_num(ii,1)=pos_in_ds;
				up=true;
			end

		end

		% Find downstream junction number
		pos_in_us1=find(us1ix==down_ix);
		pos_in_us2=find(us2ix==down_ix);
		if isempty(pos_in_us1) & isempty(pos_in_us2)
			downstream_junc_num(ii,1)=NaN;
			link_side(ii,1)='Undefined';
		elseif isempty(pos_in_us1)
			downstream_junc_num(ii,1)=pos_in_us2;
			down2=true;
		elseif isempty(pos_in_us2)
			downstream_junc_num(ii,1)=pos_in_us1;
			down1=true;
		else
			downstream_junc_num(ii,1)=NaN;
			link_side(ii,1)='Undefined';
		end

		% Determine link_type if applicable
		if up & down1
			hu=handedness(pos_in_ds);
			hd=handedness(pos_in_us1);
			if hu==hd & hu~='Undefined' & hd~='Undefined'
				link_side(ii,1)='Cis';
			elseif hu~=hd & hu~='Undefined' & hd~='Undefined'
				link_side(ii,1)='Trans';
			else 
				link_side(ii,1)='Undefined';
			end

		elseif up & down2
			hu=handedness(pos_in_ds);
			hd=handedness(pos_in_us2);

			if hu==hd & hu~='Undefined' & hd~='Undefined'
				link_side(ii,1)='Cis';
			elseif hu~=hd & hu~='Undefined' & hd~='Undefined'
				link_side(ii,1)='Trans';
			else
				link_side(ii,1)='Undefined';
			end

		end			
	end

	% Compile into link table
	links=table(link_number,downstream_x,downstream_y,downstream_IX,...
		upstream_x,upstream_y,upstream_IX,...
		link_type,downstream_junc_num,upstream_junc_num,link_side);

	table_name=fullfile(wdir,[file_name_prefix '_junction_links.txt']);
	writetable(links,table_name);

	if make_shape
		ms=makelinkshape(S,links);
		disp('Writing shapefile');
		shape_name=fullfile(wdir,[file_name_prefix '_junction_links.shp']);
		shapewrite(ms,shape_name);
	end

% Function End	
end

function [ms]=makelinkshape(S,L)
	disp('Starting shapefile construction')
	% Split stream network
	SS=split(S);

	% Make mapstructure based on split & modify
	disp('Building original mapstructure')
	ms=STREAMobj2mapstruct(SS);
	ms=rmfield(ms,{'streamorder','IX','tribtoIX'});
	ms(1,1).link_type=[]; ms(1,1).link_side=[]; ms(1,1).link_num=[];

	%Grab out data of interest from link table
	dx=L.downstream_x; dy=L.downstream_y; ux=L.upstream_x; uy=L.upstream_y;
	lty=L.link_type; lsd=L.link_side; ln=L.link_number;

	[ms,w1]=par_loop_proc(ms,dx,dy,ux,uy,lty,lsd,ln);
	close(w1);

end

function [ms,w1]=par_loop_proc(ms,dx,dy,ux,uy,lty,lsd,ln)

	DQ=parallel.pool.DataQueue;
	w1=waitbar(0,'Populating link information into mapstructure');
	afterEach(DQ,@updateBar);
	num_loop=numel(ms);
	cnt=1;

	parfor ii=1:num_loop
		xl=ms(1,ii).X;
		yl=ms(1,ii).Y;

		fun=@(DXL,UXL,DYL,UYL) any(ismember(DXL,xl)) & any(ismember(UXL,xl)) & any(ismember(DYL,yl)) & any(ismember(UYL,yl)); 
		idx=arrayfun(fun,dx,ux,dy,uy);
		ix=find(idx);

		if numel(ix)==1
			ms(1,ii).link_type=char(lty(ix));
			ms(1,ii).link_side=char(lsd(ix));
			ms(1,ii).link_num=ln(ix);
		elseif numel(ix>1)
			ms(1,ii).link_type=char(lty(ix(1)));
			ms(1,ii).link_side=char(lsd(ix(1)));
			ms(1,ii).link_num=ln(ix(1));
		end

		send(DQ,ii);
	end

		function updateBar(~)
			waitbar(cnt/num_loop,w1);
			cnt=cnt+1;
		end
end