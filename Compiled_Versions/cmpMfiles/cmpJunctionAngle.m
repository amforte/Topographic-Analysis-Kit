function cmpJunctionAngle(wdir,MatFile,fit_distance,varargin);
	% Description: 
	%	Function calculates the angle of stream junctions. The basic strategy and nomenclature
	%	used in this function is adapted from Howard, 1971, 'Optimal angles of stream junction:
	%	Geometric, stability to capture, and minimum power criteria'. The code first fits the
	%	upstream and downstream links that meet at junctions with linear fits. The lengths of 
	%	stream used to perform these fits is controlled with the 'fit_distance' parameter. Junction
	%	angles are calculated as a pair of angles (e1 and e2) between each tributary (tributary 1 and 
	%	tributary 2) and the upstream projection of the downstream link. The total junction angle 
	%	is calculated as the angle between the two upstream links, which in well behaved junctions 
	%	should be the sum of e1 and e2, but is not always the case (i.e. for some junctions both
	%	upstream links are on the same side of the upstream projection of the downstream link, in these
	%	cases the junction angle will not be the sum of the e1 and e2 angles, but instead will be the true
	%	angle between link 1 and link 2, whereas the values of e1 ane e2 will be measured from the upstream plane). 
	%	The function also calculates predicted junction angles based on the simple geometric criteria
	%	described in Howard, 1971, and the handedness of junction angle as defined by James & Krumbein,
	%	1969, 'Frequency distributions of stream link lengths'. It is important to note that angles
	%	related to junctions at which more than two upstream links meet or for which there are no 
	%	downstream nodes (i.e. the junction is an outlet) are set to NaN as these are undefined in
	%	the Howard criteria.
	%
	% Required Inputs:
	%	wdir - full path of working directory
	%	MatFile - Full path of matfile output from either 'cmpMakeStreams' or the name of a single basin mat 
	%		file from 'cmpProcessRiverBasins'
	%	fit_distance - distance in stream distance (map units) for fitting stream links, 
	%		if more than one value is provided to fit_distance, it is assumed you want to calculate
	%		junction angles fitting stream segments with the range of distances. The provided value(s)
	%		can be interpeted as the number of stream nodes over which to fit the stream links if the
	%		optional 'use_n_nodes' is set to true
	%
	% Optional Inputs:
	%	file_name_prefix ['topo'] - prefix for outputs shapefiles, textfiles, and mapfiles	
	%	use_n_nodes [false] - logical flag to interpret the value(s) provided to fit_distance as
	%		the number of stream nodes to extract up and downstream as opposed to streamwise distance
	%	previous_IX [] - Calculation of the 'IX' output is the most time consuming aspect of the
	%		calculation, so if you wish to rerun the function with different parameters (e.g. 
	%		recalculating predicted angles with a different reference concavity or different method)
	%		the presupplying the IX result from a previous run can be useful. In the compiled version
	%		this accomplished by supplying the name of a matfile generated during a previous run of
	%		the compiled version of junction angle. It is important that this IX result was produced in 
	%		a run using the same MatFile inputs and that the maximum value provided to fit distance does 
	%		not exceed the maximum fit distance that was provided when you initially ran the function to 
	%		produce the provided IX input. The code will check that the fit distances are compatible, but
	%		 will not explicity check that S, A, and the DEM are the same.
	%	predict_angle_method ['area'] - method to use to calculate predicted junction angles. Both
	%		methods are simple geometric predictions from Howard, 1971. Valid inputs are 'slope', 
	%		'area',  or 'both'. If pred_angle_method is set to 'area' (default), then you should make sure the 
	%		ref_concavity is valid for the network you are analyzing, default value is 0.5. For the 'slope'
	%		method, the slope of the two upstream and one downstream link will be the mean of the gradient
	%		of those links over the length of stream sampled (controlled by the values provided to fit_distance).
	%		Providing 'both' will calculate both the slope and area predicted angles.
	%	ref_concavity - [0.5] - reference concavity for calculating predicted junction angles using
	%		Howard, 1971 area relation. If no value is provided, default value of 0.5 will be used.
	%	verbose [false] - logical flag to report progress through function. This can be useful because
	%		certain steps of the process are very time consuming.
	% 
	% Outputs:
	%	Outputs a table as a text file. If one value is provided to fit_distance then there
	%		will be one table. If multiple values are provided to fit_distance, then there will be 
	%		a table for each fit distance, in order of increasing distances. Columns of the table are:
	%		
	%		junction_number - ID number of junction
	%		junction_x - x coordinate of junction
	%		junction_y - y coordinate of junction
	%		junction_angle - angle in degrees between tributary 1 (e1) and tributary 2 (e2)
	%		handedness - handedness of junction as defined by James & Krumbein, 1969, options
	%			are Right, Left, or Undefined, stored as a categorical
	%		split - logical indicating whether upstream projection of downstream link lies between
	%			tributary 1 and tributary 2 (true) or does not (false)
	%		e1_obs_angle - angle between tributary 1 and upstream projection of downstream link
	%		e1_Apred_angle - predicted angle between tributary 1 and upstream projection of downstream link
	%			based on Howard 1971 area relations (will appear if 'predict_angle_method' is 'area' or 'both')
	%		e1_Spred_angle - predicted angle between tributary 1 and upstream projection of downstream link
	%			based on Howard 1971 slope relations (will appear if 'predict_angle_method' is 'slope' or 'both')	
	%		e1_rotation - orientation of tributary 1 with respect to the upstream projection of the
	%			downstream link, options are CCW (counter-clockwise), CW (clockwise), or Undefined
	%			(more than two tributaries meet at junction). If split is true, one of the upstream links will be
	%			CW and one will be CCW.
	%		e1_shreve - shreve order of tributary 1
	%		e1_direction - flow azimuth of tributary 1
	%		e1_distance - stream distance used to fit tributary 1
	%		e1_num - number of upstream nodes used to fit tributary 1
	%		e1_R2 - R2 value on linear fit of tributary 1
	%		e2_obs_angle - angle between tributary 2 and upstream projection of downstream link
	%		e2_Apred_angle - predicted angle between tributary 2 and upstream projection of downstream link
	%			based on Howard 1971 area relations (will appear if 'predict_angle_method' is 'area' or 'both')
	%		e2_Spred_angle - predicted angle between tributary 2 and upstream projection of downstream link
	%			based on Howard 1971 slope relations (will appear if 'predict_angle_method' is 'slope' or 'both')	
	%		e2_rotation - orientation of tributary 2 with respect to the upstream projection of the
	%			downstream link, options are CCW (counter-clockwise), CW (clockwise), or Undefined
	%			(more than two tributaries meet at junction)
	%		e2_shreve - shreve order of tributary 2
	%		e2_direction - flow azimuth of tributary 2
	%		e2_distance - stream distance used to fit tributary 2
	%		e2_num - number of upstream nodes used to fit tributary 2
	%		e2_R2 - R2 value on linear fit of tributary 2
	%		eS_direction - flow azimuth of downstream link
	%		eS_distance - stream distance used to fit downstream link
	%		eS_num - number of downstream nodes used to fit downstream link
	%		eS_R2 - R2 value of linear fit of downstream link
	%
	%	If multiple fit_distances are provided a separate 	'*_mean_junctions.txt' file will be saved, which 
	%		will be a table containing mean and standard deviations across the angles calculated using 
	%		the different fit distances, including the mean total junction angle, standard deviation of 
	%		the total junction angle, mean angle between tributary 1 and the upstream projection of the downstream 
	%		link (e1), standard deviation of e1, mean angle between tributary 2 and the upstream projection of the
	%		downstream link (e2), standard deviation of e2. If the predict_angle_method is set to either 'slope'
	%		or 'both', this will also include the mean and standard deviations of the predicted e1 and e2 angles
	%		based on the slope method. Means and standard deviations are not calculated for the area method because
	%		these will not vary as a function of fit distance.
	%		
	% Note: 
	%	For all outputs, the tributaries are sorted so that tributary 1 is always the 'larger' of the 
	%	two tributaries. The size of the tributaries are based on the shreve order of the two tributaries.
	%	In the event that the two tributaries have the same shreve order then the drainage area is used to
	%	determine respective size.
	%
    % Examples if running for the command line, minus OS specific way of calling main TAK function:
    %	JunctionAngle /path/to/wdir Topo.mat 1000
    %	JunctionAngle /path/to/wdir Topo.mat [1000 2000 5000]
    %	JunctionAngle /path/to/wdir Topo.mat 10 use_n_nodes true
	%	
	% Related Functions:
	% 	InspectJunction JunctionLinks
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 06/26/19 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if isdeployed
		if ~isempty(varargin)
			varargin=varargin{1};
		end
	end

	% Parse Inputs
	p = inputParser;
	p.FunctionName = 'cmpJunctionAngle';	
	addRequired(p,'wdir',@(x) ischar(x));
	addRequired(p,'MatFile',@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.mat'))));
	addRequired(p,'fit_distance',@(x) isnumeric(x));

	addParameter(p,'file_name_prefix','topo',@(x) ischar(x));
	addParameter(p,'use_n_nodes',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'previous_IX',[],@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.mat'))));
	addParameter(p,'predict_angle_method','area',@(x) ischar(validatestring(x,{'slope','area','both'})));
	addParameter(p,'ref_concavity',0.50,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'verbose',false,@(x) isscalar(x) && islogical(x));

	parse(p,wdir,MatFile,fit_distance,varargin{:});
	wdir=p.Results.wdir;
	MatFile=p.Results.MatFile;
	fit_distance=p.Results.fit_distance;

	file_name_prefix=p.Results.file_name_prefix;	
	use_n_nodes=p.Results.use_n_nodes;
	PIX=p.Results.previous_IX;
	method=p.Results.predict_angle_method;
	ref_theta=p.Results.ref_concavity;
	verbose=p.Results.verbose;

	% Determine the type of input
	MatFile=fullfile(wdir,MatFile);
	D=whos('-file',MatFile);
	VL=cell(numel(D),1);
	for ii=1:numel(D);
		VL{ii}=D(ii,1).name;
	end

	if any(strcmp(VL,'DEM')) & any(strcmp(VL,'FD')) & any(strcmp(VL,'A')) & any(strcmp(VL,'S'))
		load(MatFile,'DEM','FD','A','S');
	elseif any(strcmp(VL,'DEMcc')) & any(strcmp(VL,'FDc')) & any(strcmp(VL,'Ac')) & any(strcmp(VL,'Sc'))
		load(MatFile,'DEMoc','FDc','Ac','Sc');
		DEM=DEMoc;
		FD=FDc;
		A=Ac;
		S=Sc;
	end
 
	if use_n_nodes
		if numel(fit_distance)==1
			multi=false;
			n=fit_distance;
			if n<1
				n=1;
			end
			fit_distance=n*(hypot(S.cellsize,S.cellsize));
		else
			multi=true;
			n=max(fit_distance);
			if n<1
				n=1;
			end

			nlist=fit_distance;
			nlist(nlist<1)=1;
			nlist=unique(nlist);
			if numel(nlist)==1
				n=nlist;
				multi=false;
			end			

			fit_distance=nlist.*(hypot(S.cellsize,S.cellsize));
			fit_distance=sort(fit_distance);
		end
	else
		% Convert distance to number of nodes and check if multiple distances were provided.
		if numel(fit_distance)==1
			multi=false;
			n=floor(fit_distance/hypot(S.cellsize,S.cellsize));
			if n<1
				warning('Input fit_distance is too short, will not average across more than one node, setting to minimum distance');
				n=1;
			end
		else
			n=floor(max(fit_distance)/hypot(S.cellsize,S.cellsize));
			multi=true;
			if n<1
				warning('All input fit_distance are too short, will not average across more than one node, setting to single minimum distance');
				n=1;
				mutli=false;
			end

			nlist=floor(fit_distance/hypot(S.cellsize,S.cellsize));
			if any(nlist)<1
				warning('Some input fit_distance are too short, setting these to minimum distance');
				nlist(nlist<1)=1;
				nlist=unique(nlist);
				if numel(nlist)==1
					n=nlist;
					multi=false;
				end
			end 

			fit_distance=sort(fit_distance); 
		end
	end

	% Extract stream distances
	dst=S.distance;

	% Calculate Shreve stream orders
	so=streamorder(S,'shreve');

	% Calculate gradient if pred_angle_method is slope
	switch method
	case {'slope','both'}
		if verbose
			disp('Calculating stream gradients')
		end
		z=mincosthydrocon(S,DEM,'interp',0.1);
		sg=gradient(S,z);
		G=GRIDobj(DEM);
		G.Z(S.IXgrid)=sg;	
	end

	% Find indices of nodes up and downstream
	if isempty(PIX)
		IX=nnodesupdown(S,A,n,so,verbose);
	else
		load(fullfile(wdir,PIX),'IX');
		% Check that max distance is consistent 
		UPIX=IX(:,2);
		nn=cellfun(@numel,UPIX);
		maxn=max(nn);
		if maxn<n
			error('The provided input to "previous_IX" is incompatible with the maximum value within "fit_distance"');
		end
	end

	% Convert stream x,y coordinates to lat-lon to preserve angular relationships, unless
	% there is no projection (e.g. for model results) or projection is unrecognized
	if verbose
		disp('Converting to geographic coordinates...')
	end

	if isempty(S.georef)
		xl=S.x; yl=S.y;
		warning('No projection was found, angles will be calculated based on projected coordinates')
		projcoord=false;
	else
		try
			[yl,xl]=projinv(S.georef.mstruct,S.x,S.y);
			projcoord=true;
		catch
			xl=S.x; yl=S.y;
			warning('Projection was not recognized, angles will be calculated based on projected coordinates')
			projcoord=false;
		end
	end

	% Break up output of node indicees
	num_con=size(IX,1);
	con=IX(:,1);
	ds=IX(:,2);
	us=IX(:,3:end);
	% Check and mark if there any streams with more than 2 tribs and if any
	% downstream segments of confluences don't exist to populate logical index
	if size(us,2)>2
		usvld=cellfun(@isempty,us(:,3));
	else
		usvld=true(size(con));
	end
	dsvld=~cellfun(@isempty,ds);
	vld=dsvld & usvld;

	if verbose
		disp('Finding junction angles...')
	end

	if ~multi
		% Preallocate arrays
		link_angle=zeros(num_con,11);
		strm_dir=zeros(num_con,3);
		fit_streams=zeros(num_con,9);

		switch method
		case {'area','slope'}
			pred_angle=zeros(num_con,2);
		case 'both'
			pred_angle=zeros(num_con,4);
		end

		for ii=1:num_con
			if vld(ii)
				% Extract x y node lists
				% Confluence point
				xc=xl(con{ii}); yc=yl(con{ii});
				% Upstream links
				xu1=xl(us{ii,1}); yu1=yl(us{ii,1});
				xu2=xl(us{ii,2}); yu2=yl(us{ii,2});
				% Downstream link
				xd=xl(ds{ii,1}); yd=yl(ds{ii,1});

				% Translate so all links originate from confluence 
				xu1=xu1-xc; xu2=xu2-xc;
				yu1=yu1-yc; yu2=yu2-yc;
				xd=xd-xc; yd=yd-yc;

				%% Depending on orientation of links, doing a linear fit
				%  on the link may produce spurious results. To improve 
				%  the result, each link is rotated towards horizontal, 
				%  using the mean orientation of the link to find the angle
				%  of rotation

				% Find mean angle of stream points from stream orientation
				theta1ap=median(atan2(yu1,xu1));
				theta2ap=median(atan2(yu2,xu2));
				theta3ap=median(atan2(yd,xd));

				% Extract and calculate stream distances
				e1_dst=max(dst(us{ii,1}))-dst(con{ii});
				e2_dst=max(dst(us{ii,2}))-dst(con{ii});
				ds_dst=dst(con{ii})-min(dst(ds{ii,1}));

				% Calculate predicted angles
				switch method
				case 'area'
					% Howard 1971 area method
					ix1=S.IXgrid(us{ii,1});
					ix2=S.IXgrid(us{ii,2});
					da1=max(A.Z(ix1).*A.cellsize^2);
					da2=max(A.Z(ix2).*A.cellsize^2);
					e1p=acos(((da1+da2)/da1)^-ref_theta);
					e2p=acos(((da1+da2)/da2)^-ref_theta);
					% Convert to degrees
					e1p=rad2deg(e1p);
					e2p=rad2deg(e2p);
					% Store out
					pred_angle(ii,:)=[e1p e2p];
				case 'slope'
					% Howard 1971 slope method
					ix1=S.IXgrid(us{ii,1});
					ix2=S.IXgrid(us{ii,2});
					ix3=S.IXgrid(ds{ii,1});
					sl1=mean(G.Z(ix1));
					sl2=mean(G.Z(ix2));
					sl3=mean(G.Z(ix3));
					r1=sl3/sl1;
					r2=sl3/sl2;
					if r1>1
						e1p=0;
					else
						e1p=acos(r1);
					end

					if r2>1
						e2p=0;
					else
						e2p=acos(r2);
					end
					% Convert to degrees
					e1p=rad2deg(e1p);
					e2p=rad2deg(e2p);
					% Store out
					pred_angle(ii,:)=[e1p e2p];
				case 'both'
					% Howard 1971 area method
					ix1=S.IXgrid(us{ii,1});
					ix2=S.IXgrid(us{ii,2});
					da1=max(A.Z(ix1).*A.cellsize^2);
					da2=max(A.Z(ix2).*A.cellsize^2);
					e1Ap=acos(((da1+da2)/da1)^-ref_theta);
					e2Ap=acos(((da1+da2)/da2)^-ref_theta);
					% Convert to degrees
					e1Ap=rad2deg(e1Ap);
					e2Ap=rad2deg(e2Ap);

					% Howard 1971 slope method
					ix3=S.IXgrid(ds{ii,1});
					sl1=mean(G.Z(ix1));
					sl2=mean(G.Z(ix2));
					sl3=mean(G.Z(ix3));
					r1=sl3/sl1;
					r2=sl3/sl2;
					if r1>1
						e1Sp=0;
					else
						e1Sp=acos(r1);
					end

					if r2>1
						e2Sp=0;
					else
						e2Sp=acos(r2);
					end
					% Convert to degrees
					e1Sp=rad2deg(e1Sp);
					e2Sp=rad2deg(e2Sp);
					% Store out
					pred_angle(ii,:)=[e1Ap e2Ap e1Sp e2Sp];
				end


				% Rotate all points to horizontal
				[xu1r,yu1r]=rotcoord(xu1,yu1,-theta1ap,0,0); 
				[xu2r,yu2r]=rotcoord(xu2,yu2,-theta2ap,0,0);	
				[xdr,ydr]=rotcoord(xd,yd,-theta3ap,0,0);

				% Simple approximation of linear segment on
				% rotated positions
				a1=xu1r\yu1r;
				a2=xu2r\yu2r; 
				a3=xdr\ydr;

				% Find projected y positions
				yp1r=xu1r*a1;
				yp2r=xu2r*a2;
				yp3r=xdr*a3;

				% Rotate back
				[xp1,yp1]=rotcoord(xu1r,yp1r,theta1ap,0,0);
				[xp2,yp2]=rotcoord(xu2r,yp2r,theta2ap,0,0);
				[xp3,yp3]=rotcoord(xdr,yp3r,theta3ap,0,0);

				% Calculate r-squared
				r21=rsquared(yu1,yp1);
				r22=rsquared(yu2,yp2);
				r23=rsquared(yd,yp3);

				% Convert to polar
				[theta1,rho1]=cart2pol(xp1,yp1);
				[theta2,rho2]=cart2pol(xp2,yp2);
				[theta3,rho3]=cart2pol(xp3,yp3);

				% Find maximum radii and thetas for those radii
				[mrho1,ix1]=max(rho1);
				[mrho2,ix2]=max(rho2);
				[mrho3,ix3]=max(rho3);
				theta1=theta1(ix1);
				theta2=theta2(ix2);
				theta3=theta3(ix3);

				% Calculate interlink angles
				[e1,e2,ba,split,d1,d2]=interangle(theta1,theta2,theta3);
				e1=rad2deg(e1);
				e2=rad2deg(e2);
				ba=rad2deg(ba);

				% Determine order of streams and define handedness
				% sensu James & Krumbein, 1969
				so1=max(so(us{ii,1}));
				so2=max(so(us{ii,2}));
				% so1 should be greater than or equal to so2 as output
				% from nnodesupdown function
				if so1>so2 
					if d1==1 & d2==2
						hnd=1;
					elseif d1==2 & d2==1
						hnd=2;
					elseif d1==1 & d2==1 & e1>e2 
						hnd=1;
					elseif d1==1 & d2==1 & e1<e2 
						hnd=2;
					elseif d1==2 & d2==2 & e1>e2
						hnd=2;
					elseif d1==2 & d2==2 & e1<e2
						hnd=1;
					else
						hnd=3;
					end
				else
					hnd=3;
				end

				% Convert direction angles to cardinal
				card1=mod(-90-rad2deg(theta1),360);
				card2=mod(-90-rad2deg(theta2),360);
				card3=mod(-90-rad2deg(theta3),360);

				link_angle(ii,:)=[S.x(con{ii}) S.y(con{ii}) ba e1 e2 split d1 d2 hnd so1 so2];
				strm_dir(ii,:)=[card1 card2 card3];
				fit_streams(ii,:)=[e1_dst numel(xu1) r21 e2_dst numel(xu2) r22 ds_dst numel(xd) r23];
			else
				link_angle(ii,:)=[S.x(con{ii}) S.y(con{ii}) NaN NaN NaN NaN NaN NaN NaN NaN NaN];
				pred_angle(ii,:)=[NaN NaN];
				strm_dir(ii,:)=[NaN NaN NaN];
				fit_streams(ii,:)=[NaN NaN NaN NaN NaN NaN NaN NaN NaN];

				switch method
				case {'area','slope'}
					pred_angle(ii,:)=[NaN NaN];
				case 'both'
					pred_angle(ii,:)=[NaN NaN NaN NaN];
				end
			end
		end

		% make shape file
		makejunctionshape(wdir,file_name_prefix,link_angle,pred_angle,strm_dir,fit_streams,method,1);

		% generate table
		junctions=makejunctiontable(link_angle,pred_angle,strm_dir,fit_streams,method);

		% ouptut
		out_mat=fullfile(wdir,[file_name_prefix '_junctions.mat']);
		save(out_mat,'IX','junctions','-v7.3');
		writetable(junctions,fullfile(wdir,[file_name_prefix '_junctions.txt']));
	else
		num_n=numel(nlist);

		% Preallocate cell array
		jcell=cell(2,num_n);

		for jj=1:num_n
			nOI=nlist(jj);

			% Preallocate arrays
			la=zeros(num_con,11);
			sd=zeros(num_con,3);
			fs=zeros(num_con,9);

			switch method
			case {'area','slope'}
				pa=zeros(num_con,2);
			case 'both'
				pa=zeros(num_con,4);
			end			

			for ii=1:num_con
				if vld(ii)
					% Extract ix lists and clip
					conix=con{ii};
					us1ix=us{ii,1};
					us2ix=us{ii,2};
					dsix=ds{ii,1};

					if numel(us1ix)>nOI
						us1ix=us1ix(1:nOI);
					end

					if numel(us2ix)>nOI
						us2ix=us2ix(1:nOI);
					end

					if numel(dsix)>nOI
						dsix=dsix(1:nOI);
					end

					% Extract x y node lists
					% Confluence point
					xc=xl(conix); yc=yl(conix);
					% Upstream links
					xu1=xl(us1ix); yu1=yl(us1ix);
					xu2=xl(us2ix); yu2=yl(us2ix);
					% Downstream link
					xd=xl(dsix); yd=yl(dsix);

					% Translate so all links originate from confluence 
					xu1=xu1-xc; xu2=xu2-xc;
					yu1=yu1-yc; yu2=yu2-yc;
					xd=xd-xc; yd=yd-yc;

					%% Depending on orientation of links, doing a linear fit
					%  on the link may produce spurious results. To improve 
					%  the result, each link is rotated towards horizontal, 
					%  using the mean orientation of the link to find the angle
					%  of rotation

					% Find mean angle of stream points from stream orientation
					theta1ap=median(atan2(yu1,xu1));
					theta2ap=median(atan2(yu2,xu2));
					theta3ap=median(atan2(yd,xd));

					% Extract and calculate stream distances
					e1_dst=max(dst(us1ix))-dst(conix);
					e2_dst=max(dst(us2ix))-dst(conix);
					ds_dst=dst(conix)-min(dst(dsix));

					% Calculate predicted angles 
					switch method
					case 'area'
						% Howard 1971 area method
						ix1=S.IXgrid(us1ix);
						ix2=S.IXgrid(us2ix);
						da1=max(A.Z(ix1).*A.cellsize^2);
						da2=max(A.Z(ix2).*A.cellsize^2);
						e1p=acos(((da1+da2)/da1)^-ref_theta);
						e2p=acos(((da1+da2)/da2)^-ref_theta);
						% Convert to degrees
						e1p=rad2deg(e1p);
						e2p=rad2deg(e2p);
						% Store out
						pa(ii,:)=[e1p e2p];
					case 'slope'
						% Howard 1971 slope method
						ix1=S.IXgrid(us1ix);
						ix2=S.IXgrid(us2ix);
						ix3=S.IXgrid(dsix);
						sl1=mean(G.Z(ix1));
						sl2=mean(G.Z(ix2));
						sl3=mean(G.Z(ix3));
						r1=sl3/sl1;
						r2=sl3/sl2;
						if r1>1
							e1p=0;
						else
							e1p=acos(r1);
						end

						if r2>1
							e2p=0;
						else
							e2p=acos(r2);
						end
						% Convert to degrees
						e1p=rad2deg(e1p);
						e2p=rad2deg(e2p);
						% Store out
						pa(ii,:)=[e1p e2p];
					case 'both'
						% Howard 1971 area method
						ix1=S.IXgrid(us1ix);
						ix2=S.IXgrid(us2ix);
						da1=max(A.Z(ix1).*A.cellsize^2);
						da2=max(A.Z(ix2).*A.cellsize^2);
						e1Ap=acos(((da1+da2)/da1)^-ref_theta);
						e2Ap=acos(((da1+da2)/da2)^-ref_theta);
						% Convert to degrees
						e1Ap=rad2deg(e1Ap);
						e2Ap=rad2deg(e2Ap);	

						% Howard 1971 slope method
						ix3=S.IXgrid(dsix);
						sl1=mean(G.Z(ix1));
						sl2=mean(G.Z(ix2));
						sl3=mean(G.Z(ix3));
						r1=sl3/sl1;
						r2=sl3/sl2;
						if r1>1
							e1Sp=0;
						else
							e1Sp=acos(r1);
						end

						if r2>1
							e2Sp=0;
						else
							e2Sp=acos(r2);
						end
						% Convert to degrees
						e1Sp=rad2deg(e1Sp);
						e2Sp=rad2deg(e2Sp);	
						% Store out
						pa(ii,:)=[e1Ap e2Ap e1Sp e2Sp];
					end

					% Rotate all points to horizontal
					[xu1r,yu1r]=rotcoord(xu1,yu1,-theta1ap,0,0); 
					[xu2r,yu2r]=rotcoord(xu2,yu2,-theta2ap,0,0);	
					[xdr,ydr]=rotcoord(xd,yd,-theta3ap,0,0);

					% Simple approximation of linear segment on
					% rotated positions
					a1=xu1r\yu1r;
					a2=xu2r\yu2r; 
					a3=xdr\ydr;

					% Find projected y positions
					yp1r=xu1r*a1;
					yp2r=xu2r*a2;
					yp3r=xdr*a3;

					% Rotate back
					[xp1,yp1]=rotcoord(xu1r,yp1r,theta1ap,0,0);
					[xp2,yp2]=rotcoord(xu2r,yp2r,theta2ap,0,0);
					[xp3,yp3]=rotcoord(xdr,yp3r,theta3ap,0,0);

					% Calculate r-squared
					r21=rsquared(yu1,yp1);
					r22=rsquared(yu2,yp2);
					r23=rsquared(yd,yp3);

					% Convert to polar
					[theta1,rho1]=cart2pol(xp1,yp1);
					[theta2,rho2]=cart2pol(xp2,yp2);
					[theta3,rho3]=cart2pol(xp3,yp3);

					% Find maximum radii and thetas for those radii
					[mrho1,ix1]=max(rho1);
					[mrho2,ix2]=max(rho2);
					[mrho3,ix3]=max(rho3);
					theta1=theta1(ix1);
					theta2=theta2(ix2);
					theta3=theta3(ix3);

					% Calculate interlink angles
					[e1,e2,ba,split,d1,d2]=interangle(theta1,theta2,theta3);
					e1=rad2deg(e1);
					e2=rad2deg(e2);
					ba=rad2deg(ba);

					% Determine order of streams and define handedness
					% sensu James & Krumbein, 1969
					so1=max(so(us1ix));
					so2=max(so(us2ix));
					% so1 should be greater than or equal to so2 as output
					% from nnodesupdown function
					if so1>so2 
						if d1==1 & d2==2
							hnd=1;
						elseif d1==2 & d2==1
							hnd=2;
						elseif d1==1 & d2==1 & e1>e2 
							hnd=1;
						elseif d1==1 & d2==1 & e1<e2 
							hnd=2;
						elseif d1==2 & d2==2 & e1>e2
							hnd=2;
						elseif d1==2 & d2==2 & e1<e2
							hnd=1;
						else
							hnd=3;
						end
					else
						hnd=3;
					end

					% Convert direction angles to cardinal
					card1=mod(-90-rad2deg(theta1),360);
					card2=mod(-90-rad2deg(theta2),360);
					card3=mod(-90-rad2deg(theta3),360);

					la(ii,:)=[S.x(conix) S.y(conix) ba e1 e2 split d1 d2 hnd so1 so2];
					sd(ii,:)=[card1 card2 card3];
					fs(ii,:)=[e1_dst numel(xu1) r21 e2_dst numel(xu2) r22 ds_dst numel(xd) r23];
				else
					la(ii,:)=[S.x(con{ii}) S.y(con{ii}) NaN NaN NaN NaN NaN NaN NaN NaN NaN];
					sd(ii,:)=[NaN NaN NaN];
					fs(ii,:)=[NaN NaN NaN NaN NaN NaN NaN NaN NaN];

					switch method
					case {'slope','area'}
						pa(ii,:)=[NaN NaN];
					case 'both'
						pa(ii,:)=[NaN NaN NaN NaN];
					end
				end
			end

			makejunctionshape(wdir,file_name_prefix,la,pa,sd,fs,method,jj);

			jcell{1,jj}=fit_distance(jj);
			junctions=makejunctiontable(la,pa,sd,fs,method); jcell{2,jj}=junctions;

			% ouptut
			out_mat=fullfile(wdir,[file_name_prefix '_junctions_' num2str(jj) '.mat']);
			save(out_mat,'IX','junctions','-v7.3');
			writetable(junctions,fullfile(wdir,[file_name_prefix '_junctions_' num2str(jj) '.txt']));

		end	% end nodes for

		% Calculate sensitivity to fit distance
		jamat=zeros(num_con,num_n);
		e1mat=zeros(size(jamat));
		e2mat=zeros(size(jamat));
		for jj=1:num_n
			J=jcell{2,jj};
			jamat(:,jj)=J.junction_angle;
			e1mat(:,jj)=J.e1_obs_angle;
			e2mat(:,jj)=J.e2_obs_angle;
		end

		% Calculate mean and standard deviation
		mean_junction_angle=mean(jamat,2);
		mean_e1_angle=mean(e1mat,2);
		mean_e2_angle=mean(e2mat,2);
		std_junction_angle=std(jamat,1,2);
		std_e1_angle=std(e1mat,1,2);
		std_e2_angle=std(e2mat,1,2);

		if strcmp(method,'slope') | strcmp(method,'both')
			e1pmat=zeros(size(jamat));
			e2pmat=zeros(size(jamat));
			for jj=1:num_n
				J=jcell{2,jj};
				e1pmat(:,jj)=J.e1_Spred_angle;
				e2pmat(:,jj)=J.e2_Spred_angle;
			end

			mean_e1Sp_angle=mean(e1pmat,2);
			std_e1Sp_angle=std(e1pmat,1,2);
			mean_e2Sp_angle=mean(e2pmat,2);
			std_e2Sp_angle=std(e2pmat,1,2);	

			MJ=table(mean_junction_angle,std_junction_angle,...
				mean_e1_angle,std_e1_angle,mean_e1Sp_angle,std_e1Sp_angle,...
				mean_e2_angle,std_e2_angle,mean_e2Sp_angle,std_e2Sp_angle);						
		else
			MJ=table(mean_junction_angle,std_junction_angle,...
				mean_e1_angle,std_e1_angle,...
				mean_e2_angle,std_e2_angle);	
		end

		writetable(MJ,fullfile(wdir,[file_name_prefix '_meanjunctions.txt']));
	end % end multi if
end

function [r2]=rsquared(d,pred);
	sstot=sum((d-mean(d)).^2);
	ssres=sum((d-pred).^2);
	r2=1-(ssres/sstot);
end

function [IX]=nnodesupdown(S,A,n,so,verbose)
	% Find indices of a specified number of nodes up and downstream of a confluence
	% Output IX will be a cell array containing indices. All indices refer to positions
	% in the node attributed list of the provided STREAMobj. First column is index of confluence
	% Second column are indices of n nodes downstream of confluence (not including the confluence).
	% Additional columns are n nodes upstream of confluence. Nominally, there should be two columns 
	% beyond the second column, but if there any confluences where more than 2 streams meet, then
	% there will be additional columns (which will be empty for all confluences that only have
	% two links upstream).
	%
	% In the downstream direction, if the search reaches the end of a stream (outlet) or 
	% a b-confluence point, then there will be less than n nodes downstream of this particular confluence.
	% If a confluence coincides with an outlet, then the downstream nodes cell will be an empty array.
	%
	% In the upstream direction, if the search reaches the end of a stream (channel head) or another
	% confluence, there will be less than n nodes upstream of this particular confluence (in the
	% relevant link).
	%
	% Upstream links are sorted so that tributary 1 has higher shreve order or higher drainage area
	% if shreve order of tributary 1 and 2 is the same.

	if verbose
		disp('Finding confluence points...')
	end

	% Calculate drainage area and shreve order
	da=A.Z(S.IXgrid).*A.cellsize^2;

	% Extract confluences
	cix=streampoi(S,'confluences','ix');
	clo=streampoi(S,'bconfluences','logical');

	% Find index of confluences within nal
	cixnal=find(ismember(S.IXgrid,cix));

	%% Find downstream nodes
	IXDS=cell(numel(cixnal),2);
	% Loop through confluences and find n number of nodes downstream
	if verbose
		w1=waitbar(0,'Finding nodes downstream of confluences...');
	end	

	for ii=1:numel(cixnal)
		% Position of confluence in giver list
		ix=find(S.ix==cixnal(ii));

		ixl=zeros(n,1);
		jj=1;
		while jj<=n & ~isempty(ix)
			ixl(jj)=S.ixc(ix);

			% Stop if point is b-confluence
			if clo(ixl(jj))
				break 
			end

			ix=find(S.ix==ixl(jj));
			jj=jj+1;
		end
		ixl(ixl==0)=[];

		IXDS{ii,1}=cixnal(ii);
		IXDS{ii,2}=ixl;

		if verbose
			waitbar(ii/numel(cixnal));
		end
	end

	if verbose
		close(w1);
	end

	%% Find upstream nodes

	% Preallocate cell array of most likely size (may expand columns if there are more
	%	than one tributary at a confluence)
	IXUS=cell(numel(cixnal),2);
	% Loop through confluences and find n number of nodes downstream
	if verbose
		w1=waitbar(0,'Finding nodes upstream of confluences...');
	end	

	for ii=1:numel(cixnal)
		% Position of confluence in receiver list
		lix=find(S.ixc==cixnal(ii));
		% Determine number of links
		num_links=numel(lix);
		sort_mat=zeros(num_links,2);
		for jj=1:num_links
			ixl=zeros(n,1);
			ix=lix(jj);
			kk=1;
			while kk<=n & ~isempty(ix) & numel(ix)==1;
				ixl(kk)=S.ix(ix);
				ix=find(S.ixc==ixl(kk));
				kk=kk+1;
			end
			ixl(ixl==0)=[];
			IXUS{ii,jj}=ixl;
			sort_mat(jj,1)=max(so(ixl));
			sort_mat(jj,2)=max(da(ixl));
		end
		% Sort links
		[~,six]=sortrows(sort_mat,'descend');
		% IXUS(ii,:)=IXUS(ii,six);	
		IXUS(ii,1:numel(six))=IXUS(ii,six);	

		if verbose
			waitbar(ii/numel(cixnal));
		end		
	end	

	if verbose
		close(w1);
	end

	% Join into master cell array
	IX=horzcat(IXDS,IXUS);
end

function [e1,e2,ba,split,d1,d2]=interangle(t1,t2,t3)
	% Assumes all angles are as returned from
	% cart2pol between -pi and pi
	%
	% t1 is tributary 1 theta
	% t2 is tributary 2 theta
	% t3 is downstream theta
	% 
	% e1 is angle between upstream projection and t1
	% e2 is angle between upstream projection and t2
	% ba is angle between t1 and t2
	% split is a logical indicating whether the upstream
	%	projection is between t1 and t2 (true) or not (false)
	% d1 is direction from t1 to upstream projection, CW is 1
	%	CCW is 2
	% d2 is direction from t2 to upstream direction, CW is 1
	%	CCW is 2

	% Rotate entire system so that downstream link always
	% points toward 0 and upstream projection is pi
	t1r=t1-t3;
	t2r=t2-t3;
	t3r=t3-t3;

	% More efficient way to remap rotate angles into range of
	% -pi to pi
	t1r=atan2(sin(t1r),cos(t1r));
	t2r=atan2(sin(t2r),cos(t2r));

	% Calculate angles with respect upstream projection of downstream
	% reach and interlink angle. Also determine whether the upstream
	% projection is between the two tributaries
	if t1r<0 && t2r>=0
		e1=pi-abs(t1r);
		e2=pi-t2r;
		ba=e1+e2;
		split=true;
		d1=2;
		d2=1;
	elseif t1r>=0 && t2r<0
		e1=pi-t1r;
		e2=pi-abs(t2r);
		ba=e1+e2;
		split=true;
		d1=1;
		d2=2;
	elseif t1r>=0 && t2r>=0
		e1=pi-t1r;
		e2=pi-t2r;
		ba=abs(t1r-t2r);
		split=false;
		d1=1;
		d2=1;
	elseif t1r<0 && t2r<0
		e1=pi-abs(t1r);
		e2=pi-abs(t2r);
		ba=abs(abs(t1r)-abs(t2r));
		split=false;
		d1=2;
		d2=2;
	end
end

function [n_x,n_y]=rotcoord(x,y,theta,x0,y0)
	n_x=(x-x0).*cos(theta)-(y-y0).*sin(theta);
	n_y=(x-x0).*sin(theta)+(y-y0).*cos(theta);
end

function [T]=makejunctiontable(la,pa,cr,fs,method)

	num_juncs=numel(la(:,1));

	junction_number=1:num_juncs;
	junction_number=junction_number(:);

	junction_x=la(:,1);
	junction_y=la(:,2);

	junction_angle=la(:,3);
	e1_obs_angle=la(:,4);
	e2_obs_angle=la(:,5);

	e1_shreve=la(:,10);
	e2_shreve=la(:,11);

	e1_direction=cr(:,1);
	e2_direction=cr(:,2);
	e3_direction=cr(:,3);

	e1_distance=fs(:,1);
	e2_distance=fs(:,4);
	e3_distance=fs(:,7);

	e1_num=fs(:,2);
	e2_num=fs(:,5);
	e3_num=fs(:,8);

	split=logical(zeros(num_juncs,1));
	e1_rotation=categorical(zeros(num_juncs,1));
	e2_rotation=categorical(zeros(num_juncs,1));
	handedness=categorical(zeros(num_juncs,1));
	e1_R2=zeros(num_juncs,1);
	e2_R2=zeros(num_juncs,1);
	e3_R2=zeros(num_juncs,1);
	for ii=1:num_juncs

		if la(ii,6)==1
			split(ii,1)=true;
		else
			split(ii,1)=false;
		end

		if la(ii,7)==1
			e1_rotation(ii)='CW';
		elseif la(ii,7)==2
			e1_rotation(ii)='CCW';
		else
			e1_rotation(ii)='Undefined';
		end

		if la(ii,8)==1
			e2_rotation(ii)='CW';
		elseif la(ii,8)==2
			e2_rotation(ii)='CCW';
		else
			e2_rotation(ii)='Undefined';
		end

		if la(ii,9)==1
			handedness(ii)='Right';
		elseif la(ii,9)==2
			handedness(ii)='Left';
		else
			handedness(ii)='Undefined';
		end

		E1_r2=fs(ii,3);
		E2_r2=fs(ii,6);
		E3_r2=fs(ii,9);

		if ~isnan(E1_r2) & ~isinf(E1_r2)
			e1_R2(ii,1)=E1_r2;
		else
			e1_R2(ii,1)=0;
		end

		if ~isnan(E2_r2) & ~isinf(E2_r2)
			e2_R2(ii,1)=E2_r2;
		else
			e2_R2(ii,1)=0;
		end

		if ~isnan(E3_r2) & ~isinf(E3_r2)
			e3_R2(ii,1)=E3_r2;
		else
			e3_R2(ii,1)=0;
		end
	end

	switch method
	case 'area'
		e1_Apred_angle=pa(:,1);	
		e2_Apred_angle=pa(:,2);

		T=table(junction_number,junction_x,junction_y,junction_angle,handedness,split,...
			e1_obs_angle,e1_Apred_angle,e1_rotation,e1_shreve,e1_direction,e1_distance,e1_num,e1_R2,...
			e2_obs_angle,e2_Apred_angle,e2_rotation,e2_shreve,e2_direction,e2_distance,e2_num,e2_R2,...
			e3_direction,e3_distance,e3_num,e3_R2);

	case 'slope'
		e1_Spred_angle=pa(:,1);	
		e2_Spred_angle=pa(:,2);

		T=table(junction_number,junction_x,junction_y,junction_angle,handedness,split,...
			e1_obs_angle,e1_Spred_angle,e1_rotation,e1_shreve,e1_direction,e1_distance,e1_num,e1_R2,...
			e2_obs_angle,e2_Spred_angle,e2_rotation,e2_shreve,e2_direction,e2_distance,e2_num,e2_R2,...
			e3_direction,e3_distance,e3_num,e3_R2);		

	case 'both'
		e1_Apred_angle=pa(:,1);	
		e2_Apred_angle=pa(:,2);
		e1_Spred_angle=pa(:,3);	
		e2_Spred_angle=pa(:,4);

		T=table(junction_number,junction_x,junction_y,junction_angle,handedness,split,...
			e1_obs_angle,e1_Apred_angle,e1_Spred_angle,e1_rotation,e1_shreve,e1_direction,e1_distance,e1_num,e1_R2,...
			e2_obs_angle,e2_Apred_angle,e2_Spred_angle,e2_rotation,e2_shreve,e2_direction,e2_distance,e2_num,e2_R2,...
			e3_direction,e3_distance,e3_num,e3_R2);	
	end

end

function makejunctionshape(wdir,fnp,la,pa,cr,fs,method,num_shp)
	ms=struct;

	all_nums=1:numel(pa(:,1));

	idx=~isnan(pa(:,1));


	la=la(idx,:);
	pa=pa(idx,:);
	cr=cr(idx,:);
	fs=fs(idx,:);
	all_nums=all_nums(idx);

	split=logical(la(:,6));

	num=size(la,1);
	for ii=1:num
		ms(ii,1).Geometry='Point';
		ms(ii,1).X=double(la(ii,1));
		ms(ii,1).Y=double(la(ii,2));
		ms(ii,1).Confl_ID=double(all_nums(ii));
		ms(ii,1).Junc_Ang=double(la(ii,3));
		if split(ii)
			ms(ii,1).Split='true';
		else
			ms(ii,1).Split='false';
		end		

		ms(ii,1).E1_Ang=double(la(ii,4));

		switch method 
		case 'area'
			ms(ii,1).E1_APred=double(pa(ii,1));
		case 'slope'
			ms(ii,1).E1_SPred=double(pa(ii,1));
		case 'both'
			ms(ii,1).E1_APred=double(pa(ii,1));
			ms(ii,1).E1_SPred=double(pa(ii,3));
		end			

		ms(ii,1).E2_Ang=double(la(ii,5));

		switch method 
		case 'area'
			ms(ii,1).E2_APred=double(pa(ii,2));
		case 'slope'
			ms(ii,1).E2_SPred=double(pa(ii,2));
		case 'both'
			ms(ii,1).E2_APred=double(pa(ii,2));
			ms(ii,1).E2_SPred=double(pa(ii,4));
		end	

		if la(ii,7)==1
			ms(ii,1).E1_Rot='CW';
		else
			ms(ii,1).E1_Rot='CCW';
		end

		if la(ii,8)==1
			ms(ii,1).E2_Rot='CW';
		else
			ms(ii,1).E2_Rot='CCW';
		end

		if la(ii,9)==1
			ms(ii,1).Handedness='R';
		elseif la(ii,9)==2
			ms(ii,1).Handedness='L';
		else
			ms(ii,1).Handedness='NA';
		end

		ms(ii,1).E1_Dir=double(cr(ii,1));
		ms(ii,1).E2_Dir=double(cr(ii,2));
		ms(ii,1).DS_Dir=double(cr(ii,3));

		ms(ii,1).E1_Dst=double(fs(ii,1));
		ms(ii,1).E2_Dst=double(fs(ii,4));
		ms(ii,1).DS_Dst=double(fs(ii,7));

		e1_r2=fs(ii,3);
		e2_r2=fs(ii,6);
		e3_r2=fs(ii,9);

		if ~isnan(e1_r2) & ~isinf(e1_r2)
			ms(ii,1).E1_R2=double(e1_r2);
		else
			ms(ii,1).E1_R2=double(0);
		end

		if ~isnan(e2_r2) & ~isinf(e2_r2)
			ms(ii,1).E2_R2=double(e2_r2);
		else
			ms(ii,1).E2_R2=double(0);
		end

		if ~isnan(e3_r2) & ~isinf(e3_r2)
			ms(ii,1).DS_R2=double(e3_r2);
		else
			ms(ii,1).DS_R2=double(0);
		end

	end

	shp_name=[fnp '_junction_angle_' num2str(num_shp) '.shp'];
	shp_fn_name=fullfile(wdir,shp_name);
	shapewrite(ms,shp_fn_name);
end