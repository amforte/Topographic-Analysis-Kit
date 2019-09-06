function cmpAutoKsnProfiler(wdir,MatFile,thresh_ratio,varargin)
	%	Prototype function for automated ksn fits similar to what is done with KsnProfiler. Function uses the MATLAB
	%	'ischange' function to find breaks in binned ksn values. The sensitivity of the function, i.e. how small a 
	%	change in binned ksn value is interpreted as a breakpoint, is controlled by the threshold_ratio parameter. 
	%	Users can explore the affect of different values of threshold ratio, along with reference concavity and 
	%	segment length, with the optional 'explore_param' option.
	%
	% Required Inputs:
	%	wdir - full path of working directory
	%	MatFile - Full path of matfile output from either 'cmpMakeStreams' or the name of a single basin mat file from 'cmpProcessRiverBasins'
	%	thresh_ratio - value between 0 and 1 that controls how sensitive the function is to change in the running mean of binned ksn vs distance.
	%		 Values close to 0 result in a lower threshold and more break points. 
	%
	% Optional Inputs:
	%	prev_param [] - option to use the parameters from a previous run of 'AutoKsnProfiler' stored in the output '*_params.mat' where the prefix
	%		is controlled by the 'file_name_prefix' parameter. This will override any other user provided or default values for 'thresh_ratio', 
	%		'ref_concavity', and 'segment_length'. You can still run 'explore_param' to subsequently change the values stored in the input to 
	%		'prev_param', but the idea behind this is to make it easy to run the AutoKsnProfiler on multiple datasets with the exact same options.
	% 	ref_concavity [0.50] - reference concavity (as a positive value) for calculating ksn
	%	calc_concavity [false] - logical flag to turn on calculation of a best fit concavity for each stream segment. This will slow down the function
	%		considerably so it's only recommended that you set this to true if you are explicitly interested in the concavity of individual segments.
	%	segment_length [1000] - length in map units over which to bin ksn and distance data over which the mean ksn values will be fit.
	%		Also used as the length scale over which to average values when producing the output map structure
	%	plot_example [false] - logical flag to turn on (true) a plot of results in profile view a single stream in the network. By default, the 
	%		selected stream will be the longest trunk stream in the provided network. You can override by this providing a valid input to 'channeloi'.
	%		Provided as a static way to visualize effect of the choice of 'thresh_ratio', 'ref_concavity', and 'segment_length' on the result. See
	%		'explore_param' for a more dynamic way to explore the effect of these values.
	%	explore_param [false] - logical flag to turn on (true) a plot of results in profile view a single stream in the network and the option to change
	%		the values for 'thresh_ratio', 'ref_concavity', and 'segment_length' interactively. By default, the selected stream will be the longest trunk
	%		stream in the provided network. You can override by this providing a valid input to 'channeloi'.
	%	channeloi [] - 1 x 2 array containing the x-y coordinates of a channel head of a stream to use as the selected stream if either 'plot_example' or
	%		'explore_param' is set to true. If neither of those are true, this input has no effect on the behavior of the code.
	%	conditioned_DEM [] - option to provide a hydrologically conditioned DEM for use in this function, expects the mat file as saved by 'cmpConditionDEM'
	%		See 'cmpConditionDEM' function for options for making a hydrological conditioned DEM. If no input is provided the code defaults to using the 
	%		mincosthydrocon function.
	%	interp_value [0.1] - value (between 0 and 1) used for interpolation parameter in mincosthydrocon (not used if user provides a conditioned DEM)
	%	file_name_prefix ['auto'] - name prefix for the shapefiles, textfiles, and matfiles to be export, must have no spaces to be a valid name for ArcGIS and 
	%		should NOT include any file extensions.
	%		
	%
    % Examples if running for the command line, minus OS specific way of calling main TAK function:
    %	AutoKsnProfiler /path/to/wdir Topo.mat 0.5
    %	AutoKsnProfiler /path/to/wdir Topo.mat 0.5 prev_param auto_params.mat
    %	AutoKsnProfiler /path/to/wdir Topo.mat 0.25 ref_concavity 0.4 explore_param true
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 09/06/19 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if isdeployed
		if ~isempty(varargin)
			varargin=varargin{1};
		end
	end	

	% Parse Inputs
	p = inputParser;         
	p.FunctionName = 'cmpAutoKsnProfiler';
	addRequired(p,'wdir',@(x) ischar(x));
	addRequired(p,'MatFile',@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.mat'))));
	addRequired(p,'thresh_ratio',@(x) isscalar(x) && isnumeric(x));

	addParameter(p,'segment_length',1000,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'ref_concavity',0.50,@(x) isscalar (x) && isnumeric(x));
	addParameter(p,'calc_concavity',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'plot_example',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'interp_value',0.1,@(x) isnumeric(x) && x>=0 && x<=1);
	addParameter(p,'conditioned_DEM',[],@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.mat'))));
	addParameter(p,'explore_param',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'channeloi',[],@(x) isnumeric(x) && size(x,1)==1 && size(x,2)==2);
	addParameter(p,'file_name_prefix','auto',@(x) ischar(x));
	addParameter(p,'prev_param',[],@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.mat'))));

	parse(p,wdir,MatFile,thresh_ratio,varargin{:});
	wdir=p.Results.wdir;
	MatFile=p.Results.MatFile;
	thresh_ratio=p.Results.thresh_ratio;

	segment_length=p.Results.segment_length;
	theta=p.Results.ref_concavity;
	plot_example=p.Results.plot_example;
	iv=p.Results.interp_value;
	DEMc=p.Results.conditioned_DEM;
	explore_param=p.Results.explore_param;
	channeloi=p.Results.channeloi;
	file_name_prefix=p.Results.file_name_prefix;
	ppfile=p.Results.prev_param;
	cc=p.Results.calc_concavity;

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

	% Hydrologically condition dem
	if isempty(DEMc)
		zc=mincosthydrocon(S,DEM,'interp',iv);
		DEMc=GRIDobj(DEM);
		DEMc.Z(DEMc.Z==0)=NaN;
		DEMc.Z(S.IXgrid)=zc;
	else
		load(fullfile(wdir,DEMc),'DEMc');
		if ~validatealignment(DEMc,DEM)
			error('Provided conditioned DEM does not appear to align with the original DEM')
		end
	end

	% Convert file name prefix
	file_name_prefix=fullfile(wdir,file_name_prefix);

	% Load in and overwrite other choices for important parameters if an input is provided to 'prev_param'
	if ~isempty(ppfile)
		load(fullfile(wdir,ppfile),'param')
		thresh_ratio=param.thresh_ratio;
		theta=param.reference_concavity;
		segment_length=param.segment_length;
	end

	% Initiate dialogs to change values
	if explore_param
		% Extract example channel
		if isempty(channeloi)
			ST=trunk(klargestconncomps(S,1));
		else
			chxy=streampoi(S,'channelheads','xy');
			[~,mix]=min(hypot(chxy(:,1)-channeloi(:,1),chxy(:,2)-channeloi(:,2)));
			chix=coord2ind(DEM,chxy(mix,1),chxy(mix,2));
			CHIX=GRIDobj(DEM);
			CHIX.Z(chix)=true;
			ST=modify(S,'downstreamto',CHIX);
		end


		% Run using initial values
		[fig_handle]=PlotFunc(DEMc,FD,ST,A,segment_length,theta,thresh_ratio);
		% Query user for the first time
		qa1=questdlg('Accept Values?','Fit Parameters','Change Values','Accept Values','Accept Values');

		switch qa1
		case 'Accept Values'
			% Close figure and continue with code using the original values
			close(fig_handle);			
		case 'Change Values'
			% Initiate loop state
			cont='y';
			% Start loop
			while strcmp(cont,'y');
				% Query user for new values
				prompt={'Enter new threshold ratio:','Enter new reference concavity:','Enter new segment length:'};
				definput={num2str(thresh_ratio),num2str(theta),num2str(segment_length)};
				vals=inputdlg(prompt,'Define New Values',[1 50],definput);
				% Check outputs and convert to numeric values
				if ~isempty(vals)
					if ~isempty(vals{1})
						thresh_ratio=str2num(vals{1});
					end

					if ~isempty(vals{2})
						theta=str2num(vals{2});
					end

					if ~isempty(vals{3})
						segment_length_temp=str2num(vals{3});
						if segment_length_temp<DEM.cellsize*3;
							warning('Input segment length is below 3 times the cellsize of your DEM, reverting to previous value');
						else
							segment_length=segment_length_temp;
						end
					end
				end

				% Replot with new values
				[fig_handle]=PlotFunc(DEMc,FD,ST,A,segment_length,theta,thresh_ratio);
				% Requery user
				qa2=questdlg('Accept Values?','Fit Parameters','Change Values','Accept Values','Accept Values');

				switch qa2
				case 'Accept Values'
					close(fig_handle);
					cont='n';
				case 'Change Values'
					cont='y';
				end % Sub Switch End
			end %While End
		end %Main Switch End
	end 
	
	% Plot representative example for longest trunk stream in network 
	if plot_example
		if isempty(channeloi)
			ST=trunk(klargestconncomps(S,1));
		else
			chxy=streampoi(S,'channelheads','xy');
			[~,mix]=min(hypot(chxy(:,1)-channeloi(:,1),chxy(:,2)-channeloi(:,2)));
			chix=coord2ind(DEM,chxy(mix,1),chxy(mix,2));
			CHIX=GRIDobj(DEM);
			CHIX.Z(chix)=true;
			ST=modify(S,'downstreamto',CHIX);
		end

		[fig_handle]=PlotFunc(DEMc,FD,ST,A,segment_length,theta,thresh_ratio);
	end

	% Store parameters for output incase they were changed
	param=struct; 
	param.thresh_ratio=thresh_ratio;
	param.reference_concavity=theta;
	param.segment_length=segment_length;


	w1=waitbar(0,'Calculating basic ksn values...');
	% Calculate ksn node-by-node
	g=gradient(S,DEMc);
	a=getnal(S,(A.*A.cellsize^2));
	knal=g./a.^(-theta);
	
	waitbar(0.25/5,w1,'Splitting stream network...');
	% Break stream network by channelheads
	[SC,LOCS]=STREAMobj2cell(S,'channelheads');

	% Bin and find ksn segments and breakpoints
	nalC=cell(size(SC));
	nalR=cell(size(SC));
	nalP=cell(size(SC));
	nalN=cell(size(SC));
	nalCo=cell(size(SC));
	brkX=cell(size(SC));
	brkY=cell(size(SC));
	brkD=cell(size(SC));
	brkZ=cell(size(SC));
	for ii=1:numel(SC)
		prog=(ii/numel(SC))*(3/5)+(1/5);
		waitbar(prog,w1,'Calculating auto ksn values...');
		% Bin average
		[dav,kav,idx]=BinAverage(SC{ii}.distance,knal(LOCS{ii}),segment_length);
		% Find change points
		T=mean(kav)*thresh_ratio;
		[bp,mk]=ischange(kav,'variance','Threshold',T);
		% Convert binned values back to nal
		[nal_temp,brk]=UnBinNal(idx,bp,mk);
		% Calculate residuals
		[nalC{ii},nalR{ii},brkX{ii},brkY{ii},brkD{ii},brkZ{ii},~,~,nalP{ii},nalN{ii},nalCO{ii}]=ResidNal(SC{ii},DEMc,FD,A,theta,nal_temp,brk,cc);
	end

	% Accumulate values
	waitbar(4/5,w1,'Accumulating values...');
	locl=vertcat(LOCS{:});
	nalcl=vertcat(nalC{:});
	nalrl=vertcat(nalR{:});
	nalpl=vertcat(nalP{:});
	nalnl=vertcat(nalN{:});
	if cc
		nalcol=vertcat(nalCO{:});
	end
	brkX=vertcat(brkX{:});
	brkY=vertcat(brkY{:});
	brkD=vertcat(brkD{:});
	brkZ=vertcat(brkZ{:});
	brkPnts=[brkX brkY brkD brkZ];
	brkPnts=unique(brkPnts,'rows');
	nal=accumarray(locl,nalcl,[],@mean);
	res=accumarray(locl,nalrl,[],@mean);
	pos=accumarray(locl,nalpl,[],@mean);
	neg=accumarray(locl,nalnl,[],@mean);
	if cc
		con=accumarray(locl,nalcol,[],@mean);
	end
	tta=ones(size(nal))*theta;

	waitbar(4.5/5,w1,'Generating mapstructure and shape outputs...');
	% Convert nal to a mapstruct
	if cc
		ms=STREAMobj2mapstruct(S,'seglength',segment_length,'attributes',...
			{'fit_ksn' nal @mean 'ksn_neg' neg @mean 'ksn_pos' pos @mean 'resid' res @mean 'theta' tta @mean 'seg_theta' con @mean...
			 'rough_ksn' knal @mean 'up_area' a @mean 'gradient' g @mean 'cut_fill' DEMc-DEM @mean});
	else
		ms=STREAMobj2mapstruct(S,'seglength',segment_length,'attributes',...
			{'fit_ksn' nal @mean 'ksn_neg' neg @mean 'ksn_pos' pos @mean 'resid' res @mean 'theta' tta @mean...
			 'rough_ksn' knal @mean 'up_area' a @mean 'gradient' g @mean 'cut_fill' DEMc-DEM @mean});
	end

	% generate shape
	ksn_name=[file_name_prefix '_ksn.shp'];
	bnd_name=[file_name_prefix '_bounds.shp'];

	shapewrite(ms,ksn_name);

	if ~isempty(brkPnts)
		bnds=struct;
		for jj=1:numel(brkPnts(:,1));
			bnds(jj,1).Geometry='Point';
			bnds(jj,1).X=double(brkPnts(jj,1));
			bnds(jj,1).Y=double(brkPnts(jj,2));
			bnds(jj,1).Dist=double(brkPnts(jj,3));
			bnds(jj,1).Elev=double(brkPnts(jj,4));
		end
		shapewrite(bnds,bnd_name);
	end

	% generate mat files storing the parameters
	mat_name=[file_name_prefix '_params.mat'];
	save(mat_name,'param');

	waitbar(1,w1,'Complete');
	close(w1);
end

function [f]=PlotFunc(DEMc,FD,S,A,segment_length,theta,thresh_ratio)
		gT=gradient(S,DEMc);
		aT=getnal(S,(A.*A.cellsize^2));
		knalT=gT./aT.^(-theta);

		[davT,kavT,idxT]=BinAverage(S.distance,knalT,segment_length);
		TT=mean(kavT)*thresh_ratio;
		[bpT,mkT]=ischange(kavT,'variance','Threshold',TT);
		[nalT,brkT]=UnBinNal(idxT,bpT,mkT);
		[rcnalT,resT,brkTX,brkTY,brkTD,~,zp,z,nalPT,nalNT,~]=ResidNal(S,DEMc,FD,A,theta,nalT,brkT,false);

		f=figure(1);
		clf 
		set(f,'Units','normalized','Position',[0.5 0.1 0.5 0.8],'renderer','painters');	

		subplot(3,1,1)
		hold on 
		xlim([0 max(S.distance)]);
		p1=scatter(davT,kavT,10,'k','filled');
		p2=plotdz(S,rcnalT,'color','r');
		p3=plotdz(S,rcnalT+nalPT,'color','r');
		p4=plotdz(S,rcnalT-nalNT,'color','r');
		set(p3,'LineStyle',':');
		set(p4,'LineStyle',':');
		for ii=1:numel(brkTD)
			plot([brkTD(ii) brkTD(ii)],ylim,':k','LineWidth',0.5);
		end
		xlabel('River Distance (m)');
		ylabel('Binned k_{sn}');
		title(['Threshold Ratio = ' num2str(thresh_ratio) '; Refercence Concavity = ' num2str(theta) '; Segment Length = ' num2str(segment_length) ' m']);
		legend([p1 p2 p3],{'Binned Values','Auto Fit','Uncertainty'},'location','best');
		hold off

		subplot(3,1,2)
		hold on
		xlim([0 max(S.distance)]);
		p1=scatter(S.distance,zp,5,'r');
		p2=plotdz(S,z,'color','k');
		for ii=1:numel(brkTD)
			plot([brkTD(ii) brkTD(ii)],ylim,':k','LineWidth',0.5);
		end
		xlabel('River Distance (m)');
		ylabel('Elevation (m)');
		legend([p1 p2],{'Predicted Elevation','Conditioned Elevation'},'location','best');
		hold off

		subplot(3,1,3)
		hold on
		xlim([0 max(S.distance)]);
		plot(S.distance,zeros(size(S.distance)),'-k');
		p1=scatter(S.distance,resT,5,'r','filled');
		for ii=1:numel(brkTD)
			plot([brkTD(ii) brkTD(ii)],ylim,':k','LineWidth',0.5);
		end
		xlabel('River Distance (m)');
		ylabel('Residual (m)');
		title(['Mean (ABS) Residual = ' num2str(mean(abs(resT))) ' m']);
		legend(p1,{'Predicted - Conditioned'},'location','best');
		hold off

		drawnow
end

function [Xavg,Yavg,idx]=BinAverage(X,Y,bin_size);

	ix=~isnan(X);
	X=X(ix); Y=Y(ix);

	minX=min(X);
	maxX=max(X);

	b=[minX:bin_size:maxX+bin_size];

	try
		[idx]=discretize(X,b);
	catch
		[~,idx]=histc(X,b);
	end

	Xavg=accumarray(idx(:),X,[],@mean);
	Yavg=accumarray(idx(:),Y,[],@mean);
end

function [nal,brk]=UnBinNal(idx,bp,mk)

	num_bins=numel(unique(idx));
	nal=zeros(size(idx));
	brk=zeros(size(idx));

	brk_num=1;

	for ii=1:num_bins
		bidx=idx==ii;
		% nal(bidx)=mk(ii);

		if bp(ii)
			brk(bidx)=brk_num;
			brk_num=brk_num+1;
			nal(bidx)=mk(ii);
		else
			brk(bidx)=0;
			nal(bidx)=mk(ii);
		end
	end
end

function [rcnal,res,brk_x,brk_y,brk_d,brk_z,pred_elev,elev,nalpos,nalneg,nalcon]=ResidNal(S,DEMc,FD,A,theta,nal,brk,cc_flag)
	
	if nnz(brk)>0
		brk_ix=zeros(size(S.x));
		brk_nal_ix(S.ixc)=(nal(S.ix)-nal(S.ixc))~=0;
		brk_d=S.distance(brk_nal_ix);
		brk_ix=S.IXgrid(logical(brk_nal_ix));
		[brk_x,brk_y]=ind2coord(DEMc,brk_ix);
		brk_z=DEMc.Z(brk_ix);

		SP=split(S,brk_ix);

		elev=getnal(S,DEMc);
		c=chitransform(SP,A,'a0',1,'mn',theta);
		[L,nc]=conncomps(SP);

		rcnal=zeros(size(nal));
		nalpos=zeros(size(nal));
		nalneg=zeros(size(nal));
		pred_elev=zeros(size(nal));
		nalcon=zeros(size(nal));

		for jj=1:nc
			idx=L==jj;
			out=ChiSpline(elev(idx),c(idx));
			rcnal(idx)=out.ks;
			nalpos(idx)=out.ks_pos;
			nalneg(idx)=out.ks_neg;
			pred_elev(idx)=c(idx).*rcnal(idx);
			pred_elev(idx)=pred_elev(idx)+min(elev(idx));

			if cc_flag
				W=GRIDobj(DEMc);
				W.Z(S.IXgrid(idx))=true;
				SPT=STREAMobj(FD,W);
				nalcon(idx)=BestFitTheta(S,A,DEMc);
			else 
				nalcon=[];
			end
		end

		res=pred_elev-elev;
	else
		brk_x=[]; brk_y=[]; brk_ix=[]; brk_d=[]; brk_z=[];	
		elev=getnal(S,DEMc);
		c=chitransform(S,A,'a0',1,'mn',theta);

		out=ChiSpline(elev,c);
		rcnal=ones(size(nal))*out.ks;
		nalpos=ones(size(nal))*out.ks_pos;
		nalneg=ones(size(nal))*out.ks_neg;

		pred_elev=(c.*rcnal)+min(elev);

		res=pred_elev-elev;

		if cc_flag
			nalcon=ones(size(nal));
			mnF=BestFitTheta(S,A,DEMc);
			nalcon=nalcon*mnF;
		else
			nalcon=[];
		end
	end
end

function [OUT]=ChiSpline(z,c)
	zabsF=z-min(z);
	chiF=c;

	% Spline generates lots of warnings for small segments
	warning off
	chiS=linspace(0,max(chiF),numel(chiF)).';
	try
		zS=spline(chiF,zabsF,chiS);
	catch
		% cubic spline will fail if segments is nearly straight, skip spline fit
		% in this case to avoid erroring out
		zS=zabsF;
		chiS=chiF;
	end

	OUT=struct;
	try
		ft=fittype('a*x');
		fobj=fit(chiS,zS,ft,'StartPoint',chiS\zS);
		BETA=coeffvalues(fobj);
		BETA_UNC=confint(fobj);
		OUT.ks   = BETA;
		OUT.ks_neg = (BETA)-min(BETA_UNC);
		OUT.ks_pos = max(BETA_UNC)-(BETA);
	catch
		BETA = chiS\(zS);
		OUT.ks   = BETA;
		OUT.ks_neg = 0;
		OUT.ks_pos = 0;
	end

	warning on

end

function [mnF]=BestFitTheta(S,A,DEM)

	z=getnal(S,DEM);
	mnF=fminsearch(@mnfit,0.5);

	function [sqres]=mnfit(mn)
		c=chitransform(S,A,'a0',1,'mn',mn);
		c=c./max(c);
		z=z-min(z);
		z=z./max(z);
		sqres=sum((c-z).^2);
	end

end
