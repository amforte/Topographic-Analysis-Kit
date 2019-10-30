function [ERO,varargout]=EroGrid(DEM,KSN,rel_type,varargin)
	% Usage:
	%	[ERO]=EroGrid(DEM,KSN,C,phi);
	%	[ERO]=EroGrid(DEM,KSN,C,phi,'name',value,...);
	%	[ERO,ERO_P,ERO_M]=EroGrid(DEM,KSN,C,phi,'name',value,...);	
	%
	% Description;
	%	Function to produce an erosion rate map based on an empirical relationship between normalized channel steepness (ksn)
	%	and erosion rate (E) defined in the form of ksn = C*E^phi of using a more complicated stochastic threshold model described
	%	in Lague et al, 2005. This relationship would likely come from a series of catchment averaged cosmogenic erosion rates and 
	%	basin averaged normalized channel steepness. Function can produce an erosion rate map for a unique relationship between ksn 
	%	and E or a series of relationships based on different regions defined in a separate grid provided as the optional 'VAL' input. 
	%	Function is also able to incorporate uncertainty in the interpolated KSN map and/or uncertainty on the fit parameters 
	%	(i.e. C and phi for the simple power law case) to calculate min and max erosion rate estimates based on these uncertainties. 
	%
	% Required Inputs:
	%	DEM - DEM GRIDobj
	%	KSN - GRIDobj of continous ksn (e.g. as produced by KsnChiBatch with product set to 'ksngrid') or a map structure of
	%		ksn (e.g. as produced by KsnChiBatch with product set to 'ksn'). If the KSN GRIDobj was produced alternatively, 
	%		please ensure that it is the same dimensions and pixel size as the provided DEM, i.e. the result of 
	%		validatealignment(DEM,KSN) is true.
	%	rel_type - the type of relationship to use before ksn and E, options are:
	%		'power' - simple power law relationship of the form ksn = C*E^phi
	%		'stochastic_threshold' - stochastic threshold relationship from Lague et al, 2005 as implemented by DiBiase & Whipple,
	%			2011, i.e. equation 10 of DiBiase & Whipple, 2011
	%
	% Optional Inputs Required For Power Law:
	%	C [] - coefficient of power law relationship between ksn and E (ksn = C*E^phi), can be a single value or an array of values. 
	%		If an array of values is provided, it is assumed that these are multiple coefficients that refer to different relationships
	%		based on spatial subsetting that will be defined with inputs to the optional 'edges' and 'VAL'  inputs
	%	phi [] - exponent of power law relationship between ksn and E (ksn = C*E^phi), can be a single value or an array of values. 
	%		If an array of values is provided, it is assumed that these are multiple exponents that refer to different relationships
	%		based on spatial subsetting that will be defined with inputs to the optional 'edges' and 'VAL'  inputs	
	%
	% Optional Inputs Required for Stochastic Threshold Model:
	%	k_e [1e-12] - incision efficiency constant
	%	tau_crit [45] - critical shear stress in pascals
	%	Rb	[1] - mean runoff in mm/day
	%	k [0.5] - climate variability from inverse gamma function
	%	k_w [15] - amplitude factor of the channel width/mean relationship
	%	f [0.08313] - darcy-weisbach friction factor
	%	omega_a [0.55] - downstream scaling exponent between channel width and discharge
	%	omega_s [0.25] - local (at-a-station) scaling exponent between flow width and discharge
	%	alpha_val [2/3] - friction exponent on discharge
	%	beta_val [2/3] - friction exponent on slope
	%	a [1.5] - shear stress exponent
	%	
	%
	% Other Optional Inputs:
	%	radius [5000] - radius for creating a spatially averaged ksn grid if the entry to KSN is a mapstruct
	%	KSNstd - GRIDobj of standard deviation of continous ksn (e.g. as produced by KsnChiBatch with product set to 'ksngrid') to 
	%		incoporate uncertainty in ksn smoothing into the production of erosion rate maps. If you are providing a mapstructure to 
	%		the 'KSN' argument and you wish to use the KSNstd that will be calculated to consider how uncertainty in the ksn value leads to 
	%		uncertainty in the erosion rate map, provide true for this parameter.
	%	C_std - standard deviation / uncertainty on coefficient of power law. If an entry is provided to C_std (and phi_std), then these uncertainties
	%		will be incorporated into the erosion grid outputs. There must be the same number of entires to C_std as there is to 'C'.
	%	phi_std - standard deviation / uncertainty on the exponent of power law. If an entry is provided to phi_std (and C_std), then these uncertainties
	%		will be incorporated into the erosion grid outputs. There must be the same number of entires to phi_std as there is to 'phi'.
	%	VAL [] - GRIDobj defining the regions by which to index the KSN grid to establish different empirical relationships between ksn and E.
	%		For example, VAL might be a precipitation grid if you have reason to believe there are different ksn - E relationships 
	%		depending on precipitation. If you provide an input to VAL, you must also provide an input to 'edges' that define the
	%		bounds on the differnet values within VAL that define the different regions, i.e. VAL and edges must be in the same units and
	%		edges must cover the full range of values within VAL
	%	edges [] - array that define the edges to index the GRIDobj provided to VAL. There must be n+1 entries for n entries to 'C' and 'phi' if using
	%		the 'power' relationship and n entries to 'k_e', 'tau_crit', 'Rb', and 'k' (other parameters are fixed based on provided values).
	%		It's also assumed that entries for these arrays are in the same order as the bins defined by 'edges'.
	%	resample_method ['nearest'] - method to resample the provided VAL GRIDobj if it is a different dimension or pixel size than the input
	%		DEM. Valid options are 'nearest', 'bilinear', and 'bicubic'. Default is 'nearest'.
	%	plot_result [false] - option to turn on a plot of the result
	%
	% Notes:
	%	If the relation type is set to 'power', the function does not explicitly depend on the 'KSN' arugment actually being KSN. E.g., if you 
	%	have established a relationship between erosion rate and local relief in the form rlf = C*E^phi, then this function will work equally well. This is 
	%	not the case for the 'stochastic_threshold' case, where it is required that the 'KSN' arugment is actually channel steepness
	%
	% 	Uncertainties in parameter values are not implemented for the stochastic threshold model due the complicated nature of this equation and that it must
	%	solved numerically.
	%
	%	For the 'power' type of relationship, the units of erosion rate grids will vary depending on the units of erosion rate when you fit the data. For
	%	the 'stochastic_threshold' type, outputs erosion rates are in m/Myr
	%
	% Example:
	%	% Unique relationship between ksn and E
	%	[ERO]=EroGrid(DEM,KSN,'power','C',100,'phi',0.5);
	%
	%	% Three different relationships between ksn and E depending on precipitation, where based on empirical data, ksn = 100*E^(0.5) for
	%	% preciptation between 0 and 1 m/yr, ksn = 316*E^(0.5) for precipitation between 1-2 m/yr, and ksn = 1000*E^(0.5) for precipitation
	%	% between 2-4 m/yr. PRECIP is a GRIDobj of precipitation in m/yr with mininum and maximum values greater than 0 and less than 4, 
	%	% respectively
	%	[ERO]=EroGrid(DEM,KSN,'power','C',[100 316 1000],'phi',[0.5 0.5 0.5],'VAL',PRECIP,'edges',[0 1 2 4]);
	%
	%	% Calculate upper (ERO_P) and lower (ERO_M) bounds on erosion rate map based on the standard deviation of the interpolated ksn grid
	%	[ERO,ERO_P,ERO_M]=EroGrid(DEM,KSN,'power,'C',100,'C',0.5,'KSNstd',KSNstd);
	%
	%	% Calculate upper (ERO_P) and lower (ERO_M) bounds on erosion rate map based on uncertainty in fit parameters
	%	[ERO,ERO_P,ERO_M]=EroGrid(DEM,KSN,'power','C',100,'phi',0.5,'C_std',5,'phi_std',0.01);
	%
	%	% Use the stochastic threshold model and change k_e and the tail of the gamma distribution
	%	[ERO]=EroGrid(DEM,KSN,'stochastic_threshold','k_e',1e-10,'k',0.5);
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 10/30/19 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Parse Inputs
	p = inputParser;         
	p.FunctionName = 'EroGrid';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'KSN',@(x) isa(x,'GRIDobj') | isstruct(x));
	addRequired(p,'rel_type',@(x) ischar(validatestring(x,{'power','stochastic_threshold'})));

	addParameter(p,'C',[],@(x) isnumeric(x));
	addParameter(p,'phi',[],@(x) isnumeric(x));
	addParameter(p,'k_e',1e-12,@(x) isnumeric(x));
	addParameter(p,'tau_crit',45,@(x) isnumeric(x));
	addParameter(p,'Rb',1,@(x) isnumeric(x));
	addParameter(p,'k',0.5,@(x) isnumeric(x));
	addParameter(p,'k_w',15,@(x) isnumeric(x) && isscalar(x));
	addParameter(p,'f',0.08313,@(x) isnumeric(x) && isscalar(x));	
	addParameter(p,'omega_a',0.55,@(x) isnumeric(x) && isscalar(x));
	addParameter(p,'omega_s',0.25,@(x) isnumeric(x) && isscalar(x));
	addParameter(p,'alpha_val',2/3,@(x) isnumeric(x) && isscalar(x));
	addParameter(p,'beta_val',2/3,@(x) isnumeric(x) && isscalar(x));
	addParameter(p,'a',1.5,@(x) isnumeric(x) && isscalar(x));
	addParameter(p,'radius',5000,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'KSNstd',false,@(x) isa(x,'GRIDobj') | islogical(x));
	addParameter(p,'C_std',[],@(x) isnumeric(x) || isempty(x));
	addParameter(p,'phi_std',[],@(x) isnumeric(x)|| isempty(x));
	addParameter(p,'VAL',[],@(x) isa(x,'GRIDobj') || isempty(x));
	addParameter(p,'edges',[],@(x) isnumeric(x) || isempty(x));
	addParameter(p,'resample_method','nearest',@(x) ischar(validatestring(x,{'nearest','bilinear','bicubic'})));
	addParameter(p,'error_type','std',@(x) ischar(validatestring(x,{'std','std_error'})));
	addParameter(p,'plot_result',false,@(x) islogical(x) && isscalar(x));

	parse(p,DEM,KSN,rel_type,varargin{:});
	DEM=p.Results.DEM;
	KSN=p.Results.KSN;
	rel_type=p.Results.rel_type;

	C=p.Results.C;
	phi=p.Results.phi;
	k_e=p.Results.k_e;
	tau_crit=p.Results.tau_crit;
	Rb=p.Results.Rb;
	k=p.Results.k;
	k_w=p.Results.k_w;
	f=p.Results.f;
	omega_a=p.Results.omega_a;
	omega_s=p.Results.omega_s;
	alpha_val=p.Results.alpha_val;
	beta_val=p.Results.beta_val;
	a=p.Results.a;
	radius=p.Results.radius;
	ksn_std=p.Results.KSNstd;
	C_std=p.Results.C_std;
	phi_std=p.Results.phi_std;
	VAL=p.Results.VAL;
	edges=p.Results.edges;
	resample_method=p.Results.resample_method;
	error_type=p.Results.error_type;
	plot_result=p.Results.plot_result;

	switch rel_type
	case 'power'

		% Check for inputs
		if isempty(C) | isempty(phi)
			error('Relationship type is "power", you must provide values for both "C" and "phi"');
		end

		% Basic error check for inputs
		if numel(C) ~= numel(phi)
			error('Number of coefficients provided to "C" must equal the number of exponents in "phi"')
		end

		if numel(C_std) ~= numel(phi_std)
			error('There must be the same number of entries to "C_std" and "phi_std"');
		end

		if ~isempty(C_std) & numel(C_std) ~= numel(C)
			error('There must be the same number of entries to "C_std" and "C"');
		end

		if ~isempty(phi_std) & numel(phi_std) ~= numel(phi)
			error('There must be the same number of entries to "phi_std" and "phi"');
		end

		% Additional error checker
		if ~isempty(edges)
			if isempty(VAL)
				error('Cannot produce variable erosion rate map based on edges without an input for "VAL" that corresponds to the values in edges')
			elseif numel(edges) ~= numel(C)+1
				error('Number of edges is not compatible with the number of coefficients')
			elseif min(edges)> min(VAL.Z(:),[],'omitnan')
				warning('Minimum value of "edges" is greater than minimum value in "VAL", there may be empty areas of the resulting ERO grid')
			elseif max(edges)< max(VAL.Z(:),[],'omitnan')
				warning('Maximum value of "edges" is less than maximum value in "VAL", there may be empty areas of the resulting ERO grid')			
			end
		end

		% Check for compatibility of VAL and DEM (and thus KSN grid);
		if ~isempty(VAL)
			if ~validatealignment(VAL,DEM);
				disp(['Resampling VAL to be the same resolution and dimensions as the input DEM by the ' resample_method ' method']);
				VAL=resample(VAL,DEM,resample_method);
			end
		end

		% Check for form of ksn
		if isstruct(KSN)
			disp('Generating Continuous KSN Grid')
			[KSN,KSNstd]=KsnAvg(DEM,KSN,radius,error_type);
			disp('Continuous KSN Grid is complete')
		end

		% Check whether the user wants to include uncertainty in the erosion rate calculation
		if isa(ksn_std,'GRIDobj')
			ksn_std_flag=true;
			KSNstd=ksn_std;
		elseif islogical(ksn_std) & ksn_std & isa(KSNstd,'GRIDobj')
			ksn_std_flag=true;
		else
			ksn_std_flag=false;
		end

		if ~isempty(C_std)
			fit_unc_flag=true;
		else
			fit_unc_flag=false;
		end

		if fit_unc_flag | ksn_std_flag
			std_flag=true;
		else
			std_flag=false;
		end

		% Begin calculation
		if isempty(edges)
			% Simple ksn to erosion rate conversion
			n=1/phi;
			K=C^(-n);

			ERO=K.*KSN.^n;

			if ksn_std_flag & ~fit_unc_flag
				ERO_P=K.*(KSN+KSNstd).^n;
				ERO_M=K.*(KSN-KSNstd).^n;
			elseif ~ksn_std_flag & fit_unc_flag
				% Note, increasing values of phi and C lead to lower erosion
				% rates for the same ksn, so for the 'minus' parameters
				% they should be where the uncertainties are added
				n_m=1/(phi+phi_std);
				K_m=(C+C_std)^-n_m;
				n_p=1/(phi-phi_std);
				K_p=(C-C_std)^-n_p;

				ERO_P=(K_p).*KSN.^n_p;
				ERO_M=(K_m).*KSN.^n_m;
			elseif ksn_std_flag & fit_unc_flag
				n_m=1/(phi+phi_std);
				K_m=(C+C_std)^-n_m;
				n_p=1/(phi-phi_std);
				K_p=(C-C_std)^-n_p;

				ERO_P=(K_p).*(KSN+KSNstd).^n_p;
				ERO_M=(K_m).*(KSN-KSNstd).^n_m;			
			end
		else
			% ksn to erosion rate based on regions defined by VAL
			ERO=GRIDobj(DEM);

			if ksn_std_flag
				ERO_P=GRIDobj(DEM);
				ERO_M=GRIDobj(DEM);
			end		

			n=1./phi;
			K=C.^(-n);

			num_bins=numel(edges)-1;
			for ii=1:num_bins
				IDX = VAL>=edges(ii) & VAL<edges(ii+1);
				ERO.Z(IDX.Z)=K(ii).*KSN.Z(IDX.Z).^n(ii);

				if ksn_std_flag & ~fit_unc_flag
					ERO_P.Z(IDX.Z)=K(ii).*(KSN.Z(IDX.Z)+KSNstd.Z(IDX.Z)).^n(ii);
					ERO_M.Z(IDX.Z)=K(ii).*(KSN.Z(IDX.Z)-KSNstd.Z(IDX.Z)).^n(ii);
				elseif ~ksn_std_flag & fit_unc_flag
					n_m(ii)=1/(phi(ii)+phi_std(ii));
					K_m(ii)=(C(ii)+C_stdi(ii))^-n_m;
					n_p(ii)=1/(phi(ii)-phi_std(ii));
					K_p(ii)=(C(ii)-C_std(ii))^-n_p;

					ERO_P.Z(IDX.Z)=K_p(ii).*KSN.Z(IDX.Z).^n_p(ii);
					ERO_M.Z(IDX.Z)=K_m(ii).*KSN.Z(IDX.Z).^n_m(ii);
				elseif ksn_std_flag & fit_unc_flag 
					n_m(ii)=1/(phi(ii)+phi_std(ii));
					K_m(ii)=(C(ii)+C_std(ii))^-n_m;
					n_p(ii)=1/(phi(ii)-phi_std(ii));
					K_p(ii)=(C(ii)-C_std(ii))^-n_p;

					ERO_P.Z(IDX.Z)=K_p(ii).*(KSN.Z(IDX.Z)+KSNstd.Z(IDX.Z)).^n_p(ii);
					ERO_M.Z(IDX.Z)=K_m(ii).*(KSN.Z(IDX.Z)-KSNstd.Z(IDX.Z)).^n_m(ii);
				end
			end
		end

		% Set pixels that were NaN in DEM to NaN in ERO
		IDX=GRIDobj(DEM,'logical');
		IDX.Z(isnan(DEM.Z))=true;

		ERO.Z(IDX.Z)=NaN;

		% Perform check to make sure all values are real and set imaginary to NaN
		ERO.Z(imag(ERO.Z)~=0)=NaN;

		if std_flag
			ERO_P.Z(IDX.Z)=NaN;
			ERO_M.Z(IDX.Z)=NaN;

			ERO_P.Z(imag(ERO_P.Z)~=0)=NaN;
			ERO_M.Z(imag(ERO_M.Z)~=0)=NaN;

			varargout{1}=ERO_P;
			varargout{2}=ERO_M;
		end

		if plot_result
			f1=figure(1);
			clf 
			set(f1,'unit','normalized','position',[0.1 0.1 0.8 0.8]);

			ksn_min=min(KSN.Z(:),[],'omitnan');
			ksn_max=max(KSN.Z(:),[],'omitnan');
			ksn_vec=linspace(ksn_min,ksn_max,100);


			sbplt1=subplot(3,3,[1:6]);
			hold on 
			imageschs(DEM,ERO,'colorbarlabel','Erosion Rate');
			disableDefaultInteractivity(sbplt1);
			hold off

			sbplt2=subplot(3,3,[7:9]);
			hold on 
			if fit_unc_flag & isempty(edges)
				E_vec=K.*ksn_vec.^n;
				plot(E_vec,ksn_vec,'-k','LineWidth',2);
				E_vec_m=K_m.*ksn_vec.^n_m;
				E_vec_p=K_p.*ksn_vec.^n_p;
				plot(E_vec_m,ksn_vec,':k','LineWidth',1);
				plot(E_vec_p,ksn_vec,':k','LineWidth',1);
			elseif fit_unc_flag & ~isempty(edges)
				for ii=1:num_bins
					E_vec=K(ii).*ksn_vec.^n(ii);
					plt(ii)=plot(E_vec,ksn_vec,'-','LineWidth',2);
					E_vec_m=K_m(ii).*ksn_vec.^n_m(ii);
					E_vec_p=K_p(ii).*ksn_vec.^n_p(ii);
					plot(E_vec_m,ksn_vec,':','LineWidth',1);
					plot(E_vec_p,ksn_vec,':','LineWidth',1);
					leg{ii}=['Bin ' num2str(ii)];
				end
				legend(plt,leg,'location','best');
			else
				E_vec=K.*ksn_vec.^n;
				plot(E_vec,ksn_vec,'-k','LineWidth',2);
			end
			xlabel('Erosion Rate');
			ylabel('K_{sn}');
			disableDefaultInteractivity(sbplt1);
			hold off
		end

	case 'stochastic_threshold'

		% Basic error check for inputs
		if numel(k_e) ~= numel(tau_crit) | numel(k_e) ~= numel(k) | numel(k_e) ~= numel(Rb)
			error('Number of values provided to "k_e", "tau_crit", "k", and "Rb" must be equal');
		end


		% Additional error checker
		if ~isempty(edges)
			if isempty(VAL)
				error('Cannot produce variable erosion rate map based on edges without an input for "VAL" that corresponds to the values in edges')
			elseif numel(edges) ~= numel(k_e)+1
				error('Number of edges is not compatible with the number of values provided to "k_e", "tau_crit", "k", and "Rb"')
			elseif min(edges)> min(VAL.Z(:),[],'omitnan')
				warning('Minimum value of "edges" is greater than minimum value in "VAL", there may be empty areas of the resulting ERO grid')
			elseif max(edges)< max(VAL.Z(:),[],'omitnan')
				warning('Maximum value of "edges" is less than maximum value in "VAL", there may be empty areas of the resulting ERO grid')			
			end
		end

		% Check for compatibility of VAL and DEM (and thus KSN grid);
		if ~isempty(VAL)
			if ~validatealignment(VAL,DEM);
				disp(['Resampling VAL to be the same resolution and dimensions as the input DEM by the ' resample_method ' method']);
				VAL=resample(VAL,DEM,resample_method);
			end
		end

		% Check for form of ksn
		if isstruct(KSN)
			disp('Generating Continuous KSN Grid')
			[KSN,KSNstd]=KsnAvg(DEM,KSN,radius,error_type);
			disp('Continuous KSN Grid is complete')
		end

		% Check whether the user wants to include uncertainty in the erosion rate calculation
		if isa(ksn_std,'GRIDobj')
			ksn_std_flag=true;
			KSNstd=ksn_std;
		elseif islogical(ksn_std) & ksn_std & isa(KSNstd,'GRIDobj')
			ksn_std_flag=true;
		else
			ksn_std_flag=false;
		end

		% Begin calculation
		if isempty(edges)
			% Determine ksn range
			min_ksn=min(KSN.Z(:),[],'omitnan');
			max_ksn=max(KSN.Z(:),[],'omitnan');
			ERO=GRIDobj(DEM);

			% Numerically integrate across full ksn range
			[E,Ks]=stoch_thresh(min_ksn,max_ksn,k_e,tau_crit,Rb,k,k_w,f,omega_a,omega_s,alpha_val,beta_val,a);

			% Use Ks to discretize KSN grid and populate E grid
			ix=discretize(KSN.Z,Ks);
			w1=waitbar(0,'Populating Erosion Rate Grid...');
			for ii=1:numel(Ks)-1
				E_val=mean([E(ii) E(ii+1)]);
				idx=ix==ii;
				ERO.Z(idx)=E_val;
				waitbar(ii/(numel(Ks)-1));
			end
			close(w1);

			if ksn_std_flag

				KSN_MAX=KSN+KSNstd;
				max_ksn_min=min(KSN_MAX.Z(:),[],'omitnan');
				max_ksn_max=max(KSN_MAX.Z(:),[],'omitnan');

				KSN_MIN=KSN-KSNstd;
				min_ksn_min=min(KSN_MIN.Z(:),[],'omitnan');
				min_ksn_max=max(KSN_MIN.Z(:),[],'omitnan');

				ERO_P=GRIDobj(DEM);
				ERO_M=GRIDobj(DEM);

				% Numerically integrate across full ksn range
				[EMax,KsMax]=stoch_thresh(max_ksn_min,max_ksn_max,k_e,tau_crit,Rb,k,k_w,f,omega_a,omega_s,alpha_val,beta_val,a);
				[EMin,KsMin]=stoch_thresh(min_ksn_min,min_ksn_max,k_e,tau_crit,Rb,k,k_w,f,omega_a,omega_s,alpha_val,beta_val,a);

				ix_max=discretize(KSN_MAX.Z,KsMax);
				ix_min=discretize(KSN_MIN.Z,KsMin);

				w1=waitbar(0,'Populating Erosion Rate Min and Max Grids...');
				for ii=1:numel(Ks)-1
					E_val_max=mean([EMax(ii) EMax(ii+1)]);
					idx_max=ix_max==ii;
					ERO_P.Z(idx_max)=E_val_max;

					E_val_min=mean([EMin(ii) EMin(ii+1)]);
					idx_min=ix_min==ii;
					ERO_M.Z(idx_min)=E_val_min;	
					waitbar(ii/(numel(Ks)-1));				
				end
				close(w1);
			end

		else
			% ksn to erosion rate based on regions defined by VAL
			ERO=GRIDobj(DEM);

			if ksn_std_flag
				ERO_P=GRIDobj(DEM);
				ERO_M=GRIDobj(DEM);
			end		

			num_bins=numel(edges)-1;
			for kk=1:num_bins
				IDX = VAL>=edges(kk) & VAL<edges(kk+1);

				min_ksn=min(KSN.Z(IDX.Z),[],'omitnan');
				max_ksn=max(KSN.Z(IDX.Z),[],'omitnan');	

				% Numerically integrate across full ksn range
				[E,Ks]=stoch_thresh(min_ksn,max_ksn,k_e,tau_crit,Rb,k,k_w,f,omega_a,omega_s,alpha_val,beta_val,a);

				% Use Ks to discretize KSN grid and populate E grid
				KSN_TEMP=KSN;
				KSN_TEMP.Z(~IDX.Z)=NaN;
				ix=discretize(KSN_TEMP.Z,Ks);
				w1=waitbar(0,['Populating Erosion Rate Grid For Bin ' num2str(kk) '...']);
				for ii=1:numel(Ks)-1
					E_val=mean([E(ii) E(ii+1)]);
					idx=ix==ii;
					ERO.Z(idx)=E_val;
					waitbar(ii/(numel(Ks)-1));	
				end							
				close(w1);	

				if ksn_std_flag
					KSN_MAX=KSN+KSNstd;
					max_ksn_min=min(KSN_MAX.Z(IDX.Z),[],'omitnan');
					max_ksn_max=max(KSN_MAX.Z(IDX.Z),[],'omitnan');

					KSN_MIN=KSN-KSNstd;
					min_ksn_min=min(KSN_MIN.Z(IDX.Z),[],'omitnan');
					min_ksn_max=max(KSN_MIN.Z(IDX.Z),[],'omitnan');

					ERO_P=GRIDobj(DEM);
					ERO_M=GRIDobj(DEM);

					% Numerically integrate across full ksn range
					[EMax,KsMax]=stoch_thresh(max_ksn_min,max_ksn_max,k_e,tau_crit,Rb,k,k_w,f,omega_a,omega_s,alpha_val,beta_val,a);
					[EMin,KsMin]=stoch_thresh(min_ksn_min,min_ksn_max,k_e,tau_crit,Rb,k,k_w,f,omega_a,omega_s,alpha_val,beta_val,a);

					KSN_MAX.Z(~IDX.Z)=NaN;
					KSN_MIN.Z(~IDX.Z)=NaN;
					ix_max=discretize(KSN_MAX.Z,KsMax);
					ix_min=discretize(KSN_MIN.Z,KsMin);

					w1=waitbar(0,['Populating Erosion Rate Min and Max Grids For Bin ' num2str(kk) '...']);
					for ii=1:numel(Ks)-1
						E_val_max=mean([EMax(ii) EMax(ii+1)]);
						idx_max=ix_max==ii;
						ERO_P.Z(idx_max)=E_val_max;

						E_val_min=mean([EMin(ii) EMin(ii+1)]);
						idx_min=ix_min==ii;
						ERO_M.Z(idx_min)=E_val_min;	
						waitbar(ii/(numel(Ks)-1));				
					end
					close(w1);
				end
			end
		end

		% Set pixels that were NaN in DEM to NaN in ERO
		IDX=GRIDobj(DEM,'logical');
		IDX.Z(isnan(DEM.Z))=true;

		ERO.Z(IDX.Z)=NaN;

		% Perform check to make sure all values are real and set imaginary to NaN
		ERO.Z(imag(ERO.Z)~=0)=NaN;

		if ksn_std_flag
			ERO_P.Z(IDX.Z)=NaN;
			ERO_M.Z(IDX.Z)=NaN;

			ERO_P.Z(imag(ERO_P.Z)~=0)=NaN;
			ERO_M.Z(imag(ERO_M.Z)~=0)=NaN;

			varargout{1}=ERO_P;
			varargout{2}=ERO_M;
		end

		if plot_result
			f1=figure(1);
			clf 
			set(f1,'unit','normalized','position',[0.1 0.1 0.8 0.8]);

			sbplt1=subplot(3,3,[1:6]);
			hold on 
			imageschs(DEM,ERO,'colorbarlabel','Erosion Rate [m/Myr]');
			disableDefaultInteractivity(sbplt1);
			hold off

			sbplt2=subplot(3,3,[7:9]);
			hold on 
			plot(E,Ks,'-k','LineWidth',2);
			xlabel('Erosion Rate [m/Myr]');
			ylabel('K_{sn}');
			disableDefaultInteractivity(sbplt1);
			hold off
		end



	end

end

function [E,Ks]=stoch_thresh(ksn_min,ksn_max,k_e,tau_crit,Rb,k,k_w,f,omega_a,omega_s,alpha_val,beta_val,a)
	% Numerical integration of Lague et al 2005, equation 16. Adapted from a code originally written
	% by Roman DiBiase

	% Convert Rb to k_q
	k_q=Rb/(24*60*60*10*100);

	% Derived Parameters
	k_t = 0.5*1000*(9.81^(2/3))*(f^(1/3));      % set to 1000 a la Tucker 2004
	y = a*alpha_val*(1-omega_s);                    % gamma exponent
	m = a*alpha_val*(1-omega_a);                    % m in erosion law
	n = a*beta_val;                                 % n in erosion law
	psi_crit = k_e*tau_crit^a;                   % threshold term in erosion law
	K = k_e*(k_t^a)*(k_w^(-a*alpha_val));           % erosional efficiency

	Ks=linspace(floor(ksn_min),ceil(ksn_max),1000);
	E = zeros(size(Ks));        
	Q_starc = zeros(size(Ks));

	% Set Integration Parameters
	q_min = 0.00368*k;           % minimum q needed to have frequency > 1e-8
	q_max = 1000000*exp(-k);      % maximum q, above which frequency is < 1e-8	

	% Numerical integration across Ks range
	w1=waitbar(0,'Numerically Integrating...');
	for ii = 1:length(Ks)
	    % calculate critical discharge for each value of Ks (equation 27)
	    Q_starc(ii) = ((K./psi_crit).*(Ks(ii).^(n))*(k_q^m)).^(-1./y);
	    if Q_starc(ii) < q_min
	        Q_starc(ii) = q_min;
	    elseif Q_starc(ii) > q_max
	        Q_starc(ii) = q_max - 1;
	    end
	    
	    Er =    @erosion_law;                       % incision law (equation 13)
	    PDF =   @inv_gamma;                              % discharge PDF (equation 3)
	    ErPow = @(ks,q,k,kq,kw,ke,tc,f,omega_a,omega_s,alpha_val,beta_val,a) Er(ks,q,kq,kw,ke,tc,f,omega_a,omega_s,alpha_val,beta_val,a).*PDF(q,k);       % integrand of equation 16
	    
	    E(ii) = integral(@(x) ErPow(Ks(ii),x,k,k_q,k_w,k_e,tau_crit,f,omega_a,omega_s,alpha_val,beta_val,a),Q_starc(ii),q_max);
	    waitbar(ii/length(Ks));
	end
	close(w1);
end

function [I] = erosion_law(Ks,Q_star,k_q,k_w,k_e,tau_crit,f,omega_a,omega_s,alpha_val,beta_val,a)
%EROSION_LAW calculates "daily" incision law as either eq. 13 in Lague 2005
%   or eq. 10 in Tucker 2004

sec_per_yr = 31556926;     % seconds per year
Ma =  1000000;      % years per Ma

%derived parameters
k_t = 0.5*1000*(9.81^(2/3))*(f^(1/3));      % set to 1000 a la Tucker 2004
y = a*alpha_val*(1-omega_s);                    % gamma exponent
m = a*alpha_val*(1-omega_a);                    % m in erosion law
n = a*beta_val;                                 % n in erosion law
psi_crit = k_e*tau_crit^a;                   % threshold term in erosion law
K = k_e*(k_q^m)*(k_t^a)*(k_w^(-a*alpha_val));   % erosional efficiency

%main calculation

I = sec_per_yr*Ma*(K*(Ks.^n).*(Q_star.^y)-psi_crit); %incision in m/Ma

end


function [f] = inv_gamma(Q_star,k)
%inv_gamma calculates the inverse gamma frequency distribution
%	used in Lague et al. 2005 eq. 3

%   Last looked at by Roman DiBiase 10/23/2019

f=exp(-k./Q_star).*((k^(k+1))*Q_star.^(-(2+k)))/(gamma(k+1));

end

function [KSNGrid,KSNstdGrid] = KsnAvg(DEM,ksn_ms,radius,er_type)

	% Calculate radius
	radiuspx = ceil(radius/DEM.cellsize);
	SE = strel('disk',radiuspx,0);

	% Record mask of current NaNs
	MASK=isnan(DEM.Z);

	% Make grid with values along channels
	KSNGrid=GRIDobj(DEM);
	KSNGrid.Z(:,:)=NaN;
	for ii=1:numel(ksn_ms)
		ix=coord2ind(DEM,ksn_ms(ii).X,ksn_ms(ii).Y);
		ix(isnan(ix))=[];
		KSNGrid.Z(ix)=ksn_ms(ii).ksn;
	end

	% Local mean based on radius
	ISNAN=isnan(KSNGrid.Z);
    [~,L] = bwdist(~ISNAN,'e');
    ksng = KSNGrid.Z(L);           
    FLT   = fspecial('disk',radiuspx);
    ksng   = imfilter(ksng,FLT,'symmetric','same','conv');

    nhood   = getnhood(SE);
    ksnstd   = stdfilt(ksng,nhood); 

    switch er_type
    case 'std_error'
    	II=~MASK; II=single(II);
    	avg_num=imfilter(II,FLT,'symmetric','same','conv');
    	num_nhood_pix=sum(SE.Neighborhood(:));
    	num_pix=avg_num.*num_nhood_pix;
    	ksnstder=ksnstd./sqrt(num_pix);
    	ksnstder(MASK)=NaN;
    end

    % Set original NaN cells back to NaN
    ksng(MASK)=NaN;
    ksnstd(MASK)=NaN;

    % Output
    KSNGrid.Z=ksng;

    switch er_type
    case 'std'
	    KSNstdGrid=GRIDobj(DEM);
	    KSNstdGrid.Z=ksnstd;
	case 'std_error'
	    KSNstdGrid=GRIDobj(DEM);
	    KSNstdGrid.Z=ksnstder;	
    end	
end
