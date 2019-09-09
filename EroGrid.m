function [ERO,varargout]=EroGrid(DEM,KSN,C,phi,varargin)
	% Usage:
	%	[ERO]=EroGrid(DEM,KSN,C,phi);
	%	[ERO]=EroGrid(DEM,KSN,C,phi,'name',value,...);
	%	[ERO,ERO_P,ERO_M]=EroGrid(DEM,KSN,C,phi,'name',value,...);	
	%
	% Description;
	%	Function to produce an erosion rate map based on an empirical relationship between normalized channel steepness (ksn)
	%	and erosion rate (E) defined in the form of ksn = C*E^phi. This relationship would likely come from a series of 
	%	catchment averaged cosmogenic erosion rates and basin averaged normalized channel steepness. Function can produce an 
	%	erosion rate map for a unique relationship between ksn and E or a series of relationships based on different regions 
	%	defined in a separate grid provided as the optional 'VAL' input. Function is also able to incorporate uncertainty in
	%	the interpolated KSN map and/or uncertainty on the fit parameters (i.e. C and phi) to calculate min and max erosion rate
	%	estimates based on these uncertainties. 
	%
	% Required Inputs:
	%	DEM - DEM GRIDobj
	%	KSN - GRIDobj of continous ksn (e.g. as produced by KsnChiBatch with product set to 'ksngrid') or a map structure of
	%		ksn (e.g. as produced by KsnChiBatch with product set to 'ksn'). If the KSN GRIDobj was produced alternatively, 
	%		please ensure that it is the same dimensions and pixel size as the provided DEM, i.e. the result of 
	%		validatealignment(DEM,KSN) is true.
	%	C - coefficient of power law relationship between ksn and E (ksn = C*E^phi), can be a single value or an array of values. 
	%		If an array of values is provided, it is assumed that these are multiple coefficients that refer to different relationships
	%		based on spatial subsetting that will be defined with inputs to the optional 'edges' and 'VAL'  inputs
	%	phi - exponent of power law relationship between ksn and E (ksn = C*E^phi), can be a single value or an array of values. 
	%		If an array of values is provided, it is assumed that these are multiple exponents that refer to different relationships
	%		based on spatial subsetting that will be defined with inputs to the optional 'edges' and 'VAL'  inputs
	%
	% Optional Inputs:
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
	%	edges [] - array that define the edges to index the GRIDobj provided to VAL. There must be n+1 entries for n entries to 'C' (and 'phi').
	%		It's also assumed that entries for 'C' and 'phi' are in the same order as the bins defined by 'edges'.
	%	resample_method ['nearest'] - method to resample the provided VAL GRIDobj if it is a different dimension or pixel size than the input
	%		DEM. Valid options are 'nearest', 'bilinear', and 'bicubic'. Default is 'nearest'.
	%
	% Note:
	%	The function does not explicitly depend on the 'KSN' arugment actually being KSN. E.g., if you have established a relationship between 
	%	erosion rate and local relief in the form rlf = C*E^phi, then this function will work equally well.
	%
	% Example:
	%	% Unique relationship between ksn and E
	%	[ERO]=EroGrid(DEM,KSN,100,0.5);
	%
	%	% Three different relationships between ksn and E depending on precipitation, where based on empirical data, ksn = 100*E^(0.5) for
	%	% preciptation between 0 and 1 m/yr, ksn = 316*E^(0.5) for precipitation between 1-2 m/yr, and ksn = 1000*E^(0.5) for precipitation
	%	% between 2-4 m/yr. PRECIP is a GRIDobj of precipitation in m/yr with mininum and maximum values greater than 0 and less than 4, 
	%	% respectively
	%	[ERO]=EroGrid(DEM,KSN,[100 316 1000],[0.5 0.5 0.5],'VAL',PRECIP,'edges',[0 1 2 4]);
	%
	%	% Calculate upper (ERO_P) and lower (ERO_M) bounds on erosion rate map based on the standard deviation of the interpolated ksn grid
	%	[ERO,ERO_P,ERO_M]=EroGrid(DEM,KSN,100,0.5,'KSNstd',KSNstd);
	%
	%	% Calculate upper (ERO_P) and lower (ERO_M) bounds on erosion rate map based on uncertainty in fit parameters
	%	[ERO,ERO_P,ERO_M]=EroGrid(DEM,KSN,100,0.5,'C_std',5,'phi_std',0.01);
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 06/16/19 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Parse Inputs
	p = inputParser;         
	p.FunctionName = 'EroGrid';
	addRequired(p,'DEM',@(x) isa(x,'GRIDobj'));
	addRequired(p,'KSN',@(x) isa(x,'GRIDobj') | isstruct(x));
	addRequired(p,'C',@(x) isnumeric(x));
	addRequired(p,'phi',@(x) isnumeric(x));

	addParameter(p,'radius',5000,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'KSNstd',false,@(x) isa(x,'GRIDobj') | islogical(x));
	addParameter(p,'C_std',[],@(x) isnumeric(x));
	addParameter(p,'phi_std',[],@(x) isnumeric(x));
	addParameter(p,'VAL',[],@(x) isa(x,'GRIDobj'));
	addParameter(p,'edges',[],@(x) isnumeric(x));
	addParameter(p,'resample_method','nearest',@(x) ischar(validatestring(x,{'nearest','bilinear','bicubic'})));
	addParameter(p,'error_type','std',@(x) ischar(validatestring(x,{'std','std_error'})));

	parse(p,DEM,KSN,C,phi,varargin{:});
	DEM=p.Results.DEM;
	KSN=p.Results.KSN;
	C=p.Results.C;
	phi=p.Results.phi;

	radius=p.Results.radius;
	ksn_std=p.Results.KSNstd;
	C_std=p.Results.C_std;
	phi_std=p.Results.phi_std;
	VAL=p.Results.VAL;
	edges=p.Results.edges;
	resample_method=p.Results.resample_method;
	error_type=p.Results.error_type;


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
				n_m=1/(phi(ii)+phi_std(ii));
				K_m=(C(ii)+C_stdi(ii))^-n_m;
				n_p=1/(phi(ii)-phi_std(ii));
				K_p=(C(ii)-C_std(ii))^-n_p;

				ERO_P.Z(IDX.Z)=K_p.*KSN.Z(IDX.Z).^n_p;
				ERO_M.Z(IDX.Z)=K_m.*KSN.Z(IDX.Z).^n_m;
			elseif ksn_std_flag & fit_unc_flag 
				n_m=1/(phi(ii)+phi_std(ii));
				K_m=(C(ii)+C_std(ii))^-n_m;
				n_p=1/(phi(ii)-phi_std(ii));
				K_p=(C(ii)-C_std(ii))^-n_p;

				ERO_P.Z(IDX.Z)=K_p.*(KSN.Z(IDX.Z)+KSNstd.Z(IDX.Z)).^n_p;
				ERO_M.Z(IDX.Z)=K_m.*(KSN.Z(IDX.Z)-KSNstd.Z(IDX.Z)).^n_m;
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
