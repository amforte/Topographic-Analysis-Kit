function cmpEroGrid(wdir,MatFile,KsnFile,C,phi,varargin)
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
	%	wdir - full path of working directory
	%	MatFile - Full path of matfile output from either 'cmpMakeStreams' or the name of a single basin mat file from 'cmpProcessRiverBasins'
	%	KsnFile - ascii grid of continous ksn (e.g. as produced by KsnChiBatch with product set to 'ksngrid') or a shapefile of
	%		ksn (e.g. as produced by KsnChiBatch with product set to 'ksn'). Input to KsnFile must have either a '.txt' or '.shp' file extension
	%	C - coefficient of power law relationship between ksn and E (ksn = C*E^phi), can be a single value or an array of values. 
	%		If an array of values is provided, it is assumed that these are multiple coefficients that refer to different relationships
	%		based on spatial subsetting that will be defined with inputs to the optional 'edges' and 'VAL'  inputs
	%	phi - exponent of power law relationship between ksn and E (ksn = C*E^phi), can be a single value or an array of values. 
	%		If an array of values is provided, it is assumed that these are multiple exponents that refer to different relationships
	%		based on spatial subsetting that will be defined with inputs to the optional 'edges' and 'VAL'  inputs
	%
	% Optional Inputs:
	%	radius [5000] - radius for creating a spatially averaged ksn grid if the entry to KSN is a mapstruct
	%	KSNstd - ascii grid of standard deviation of continous ksn (e.g. as produced by KsnChiBatch with product set to 'ksngrid') to 
	%		incoporate uncertainty in ksn smoothing into the production of erosion rate maps. If you are providing a mapstructure to 
	%		the 'KSN' argument and you wish to use the KSNstd that will be calculated to consider how uncertainty in the ksn value leads to 
	%		uncertainty in the erosion rate map, use the 'use_ksnstd' logical flag
	%	use_ksnstd [false] - logical flag to calculate and use the standard deviation of ksn in estimating erosion rates, in the event that
	%		you are providing a shapefile of Ksn values to the 'KsnFile' input.
	%	C_std - standard deviation / uncertainty on coefficient of power law. If an entry is provided to C_std (and phi_std), then these uncertainties
	%		will be incorporated into the erosion grid outputs. There must be the same number of entires to C_std as there is to 'C'.
	%	phi_std - standard deviation / uncertainty on the exponent of power law. If an entry is provided to phi_std (and C_std), then these uncertainties
	%		will be incorporated into the erosion grid outputs. There must be the same number of entires to phi_std as there is to 'phi'.
	%	VAL [] - ascii grid (.txt) defining the regions by which to index the KSN grid to establish different empirical relationships between ksn and E.
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
    % Examples if running for the command line, minus OS specific way of calling main TAK function:
	%	EroGrid /path/to/wdir Topo.mat batch_ksngrid.txt 100 0.5
	%	EroGrid /path/to/wdir Topo.mat batch_ksngrid.txt 100 0.5 KSNstd batch_ksngridstd.txt
	%	EroGrid /path/to/wdir Topo.mat batch_ksngrid.txt [100 316 1000] [0.5 0.5 0.5] VAL precip.txt edges [0 1 2 4] 		
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Function Written by Adam M. Forte - Updated : 06/16/19 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if isdeployed
		if ~isempty(varargin)
			varargin=varargin{1};
		end
	end

	% Parse Inputs
	p = inputParser;         
	p.FunctionName = 'cmpEroGrid';
	addRequired(p,'wdir',@(x) ischar(x));
	addRequired(p,'MatFile',@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.mat'))));
	addRequired(p,'KsnFile',@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.txt'))) | ~isempty(regexp(x,regexptranslate('wildcard','*.shp'))));
	addRequired(p,'C',@(x) isnumeric(x));
	addRequired(p,'phi',@(x) isnumeric(x));

	addParameter(p,'file_name_prefix','batch',@(x) ischar(x));
	addParameter(p,'radius',5000,@(x) isscalar(x) && isnumeric(x));
	addParameter(p,'KSNstd',[],@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.txt'))));
	addParameter(p,'use_ksnstd',false,@(x) isscalar(x) && islogical(x));
	addParameter(p,'C_std',[],@(x) isnumeric(x));
	addParameter(p,'phi_std',[],@(x) isnumeric(x));
	addParameter(p,'VAL',[],@(x) ~isempty(regexp(x,regexptranslate('wildcard','*.txt'))));
	addParameter(p,'edges',[],@(x) isnumeric(x));
	addParameter(p,'resample_method','nearest',@(x) ischar(validatestring(x,{'nearest','bilinear','bicubic'})));
	addParameter(p,'error_type','std',@(x) ischar(validatestring(x,{'std','std_error'})));

	parse(p,wdir,MatFile,KsnFile,C,phi,varargin{:});
	wdir=p.Results.wdir;
	MatFile=p.Results.MatFile;
	KsnFile=p.Results.KsnFile;
	C=p.Results.C;
	phi=p.Results.phi;

	file_name_prefix=p.Results.file_name_prefix;
	radius=p.Results.radius;
	ksn_std=p.Results.KSNstd;
	use_ksnstd=p.Results.use_ksnstd;
	C_std=p.Results.C_std;
	phi_std=p.Results.phi_std;
	VAL=p.Results.VAL;
	edges=p.Results.edges;
	resample_method=p.Results.resample_method;
	error_type=p.Results.error_type;

	% Determine the type of input for the DEM
	MatFile=fullfile(wdir,MatFile);
	D=whos('-file',MatFile);
	VL=cell(numel(D),1);
	for ii=1:numel(D);
		VL{ii}=D(ii,1).name;
	end

	if any(strcmp(VL,'DEM'))
		load(MatFile,'DEM');
	elseif any(strcmp(VL,'DEMcc'))
		load(MatFile,'DEMoc','FDc','Ac','Sc');
		DEM=DEMoc;
	end

	% Determine the type of input for the KSN and load
	KsnFile=fullfile(wdir,KsnFile);
	if ~isempty(regexp(KsnFile,regexptranslate('wildcard','*.txt')))
		KSN=GRIDobj(KsnFile);
		if ~validatealignment(KSN,DEM)
			KSN=resample(KSN,DEM,resample_method);
		end
	elseif ~isempty(regexp(KsnFile,regexptranslate('wildcard','*.shp')))
		KSN=shaperead(KsnFile);
	end

	% Determine if there is an input for KSNstd
	if ~isempty(ksn_std)
		if ~isempty(regexp(ksn_std,regexptranslate('wildcard','*.txt')))
			ksn_std=fullfile(wdir,ksn_std);
			ksn_std=GRIDobj(ksn_std);
			if ~validatealignment(ksn_std,DEM)
				ksn_std=resample(ksn_std,DEM,resample_method);
			end
		end
	end

	% Determine if there is an input for VAL and load
	if ~isempty(VAL)
		ValFile=fullfile(wdir,VAL);
		VAL=GRIDobj(ValFile);
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
		struct_flag=true;
	else
		struct_flag=false;
	end

	% Check whether the user wants to include uncertainty in the erosion rate calculation
	if isa(ksn_std,'GRIDobj')
		ksn_std_flag=true;
		KSNstd=ksn_std;
	elseif use_ksnstd & struct_flag & isa(KSNstd,'GRIDobj')
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

	% Write out
	ero_name=fullfile(wdir,[file_name_prefix '_erogrid.txt']);
	GRIDobj2ascii(ERO,ero_name);

	if std_flag
		ERO_P.Z(IDX.Z)=NaN;
		ERO_M.Z(IDX.Z)=NaN;

		ERO_P.Z(imag(ERO_P.Z)~=0)=NaN;
		ERO_M.Z(imag(ERO_M.Z)~=0)=NaN;

		% Write out
		ero_p_name=fullfile(wdir,[file_name_prefix '_ero_max_grid.txt']);
		ero_m_name=fullfile(wdir,[file_name_prefix '_ero_min_grid.txt']);		
		GRIDobj2ascii(ERO_P,ero_p_name);
		GRIDobj2ascii(ERO_M,ero_m_name);
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
