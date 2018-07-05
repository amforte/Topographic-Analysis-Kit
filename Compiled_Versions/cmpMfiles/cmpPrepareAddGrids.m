function cmpPrepareAddGrids(wdir,out_file_name,varargin)
	% Function to prepare additional grids for use in 'cmpProcessRiverBasins'
	%
	% Inputs:
	%	out_file_name - name for the mat file to be produced, do not include the '.mat'
	%
	%	Additional inputs must be given in groups of twos, and be in the order: 
	%		1) name of the ascii or geotiff of the extra raster data you want to include (must be in the same projection
	%			as the original grid you provided to 'cmpMakeStreams' and will use in 'cmpProcessRiverBasins')
	%		2) a reference name for the produced grid
	%
	%	
	if isdeployed
		if ~isempty(varargin)
			varargin=varargin{1};
		end
	end
	
	if mod(numel(varargin),2)~=0
		error('You must provide two entries for each additional grid you wish to prepare');
	end

	num_grids=numel(varargin)/2;

	% Generate indices
	ix1=[1:2:numel(varargin)];
	ix2=[2:2:numel(varargin)];

	% Generate empty cell
	AG=cell(num_grids,2);

	for ii=1:num_grids
		try
			AG{ii,1}=GRIDobj(fullfile(wdir,varargin{ix1}));
		catch
			error('First entry not recognized as either a geotiff or ascii');
		end

		if ischar(varargin{ix2})
			AG{ii,2}=varargin{ix2};
		else
			error('Second entry not recognized as a valid character entry');
		end
	end

	fn=fullfile(wdir,[out_file_name '.mat']);
	save(fn,'AG');
end