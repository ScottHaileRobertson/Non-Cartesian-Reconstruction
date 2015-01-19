classdef ReconModel
	properties
		model;			% System Matrix
		verbose;
		unique_string;
	end
	
	methods
		function obj = ReconModel(system_model, is_verbose)
			% Save the System Matrix
			obj.model = system_model;
			obj.verbose = is_verbose;
		end
	end
	
	methods (Abstract)
		% Reconstructs an image volume using the given data
		reconVol = reconstruct(obj,data, cropVolume, kspace_recon);
	end
end
