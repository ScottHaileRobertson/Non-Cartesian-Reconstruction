classdef Unity < Recon.DCF.DCF
	methods
		% Constructor
		function obj = Unity(traj, verbosity)
			% Store properties of DCF
			obj.verbose = verbosity;
			obj.dcf_type = 'unity';
			obj.unique_string = 'unityDcf';
			obj.dcf_style = 'dataspace';
			
			obj.dcf = ones(size(traj,1),1);
		end
	end
end
