classdef Analytical < Recon.DCF.DCF
	methods
		% Constructor
		function obj = Analytical(traj, verbosity)
			% Store properties of DCF
			obj.verbose = verbosity;
			obj.dcf_type = 'analytical';
			obj.unique_string = 'analyticalDcf';
			obj.dcf_style = 'dataspace';
			
			
			r = sqrt(sum(traj.^2,2));
			
			obj.dcf = r.^2;
			
% 			nonzero_vals = (obj.dcf ~= 0);
% 			obj.dcf(nonzero_vals) = 1./obj.dcf(nonzero_vals);
		end
	end
end
