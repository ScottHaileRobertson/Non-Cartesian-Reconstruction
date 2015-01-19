%   Pipe, J. G., & Menon, P. (1999). Sampling density compensation in MRI:
%   rationale and an iterative numerical solution. Magnetic resonance in
%   medicine : official journal of the Society of Magnetic Resonance in
%   Medicine / Society of Magnetic Resonance in Medicine, 41(1), 179–86.
%   Retrieved from http://www.ncbi.nlm.nih.gov/pubmed/10025627
classdef Iterative < Recon.DCF.DCF
		properties
			dcf_iterations;
		end
	methods
		% Constructor
		function obj = Iterative(model, iterations, verbosity)
			% Store properties of DCF
			obj.verbose = verbosity;
			obj.dcf_iterations = iterations;
			obj.dcf_type = 'iterative';
			obj.unique_string = ['iterativeDcf' num2str(obj.dcf_iterations) 'iter'];
			obj.dcf_style = 'dataspace';
						
% 			dcf_model = double(model.A~=0);
% 			spmd % Take advantage of parallel resources
			obj.dcf = 1./(model.A * ones(size(model.A,2),1)); % Reasonable first guess
			for iter = 1:obj.dcf_iterations
				if(obj.verbose)
					disp(['   DCF Iteration:' num2str(iter)]);
				end
				obj.dcf = obj.dcf ./ (model.A * (model.A'*obj.dcf));
				
			end
% 			end
		end
	end
end