classdef Kernel
	properties
		verbose;
		unique_string;
	end
	
	methods
		function obj = Kernel(verbosity)
			% Save properties of object
			obj.verbose = verbosity;
		end
	end
	
	methods (Abstract)
		 kernel_vals = kernelValues(obj, distances);	
	end
end