classdef Proximity
	properties
		kernelObj;
		unique_string;
	end
	
	methods
		function obj = Proximity(kernObj)
			% Save properties to object
			obj.kernelObj = kernObj;
		end
	end
	
	methods (Abstract)
		[kernel_vals idxOutOfBounds] = kernelValues(obj, distances);
	end
end