classdef Sinc < Recon.SysModel.Kernel.Kernel
	properties
		kernel_extent_oversamp;
		overgrid_factor;
	end
	
	methods
		% Constructor
		function obj = Sinc(kernelExtentOversamp, overgridFactor, verbose)
			% Call super constructor to build obj
			obj = obj@Recon.SysModel.Kernel.Kernel(verbose);
			
			% Store properties
			obj.kernel_extent_oversamp = kernelExtentOversamp;
			obj.overgrid_factor = overgridFactor;
						
			% Fill in unique string
			obj.unique_string = ['sinc_width' num2str(obj.kernel_extent_oversamp) ...
				'_overgrid' num2str(obj.overgrid_factor)];
		end
		
		function [kernel_vals] = kernelValues(obj, distances)
		
			% Calculate sinc Function
			kernel_vals = sin(2*pi*distances)./(2*pi*distances);
			kernel_vals(distances==0)=1;
			
			%Normalize
			kernel_vals = kernel_vals/max(kernel_vals(:));
% 			kernel_vals = kernel_vals/sum(kernel_vals(:)); % seems better than max...  
		end
	end
end