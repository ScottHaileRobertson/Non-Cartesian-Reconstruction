classdef Gaussian < Recon.SysModel.Kernel.Kernel
	properties
		kernel_extent_oversamp;
		overgrid_factor;
		sigma;
	end
	
	methods
		% Constructor
		function obj = Gaussian(kernelExtentOversamp, overgridFactor, sigma, verbose)
			% Call super constructor to build obj
			obj = obj@Recon.SysModel.Kernel.Kernel(verbose);
			
			% Store properties
			obj.kernel_extent_oversamp = kernelExtentOversamp;
			obj.overgrid_factor = overgridFactor;
			obj.unique_string = 'GaussKern';
			
			obj.sigma = sigma;
						
			% Fill in unique string
			obj.unique_string = ['gaussian_width' num2str(obj.kernel_extent_oversamp) ...
				'_overgrid' num2str(obj.overgrid_factor) '_sigma' num2str(obj.sigma)];
		end
		
		function [kernel_vals] = kernelValues(obj, distances)
			% Calculate Normalized Gaussian Function
			kernel_vals = normpdf(distances,0,obj.sigma);
			kernel_vals = kernel_vals/max(kernel_vals(:));
		end
	end
end