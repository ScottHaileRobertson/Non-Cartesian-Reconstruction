%% GAUSSIAN
%
%   A Gaussian class of gridding kernel. 
%   inputs: sigma           - The sharpness of the gaussian function
%           kernelExtent    - The nonzero range of the kernel in units of 
%                             pre-overgridded k-space voxels
%           verbose         - If 1, it will verbosely print information
%
%   Author: Scott Haile Robertson
%   Website: www.ScottHaileRobertson.com
%
classdef Gaussian < Recon.SysModel.Kernel.Kernel
	properties
		sigma;
	end
	
	methods
		% Constructor
		function obj = Gaussian(kernelSigma, kernelExtent, verbose)
			% Call super constructor to build obj
			obj = obj@Recon.SysModel.Kernel.Kernel(kernelExtent, verbose);
			
			% Store properties
			obj.sigma = kernelSigma;

			% Fill in unique string
			obj.unique_string = ['Gauss_e' num2str(obj.extent) ...
				'_s' num2str(obj.sigma)];
		end
		
		function [kernel_vals] = evaluate(obj, kdistance_preovergrid)
			% Calculate Normalized Gaussian Function
			kernel_vals = normpdf(kdistance_preovergrid,0,obj.sigma)/normpdf(0,0,obj.sigma);
		end
	end
end