%KAISERBESSEL   Creates a Kaiser Bessel kernel.
%   KAISERBESSEL(filter_dims, fermi_scale, fermi_width)
%   creates a Kaiser Bessel function from 0-kernel_extent_oversamp. Jackson et al
%   showed that the Kaiser Bessel function is a good approximation to the
%   prolate spheroidal wave function. (Selection of a Convolution Function
%   for Fourier Inversion using Gridding. Jackson et al. 1991.
%
%   Choice of Beta from "Rapid Gridding Reconstruction with a minimal
%   Oversampling ratio" Beaty et all IEEE 2005.
%
%   w = kernel width (in oversampled units)
%   a = overgridfactor
%
%   Copyright: 2012 Scott Haile Robertson.
%   Website: www.ScottHaileRobertson.com
%   $Revision: 1.0 $  $Date: 2012/07/19 $
classdef KaiserBessel < Recon.SysModel.Kernel.Kernel
	properties
		recon_size;
		kernel_extent_oversamp;
		overgrid_factor;
		beta;
	end
	
	methods
		% Constructor
		function obj = KaiserBessel(kernelExtentOversamp, overgridFactor, reconSize, betaOverride, verbose)
			% Call super constructor to build obj
			obj = obj@Recon.SysModel.Kernel.Kernel(verbose);
			
			% Store properties
			obj.kernel_extent_oversamp = kernelExtentOversamp;
			obj.overgrid_factor = overgridFactor;
			obj.recon_size = reconSize;
			obj.unique_string = 'KB';
			
			%Calculate Beta value (Rapid Gridding Reconstruction With a Minimal Oversampling Ratio. Beatty et al. 2005.)
			if(isempty(betaOverride))
				obj.beta = pi*sqrt( (0.5*obj.kernel_extent_oversamp)^2-0.8 );
			else
				obj.beta = betaOverride;
			end
			
			% Fill in unique string
			obj.unique_string = ['kaiser_width' num2str(obj.kernel_extent_oversamp) ...
				'_overgrid' num2str(obj.overgrid_factor) '_beta' num2str(obj.beta)];
		end
		
		function [kernel_vals] = kernelValues(obj, distances)	
			% Calculate Kaiser Bessel Function
			kernel_vals = besseli(0,obj.beta*...
				sqrt(1 - (2*distances/obj.kernel_extent_oversamp).^2))./obj.kernel_extent_oversamp;
			
			% Its tempting to just use some fixed size (like the matrix
			% size), however the kaiser bessel function is created to be
			% zero at the kernel size, thus its best to just use the kernel
			% size. This makes optimization a pain, but theres not much we
			% can do...
			kernel_vals = besseli(0,obj.beta*...
				sqrt(1 - (2*distances/obj.recon_size(1)).^2))./obj.recon_size(1);

			%Normalize - NOTE I SHOULD BE NORMALIZING AREA, but this is
			%just a global scale, so who cares...
			kernel_vals = kernel_vals/max(kernel_vals(:));
		end
	end
end