classdef Optimized < Recon.SysModel.Kernel.Kernel
	properties
		kernel_extent_oversamp;
		overgrid_factor;
		norm_val;
		nIter;
		interp_values;
		interp_dist;
		ikern;
	end
	
	
	%%
	%
	%  IDEA - make sure k-space and image space are well oversampled in
	%  this example so that the optimality is good in both domains
	%
	%
	%
	%
	%
	%
	%%
	
	methods
		% Constructor
		function obj = Optimized(kernelExtentOversamp, imageSize, ...
				overgridFactor, kernOversampFactor, iter, actFOV, verbose)
			% Call super constructor to build obj
			obj = obj@Recon.SysModel.Kernel.Kernel(verbose);
			
			% Store properties
			obj.kernel_extent_oversamp = kernelExtentOversamp;
			obj.overgrid_factor = overgridFactor;
			obj.nIter = iter;
			
			if(obj.verbose)
				disp('Optimizing kernel...');
			end
			
			lut_size = imageSize*kernOversampFactor^2;
			
			% Calculate k-space in units of overgrid units
			k_lsp_fine = imageSize*kernOversampFactor*...
				linspace(-0.5, 0.5, lut_size+1); % in overgrid units
			k_lsp_fine = k_lsp_fine(1:(end-1));
			
			% Calculate image space in units inverse overgrid units
			delta_k = k_lsp_fine(2)-k_lsp_fine(1);
			i_lsp_fine = (imageSize/delta_k)*...
				linspace(-0.5,0.5,lut_size+1);
			i_lsp_fine = i_lsp_fine(1:(end-1));
			
			% Create binary bounding functions
			k_bound = (abs(k_lsp_fine) <= (0.5*obj.kernel_extent_oversamp));
			i_bound = (abs(i_lsp_fine) <= 0.5*(imageSize*actFOV));
			
			% Create plot bounds
			k_plot = (abs(k_lsp_fine) <= obj.kernel_extent_oversamp);
			i_plot = (abs(i_lsp_fine) <= imageSize*obj.overgrid_factor);
			
			% Start with Spatially bounded signal
			i_kern = i_bound;
			k_kern = zeros(size(i_lsp_fine)); % dummy initialization
			
			for iBound = 1:obj.nIter
				if(obj.verbose)
					disp(['Kernel optimization ' num2str(iBound) '/' num2str(obj.nIter)]);
				end
				% Calculate kernel in freq domain
				k_kern = fftshift(fftn(i_kern));
				
				if(obj.verbose)
					figure(1);
					subplot(1,2,1);
					plot(k_lsp_fine,real(k_kern),'-r');
					hold on
					plot(k_lsp_fine,imag(k_kern),'-g');
					plot(k_lsp_fine,abs(k_kern),'-b');
					legend('Real','Imaginary','Magnitude');
					hold off;
					% 					axis([-kernelExtentOversamp*overgridFactor kernelExtentOversamp*overgridFactor -max(abs(k_kern(:))) max(abs(k_kern(:)))]);
					title('Kspace kernel');
					subplot(1,2,2);
					plot(i_lsp_fine,real(i_kern));
					axis([min(i_lsp_fine(i_plot(:)))...
						max(i_lsp_fine(i_plot(:)))...
						-max(real(i_kern(i_plot(:))))...
						max(real(i_kern(i_plot(:))))]);
										set(gca,'XTick',0.5*[-imageSize*overgridFactor -imageSize imageSize imageSize*overgridFactor],...
											'XTickLabel',{'-eFOV/2', '-FOV/2','FOV/2', 'eFOV/2'});
					title('Image kernel');
				end
				
				% Enforce that its all real
				% 				k_kern = abs(k_kern);
				
				% Bound kernel in freq domain
				k_kern = k_bound.*k_kern;
				
				% Normalize kernel
				k_kern = k_kern/max(k_kern(:));
				
				% Calculate spatial kernel
				i_kern = ifftn(ifftshift(k_kern));
				
				if(obj.verbose)
					figure(1);
					subplot(1,2,1);
					plot(k_lsp_fine,real(k_kern),'-r');
					hold on
					plot(k_lsp_fine,imag(k_kern),'-g');
					plot(k_lsp_fine,abs(k_kern),'-b');
					hold off;
					legend('Real','Imaginary','Magnitude');
					% 					axis([-kernelExtentOversamp*overgridFactor kernelExtentOversamp*overgridFactor -max(abs(k_kern(:))) max(abs(k_kern(:)))]);
					title('Kspace kernel');
					subplot(1,2,2);
					plot(i_lsp_fine,real(i_kern));
					axis([min(i_lsp_fine(i_plot(:)))...
						max(i_lsp_fine(i_plot(:)))...
						-max(real(i_kern(i_plot(:))))...
						max(real(i_kern(i_plot(:))))]);
										set(gca,'XTick',0.5*[-imageSize*overgridFactor -imageSize imageSize imageSize*overgridFactor],...
											'XTickLabel',{'-eFOV/2', '-FOV/2','FOV/2', 'eFOV/2'});
					title('Image kernel');
				end
				
				% Bound kernel in spatial domain
				i_kern = i_bound.*i_kern;
			end
			
			
			
			% Calculate kernel in freq domain
			k_kern = fftshift(fftn(i_kern));
			
			% Enforce that its all real
			% 				k_kern = abs(k_kern);
			
			% Bound kernel in freq domain
			k_kern_last = k_kern/max(k_kern);
			k_kern = k_bound.*k_kern;
			
			% Calculate spatial kernel
			i_kern = ifftn(ifftshift(k_kern));
			
			% Normalize kernel
			k_kern = k_kern/max(k_kern(:));
			k_kern_last = k_kern;
			
			% Calculate normalization value
			obj.interp_values = k_kern/max(k_kern(k_bound));
			obj.interp_dist = k_lsp_fine;
			obj.ikern = i_kern;
			
			if(obj.verbose)
				
				figure(1);
				subplot(1,2,1);
				plot(obj.interp_dist,real(obj.interp_values),'-r');
				hold on;
				plot(k_lsp_fine,imag(obj.ikern),'-g');
				plot(obj.interp_dist,abs(obj.interp_values),'-b');
				hold off;
				legend('Real','Imaginary','Magnitude');
% 				axis([-kernelExtentOversamp*overgridFactor kernelExtentOversamp*overgridFactor -max(abs(k_kern(:))) max(abs(obj.interp_values))]);
				title('Kspace kernel');
				subplot(1,2,2);
						plot(i_lsp_fine,real(i_kern));
					axis([min(i_lsp_fine(i_plot(:)))...
						max(i_lsp_fine(i_plot(:)))...
						-max(real(i_kern(i_plot(:))))...
						max(real(i_kern(i_plot(:))))]);
						set(gca,'XTick',0.5*[-imageSize*overgridFactor -imageSize imageSize imageSize*overgridFactor],...
									'XTickLabel',{'-eFOV/2', '-FOV/2','FOV/2', 'eFOV/2'});
				title('Image kernel');
			end
			
			x_test = linspace(-4,4,length(obj.interp_dist));
			test_vals1 = obj.kernelValues(x_test);
			test_vals2 = obj.kernelValues2(x_test);
			
			figure();
			subplot(1,3,1)
			plot(x_test,abs(test_vals1),'x-r','MarkerSize',5);
			hold on;
			plot(x_test,abs(test_vals2),'+-g','MarkerSize',5);
			plot(obj.interp_dist,abs(k_kern_last),'.b','MarkerSize',5);		
			title('Magnitude');
			
			subplot(1,3,2)
			plot(x_test,real(test_vals1),'x-r','MarkerSize',5);
			hold on;
			plot(x_test,real(test_vals2),'+-g','MarkerSize',5);
			plot(obj.interp_dist,real(k_kern_last),'.b','MarkerSize',5);	
			title('Real');
			
			subplot(1,3,3)
			plot(x_test,imag(test_vals1),'x-r','MarkerSize',5);
			hold on;
			plot(x_test,imag(test_vals2),'+-g','MarkerSize',5);
			plot(obj.interp_dist,imag(k_kern_last),'.b','MarkerSize',5);	
			title('Imaginary');
			legend('Fourier series','interpolated','Measured')
			
			if(obj.verbose)
				disp('Finished optimizing kernel.');
			end
% 			
% 			error('Stopping');
			
			% Fill in unique string
			obj.unique_string = ['optimal_width' num2str(obj.kernel_extent_oversamp) ...
				'_obj.overgrid_factor' num2str(obj.overgrid_factor)...
				'_nIter' num2str(obj.nIter)];
		end
		
		function [kv2] = kernelValues2(obj, distances)

			delta_k = obj.interp_dist(2)-obj.interp_dist(1);
			fov_i = 1/delta_k;
			lsp_i = linspace(-0.5*fov_i,0.5*fov_i,length(obj.interp_dist)+1)+0.5*fov_i;
			lsp_i = lsp_i(1:(end-1));
			
			npts = length(distances);
			kv2 = zeros(size(distances));
			parfor ik = 1:npts
				kv2(ik) = sum(obj.ikern.*exp(-2*pi*1i*(lsp_i)*distances(ik)));
			end
		end
		
		function [kernel_vals] = kernelValues(obj, distances)
			if(obj.verbose)
				disp('Calculating kernel values...');
			end
	
			% Calculate magnitude
			mag_vals = interp1(obj.interp_dist,abs(obj.interp_values),distances,'linear',0);
			
			% Calculate phase values
			phase_vals = interp1(obj.interp_dist,unwrap(angle(obj.interp_values)),distances);
			if(any(isnan(phase_vals)))
				if(any(mag_vals(isnan(phase_vals))~=0))
					error('Impossible state');
				else
					% Doesnt matter since mag is zero
					phase_vals(isnan(phase_vals)) = 0;
				end
			end
			
			% Compine magnitude and phase
			% 			kernel_vals =  mag_vals;
			kernel_vals =  mag_vals.*exp(phase_vals*1i);
			
			if(obj.verbose)
				disp('Finished calculating kernel values.');
			end
		end
	end
end