classdef L1Proximity < Recon.SysModel.Proximity.Proximity
	properties
		verbose;
	end
	
	methods
		function obj = L1Proximity(kernelObj, verbosity)
			% Call super constructor to build obj
			obj = obj@Recon.SysModel.Proximity.Proximity(kernelObj);
			
			% Save properties of object
			obj.verbose = verbosity;
			obj.unique_string = ['L1_' kernelObj.unique_string];
		end
		
		function [sample_idx,voxel_idx,kernel_vals] = kernelValues(obj, traj,...
				kernel_extent,...
				reconMatrixSize)
			
			if(obj.verbose)
				disp('Calculating L1 distances...');
			end
			
			[nSamps nDims] = size(traj);
			keepValues = true(1,1);
			kernel_vals = 1;
			for iDim = 1:nDims
				if(obj.verbose)
					disp(['Dim ' num2str(iDim) '...']);
				end
							
				[sample_idx,voxel_idx,distDim] = ...
					Recon.SysModel.Proximity.sparse_gridding_distance_mex(traj',...
					kernel_extent,...
					uint32(reconMatrixSize'), iDim-1);
				
				% Look for any values that are still out of bounds
				keepValues = keepValues & (sample_idx > 0) & (voxel_idx > 0);
				
				% Accumulate kernel values
				kernel_vals = kernel_vals.* ...
					obj.kernelObj.kernelValues(abs(distDim));
			end
			
			if(obj.verbose)
				disp('Finished calculating L1 distances.');
			end
			
			if(obj.verbose)
				disp('Applying L1 Bounds...');
			end
			

			sample_idx = sample_idx(keepValues);
			voxel_idx = voxel_idx(keepValues);
			kernel_vals = kernel_vals(keepValues);
			clear keepValues;
			
			if(obj.verbose)
				disp('Finished applying L1 Bounds...');
			end
		end
	end
end