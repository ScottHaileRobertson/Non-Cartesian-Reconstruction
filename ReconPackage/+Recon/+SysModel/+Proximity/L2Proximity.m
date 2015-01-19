classdef L2Proximity < Recon.SysModel.Proximity.Proximity
	properties
		verbose;
	end
	
	methods
		function obj = L2Proximity(kernelObj, verbosity)
			% Call super constructor to build obj
			obj = obj@Recon.SysModel.Proximity.Proximity(kernelObj);
			
			% Save properties of object
			obj.unique_string = ['L2_' kernelObj.unique_string];
			obj.verbose = verbosity;
		end

		function [sample_idx,voxel_idx,kernel_vals] = kernelValues(obj, traj,...
				kernel_extent,...
				reconMatrixSize)
						% Calculate gridding oversampled_distances
			if(obj.verbose)
				disp('Calculating L2 oversampled_distances...');
			end
			
			[sample_idx,voxel_idx,oversampled_distances] = ...
				Recon.SysModel.Proximity.sparse_gridding_distance_mex(traj',...
				kernel_extent,...
				uint32(reconMatrixSize'), -1);
			if(obj.verbose)
				disp('Finished calculating L2 oversampled_distances.');
			end
			
			if(obj.verbose)
				disp('Applying L2 Bounds...');
			end
			
			% Look for any values that are still out of bounds
			keepValues =  (sample_idx > 0) & (voxel_idx > 0);
			sample_idx = sample_idx(keepValues);
			voxel_idx = voxel_idx(keepValues);
			oversampled_distances = oversampled_distances(keepValues);
			clear keepValues;
			
			if(obj.verbose)
				disp('Finished applying L2 Bounds...');
			end
			
			if(obj.verbose)
				disp('Applying Kernel...');
			end
			% Calculate kernel values
			kernel_vals = obj.kernelObj.kernelValues(oversampled_distances);
			clear oversampled_distances;
			if(obj.verbose)
				disp('Finished applying Kernel...');
			end
		end
	end
end