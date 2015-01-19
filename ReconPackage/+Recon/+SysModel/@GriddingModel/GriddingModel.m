classdef GriddingModel < Recon.SysModel.SystemModel
	properties
		outputImageSize;
		overgridFactor;
		kernel_extent_overgrid;
		kaiser_b;
		proximetyObj;
		verbose;
        nonEmptyCols
	end
	methods
		function obj = GriddingModel(traj, output_image_size, ...
				overgrid_factor, kernel_extent_overgrid, proxObj, verbosity)
			
			% Initialize properties
			obj.outputImageSize = output_image_size;
			obj.overgridFactor = overgrid_factor;
			obj.kernel_extent_overgrid = kernel_extent_overgrid;
			obj.reconMatrixSize = ceil(output_image_size * overgrid_factor);
			obj.proximetyObj = proxObj;
			obj.unique_string = ['griddingModel_' proxObj.unique_string];
			obj.verbose = verbosity;
						
			%Apply proximity info and create sparse system matrix;			
			if(obj.verbose)
				disp('Calculating kernel values...');
			end
			
            [sample_idx,voxel_idx,kernel_vals] = obj.proximetyObj.kernelValues(traj,...
				obj.kernel_extent_overgrid,...
				obj.reconMatrixSize);
            
			if(obj.verbose)
				disp('Finished calculating kernel values.');
            end
			
            	kernel_vals = gather(kernel_vals);
			
			if(obj.verbose)
				disp('Finished calculating kernel values.');
			end
			
			if(obj.verbose)
				disp('Making super sparse system matrix...');
			end
			
			% Make the system matrix more sparse
			[uniqueVox,iVox,iUniqueVox] = unique(voxel_idx);
			clear voxel_idx;
			obj.nonEmptyCols = uniqueVox;
			
			% Make sparse matrix
			A = sparse(sample_idx,iUniqueVox,kernel_vals,size(traj,1),...
				max(iUniqueVox),length(sample_idx));
			
			% Distribute sparse matrix
% 			spmd(6)
% 			codist = codistributor2dbc([2 3],64);
% 			A = codistributed(A,codist);
% 			end
			% Store in object
			obj.A = A;
        end
		
        function reconVol = crop(obj,reconVol)
			%  Shift data by fov/4
			% 			reconVol = circshift(reconVol,round(0.5*(obj.reconMatrixSize-obj.outputImageSize)));
			reconVol = circshift(reconVol,0.5*obj.outputImageSize);
			
			% Crop BL corner
			reconVol = reconVol(1:obj.outputImageSize(1),...
				1:obj.outputImageSize(2), ...
				1:obj.outputImageSize(3));
		end
        
		function obj = imageSpace(obj)
			obj.reconVol = imageSpace@SystemModel(obj,obj.reconVol);
			
% 			%  Shift data by fov/4
% 			% 			reconVol = circshift(reconVol,round(0.5*(obj.reconMatrixSize-obj.outputImageSize)));
% 			obj.reconVol = circshift(obj.reconVol,0.5*obj.outputImageSize);
% 			
% 			% Crop BL corner
% 			obj.reconVol = obj.reconVol(1:obj.outputImageSize(1),...
% 				1:obj.outputImageSize(2), ...
% 				1:obj.outputImageSize(3));
% 			
% 						% Deapodize  at lower res
% 			obj.reconVol = obj.proximetyObj.deapodize(obj.neighborhoodSize,1,obj.outputImageSize,obj.reconVol);
			
			% center the DC freq
			obj.reconVol = fftshift(obj.reconVol);
			
			% Deapodize - consider making at lower res
			obj.reconVol = obj.proximetyObj.deapodize(obj.neighborhoodSize,obj.overgridFactor,obj.reconMatrixSize,obj.reconVol);
			
			% Uncenter the DC freq
			obj.reconVol = fftshift(obj.reconVol);
			
			%  Shift data by fov/4
			% 			reconVol = circshift(reconVol,round(0.5*(obj.reconMatrixSize-obj.outputImageSize)));
			obj.reconVol = circshift(obj.reconVol,0.5*obj.outputImageSize);
			
			% Crop BL corner
			obj.reconVol = obj.reconVol(1:obj.outputImageSize(1),...
				1:obj.outputImageSize(2), ...
				1:obj.outputImageSize(3));
		end
		
		function true_false = isCompatible(obj,traj,header,overgridfactor,nNeighbors);
			true_false = false;
		end
	end
end