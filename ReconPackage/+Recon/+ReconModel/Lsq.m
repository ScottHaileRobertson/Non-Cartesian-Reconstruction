classdef Lsq < Recon.ReconModel.ReconModel
	properties
		dcfObj;
	end
	methods
		% Constructor
		function obj = Lsq(system_model, dcf_obj, verbose)
			% Call super constructor to build recon obj
			obj = obj@Recon.ReconModel.ReconModel(system_model,verbose);
			
			% Store properties
			obj.dcfObj = dcf_obj;
			
			obj.unique_string = ['lsqRecon_' obj.model.unique_string '_' obj.dcfObj.unique_string];
		end
		
		% Reconstructs an image volume using the given data
		function reconVol = reconstruct(obj,data,cropVolume,kspace_recon)
			if(obj.verbose)
				disp('Reconstructing...');
			end
			switch(obj.dcfObj.dcf_style)
				case 'gridspace'
					nonzero_dcf = (obj.dcfObj.dcf~=0);
					reconVol = (obj.model.A' * data);
					reconVol(nonzero_dcf) = reconVol(nonzero_dcf) .* obj.dcfObj.dcf(nonzero_dcf);
				case 'dataspace'
					reconVol = obj.model.A' * (obj.dcfObj.dcf .* data);
				otherwise
					error('DCF style not recognized');
			end
			
			% 			obj.model.A = [];
			
			% Make sparse k-space representation
% 			if(isa(obj.model,'GriddingSystemModel'))
				reconVol = sparse(ones(size(obj.model.nonEmptyCols)),...
					obj.model.nonEmptyCols,reconVol,1,prod(obj.model.reconMatrixSize),length(reconVol));
				
				obj.model.nonEmptyCols = [];
% 			end
			
			% 			obj.dcfObj = [];
			
			% Reshape from vector to matrix
			reconVol = reshape(full(reconVol),obj.model.reconMatrixSize);
			
			% Crop volume
			if(cropVolume)
				if(obj.verbose)
					disp('Calculating IFFTN...');
				end
				
				% Calculate image space
				reconVol = ifftn(reconVol);
				
				if(obj.verbose)
					disp('Finished calculating IFFTN.');
					disp('Cropping...');
				end
				
				% Crop volume
				reconVol = obj.model.crop(reconVol);
				
				if(obj.verbose)
					disp('Finished cropping.');
				end
				
				% if kspace is wanted, take fft of cropped volume
				if(kspace_recon)
					if(obj.verbose)
						disp('Calculating FFTN...');
					end
					
					reconVol = ifftn(reconVol);
					
					if(obj.verbose)
						disp('Finished calculating FFTN.');
					end
				end
			else
				% Put into spatial domain
				if(~kspace_recon)
					if(obj.verbose)
						disp('Calculating IFFTN...');
					end
					
					reconVol = ifftn(reconVol);
					
					% Shift by 
					reconVol = fftshift(reconVol);
% 					circshift(reconVol,-0.5*obj.model.outputImageSize);
					
					if(obj.verbose)
						disp('Finished calculating IFFTN.');
					end
				end
			end
			
			if(obj.verbose)
				disp('Finished Reconstructing.');
			end
		end
	end
end
