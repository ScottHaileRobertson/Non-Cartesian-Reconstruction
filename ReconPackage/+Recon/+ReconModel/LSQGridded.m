classdef LSQGridded < Recon.ReconModel.GriddedReconModel
	properties
		D;
	end
	methods
		% Constructor
		function [obj] = LSQGridded(systemModel, D, verbose)
			% Call super constructor to build recon obj
			obj = obj@Recon.ReconModel.GriddedReconModel(systemModel, verbose);
		
            % Store properties
            obj.D = D;
			obj.unique_string = ['grid_' obj.system.unique_string '_' obj.D.unique_string];
        end
        
        function gridVol = grid(obj, data)
			switch(obj.D.space)
				case 'gridspace'
					gridVol = (obj.system' * data) .* obj.D;
				case 'dataspace'
					gridVol = obj.system' * (obj.D .* data);
				otherwise
					error('DCF style not recognized');
            end
		end
	end
end
