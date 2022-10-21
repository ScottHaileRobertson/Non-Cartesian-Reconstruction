classdef GriddedReconModel < Recon.ReconModel.ReconModel
    methods
        % Constructor
        function [obj] = GriddedReconModel(systemModel, verbose)
            % Call super constructor to build recon obj
            obj = obj@Recon.ReconModel.ReconModel(systemModel, verbose);
        end
        
        function reconVol = reconstruct(obj, data, traj, varargin)
            if(obj.verbose)
                disp('Reconstructing...');
            end
            
            if(nargin<=3 || isempty(varargin{1}))
                % If initial guess isnt given, assume zeros
                gridVol = zeros(obj.system.fullsize);
            else
                % first input is initial guess
                gridVol = varargin{1};
            end

            % Grid data
            if(obj.verbose)
                disp('Gridding Data...');
            end
            reconVol = obj.grid(data,gridVol);
            if(obj.verbose)
                disp('Finished gridding Data...');
            end
            
            % Reshape from vector to matrix
            reconVol = reshape(full(reconVol),ceil(obj.system.fullsize));
            
            if(obj.verbose)
                disp('Calculating IFFTN...');
            end
            % Calculate image space
            reconVol = ifftn(reconVol);
            reconVol = fftshift(reconVol);
            if(obj.verbose)
                disp('Finished calculating IFFTN.');
            end
            
            if(obj.crop)
                reconVol = obj.system.crop(reconVol);
            end
            
            if(obj.deapodize)
                % Calculate deapodization volume and deapodize
                if(obj.verbose)
                    disp('Calculating k-space deapodization function...');
                end
                deapVol = obj.grid(double(~any(traj,2)));
                
                % Reshape from vector to matrix
                deapVol = reshape(full(deapVol),ceil(obj.system.fullsize));
                if(obj.verbose)
                    disp('Finished calculating k-space deapodization function.');
                end
                
                % Calculate image domain representation
                if(obj.verbose)
                    disp('Calculating Image domain deapodization function...');
                end
                deapVol = ifftn(deapVol);
                deapVol = fftshift(deapVol);
                if(obj.verbose)
                    disp('Calculating Image domain deapodization function...');
                end
                
                if(obj.crop)
                    deapVol = obj.system.crop(deapVol);
                end
                
                if(obj.verbose)
                    disp('deapodizing...');
                end
                reconVol = reconVol./deapVol;
                clear deapVol;
                if(obj.verbose)
                    disp('Finished deapodizing.');
                end
            end
            
            if(obj.verbose)
                disp('Finished Reconstructing.');
            end
        end
    end
    methods (Abstract)
        gridVol = grid(obj, data);
    end
end
