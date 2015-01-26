%% ITERATIVE
%
%   An iterative DCF class based off:
%   Pipe, J. G., & Menon, P. (1999). Sampling density compensation in MRI:
%   rationale and an iterative numerical solution. Magnetic resonance in
%   medicine : official journal of the Society of Magnetic Resonance in
%   Medicine / Society of Magnetic Resonance in Medicine, 41(1), 179–86.
%   Retrieved from http://www.ncbi.nlm.nih.gov/pubmed/10025627
%
%   Author: Scott Haile Robertson
%   Website: www.ScottHaileRobertson.com
%
classdef Iterative < Recon.DCF.DCF
    properties
        iterations;
    end
    methods
        % Constructor
        function obj = Iterative(modelObj, dcfIiterations, verbosity)
            % Call super constructor to build obj
            obj = obj@Recon.DCF.DCF(verbosity);
            
            % Store properties of DCF
            obj.iterations = dcfIiterations;
            obj.unique_string = ['iter' num2str(obj.iterations)];
            obj.space = 'dataspace';
            
            % 			obj.model.systemModel = [];
            
            % Be clever and only calculate non-sparse k-space voxels by
            % using a sparser system model. This way you dont grid or
            % ungrid zeros unnecessarily
            if(isa(modelObj,'Recon.SysModel.MatrixSystemModel'))
                modelObj = modelObj.makeSuperSparse();
            end
            
            obj.dcf = 1./(modelObj.A * ones(size(modelObj.A,2),1)); % Reasonable first guess
            for iter = 1:obj.iterations
                if(obj.verbose)
                    disp(['   DCF Iteration:' num2str(iter)]);
                end
                obj.dcf = obj.dcf ./ (modelObj.A * (modelObj.A'*obj.dcf));
            end
            
            % Undo sparser system model so we dont break anything
            if(isa(modelObj,'Recon.SysModel.MatrixSystemModel'))
                modelObj = modelObj.revertSparseness();
            end
        end
    end
end