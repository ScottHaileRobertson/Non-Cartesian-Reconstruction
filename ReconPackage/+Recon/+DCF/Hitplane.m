classdef Hitplane < Recon.DCF.DCF
	methods
		% Constructor
		function obj = Hitplane(model, verbosity)
			% Store properties of DCF
			obj.verbose = verbosity;
			obj.dcf_type = 'hitplane';
			obj.unique_string = 'hitplaneDcf';
			obj.dcf_style = 'gridspace';
			
			% Note - hitplane DCF is on the overall kspace image,
			% not the data itself (dimension of kspace not data)
			if(isa(model.A,'Fatrix'))
				%Handle silly Fessler object
				obj.dcf = sum(model.A.arg.G~=0,1)';
				
				% 				Doesnt seem to work - perhaps complex data
				% 				cancels itself out, resulting in less
				% 				density compensation?
				% 				obj.dcf = max(sum(model.A.arg.G,1)',1); %
				
				% This doesn't work either, I guess we just have to
				% count...
				% obj.dcf = max(sum(abs(model.A.arg.G),1)',1);
			else
				obj.dcf = sum(model.A~=0,1)';
				
				% 				Doesnt seem to work - perhaps complex data
				% 				cancels itself out, resulting in less
				% 				density compensation?
				% 				obj.dcf = max(sum(model.A,1)',1);
				
				% This doesn't work either, I guess we just have to
				% count...
				% obj.dcf = max(sum(abs(model.A),1)',1);
			end

			obj.dcf = full(obj.dcf);
			nonzero_vals = (obj.dcf~=0);
			obj.dcf(nonzero_vals) = 1./(model.kernel_extent_oversamp(1)*obj.dcf(nonzero_vals));
			
% % 			Im not sure why this works better on the grid than the data,
% % 			but I tested it and its terrible with dataspace
% 						obj.dcf_style = 'dataspace';
% 			
% 						% Note - hitplane DCF is on the overall kspace image,
% 						% not the data itself (dimension of kspace not data)
% 						if(isa(model.A,'Fatrix'))
% 							%Handle silly Fessler object
% 							obj.dcf = sum(model.A.arg.G~=0,2);
% 							% 				obj.dcf = max(sum(model.A.arg.G,1)',1);
% 						else
% 							obj.dcf = sum(model.A~=0,2);
% 							% 				obj.dcf = max(sum(model.A,1)',1);
% 						end
			
% 			obj.dcf = 1./(model.neighborhoodSize(1)*obj.dcf);
		end
	end
end