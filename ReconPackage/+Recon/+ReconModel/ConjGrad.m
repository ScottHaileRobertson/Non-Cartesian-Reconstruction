classdef ConjGrad < Recon.ReconModel.ReconModel
	properties
		iterations;
		saveIterations;
		dcfObj;
	end
	methods
		function obj = ConjGrad(system_model, dcf_obj, iter, saveIter, verbose)
			% Call super constructor to build recon obj
			obj = obj@Recon.ReconModel.ReconModel(system_model, verbose);
			
			% Save properties
			obj.iterations = iter;
			obj.saveIterations = saveIter;
			obj.dcfObj = dcf_obj;
			obj.unique_string = ['cgRecon_iter' num2str(obj.iterations) obj.model.unique_string];
		end
		
		% Reconstructs an image volume using the given data
		function reconVol = reconstruct(obj,data, weights, x,userFun,details)
			% 			C = 0;
			%
			% 			reconVol = qpwls_pcg1(startingGuess, obj.model.A, Gdiag((1/mean(abs(data(:))))*weights), data, C, ...
			% 				'niter', obj.iterations, 'isave',obj.saveIterations, 'userfun', userFun);
			%
			% 			reconVol  = reshape(reconVol, [obj.model.reconMatrixSize ...
			% 				length(obj.saveIterations)]);
			
				
			
			
			% 								We want to solve A'*W*A*x = A'*W*y
			if(obj.verbose)
				disp('Starting itterative SNR compensation');
			end
			
			x = x(:);
			
			switch(obj.dcfObj.dcf_style)
				case 'gridspace'
					error('Not implemented');
% 					Afun = @(x,transp_flag) (obj.model.A'*(weights.*obj.dcfObj.dcf.*(obj.model.A*x)));
% 					b = obj.model.A'*weights.*obj.dcfObj.dcf.*data(:);
% 					for(iIter=1:obj.iterations)
% 						disp(['Iteration ' num2str(iIter)]);
% 						
% 						x = lsqr(@aFunNoSNR, b,0,1,[],[],x);
% 						
% 						userFun(reshape(x,[obj.model.reconMatrixSize]), iIter, details);
% 					end
% 					% 					nonzero_dcf = (obj.dcfObj.dcf~=0);
% 					% 					reconVol = (obj.model.A' * data);
% 					% 					reconVol(nonzero_dcf) = reconVol(nonzero_dcf) .* obj.dcfObj.dcf(nonzero_dcf);
				case 'dataspace'
					% 					reconVol = obj.model.A' * (obj.dcfObj.dcf .* data);
					b = data.*obj.dcfObj.dcf;
					for(iIter=1:obj.iterations)
						disp(['Iteration ' num2str(iIter)]);
						
						% 						[x,flag,relres,iter,resvec] = lsqr(@aFunSNR, b,0,obj.iterations,[],[],x);
						[x,flag,relres,iter,resvec] = lsqr(@aFunNoSNR, b,0,1,[],[],x);
						
						% 						userFun(reshape(x,[obj.model.reconMatrixSize]), iIter, details);
					end
				otherwise
					error('DCF style is not supported');
			end
			
			% Make sparse k-space representation
			if(isa(obj.model,'GriddingSystemModel'))
				x = sparse(ones(size(obj.model.nonEmptyCols)),...
					obj.model.nonEmptyCols,x,1,prod(obj.model.reconMatrixSize),length(x));
				obj.model.nonEmptyCols = [];
			end
			
			reconVol = reshape(full(x),obj.model.reconMatrixSize);
			reconVol = ifftshift(reconVol);
			reconVol = ifftn(reconVol);
			reconVol = fftshift(reconVol);
			
			
			function y = aFunNoSNR(x_,transp_flag)
				if strcmp(transp_flag,'transp')      % x = A'*y
					% 					disp(['Transp, size=' num2str(length(x_(:)))]);
					% 					y = (obj.model.A'*(weights.*obj.dcfObj.dcf.*x_));
					y = (obj.model.A'*(x_.*obj.dcfObj.dcf));
				elseif strcmp(transp_flag,'notransp') % y = A*x
					% Make full matrix
					xFull = sparse(ones(size(obj.model.nonEmptyCols)),...
						obj.model.nonEmptyCols,x_,1,...
						prod(obj.model.reconMatrixSize),length(x_));
					xFull = reshape(full(xFull),obj.model.reconMatrixSize);
					
					% Wavelet filter for compressed Sensing
					x_ = waveletFilter(xFull,0.01,'db1');
					
					% Make sparse and Vectorize
					x_ = x_(obj.model.nonEmptyCols);
					
					% 					disp(['NO Transp, size=' num2str(length(x_(:)))]);
					y = (obj.model.A*x_);
					
				end
			end
			
		end
	end
end
