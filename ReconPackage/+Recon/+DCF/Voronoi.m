classdef Voronoi < Recon.DCF.DCF
	methods
		% Constructor
		function obj = VoronoiDcf(traj, header, verbosity)
			% Store properties of DCF
			obj.verbose = verbosity;
			obj.dcf_type = 'voronoi';
			obj.unique_string = 'voronoiDcf';
			obj.dcf_style = 'dataspace';
			
			[unique_traj iTraj iUnique] = unique(traj,'rows');
			num_unique = size(unique_traj,1);
			count_data = histc(iUnique,1:max(iUnique));
			
% 			[rows cols] = size(traj);
% 			[sorted_traj sort_idx] = sortrows(traj);
% 			first_unique = [true; any((sorted_traj(1:rows-1,:) ~= sorted_traj(2:rows,:)),2)];
% 			unique_cumsum = cumsum(first_unique);
% 			count_data = histc(unique_cumsum,cumsum(first_unique(first_unique)));
% 			unique_data = sorted_traj(first_unique,:);
% 			unique_idx = sort_idx(first_unique);
% 			num_unique = size(unique_traj,1);	
			
			% Do voronoi calculation
			if(obj.verbose)
				disp('Calculating Voronoi cells...');
			end
			[v c] = voronoin(unique_traj,{'Qbb '});
			if(obj.verbose)
				disp('Finished calculating Voronoi cells.');
			end
			
			%use cellfun
			if(obj.verbose)
				disp('Calculating area of Voronoi cells...');
			end
			obj.dcf = -1*ones(num_unique,1);
			for i = 1:num_unique
				%We Dont know what to do with outermost ones...
				if all(c{i}~=1) 
					% Get all the vertices in the current cell, then reate convex hull
					% and calculate area
					[K, obj.dcf(i)] = convhulln(v(c{i},:));
				else
					% for now, just make them -1
					obj.dcf(i) = -1;
				end
			end
			if(obj.verbose)
				disp('Finished calculating ara of Voronoi cells.');
			end
			
			%Compensate for repeated measurements
			obj.dcf = obj.dcf./count_data;
			
			%Apply dcf to each point
			obj.dcf = obj.dcf(iUnique);
			
			% Set outermost points to be equal to second outermost
			nFrames = header.rdb.rdb_hdr_user20;
			nPts = header.rdb.rdb_hdr_frame_size;

			obj.dcf = reshape(obj.dcf,[nPts nFrames]);
			for iFrame = 1:nFrames
				firstBad = find(obj.dcf(:,iFrame)==-1,1,'first');
				obj.dcf(firstBad:end,iFrame) = obj.dcf(firstBad-1,iFrame);
			end
			obj.dcf = obj.dcf(:);
			
% 			outermost_points = (obj.dcf < 0)';
% 			obj.dcf(outermost_points) = obj.dcf([outermost_points(2:end) false]);
% 			
			
			r = sqrt(traj(:,1).^2+traj(:,2).^2+traj(:,3).^2);
			dc_mean = mean(obj.dcf(r==min(r(:))));
			
			%Normalize DCF area to one
			obj.dcf = obj.dcf./dc_mean;
		end
	end
end
