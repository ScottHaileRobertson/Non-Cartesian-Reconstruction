output_image_size = 128*[1 1 1];
overgrid_factor = 2;
kernel.sharpness = 1/3;
kernel.sharpness_overgrid = kernel.sharpness*overgrid_factor;
kernel.extent = 6*kernel.sharpness;
kernel.extent_overgrid = kernel.extent*overgrid_factor;
verbose = 1;
dcf_iter = 10;

% Load demo data
load demo_radial_mri_data
load demo_radial_mri_traj

% Create kernel
kernelObj = Recon.SysModel.Kernel.Gaussian(...
    kernel.extent_overgrid , overgrid_factor, kernel.sharpness_overgrid, verbose);

% Create Proximity object
proxObj = Recon.SysModel.Proximity.L2Proximity(kernelObj, verbose);
clear kernelObj;

% Create System model
modelObj = Recon.SysModel.GriddingModel(traj, output_image_size, ...
    overgrid_factor, kernel.extent_overgrid, proxObj, verbose);

% Calculate DCF
dcfObj = Recon.DCF.Iterative(modelObj, dcf_iter, verbose);

% Create ReconModel
reconObj = Recon.ReconModel.Lsq(modelObj, dcfObj, verbose);

% Reconstruct
reconVol = reconObj.reconstruct(data,1, 0);

% Deapodize
% bak_dcf = reconObj.dcfObj;
% reconObj.dcfObj = UnityDcf(traj,verbose);
% deapVol = reconObj.reconstruct(dc_data, cropVolume, kspace_recon);
% reconObj.dcfObj = bak_dcf;
% clear dc_data bak_dcf;
%     % 	isSignificant = (abs(deapVol) > 1E-3*max(abs(deapVol(:))));
%     % 	reconVol(isSignificant) = reconVol(isSignificant)./deapVol(isSignificant);
%     reconVol = reconVol./deapVol;

figure()
imagesc(abs(reconVol(:,:,70)));
colormap(gray)
