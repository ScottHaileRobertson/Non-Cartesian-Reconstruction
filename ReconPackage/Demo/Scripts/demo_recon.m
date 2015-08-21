%% DEMO_RECON
%
% A demonstration of the capabilities of the Recon toolbox. Note that there
% are many additional kernels that can be uncommented to show how the
% various components can "plug and play" with eachother. An alternative to
% this "manual" method, is to just use "recontool" which uses a GUI to
% setup the objects.
%
% Author: Scott Haile Robertson
% Website: www.ScottHaileRobertson.com
%

%% 1. Define reconstruction parameters
output_image_size = 128*[1 1 1];
overgrid_factor = 3;
kernel.sharpness = 1/3;
kernel.extent = 6*kernel.sharpness;
verbose = 1;

%% 2. Load demo data
load demo_radial_mri_data
load demo_radial_mri_traj

%% 3. Choose kernel
kernelObj = Recon.SysModel.Kernel.Gaussian(kernel.sharpness, kernel.extent, verbose);
% kernelObj = Recon.SysModel.Kernel.KaiserBessel(kernel.sharpness, kernel.extent, verbose);
% kernelObj = Recon.SysModel.Kernel.Sinc(kernel.sharpness, kernel.extent, verbose);

%% 4. Choose Proximity object
proxObj = Recon.SysModel.Proximity.L2Proximity(kernelObj, verbose);
% proxObj = Recon.SysModel.Proximity.L1Proximity(kernelObj, verbose);
clear kernelObj;

%% 5. Create System model
systemObj = Recon.SysModel.MatrixSystemModel(traj, overgrid_factor, ...
    output_image_size, proxObj, verbose);

%% 6. Choose density compensation function (DCF)
dcfObj = Recon.DCF.Analytical3dRadial(traj, verbose);
dcfObj = Recon.DCF.Iterative(systemObj, 10, verbose);
% dcfObj = Recon.DCF.Voronoi(traj, header, verbose);
% dcfObj = Recon.DCF.Unity(traj, verbose);

%% 7. Choose Reconstruction Model
reconObj = Recon.ReconModel.LSQGridded(systemObj, dcfObj, verbose);
% reconObj = Recon.ReconModel.ConjugateGradient(systemObj, 10, verbose);
clear modelObj;
clear dcfObj;

%% 8. Reconstruct Data
reconVol = reconObj.reconstruct(data, traj);

%% 10. Display the result
figure()
imagesc(abs(reconVol(:,:,round(0.5*size(reconVol,3)))));
colormap(gray)