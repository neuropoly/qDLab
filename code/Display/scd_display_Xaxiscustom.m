function absc = scd_display_Xaxiscustom(scheme, data)
% Input :  protocol Nx9 : Gx Gy Gz |G|(mT/um) Delta(ms) delta(ms) TE(ms) q(um-1) identifier
% data : MRI Signal in the selected voxel
% absc : Custom Xaxis along which will be plotted the functions

dbstop if error
error('Edit this file to insert your own xaxis')

bvec = scheme(:,[1 2 3]);

% Xaxis
absc=bvec(:,3);