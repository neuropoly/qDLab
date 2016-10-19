function absc = scd_display_Xaxiscustom(scheme)
% Input :  protocol Nx9 : Gx Gy Gz |G|(mT/um) Delta(ms) delta(ms) TE(ms) q(um-1) identifier
%Custom Xaxis

error('Edit this file to insert your own xaxis')
Gnorm = scheme(:,4);
Gz=scheme(:,3);
absc = Gz./Gnorm;