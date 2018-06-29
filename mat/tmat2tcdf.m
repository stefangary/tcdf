function tmat2tcdf(lam,phi,dep,tcdf_out)
%-------------------------------------------
%
% tmat2tcdf(lam,phi,dep,tcdf_out)
%
%-------------------------------------------
% This script will write the lam, phi, dep
% variables defined in a trajectory data
% set loaded in MATLAB memory to a netcdf
% file.
%
% The user must verify that the MATLAB
% traj arrays are in the correct orientation
% for tcdf (FLAME trajectory format) output:
%
%    +----ntraj----+
%    |             |
%    |             |
%   npts(time)     |
%    |             |
%    |             |
%    +-------------+
%
% .i.e. columns are single trajectories,
%       rows are an instant in time.
%
% This software is distributed under the
% terms of the GNU GPL v3 and any later version.
% Copyright Stefan Gary 2018
%-------------------------------------------

% Get information about the inputs
[npts,ntraj] = size(lam);

% Preliminary error checking
if( (size(lam) == size(phi)) & (size(lam) == size(dep)) )
    % Input variables are the same sizes.
else
    disp('Error: Input variables are not the same size.')
    return
end

% Filter out NaNs
lam(isnan(lam)) = -999;
phi(isnan(phi)) = -999;
dep(isnan(dep)) = -999;

% Filter out any data points that make no sense
% (i.e. out of bounds of global lon,lat,dep
% coordinates).
dep( (dep<0) & (dep>10000) ) = -999;
lam( (lam<-360) & (lam>360) ) = -999;
phi( (phi<-90) & (phi>90) ) = -999;

% Synchonize lam and phi only.
lam( (lam == -999) | (phi == -999) ) = -999;
phi( (lam == -999) | (phi == -999) ) = -999;

%------------Write to output file------------

disp(strcat('Creating ',tcdf_out,'...'))
ncfid = netcdf.create(tcdf_out,'64BIT_OFFSET');

disp('Create the dimensions...')
timedid = netcdf.defDim(ncfid,'time',npts);
trajdid = netcdf.defDim(ncfid,'traj',ntraj);

vdims(1) = timedid;
vdims(2) = trajdid;

disp('Create the variables...')
lamid = netcdf.defVar(ncfid,'lam','float',vdims);
phiid = netcdf.defVar(ncfid,'phi','float',vdims);
depid = netcdf.defVar(ncfid,'dep','float',vdims);

disp('Exit define mode...')
netcdf.endDef(ncfid);

disp('Write lam to output...')
netcdf.putVar(ncfid,lamid,lam);

disp('Write phi to output...')
netcdf.putVar(ncfid,phiid,phi);

disp('Write dep to output...')
netcdf.putVar(ncfid,depid,dep);

disp('Close output file...')
netcdf.close(ncfid);

disp('Done converting tmat to tcdf.')
return

end
