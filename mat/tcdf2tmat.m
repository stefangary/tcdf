function tcdf2tmat(tcdf_in,tmat_out,lt,ls,lr,luvw,ll,lymd,lbot)
%-------------------------------------------
%
% tcdf2tmat(tcdf_in,tmat_out,lt,ls,lr,luvw,ll,lymd,lbot)
%
%-------------------------------------------
% This script will open the cdf file
% containing trajectories with variable
% names lam, phi, dep, etc.  The traj
% will then be saved in .mat format for
% use in other programs.
%
% If lt, ls, lr, ll, lymd, lbot, or luvw=1, then this script
% will attempt to load temp, salt, density,
% launch date, date along traj, bottom depth, or velocity data, respectively,
% from the file.
%
% This software is distributed under the
% terms of the GNU GPL v3 and any subsequent
% version.
% Copyright Stefan Gary 2018
%-------------------------------------------

disp(strcat('Opening ',tcdf_in,'...'))
ncfid = netcdf.open(tcdf_in,'NOWRITE');

% Everything is loaded and everything is written.  For large ensembles of
% trajectory data, this will break on low RAM systems or old versions of
% MATLAB.
disp('Loading variables')
lam = get_lag_var(ncfid,'lam');
phi = get_lag_var(ncfid,'phi');
dep = get_lag_var(ncfid,'dep');
if(lt)
    temp = get_lag_var(ncfid,'temp');
end
if (ls)
    salt = get_lag_var(ncfid,'salt');
end
if(lr)
    rho = get_lag_var(ncfid,'rho');
end
if(luvw)
    u = get_lag_var(ncfid,'u');
    v = get_lag_var(ncfid,'v');
    w = get_lag_var(ncfid,'w');
end
if (ll)
	launch_year = get_lag_var(ncfid,'l_y');
	launch_month = get_lag_var(ncfid,'l_m');
	launch_day = get_lag_var(ncfid,'l_d');
end
if (lymd)
    year = get_lag_var(ncfid,'year');
    mon = get_lag_var(ncfid,'mon');
    day = get_lag_var(ncfid,'day');
end
if ( lbot)
    bdep = get_lag_var(ncfid,'bdep');
end
disp(strcat('Saving trajectories to ',tmat_out,'...'))
save(tmat_out,'lam','phi','dep','-v7.3')
if(lt)
    save(tmat_out,'temp','-v7.3','-APPEND')
end
if(ls)
    save(tmat_out,'salt','-v7.3','-APPEND')
end
if(lr)
    save(tmat_out,'rho','-v7.3','-APPEND')
end
if(luvw)
    save(tmat_out,'u','v','w','-v7.3','-APPEND')
end
if(ll)
	save(tmat_out,'launch_year','launch_month','launch_day','-v7.3','-APPEND')
end
if(lymd)
	save(tmat_out,'year','mon','day','-v7.3','-APPEND')
end
if(lbot)
	save(tmat_out,'bdep','-v7.3','-APPEND')
end

end

%------------------------------------------------
% This function is the MATLAB equivalent of
% get_lag_var in the fortran source.  Given
% a file id and a variable name, extract the
% variable (taking into account the packing)
% to a real valued array real_out.
%------------------------------------------------
function real_out = get_lag_var(file_id,var_name)

    % Get variable ID
    var_id = netcdf.inqVarID(file_id,var_name);
    
    % Get variable information
    [~,~,~,natts] = netcdf.inqVar(file_id,var_id);
    
    % Number of attributes determines the
    % way files are unpacked.
    if(natts == 0)
        % Data are unpacked so read directly
        real_out = double(netcdf.getVar(file_id,var_id));
        
    elseif (natts == 1)
        % Data have only an add_offset,
        % must be years stored in bytes
        add_offset = netcdf.getAtt(file_id,var_id,'add_offset');
        real_out = double(netcdf.getVar(file_id,var_id)) + add_offset;

    elseif (natts == 2)
        % Data have both a scale_factor and an add_offset
        % and are short integers.
        add_offset = netcdf.getAtt(file_id,var_id,'add_offset');
        scale_factor = netcdf.getAtt(file_id,var_id,'scale_factor');
        real_out = double(netcdf.getVar(file_id,var_id))*scale_factor + add_offset;
    end

    return

end
