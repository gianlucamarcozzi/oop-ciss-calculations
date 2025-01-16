function [Sys, Exp] = importparam(pathToParam)
% Sys = addnucspin(Sys, nucStr, A, nNuc)
%
% Input
%   Sys
%   nucStr
%   A0      (1, nElectrons x 3)
%   A0Frame (1, nElectrons x 3)
%   nNuc
%
% Output
%   Sys
%     .Nucs         of the form [nucStr, nucStr, nucStr, ...]
%     .Sys.A        (nNuc, nElectrons x 3)
%     .Sys.AFrame   (nNuc, nElectrons x 3)

if ~isfile(pathToParam)
    error("File not found.")
end

M = readlines(pathToParam);

sysPar = ["S", "g", "gFrame", "eeFrame", "A", "AFrame", "Nucs", "nNuc", ...
    "J", "dip", "lwpp", "initState"];
expPar = ["x", "mwFreq", "CenterSweep", "Range", "nPoints", "Harmonic", ...
    "nTheta", "nPhi"];

Sys = struct();
Exp = struct();
for ii = 1:numel(M)
    m = M(ii);  % Current line
    % Get parName from the file row
    temp = strsplit(M{ii}, '=');
    parName = strrep(temp(1), ' ', '');
    % Compare with expected parameters
    [sysFlag, iSys] = ismember(parName, sysPar);
    [expFlag, iExp] = ismember(parName, expPar);
    if sysFlag
        par = sysPar(iSys);
    elseif expFlag
        par = expPar(iExp);
    elseif isscalar(temp)  % no '=' in this line
        continue
    else
        disp(parName)
        error("Parameter %s not found in the list " + ...
            "of expected parameters.\n", parName{:})
    end
    zz = ii;
    while contains(M{zz}, '...')
        m(end + 1) = M(zz + 1);
        zz = zz + 1;
    end
    % fprintf("par = %s:\n\t", par)
    % disp(m)
    if sysFlag
        Sys = addpar(Sys, m, par);
    else
        Exp = addpar(Exp, m, par);
    end

end
Sys = getrhofrominitstate(Sys);
Exp = getxfromrange(Exp);
end

function Sys = addpar(Sys, m, par)
% par should be a string vector
    temp = strsplit(m{1}, '=');
    m(1) = temp(end);
    for ii = 1:numel(m)
        temp2 = strsplit(m{ii}, '%');
        tt{ii} = temp2(1);
    end

    tfin = "";
    for ii = 1:numel(m)
        tt(ii) = strrep(tt{ii}, '...', '');
        tfin = append(tfin, tt{ii});
    end
    Sys.(par) = eval(tfin);
end

function Sys = getrhofrominitstate(Sys)
    if ~isfield(Sys, 'initState')
        warning("The initial density matrix is set to 'singlet'" + ...
            " since no initState was declared in the parameters.")
        initState = 'singlet';
    else
        initState = lower(Sys.initState);
    end

    Sys.rho = zeros(4, 4);
    if strcmp(initState, 'singlet')
        % Singlet
        Sys.rho(2, 2) = 1;
        Sys.rho(2, 3) = -1;
        Sys.rho(3, 2) = -1;
        Sys.rho(3, 3) = 1;
        Sys.rho = 1/2*Sys.rho;
    elseif strcmp(initState, 'ciss')
        % Up-down and down-up
        Sys.rho(2, 2) = 1;
        Sys.rho(3, 3) = 1;
        Sys.rho = 1/2*Sys.rho;
    end

end

function Exp = getxfromrange(Exp)
    if isfield(Exp, 'Range')
        Exp.x = linspace(Exp.Range(1), Exp.Range(2), Exp.nPoints);
    elseif isfield(Exp, 'CenterSweep')
        Exp.x = linspace(Exp.CenterSweep(1) - Exp.CenterSweep(2)/2, ...
                         Exp.CenterSweep(1) + Exp.CenterSweep(2)/2, ...
                         Exp.nPoints);
    else
        error("Specify either Range or CenterSweep in the parameters.")
    end

    if isfield(Exp, 'Range') && isfield(Exp, 'CenterSweep')
        warning("Both Range and CenterSweep were specified in the " + ...
            "parameters. Range was used to calculate x.")
    end
end