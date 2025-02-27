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

sysPar = ["S", "g", "gFrame", "eeFrame", "A", "AFrame", "Nucs", "nNucs", ...
    "J", "dip", "lwpp", "initState", "chi", "lw1", "lw2", "trlwpp"];
expPar = ["x", "mwFreq", "CenterSweep", "Range", "nPoints", "Harmonic", ...
    "nThetas", "nPhis"];

Sys = struct();
Exp = struct();
M = deletecomments(M);
M = attachlineswithdots(M);

for ii = 1:numel(M)
    % Get parName
    temp = strsplit(M(ii), '=');
    if isscalar(temp)  % no '=' in this line
        continue    
    end
    parName = strrep(temp(1), ' ', '');
    % Compare with expected parameters
    [sysFlag, iSys] = ismember(parName, sysPar);
    [expFlag, iExp] = ismember(parName, expPar);
    if sysFlag
        par = sysPar(iSys);
    elseif expFlag
        par = expPar(iExp);
    else
        error("Parameter %s not found in the list " + ...
            "of expected parameters.\n", parName{:})
    end
    
    if sysFlag
        Sys = addpar(Sys, M(ii), par);
    else
        Exp = addpar(Exp, M(ii), par);
    end

end
Sys = getrhofromparamfile(Sys);
Exp = getxfromrange(Exp);
end

function M = deletecomments(M)
    for ii = 1:numel(M)
        temp = strsplit(M(ii), "%");
        M(ii) = temp(1);
    end
end

function M = attachlineswithdots(M)
    for ii = numel(M):-1:1
        if contains(M{ii}, '...')
            if ii == numel(M)
                error("The last line of the parameter file contains " + ...
                    "'...'. Import failed.")
            end
            % Remove dots
            temp = strsplit(M{ii}, '...');
            part1 = temp(1);
            % Append the rows
            M(ii) = append(part1, M(ii + 1));
            % Clean the row that was appended
            M(ii + 1) = "";
        end
    end
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