function Sys = getrhofromparamfile(Sys)
    
    Sys.rho = zeros(4, 4);
    if ~isfield(Sys, 'initState') && ~isfield(Sys, 'chi')
        warning("The initial density matrix is set to 'singlet'" + ...
            " since no initState nor chi was declared in the parameters.")
        chi = 0;
    elseif isfield(Sys, 'initState') && isfield(Sys, 'chi')
        warning("Both initState and chi are declared in the " + ...
            "parameters. initState is used to construct rho. chi is " + ...
            "ingored.")
    end
    if isfield(Sys, 'initState')  % initState used for rho
        initState = lower(Sys.initState);
        if strcmp(initState, 'singlet')
            chi = 0;
        elseif strcmp(initState, 'ud')
            chi = 90;
        elseif strcmp(initState, 'du')
            chi = -90;
        elseif strcmp(initState, 'ciss')
            % Up-down and down-up
            Sys.rho(2, 2) = 1;
            Sys.rho(3, 3) = 1;
            Sys.rho = 1/2*Sys.rho;
            return
        else
            error("initState not implemented.")
        end
    else  % only chi is present
        chi = Sys.chi;
    end

    rhoS = zeros(4, 4);
    rhoS(2, 2) = 1;
    rhoS(2, 3) = -1;
    rhoS(3, 2) = -1;
    rhoS(3, 3) = 1;
    rhoS = 1/2*rhoS;
    rhoT0 = abs(rhoS);
    chi = chi*pi/180;  % rad
    Sys.rho = cos(chi/2)*rhoS + sin(chi/2)*rhoT0; 

end

