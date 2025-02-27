function energ = myeigenenergies(w0, deltaw, J, d)
    
    signJmd(1, 1, :) = [-1, 1, -1, 1];  % (1, 1, 4)
    signOmega(1, 1, :) = [-1, -1, 1, 1];  % (1, 1, 4)
    
    Omega = hypot(deltaw, J + d/2);
    energ = w0 + pagemtimes(signJmd, J - d) + pagemtimes(signOmega, Omega);
end

