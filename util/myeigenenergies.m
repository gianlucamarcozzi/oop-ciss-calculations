function energ = myeigenenergies(w0, deltaw, J, d, itrans)
    
    signJmd = [-1, 1, -1, 1];
    signOmega = [-1, -1, 1, 1];
    
    Omega = hypot(deltaw, J + d/2);
    energ = w0 + signJmd(itrans)*(J - d) + signOmega(itrans)*Omega;
end

