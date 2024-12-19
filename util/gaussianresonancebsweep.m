function y1 = gaussianresonancebsweep(wReson, mwFreq, c, mode)
    arguments
        wReson
        mwFreq  
        c         
        mode string = "lwpp"
    end    

    normTerm = sqrt(2/pi)/c;
    y1 = normTerm*exp(-2 * (wReson - mwFreq).^2 ./ c^2);

end
