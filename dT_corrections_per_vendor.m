function [dT] = dT_corrections_per_vendor(dT,PARAMS)
    
    switch PARAMS.scanner_vendor
        case 'Philips'
            % here we make some simple corrections to the data due to the built-in
            % setup of Philips scanners, which invert the phase and perform fftshift
            % along a single dimension (not along both dimensions!).
            dT = -1*dT;
            dT = fftshift(dT,2);
    end