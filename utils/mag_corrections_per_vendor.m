function [im_mag] = mag_corrections_per_vendor(im_mag,PARAMS)
    
    switch PARAMS.scanner_vendor
        case 'Philips'
            % here we make some simple corrections to the data due to the built-in
            % setup of Philips scanners, which invert the phase and perform fftshift
            % along a single dimension (not along both dimensions!).
            im_mag = fftshift(im_mag,2);
    end