function [Fp, Fc] = det_response_wrapper(ra, dec, det, gpstime, psi)

% [Fp, Fc] = det_response_wrapper(ra, dec, det, gpstime, psi)
%
% This function is a wrapper around Patrick Sutton's ComputeAntennaResponse
% function. The inputs are:
%   ra - the right ascension given as either a string format 'hh:mm:ss.s',
%       or as an angle in radians (this can be a vector)
%   dec - the declination given as either a string format 'dd:mm:ss.s', or
%       as an angle in radians (this can be a vector)
%   det - a string containing the required detector name. Valid detector
%       names are: 'H1' for the LIGO Hanford detectors, 'G1' for the 
%       GEO600 detector, 'V1' for the Virgo detector. For other available
%       detectors see LoadDetectorData.m
%   gpstime - the required time in GPS seconds (this can be a vector of
%       times)
%   psi - the gravitational wave polarisation angle in radians. This should
%       be between -pi/4 and pi/4
%
% The function will return the plus and cross polarisation antenna response
% functions (Fp and Fc) for the given sky location, detector and
% polarisation angle.

% if RA and Dec are strings convert into radians;
if ischar(ra) && ischar(dec)
    % check strings of the form 'XX:XX:XX'
    dr = strfind(ra, ':');
    dd = strfind(dec, ':');
    
    if length(dr) == 2 && length(dd) == 2
        [rar, decr] = radec_to_rads(ra, dec);
    else
        error('Error... RA and dec strings are not in the form XX:XX:XX');
    end
elseif isnumeric(ra) && isnumeric(dec)
    rar = ra;
    decr = dec;
end
    
% convert gps time to GMST
gmst = gps2sidereal(gpstime);

% convert gmst from hours into rads
gmst = (gmst/24)*2*pi;

% convert ra to sidereal hour angle for given time (at GMST) - runs in
% opposite east to west direction than regular hour angle
hr = rar - gmst;

% convert dec to polar angle needed by ComputeAntennaResponse (in which the
% north pole is 0 rads and the south pole is pi rads
theta = -decr + pi/2;

% run detector repsonse function
[Fp, Fc] = ComputeAntennaResponse(hr,theta,psi,det);