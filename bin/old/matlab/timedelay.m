function td=timedelay(ra,source_declination_radians, gpstime, detector)

if (ischar(detector))
   %----- Load data on detector.
   DetData = LoadDetectorData(detector);
   d = DetData.V;
elseif (isnumeric(detector) && isequal(size(detector),[3 3]))
   d = detector;
else
   error(['Detector not recognized. 4th argument should be a ' ... 
       'detector/site name or a 3x3 array.']);  
end

gmst = gps2sidereal(gpstime);

greenwich_hour_angle=gmst-ra;

earth_center = [0.0, 0.0, 0.0];

ehat_src(1) = cos(source_declination_radians) * cos(greenwich_hour_angle);
ehat_src(2) = cos(source_declination_radians) * -sin(greenwich_hour_angle);
ehat_src(3) = sin(source_declination_radians);

%position of detector 2 with respect to detector 1
	

delta_xyz(1) = earth_center(1) - d(1);
delta_xyz(2) =  earth_center(2) - d(2);
delta_xyz(3) = earth_center(3) - d(3);
	
% Arrival time at detector 1 - arrival time at detector 2.  This
% is positive when the wavefront arrives at detector 1 after
% detector 2 (and so t at detector 1 is greater than t at detector
% 2).
	 
dp=dot(ehat_src, delta_xyz);
td= dp / 3e8;