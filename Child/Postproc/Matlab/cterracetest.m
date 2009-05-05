function result = cterracetest( fname, ts, x, y, min_drarea, ...
                                min_age, max_age, min_elev, max_elev )
% cterracetest: Determines whether the closest point to (x,y) in a CHILD
% run passes through the elevation range (min_elev,max_elev) during the
% time period (min_age,max_age). This is used to test whether a main stream
% (with drainage area > min_drarea) passes through the known altitude of a
% strath terrace within the right time span.
%
% GT, Apr 09

% calculate the elevation time series
[t,z]=cpointelevtseries(fname,x,y,ts,min_drarea);

% convert ages into simulation time
start_age = t(end);
min_sim_time = start_age - max_age;
max_sim_time = start_age - min_age;

% interpolate the time series to 1-year level
ti = 0:start_age;
zi = interp1( t, z, ti, 'linear' );

% set up time-elevation "polygon" (rectangle)
time_box = [min_sim_time max_sim_time max_sim_time min_sim_time];
elev_box = [min_elev min_elev max_elev max_elev];

% test whether any of the points lie within the polygon (rectangle) defined
% by ages and elevations
in_time_elev_box = inpolygon(ti,zi,time_box,elev_box);
result = max(in_time_elev_box)>0;
