%  
% the procedure for running the mht (multiple hypothesis testing) code

1. prepare the input file (the main code is 'mht_test.m')

a). temperature data for all acquisitions: download from 
KNMI using get_meteo.m, save as .mat file (variable 'temp': temperature difference between slaves and master)

b). project .mat file: This should be the standard output file from the DePSI PSI toolbox. But the enduser can also manufacture it, and manually create the matrix of perpendicular baseline 'Bperp', temporal baseline 'Btemp', acquisition dates 'dates',wavelength 'lambda' and project name 'project_id'.
See an example 'limburg_tsx_project.mat'

c). critical value and the level of significance for all alternative hypotheses --> K_a.mat: run the code K_alpha.m
The project file needs to be defined in this code. 

d). deformation time series of all points, variable 'defo' in mht_test.m

e). variance of unit weight (vuw): e.g. 3^2 mm^2 for terrasar-x data, or it can be defined in terms of Eq.3.28, Page 43 in my PhD thesis 'Monitoring civil infrastructure using satellite radar interferometry' https://repository.tudelft.nl/islandora/object/uuid%3Af4c6a3a2-73a8-4250-a34f-bc67d1e34516.

2. MHT Output
   the parameter estimation for each point is saved in '..defo_model_param.csv' file
CSV Format: col = column, H0 is default, linear model, Ha is the best alternative model
----------------------------------------------
col_1  col_2     col_3        col_4        
No.    ModelNo. velocity(H0)  velocity(Ha)
----------------------------------------------

ModelNo.
----------------------------------------------
ModelNo.1 is H0, linear model
col_4 col_5 col_6 col_7 col_8 col_9
 0      0    0     0      0    0
----------------------------------------------
ModelNo.2 is velocity + temperature-related
 col_5              col_6 col_7 col_8 col_9
 temperature_param   0    0      0     0
----------------------------------------------
ModelNo.3 is velocity + temperature-related + heaviside (1 offset)
col_5              col_6    col_7     col_8             col_9
temperature_param   offset    intercept 0      offset_acquisition
----------------------------------------------
ModelNo.4 is velocity + heaviside (1 offset)
----------------------------------------------
col_5     col_6    col_7     col_8          col_9
0         offset   intercept 0      offset_acquisition
----------------------------------------------
ModelNo.5 is two velocities + breakpoint
----------------------------------------------
col_4       col_5         col_9
velocity1 velocity2-velocity1 the breakpoint acquisition

The degree of freedom is always recommended to be double checked. 

%%%%%%%%%%%%
% 1 create Btemp 
fileIn = 'limburg_tsx_dates.txt';
dates = load(fileIn);
dates = datenum(num2str(dates),'yyyymmddHH');
Btemp = (dates(2:end) - dates(1))./365;
% 2 define Bperp
 Bperp = ones(length(Btemp),1);
%%%% since I don't know the Bperp, I assumed they are 1
% 3. lambda, wavelength [m]
lambda = 0.031
% 4. project_id 
 project_id = 'limburg_tsx';
% 5. dates, string

dates = datestr(dates,'dd-mmm-yyyy');
dates = [dates(2:end,:); dates(1,:)];

save([fileIn(1:12) 'project.mat'], 'Btemp','Bperp','dates','lambda','project_id');






