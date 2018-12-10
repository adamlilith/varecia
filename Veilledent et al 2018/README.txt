Copyright (C) March 2012:
Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>
Clovis Grinand <clovis@goodplanet.org>
Oliver Gardi <oliver.gardi@helvetas.org>

Reference: Vieilledent G., C. Grinand and R. Vaudry. 2013. Forecasting
deforestation and carbon emissions in tropical developing countries
facing demographic expansion: a case study in Madagascar. Ecology and
Evolution. DOI: 10.1002/ece3.550

For MS-Windows users
====================

Scripts for MS-Windows are available in the folder "scripts_windows". 

To use these scripts, you must install:
1. R for MS-Windows (http://cran.r-project.org/bin/windows/base/)
2. GRASS for MS-Windows (http://grass.osgeo.org/download/software/ms-windows/)
3. a GUI for R such as RStudio (http://www.rstudio.com/ide/download/desktop)
4. MSYS (a collection of GNU utilities such as bash, make, gawk and grep, http://www.mingw.org/wiki/MSYS)
5. ImageMagick for Windows (http://www.imagemagick.org/script/binary-releases.php#windows)

Don't forget to indicate the path to these programs in MS-Windows
environment variables (http://en.wikipedia.org/wiki/Path_(variable)):
ex. C:\Program Files\RStudio\bin;C:\msys\1.0\bin;C:\Program Files
(x86)\ImageMagick-6.8.3-Q16

Then, to use R/GRASS scripts
1. Start GRASS in command line (GRASS Command Line)
2. Select the appropriate location (ex. phcfM_SM) and mapset (ex. study_area_4)
3. Start your R GUI in the GRASS terminal
4. Run the scripts

For Linux users
===============

Scripts for Linux are available in the folder "scripts_linux". 

To use these scripts, you must install:
1. R for your Linux distribution (http://cran.r-project.org/bin/linux)
2. GRASS for your Linux distribution (http://grass.osgeo.org/download/software/linux/)
3. a GUI for R such as RStudio (http://www.rstudio.com/ide/download/desktop)
4. ImageMagick for your Linux distribution

Then, to use R/GRASS scripts
1. Start GRASS in command line
2. Select the appropriate location (ex. phcfM_SM) and mapset (ex. study_area_4)
3. Start your R GUI in the GRASS terminal
4. Run the scripts

For developers: script changes from Linux to Windows
====================================================

Scripts have been initially written under Linux. Some minor changes
were made in the scripts so that they can be run under Windows:
- The "system()" command must be replaced by "shell()"
- The apostrophe "'" must be replaced by quotation marks "\""
- When using mapcalc, a same raster cannot be used for input and output:
example: 
system("r.mapcalc 'a=if(isnull(a),0,1)'") 
must be replaced by:
shell("r.mapcalc \"atemp=if(isnull(a),0,1)\"")
shell("r.mapcalc \"a=atemp\"")
- Call to awk is done via the command "gawk"



