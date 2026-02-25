/*********************************************************
# Stata do file:    0.globals.do
#
# Author:      Sharon L Cadogan
#
# Date:        11/04/25

# Description: Globals for Infections and SMI study
#
***********************************************************/

/*Set Stata version number*/
version 17


*raw data
global rawdrivedir "Z:/GPRD_GOLD"
global rawdatadir"Z:\GPRD_GOLD\Georgia\Georgia_extract_2024_12_13\processed_data\smi_cohort\20250715\"

/*Analysis data and other files*/
global maindir "J:"
global ehrdir "$maindir/EHR-Working"
global projectdir "$ehrdir/Charlotte_Sharon\Infections and SMI"

/*do files*/
/**********************************************************
# PROJECT FOLDER FILEPATHS
********************************************************* */

global dodir "$projectdir/dofiles" /*Stata do files*/
global logdir "$projectdir/logfiles" /*Stata log files*/

global cleaning 

include "$dofiles/1.cleaning.do"
