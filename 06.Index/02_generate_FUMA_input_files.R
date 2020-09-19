###################################
### FUMA - Generate input files ###
###################################
# Invoke the script in the commandline as follows: R CMD BATCH --slave 06_generate_FUMA_input_files.R &>> generate_FUMA_input_files.log
# or with Rscript: Rscript --slave 06_generate_FUMA_input_files.R &>> generate_FUMA_input_files.log

# add installed packages from home path to the library search path
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.6/")

# ...