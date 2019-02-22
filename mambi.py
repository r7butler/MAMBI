


#
#       Here is a script where we will attempt to run David Gillet's MAMBI.DJG function in the python environment
#
#           MAMBI.DJG is an R function, so likely it will be translated as MAMBI_DJG
#
#   Through the comments in this script, I am also trying to lay out a sort of documentation for a basic workflow for this idea, 
#       since it seems that we will be doing this same sort of thing quite often.



# Imports
import pandas as pd
import numpy as np
import sqlalchemy
import rpy2
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
import rpy2.robjects.packages as rpackages
from rpy2.robjects.packages import importr
import rpy2.rinterface as rinterface
r = robjects.r

# We need these functions to convert datatypes
as_character = r['as.character']
as_double = r['as.double']


# So importr did NOT work as is did with Marcus' code
# importr will only work if there is an installed R package.
# in this case, we need to import the function a bit differently

# What we do is read the file that contains the R function, and then pass it in as an argument to the rpy2.robjects.r function
# Here, the rpy2.robjects.r function is simply listed as r, as you can see what we did in the imports section
# When the string of r code is passed into the r function, it interprets the string as R code, and returns the function that was in the R code
#   and it enables us to use it
mambi_file = open("MAMBI.R", 'r')
MAMBI_DJG_string = mambi_file.read()
MAMBI_DJG = r(MAMBI_DJG_string)
mambi_file.close()

# Get the pandas dataframe that will be later passed into the MAMBI function
benthicdata = pd.read_excel("testbenthic.xlsx")

print benthicdata
print '\n'

# Here we need to convert to pandas dataframe to an R object that can be passed into the MAMBI function
r_benthicdata = pandas2ri.py2ri(benthicdata)


# It is ABSOLUTELY NECESSARY to convert the datatypes of the columns after converting a dataframe from pandas world to R world.
# Without converting the datatypes to what they need to be, the function will break with 100% certainty.

# Convert StationID to a character
r_benthicdata[0] = as_character(r_benthicdata[0])

# Convert Latitude to double
r_benthicdata[1] = as_double(r_benthicdata[1])

# Convert Longitude to double
r_benthicdata[2] = as_double(r_benthicdata[2])

# Convert Replicate to double
r_benthicdata[3] = as_double(r_benthicdata[3])

# Convert SampleDate to double
r_benthicdata[4] = as_double(r_benthicdata[4])

# Convert HabClassName to character
r_benthicdata[5] = as_character(r_benthicdata[5])

# Convert Salinity to double
r_benthicdata[6] = as_double(r_benthicdata[6])

# Convert Species to character
r_benthicdata[7] = as_character(r_benthicdata[7])

# Convert Abundance to double
r_benthicdata[8] = as_double(r_benthicdata[8])



print r_benthicdata
print '\n'


# And finally, we call the function
MAMBI_DJG(r_benthicdata)








# ---------------------------- #
# ------- CODE ARCHIVE ------- #
# ---------------------------- #

# This was me messing around with importing R functions into python
# I saved it because I may use this for documentation later
'''
testr = open("testr.R", 'r')
average = r(testr.read())
testr.close()
'''





