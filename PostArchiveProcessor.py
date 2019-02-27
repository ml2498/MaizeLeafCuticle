"""Apply initial processing operations to data tables in the raw archive. 

Use raw master tables in archives as starting point, and apply a set of changes.
Split data into separate streams for weight and leaf curling data, convert long-
form timestamps to Unix ones, and do genotype name-lookup and replacement for
the plot barcodes.

Like PhenotypeArchiver, this doesn't modify any data, but it does perform one QC
step: it verifies that all the readings from one sample have the same metadata
(applicable only to wet weights).

NOTE: I've eliminated the code for writing a new leaf curling table, as we've
decided we're probably not using that data. The processing code remains in place
in case it's decided that we still want to.

Author:         James Chamness
Last modified:  June 5, 2017
"""

"""Dependencies"""
import datetime
import logging
import os
import time
from importlib.machinery import SourceFileLoader
Common = SourceFileLoader("Common","../../pipeline/Common.py").load_module()

""""Set up logger routing to console and to file"""
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
logger.addHandler(logging.StreamHandler())
#logger.addHandler(logging.FileHandler('PostArchiveProcessor.log', mode='w'))

"""Apply processing to the master ww/lc table for the given environment.

    Arguments:
    environment -- e.g. AZ16, SD16, etc.
    dict -- a dictionary mapping plot barcodes to genotype names for the given environment
    wwMasterFilename -- filepath to input file for given environment
    wwOutputFilename -- filepath to output file for weight data
"""
def processMasterWWLCTable(environment, dict, wwMasterFilename, wwOutputFilename):
    
    """Initialize variables"""
    masterTable = Common.readTableFromFile(wwMasterFilename)
    masterWWTable = []
    #masterLCTable = []
    wwHeader = "Genotype\tLeaf_Number\tWeight\tTimestamp\tTimepoint\tClip_Weight\tScorer\tComputer\tScanner\tBalance\tSensor\tDate"
    #lcHeader = "Genotype\tLeaf_Number\tLeaf_Curl\tTimestamp\tTimepoint\tScorer\tComputer\tScanner\tBalance\tSensor"
    
    """Process table entries"""
    for entry in masterTable:
        if not entry[0] == "Chk":
            genotype = dict[entry[0]]
        else:
            genotype = "Chk_" + entry[2].replace("/","-") # append the date so that the same check dish used on different days does not appear as one
        wwTimestamp = str(getTimestamp2(entry[5],entry[6],environment))
        #lcTimestamp = str(getTimestamp2(entry[8],entry[9],environment))
        wwEntry = [genotype,entry[1],entry[4],wwTimestamp] + entry[10:] + [entry[2]]
        #lcEntry = [genotype,entry[1],entry[7],lcTimestamp,entry[10]] + entry[12:]
        masterWWTable.append(wwEntry)
        #masterLCTable.append(lcEntry)
    
    """Verify common metadata by sample"""
    logger.debug("Verifying metadata...")
    samples = set(list(map(lambda x: (x[0],x[1]), masterWWTable)))
    for sample in samples:
        found = False
        metadata = []
        for entry in masterWWTable:
            if sample[0] == entry[0] and sample[1] == entry[1]:
                if not found:
                    metadata = entry[5:]
                    found = True
                    continue
                else:
                    if not metadata == entry[5:]:
                        logger.debug("WARNING: conflicting metadata between sample readings!")
                        logger.debug("Sample: " + str(sample))
                        logger.debug("Terminating PostArchiveProcessor")
                        quit()
    logger.debug("...done.")
    
    """Write new tables to file"""
    logger.debug("Writing tables to file")
    Common.writeTableToFile(masterWWTable,wwOutputFilename,wwHeader)
    #Common.writeTableToFile(masterLCTable,lcOutputFilename,lcHeader)
    
    return

"""Apply processing to the master dw table for the given environment.

    Arguments:
    environment -- e.g. AZ16, SD16, etc.
    dict -- a dictionary mapping plot barcodes to genotype names for the given environment
    dwMasterFilename -- filepath to input file for given environment
    dwOutputFilename -- filepath to output file
"""
def processMasterDWTable(environment,dict,dwMasterFilename,dwOutputFilename):
    
    """Initialize variables"""
    masterTable = Common.readTableFromFile(dwMasterFilename)
    masterDWTable = []
    dwHeader = "Genotype\tLeaf_Number\tDry_Weight\tTimestamp\tScorer\tComputer\tScanner\tBalance\tClipWeight"
    
    """Process table entries"""
    for entry in masterTable:
        genotype = dict[entry[0]]
        dwTimestamp = str(getTimestamp1(entry[2],environment))
        dwEntry = [genotype,entry[1],entry[3],dwTimestamp] + entry[4:]
        masterDWTable.append(dwEntry)
    
    """Write new tables to file"""
    Common.writeTableToFile(masterDWTable,dwOutputFilename,dwHeader)
    
    return

"""Convert datetime from a given environment to seconds since Unix epoch.

    For use with dw spreadsheets.

    Arguments:
    dt -- a string of the form "6/28/2016 13:56:19"
    environment -- e.g. AZ16, SD16, etc.
"""
def getTimestamp1(dt, environment):
    
    # UTC offset -- number of hours behind UTC
    # note I think both Arizona and San Diego have same offset because Arizona doesn't observe daylight savings
    if environment[:2] == "AZ":
        utcOffset = 7
    elif environment[:2] == "SD":
        utcOffset = 7
    else:
        raise Exception
    
    date, time = dt.split(" ")
    
    # Construct a naive datetime object (NOT timezone aware)
    month,day,year = list(map(lambda x: int(x),date.split("/")))
    hour, minute, second = list(map(lambda x: int(x),time.split(":")))
    dateobj = datetime.datetime(year,month,day,hour,minute,second)
    
    # Construct and apply offset to convert naive datetime object to UTC
    delta = datetime.timedelta(hours=utcOffset)
    dateobj = dateobj + delta
    
    # Convert the UTC datetime object to the number of seconds since the Unix epoch
    epoch = datetime.datetime.utcfromtimestamp(0) # a UTC datetime object corresponding to the Unix epoch
    unix_time = int((dateobj - epoch).total_seconds()) # create a timedelta object and measure in seconds
    
    return unix_time

"""Convert datetime from a given environment to seconds since Unix epoch.

    For use with ww/lc spreadsheets.

    Arguments:
    date -- a string of the form "6/28/2016"
    time -- a string of the form "1:56:19 PM"
    environment -- e.g. AZ16, SD16, etc.
"""
def getTimestamp2(date, time, environment):
    
    # UTC offset -- number of hours behind UTC
    # note I think both Arizona and San Diego have same offset because Arizona doesn't observe daylight savings
    if environment[:2] == "AZ":
        utcOffset = 7
    elif environment[:2] == "SD":
        utcOffset = 7
    elif environment[:2] == "PT":
        utcOffset = 7
    else:
        raise Exception
    
    # Construct a naive datetime object (NOT timezone aware)
    month,day,year = list(map(lambda x: int(x),date.split("/")))
    time, meridian = time.split(" ")
    hour, minute, second = list(map(lambda x: int(x),time.split(":")))
    if meridian == "PM" and not hour == 12:
        hour += 12
    dateobj = datetime.datetime(year,month,day,hour,minute,second)
    
    # Construct and apply offset to convert naive datetime object to UTC
    delta = datetime.timedelta(hours=utcOffset)
    dateobj = dateobj + delta
    
    # Convert the UTC datetime object to the number of seconds since the Unix epoch
    epoch = datetime.datetime.utcfromtimestamp(0) # a UTC datetime object corresponding to the Unix epoch
    unix_time = int((dateobj - epoch).total_seconds()) # create a timedelta object and measure in seconds
    
    return unix_time

"""Executable"""
if __name__ == "__main__":
    
    """
    ============================================================================
    ==== CONFIGURATION
    ============================================================================
    """
    
    environment = "AZ16"
    environment = "SD16"
    environment = "PT17"
    environment = "AZ17"
    
    """
    ============================================================================
    ==== SCRIPT
    ============================================================================
    """
    
    logger.debug('PostArchiveProcessor launched ' + str(time.asctime()))
    
    """Path to the parent directory for the archive for the environment"""
    envArchivePath = Common.archivePath + os.sep + "raw" + os.sep + environment
    
    """Input paths for the environment-specific master tables"""
    wwMasterFilename = envArchivePath + os.sep + "MLC_" + environment + "_Master_WetWeights_LeafCurl.csv"
    dwMasterFilename = envArchivePath + os.sep + "MLC_" + environment + "_Master_DryWeights.csv"
    
    """Output paths for the environment-specific master tables"""
    wwOutputFilename = Common.archivePath + os.sep + "processed" + os.sep + "MLC_" + environment + "_Master_WetWeights_Processed.csv"
    dwOutputFilename = Common.archivePath + os.sep + "processed" + os.sep + "MLC_" + environment + "_Master_DryWeights_Processed.csv"
    #lcOutputFilename = Common.lcPath + os.sep + "MLC_" + environment + "_Master_LeafCurl_Processed.csv"
    
    """Look up the appropriate barcode dictionary"""
    dict = Common.getPlotBarcodeDict(environment)
    
    logger.debug('Processing master wet weight & leaf curling table...')
    processMasterWWLCTable(environment,dict,wwMasterFilename,wwOutputFilename)
    logger.debug("...done.")
    
    logger.debug('Processing master dry weight table...')
    processMasterDWTable(environment,dict,dwMasterFilename,dwOutputFilename)
    logger.debug("...done.")
    
    logger.debug('PostArchiveProcessor has finished running!')
    
