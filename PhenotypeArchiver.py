"""Archives raw CE data from the original scoring sheets on the Google Drive.

This script does NOT modify any data or perform any QC. It serves only to
archive raw data from the cloud to local storage, e.g. the GoreLab server.  

Depending on the selected function, this script may perform a few different
operations. All of these require selecting the "environment" (e.g. AZ16) to work
with. This is done by setting the only variable shown in the section of the
__main__ executable entitled "Configuration".

For every environment, the script provides a way to archive data by day (or for
a list of days). For environments where we're done taking data, I also wrote a
few convenience functions to archive every sheet. There are also functions to
recompile the master tables from the local summary tables. See the descriptions
at the bottom of this script.

Author:         James Chamness
Last Modified:  06/05/2017
"""

"""Dependencies"""
import gspread
import logging
from oauth2client.service_account import ServiceAccountCredentials
import sys
import time
import datetime
from datetime import datetime
import os
from importlib.machinery import SourceFileLoader
Common = SourceFileLoader("Common","../../pipeline/Common.py").load_module()

""""Set up logger routing to console and to file"""
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
logger.addHandler(logging.StreamHandler())
logger.addHandler(logging.FileHandler('PhenotypeArchiver.log', mode='w'))

"""Update the master tables with all data on the Drive collected on the given date.
    
    Scoring sheets should follow these naming conventions:
        MLC_AZ_UserInitials_mmddyyyy_computerNumber (for wet weight/curling)
        MLC_AZ_UserInitials_mmddyyyy_computerNumber_DW (for dry weight)
    
    If any mistakes are found in the scoring sheets from the given date, the
    update attempt will be abandoned and no changes made to local archive.
    
    Arguments:
    oauthKeyFilename -- the filename of the .json authentication file for oauth2 client validation
    environment -- the 4-char string identifying the environment: "AZ16", "SD17", etc.
    archivePath -- the path of the directory in which output files will be written
    wwMasterFilename -- the filename of the file containing the master wet weight table
    dwMasterFilename -- the filename of the file containing the master dry weight table
    date -- date for which to fetch data; mmddyyyy format, including leading 0s
    includeDailySummaryTables -- if true, write a local summary file for each available spreadsheet
"""
def updateFromCloudForDate(oAuthKeyFilename, environment, archivePath, wwMasterFilename, dwMasterFilename, date, includeDailySummaryTables):
    
    logger.debug("Archiving phenotypes from Google Drive for " + date)
    
    ### Get the ww and dw tables from the Drive
    client = getDriveClient(oAuthKeyFilename)
    genotypeBarcodeDictionary = Common.getPlotBarcodeDict(environment)
    wwlcSpreadsheets, dwSpreadsheets = getAllScoringSpreadsheetsFromDate(client, environment, date)
    if len(wwlcSpreadsheets) == 0 and len(dwSpreadsheets) == 0:
        logger.debug("No spreadsheets found for " + date)
        return
    
    logger.debug("Parsing data...")
    wwTable = parseWetWeightTablesFromSpreadsheets(wwlcSpreadsheets, genotypeBarcodeDictionary, date)
    dwTable = parseDryWeightTablesFromSpreadsheets(dwSpreadsheets, genotypeBarcodeDictionary, date)
    if not wwTable == None and not dwTable == None:
        logger.debug("...done.")
    else:
        logger.debug("Parsing failed for at least one sheet! Archive attempt for " + date + " terminated.")
        return
     
    ### Sync these with the master tables on file (and write day-separated table
    ### files if applicable)
    logger.debug("Syncing data with local archive...")
    if wwTable != []:
        header = "Plot\tLeaf_Number\tDate\tTime\tWeight\tDate\tTime\tLeaf_Curl\tDate\tTime\tTime_Point\tClip_Weight\tScorer\tComputer\tScanner\tBalance\tSensor"
        if not os.path.isfile(wwMasterFilename):
            logger.debug("No master ww table file found for " + environment + " archive: writing as new master table file...")
            Common.writeTableToFile(wwTable,wwMasterFilename,header)
            logger.debug("...done.")
        else:
            logger.debug("Adding " + date + " data to master table for " + environment + "...")
            masterWWTable = Common.readTableFromFile(wwMasterFilename, header=True)
            newMasterWWTable = updateMasterWWTableWithNewTable(masterWWTable, wwTable)
            Common.writeTableToFile(newMasterWWTable, wwMasterFilename, header)
            logger.debug("...done.")
        if includeDailySummaryTables:
            logger.debug("Writing single-day ww table file for " + date + "...")
            wwTable = sorted(wwTable, key=lambda x: (x[0], x[1]))
            wwFilename = archivePath + os.sep + "MLC_" + environment + "_" + date + "_WetWeights_LeafCurl_Summary.csv"
            Common.writeTableToFile(wwTable, wwFilename, header)
            logger.debug("...done.")
        logger.debug("All wet weight/leaf curling data synced for " + date)
    else:
        logger.debug("No wet weight/leaf curling data to sync for " + date)
         
    if dwTable != []:
        header = "Plot\tLeaf_Number\tTimestamp\tDry_Weight\tScorer\tComputer\tScanner\tBalance\tClipWeight"
        if not os.path.isfile(dwMasterFilename):
            logger.debug("No master dw table file found for " + environment + " archive: writing as new master table file...")
            Common.writeTableToFile(dwTable, dwMasterFilename, header)
            logger.debug("...done.")
        else:
            logger.debug("Adding " + date + " data to master table for " + environment + "...")
            masterDWTable = Common.readTableFromFile(dwMasterFilename, header=True)
            newMasterDWTable = updateMasterDWTableWithNewTable(masterDWTable, dwTable)
            Common.writeTableToFile(newMasterDWTable, dwMasterFilename, header)
            logger.debug("...done.")
        if includeDailySummaryTables:
            logger.debug("Writing single-day dw table file for " + date + "...")
            dwTable = sorted(dwTable, key=lambda x: (x[0], x[1]))
            dwFilename = archivePath + os.sep + "MLC_" + environment + "_" + date + "_DryWeights.csv"
            Common.writeTableToFile(dwTable, dwFilename, header)
    else:
        logger.debug("No dry weight data to sync for " + date)    
     
    logger.debug("Phenotypes successfully archived from Google Drive for " + date)
    
"""Update the master tables with all data in a list of local filenames.

    This is much faster than pulling them from the cloud.

    Arguments:
    localFiles -- a list of filenames containing (preformatted) data tables
        follow naming convention: "MLC_AZ16_06132016_JC_1", MLC_AZ16_06102016_JC_1_DW", e.g.
    wwMasterFilename -- the filename of the file containing the master wet weight table
    dwMasterFilename -- the filename of the file containing the master dry weight table
"""
def updateFromLocalSummaryTableFiles(localFiles, wwMasterFilename, dwMasterFilename):
    
    logger.debug("Updating master tables from the following local day-summary files:")
    for filename in localFiles:
        logger.debug("    " + filename)
    
    wwTable = []
    dwTable = []
    for filename in localFiles:
        if filename.split("_")[len(filename.split("_"))-1] != "DryWeights.csv":
            if not os.path.isfile(wwMasterFilename):
                logger.debug("No master ww table file found for " + environment + " archive: writing as new master table file...")
                masterWWTable = []
            else:
                masterWWTable = Common.readTableFromFile(wwMasterFilename)
            logger.debug("Adding " + filename + " data to master table for " + environment + "...")
            header = "Plot\tLeaf_Number\tDate\tTime\tWeight\tDate\tTime\tLeaf_Curl\tDate\tTime\tTime_Point\tClip_Weight\tScorer\tComputer\tScanner\tBalance\tSensor"
            newMasterWWTable = updateMasterWWTableWithNewTable(masterWWTable, Common.readTableFromFile(filename))
            Common.writeTableToFile(newMasterWWTable, wwMasterFilename, header)
            logger.debug("...done.")
        else:
            if not os.path.isfile(dwMasterFilename):
                logger.debug("No master dw table file found for " + environment + " archive: writing as new master table file...")
                masterDWTable = []
            else:
                masterDWTable = Common.readTableFromFile(dwMasterFilename)
            logger.debug("Adding " + filename + " date to master table for " + environment + "...")
            header = "Plot\tLeaf_Number\tTimestamp\tDry_Weight\tScorer\tComputer\tScanner\tBalance\tClipWeight"
            newMasterDWTable = updateMasterDWTableWithNewTable(masterDWTable, Common.readTableFromFile(filename))
            Common.writeTableToFile(newMasterDWTable, dwMasterFilename, header)
            logger.debug("...done.")
    
    logger.debug("Master tables updated.")

"""Return two lists: all wet weight/leaf curling, and all dry weight scoring spreadsheets from the given date.

    Arguments:
    client -- a gspread client for interacting with the Google Drive
    environment -- the 4-char string identifying the environment: "AZ16", "SD17", etc.
    date -- date for which to fetch data; mmddyyyy format, including leading 0s
"""
def getAllScoringSpreadsheetsFromDate(client, environment, date):
    
    logger.debug("Gathering scoring spreadsheets from " + date + "...")
    
    wwlcSpreadsheets = []
    dwSpreadsheets = []
    all_available_spreadsheets = client.openall()
    if not all_available_spreadsheets:
        logger.debug("No spreadsheets found.")
    else:
        logger.debug("    Spreadsheets found:")
    
    loc = environment[0:2]
    
    for spreadsheet in all_available_spreadsheets:
        title_els = spreadsheet.title.strip().split("_")
        if title_els[0] == "MLC" and title_els[1] == loc:
            spreadsheetDate = title_els[2]    
            if spreadsheetDate == date:
                if title_els[len(title_els)-1] != "DW":
                    wwlcSpreadsheets.append(spreadsheet)
                    #pass
                else:
                    dwSpreadsheets.append(spreadsheet)
                    #pass
                logger.debug("        " + spreadsheet.title)
    
    return(wwlcSpreadsheets, dwSpreadsheets)

"""Return a table of wet weight and leaf curling data parsed from a list of spreadsheets.

    Consolidate wet weight and leaf curling data from all the scoring sheets
    into a single table. If a mistake is found in any of these sheets, abandon
    consolidation attempt and return None (the function called to consolidate
    that sheet will report the mistake).
    
    If any expected worksheets are missing from any of the spreadsheets, report
    so, abandon consolidation attempt, and return None.
    
    Arguments:
    wwlcSpreadsheets -- a list of gpsread spreadsheets of wet weight & leaf curling data
    genotypeBarcodeDictionary -- a dictionary mapping line names to plot barcodes
    date -- date for which to fetch data; mmddyyyy format, including leading 0s
"""
def parseWetWeightTablesFromSpreadsheets(wwlcSpreadsheets, genotypeBarcodeDictionary, date):
    
    wwlcTable = [] # summary table containing data from all sheets that scan without problems
    clipWeightsUsedWithBarcodes = {} # structure common between all sheets
    for spreadsheet in wwlcSpreadsheets:
        
        wwlcTableForSheet = parseWetWeightTableFromSpreadsheet(spreadsheet, genotypeBarcodeDictionary, date)
        
        ### if there is a mistake with any sheet, abandon entire update attempt
        if wwlcTableForSheet == None:
            return None
        
        ### update the dictionary mapping barcodes to clip weights
        clipWeightSheet = spreadsheet.worksheet("ClipWeights")
        clipWeightsUsedWithBarcodes = parseClipWeightsFromWorksheet(clipWeightsUsedWithBarcodes, clipWeightSheet, genotypeBarcodeDictionary, date, spreadsheet)
        
        ### if there is a mistake with the ClipWeights sheet, abandon entire update attempt
        if clipWeightsUsedWithBarcodes is None:
            logger.debug("Exception occurred while attempting to parse ClipWeights sheet!")
            return None
        
        ### if there were no errors with the table for the current sheet, add it to the summary table
        if wwlcTableForSheet is not None:        
            wwlcTable = wwlcTable + wwlcTableForSheet
            logger.debug(spreadsheet.title + " added to summary for " + date)
            
    ### now that all the data has been read, add a field for clip weight to the summary table
    ### if unknown, use NA
    updatedWWLCTable = []
    for entry in wwlcTable:
        if len(entry) < 16: # only need to add the clip weight if not done already
            #if clipWeightsUsedWithBarcodes.has_key(entry[0]): # grammar for python2
            if entry[0] in clipWeightsUsedWithBarcodes.keys(): # modified for python3 grammar
               #updatedWWLCTable.append(entry+[clipWeightsUsedWithBarcodes[entry[0]]])
               updatedWWLCTable.append(entry[:10] + [clipWeightsUsedWithBarcodes[entry[0]]] + entry[10:] )
            else:
                #updatedWWLCTable.append(entry+["NA"])
                updatedWWLCTable.append(entry[:10] + ["NA"] + entry[10:] )
    wwlcTable = updatedWWLCTable
    
    wwlcTable = list(map(lambda x: x[0].split("_") + x[1:], wwlcTable)) # split plot barcode from leaf number
    
    return wwlcTable

"""Return a table of wet weight and leaf curling data parsed from a spreadsheet.

    Consolidate wet weight and leaf curling data from the given scoring sheet
    into a single table. Also, parse the "User" sheet and associate a username,
    computer, scanner, balance, and sensor with every reading. If any error is
    found in the given spreadsheet, report the error, abandon consolidation
    attempt, and return None.
    
    Does NOT parse clip weight data (this is done separately).
    
    Checks the titles of the worksheets in each spreadsheet against those that
    are expected. If an expected worksheet is missing, report so, abandon
    consolidation attempt, and return None. 

    Arguments:
    wwlcSpreadsheet -- a gspread spreadsheet with wet weight & leaf curling data
    genotypeBarcodeDictionary -- a dictionary mapping line names to plot barcodes
    date -- date for which to fetch data; mmddyyyy format, including leading 0s
"""
def parseWetWeightTableFromSpreadsheet(spreadsheet, genotypeBarcodeDictionary, date):
    
    logger.debug("    Parsing data from: " + spreadsheet.title)
    
    wwlcTableForSheet = []
    abandonSheet = False # if any problems are found with the current sheet, abandon attempt to add it
    
    expectedWorksheetTitles = ["User","Time1","Time2","Time3","Time4","Time5","ClipWeights"]
    wwSheetTitles = expectedWorksheetTitles[1:-1]
    
    ### Make sure all expected worksheets appear
    for worksheetTitle in expectedWorksheetTitles:
        if not worksheetTitle in list(map(lambda x: x.title, spreadsheet.worksheets())):
#             """Hardcoding the following dates, the first few experiment days when we weren't taking clipweights"""
#             if date in ["06062016", "06072016", "06092016"]:
#                 continue
#             """Otherwise, assume there must be a clipweight sheet"""
            logger.debug("Unable to find " + worksheetTitle + " in " + spreadsheet.title + ": abandoning attempt to update for " + date)
            return None
        
    ### Parse each worksheet based on the title
    for worksheetTitle in wwSheetTitles:
        worksheet = spreadsheet.worksheet(worksheetTitle)
        worksheetVals = worksheet.get_all_values()[1:]
        dataEntryColumn = list(map(lambda x: x[1], worksheetVals[1:])) # selects B2:BXXXX
        ### figure out index of the row of the last data point
        indexOfLastDatapoint = 0
        for i in range(0,len(dataEntryColumn)):
            if dataEntryColumn[i] == "":
                indexOfLastDatapoint = i-1
                break
        ### screen for empty sheets
        if indexOfLastDatapoint == -1:
            logger.debug("    " + spreadsheet.title + ":" + worksheet.title + " is an empty sheet, skipping")
            break
        ### iterate over the "data" rows and add each datapoint to the table
        for i in range(0,((indexOfLastDatapoint+1)//3)+1): #modified for python3 grammar
            indexOfSummaryLine = i*3
            newEntry = worksheetVals[indexOfSummaryLine][6:16]
             
            """Quality-check each component of the datapoint entry. If a mistake
            is found, abandon attempt to add the current spreadsheet to the day
            summary"""
             
            ### Check for missing data
            for val in newEntry:
                if val == "":
                    logger.debug("Missing value found in "+spreadsheet.title+", "+worksheet.title+", line "+str(indexOfSummaryLine+2))
                    abandonSheet = True
                    break
                     
            ### Check for valid barcode (excluding water dish checks)
            if not newEntry[0].split("_")[0] == "Chk":
                if not newEntry[0].split("_")[0] in genotypeBarcodeDictionary.keys():
                    logger.debug("Unknown or improperly formatted barcode found in "+spreadsheet.title+", "+worksheet.title+", line "+str(indexOfSummaryLine+2) + ": " + newEntry[0].split("_")[0])
                    abandonSheet = True
                    break
             
            ### Check for valid dates
            sheetDate = str(int(date[0:2])) + "/" + str(int(date[2:4])) + "/" + date[4:]
            if not newEntry[1] == sheetDate or not newEntry[4] == sheetDate or not newEntry[7] == sheetDate:
                logger.debug("Date for datapoint does not match sheet date, or is improperly formatted, in "+spreadsheet.title+", "+worksheet.title+", line "+str(indexOfSummaryLine+2))
                abandonSheet = True
                break
             
            ### Check for valid times
            if not validateTimeString(newEntry[2], True) or not validateTimeString(newEntry[5], True) or not validateTimeString(newEntry[8], True):
                logger.debug("Time for datapoint is improperly formatted, in "+spreadsheet.title+", "+worksheet.title+", line "+str(indexOfSummaryLine+2))
                abandonSheet = True
                break
             
            ### This step will not be reached unless datapoint values check out
            wwlcTableForSheet.append(newEntry)
     
        if abandonSheet:
            logger.debug("Abandoning attempt to consolidate wet weight data for " + date)
            logger.debug("Please correct mistaken or missing value(s) and try again")
            return None
    
    ### Parse metadata from "User" sheet and add to every reading from the sheet
    worksheet = spreadsheet.worksheet("User")
    try:
        vals = list(map(lambda x: x[1], worksheet.get_all_values())) # should have 5 values
        metadata = []
        for i in range(0,5):
            if not i+1 > len(vals) and not vals[i] == []:
                metadata.append(vals[i])
            else:
                logger.debug("WARNING: " + spreadsheet.title + " is missing some metadata!")
                logger.debug("    Substituting NA")
                metadata.append("NA")
                if i == 0:
                    missing = "name"
                elif i == 1:
                    missing = "computer"
                elif i == 2:
                    missing = "scanner"
                elif i == 3:
                    missing = "balance"
                else:
                    missing = "sensor"
    except IndexError:
        logger.debug("WARNING: no metadata recorded for " + spreadsheet.title + "!")
        logger.debug("    Substituting NA")
        metadata = ["NA","NA","NA","NA","NA"]
    
    wwlcTableForSheet = list(map(lambda x: x + metadata, wwlcTableForSheet))
    
    return wwlcTableForSheet

"""Return true if the given timestring is formatted correctly.

    Should be formatted as '16:15:47', e.g. This function attempts to parse
    the timestring as it would be to generate a timestamp; if this fails for any
    reason, return False.

    Arguments:
    timeString -- a timestamp string as recorded in the scoring spreadsheets
    asTwelveHour -- if True, expected format is instead '1:15:47 PM', e.g.
"""
def validateTimeString(timeString, asTwelveHour=False):
    
    els = timeString.split(":")
    if not len(els) == 3:
        return False
    
    # hour must be zero-padded for parsing
    if len(timeString.split(":")[0]) == 1:
        timeString = "0" + timeString    
    try:
        if asTwelveHour:
            t = datetime.strptime(timeString, "%I:%M:%S %p")
        else:
            t = datetime.strptime(timeString, "%H:%M:%S")
    except ValueError:
            return False
    return True

"""Update and return the barcode->clip weight dictionary with the data from the given worksheet.

    If any error is found in the clip weight sheet, report the error, and return
    None.
    
    The same dictionary object must be shared between spreadsheets, because it
    may be that the clip weight for a given sample was recorded on a different
    sheet than its wet weight/leaf curling data.
    
    Arguments:
    clipWeightsUsedWithBarcodes -- a dictionary mapping sample barcodes to clip weights
    clipweightSheet -- the worksheet (not spreadsheet) containing the clip weights
    genotypeBarcodeDictionary -- a dictionary mapping line names to plot barcodes
    date -- date for which to fetch data; mmddyyyy format, including leading 0s
    spreadsheet -- the spreadsheet document containing the clipWeightSheet worksheet
"""
def parseClipWeightsFromWorksheet(clipWeightsUsedWithBarcodes, clipWeightSheet, genotypeBarcodeDictionary, date, spreadsheet):
    
    #logger.debug("Scanning " + spreadsheet.title + ":" + clipWeightSheet.title)
    abandonSheet = False
    worksheetVals = clipWeightSheet.get_all_values()[1:]
    dataEntryColumn = list(map(lambda x: x[1], worksheetVals[1:])) # selects B2:BXXXX
    ### figure out row of the last data point
    indexOfLastDatapoint = 0
    for i in range(0,len(dataEntryColumn)):
        if dataEntryColumn[i] == "":
            indexOfLastDatapoint = i-1
            break
    ### call just the rows needed
    #for i in range(0,((indexOfLastDatapoint+1)/2)+1):
    for i in range(0,((indexOfLastDatapoint+1)//2)+1): # modified for python3 grammar
        indexOfSummaryLine = i*2
        vals = worksheetVals[indexOfSummaryLine]
        
        """Quality-check each component of the datapoint entry. If a mistake
        is found, abandon attempt to update clipweight dictionary."""
          
        ### Check for missing data
        for val in vals:
            if val == "":
                logger.debug("Missing value found in "+spreadsheet.title+", "+clipWeightSheet.title+", line "+str(indexOfSummaryLine+2))
                abandonSheet = True
                break
                    
        ### Check for valid barcode
        if not vals[3].split("_")[0] in genotypeBarcodeDictionary.keys():
            logger.debug("Unknown or improperly formatted barcode found in "+spreadsheet.title+", "+clipWeightSheet.title+", line "+str(indexOfSummaryLine+2) + ": " + vals[3].split("_")[0])
            abandonSheet = True
            break
        """Commenting this out because we don't really care about timestamps for clip weights.    
        ### Check for valid dates
        sheetDate = str(int(date[0:2])) + "/" + str(int(date[2:4])) + "/" + date[4:]
        timestamp = vals[2].split(" ")
        if not len(timestamp) == 2:
            logger.debug("Timestamp for datapoint is improperly formatted, in "+spreadsheet.title+", "+clipWeightSheet.title+", line "+str(indexOfSummaryLine+2))
            abandonSheet = True
            break
        if not timestamp[0] == sheetDate:
            logger.debug("Date for datapoint does not match sheet date, in "+spreadsheet.title+", "+clipWeightSheet.title+", line "+str(indexOfSummaryLine+2))
            abandonSheet = True
            break
            
        ### Check for valid times
        if not validateTimeString(timestamp[1]):
            logger.debug("Time for datapoint is improperly formatted, in "+spreadsheet.title+", "+clipWeightSheet.title+", line "+str(indexOfSummaryLine+2))
            abandonSheet = True
            break
        """
        ### This step will not be reached unless datapoint values check out
        barcode = vals[3]
        clipWeight = vals[4]
        clipWeightsUsedWithBarcodes[barcode] = clipWeight
    
    if abandonSheet:
        logger.debug("Abandoning attempt to consolidate clip weight data for " + date)
        logger.debug("Please correct mistaken or missing value(s) and try again")
        return None
    
    return clipWeightsUsedWithBarcodes
        
"""Return a table of all dry weight data parsed from a list of spreadsheets.

    Consolidate dry weight data from all the scoring sheets into a single table.
    If a mistake is found in any of these sheets, the consolidation attempt will
    be abandoned (the function called to consolidate that sheet will report the
    mistake), and None will be returned.

    Arguments:
    dwSpreadsheets -- a list of gspread spreadsheets with dry weight data
    genotypeBarcodeDictionary -- a dictionary mapping line names to plot barcodes
    date -- date for which to fetch data; mmddyyyy format, including leading 0s
"""
def parseDryWeightTablesFromSpreadsheets(dwSpreadsheets, genotypeBarcodeDictionary, date):
    dwTable = []
    abandonSheet = False
    for spreadsheet in dwSpreadsheets:
        
        dwTableForSheet = []
        
        ### Parse metadata from "User" sheet and add to every reading from the sheet
        worksheet = spreadsheet.worksheet("User")
        vals = list(map(lambda x: x[0:2], worksheet.get_all_values())) # should be a 5x2 table
        name = vals[0][1]
        computer = vals[1][1]
        scanner = vals[2][1]
        balance = vals[3][1]
        clipweight = vals[4][1]
        readingMetadata = [name,computer,scanner,balance,clipweight]
        
        ### Parse values from scoring spreadsheet
        worksheet = spreadsheet.worksheet("Scoring")
        worksheetVals = worksheet.get_all_values()[1:]
        dataEntryColumn = list(map(lambda x: x[1], worksheetVals[1:])) # selects B2:BXXXX
        ### figure out row of the last data point
        indexOfLastDatapoint = 0
        for i in range(0,len(dataEntryColumn)):
            if dataEntryColumn[i] == "":
                indexOfLastDatapoint = i-1
                break
        ### call just the rows needed
        #for i in range(0,((indexOfLastDatapoint+1)/2)+1):
        for i in range(0,((indexOfLastDatapoint+1)//2)+1): # modified for python3 grammar
            indexOfSummaryLine = i*2
            vals = worksheetVals[indexOfSummaryLine]
            
            """Quality-check each component of the datapoint entry. If a mistake
            is found, abandon attempt to add the current spreadsheet to the day
            summary."""
            
            ### Check for missing data
            for val in vals:
                if val == "":
                    logger.debug("Missing value found in "+spreadsheet.title+", "+worksheet.title+", line "+str(indexOfSummaryLine+2))
                    abandonSheet = True
                    break
                    
            ### Check for valid barcode
            if not vals[3].split("_")[0] in genotypeBarcodeDictionary.keys():
                logger.debug("Unknown or improperly formatted barcode found in "+spreadsheet.title+", "+worksheet.title+", line "+str(indexOfSummaryLine+2) + ": " + vals[3].split("_")[0])
                abandonSheet = True
                break
            
            ### Check for valid dates
            sheetDate = str(int(date[0:2])) + "/" + str(int(date[2:4])) + "/" + date[4:]
            timestamp = vals[2].split(" ")
            if not len(timestamp) == 2:
                logger.debug("Timestamp for datapoint is improperly formatted, in "+spreadsheet.title+", "+worksheet.title+", line "+str(indexOfSummaryLine+2))
                abandonSheet = True
                break
            if not timestamp[0] == sheetDate:
                logger.debug("Date for datapoint does not match sheet date, in "+spreadsheet.title+", "+worksheet.title+", line "+str(indexOfSummaryLine+2))
                abandonSheet = True
                break
            
            ### Check for valid times
            if not validateTimeString(timestamp[1]):
                logger.debug("Time for datapoint is improperly formatted, in "+spreadsheet.title+", "+worksheet.title+", line "+str(indexOfSummaryLine+2))
                abandonSheet = True
                break
            
            ### This step will not be reached unless datapoint values check out
            entry = [vals[3], vals[2], vals[5]]
            dwTableForSheet.append(entry)
            
        ### Add metadata to table entries
        dwTableForSheet = list(map(lambda x: x + readingMetadata, dwTableForSheet))
        dwTable += dwTableForSheet
            
        if abandonSheet:
            logger.debug("Abandoning attempt to consolidate dry weight data for " + date)
            logger.debug("Please correct mistaken or missing value(s) and try again")
            return None
    
    dwTable = list(map(lambda x: x[0].split("_") + x[1:], dwTable))
    
    return dwTable

"""Return master ww/lc table integrating novel data from new table.

    This function takes an existing master table and a table of new data, and
    integrates them into an "updated" master table. When integrating, it
    performs a basic redundancy check for duplicated samples; if a measurement
    is found in the new table with the same genotype, leaf number and timepoint
    as a measurement already in the master table, it will not be added. There is
    a line in this function to throw a warning, but I've commented it out, for
    the reason that one may want to update the master table with the same sheet
    as part of a bulk re-build from local files.
    
    EDIT: this redundancy check is skipped for any entries with plot = "Chk", as
    there will be repeat entries of the same check dish barcodes from different
    days.
        
    Arguments:
    masterWWTable -- the master wet weight/leaf curling table, in the state present on file locally
    wwTable -- the table of new data to update with
"""
def updateMasterWWTableWithNewTable(masterWWTable, wwTable):

    newMasterWWTable = masterWWTable
    if not masterWWTable == []:
        
        archivedSamples = list(map(lambda x: (x[0],x[1]),masterWWTable))
        duplicateSampleFound = False
        for newEntry in wwTable:
            if (newEntry[0],newEntry[1]) in archivedSamples:
                duplicateSampleFound = True
        
        archivedDataPoints = list(map(lambda x: (x[0],x[1],x[10]),masterWWTable))
        for newEntry in wwTable:
            if newEntry[0] == "Chk":
                newMasterWWTable.append(newEntry)
                continue
            if not (newEntry[0],newEntry[1],newEntry[10]) in archivedDataPoints:
                newMasterWWTable.append(newEntry)
            else:
                #logger.debug("Warning: duplicate datapoint detected!")
                #logger.debug(newEntry[0] + ", leaf# " + str(newEntry[1]) + " already has a datapoint entered for timepoint " + str(newEntry[10]))
                pass

    else:
        newMasterWWTable = wwTable

    newMasterWWTable = sorted(newMasterWWTable, key=lambda x: (x[0], x[1]))
    return newMasterWWTable

"""Return a master dry weight table updated with data from the new table.

    Checks for duplicate sample while adding. If a sample in the master is
    repeated in the new table, it will be ignored. If the new table is empty the
    old master table will be returned unchanged.
    
    Arguments:
    masterDWTable -- the master dry weight table, in the state present on file locally
    dwTable -- the table of new data to update with
"""
def updateMasterDWTableWithNewTable(masterDWTable, dwTable):
    
    newMasterDWTable = masterDWTable
    if not masterDWTable == []:
        archivedSamples = list(map(lambda x: (x[0],x[1]),masterDWTable))
        for newEntry in dwTable:
            if not (newEntry[0],newEntry[1]) in archivedSamples:
                newMasterDWTable.append(newEntry)
    else:
        newMasterDWTable = dwTable
        
    newMasterDWTable = sorted(newMasterDWTable, key=lambda x: (x[0], x[1]))
    return newMasterDWTable

"""Return a Google Drive API client with Oauth2 credentials.

    NOTE: to be accessible via the client, a sheet on the drive must be shared
    with maizecuticledatamanagement@appspot.gserviceaccount.com
    
    Arguments:
    jsonKeyFilename -- filename of the .json key file with credentials for accessing the drive account
"""
def getDriveClient(jsonKeyFilename):
    scope = ['https://spreadsheets.google.com/feeds']
    credentials = ServiceAccountCredentials.from_json_keyfile_name(jsonKeyFilename, scope)
    client = gspread.authorize(credentials)
    return client

"""Recompile all the data from AZ16 from the Google drive."""
def recompileAZ16FromDrive(oAuthKeyFilename, environment, archivePath, wwMasterFilename, dwMasterFilename, includeDailySummaryTables):
    datesToInclude = [#"06062016",
                      #"06072016",
                      #"06092016",
                      "06102016",
                      "06112016",
                      "06122016",
                      "06132016",
                      "06142016",
                      "06152016",
                      "06172016",
                      "06182016",
                      "06192016",
                      "06202016",
                      "06212016",
                      "06222016",
                      "06232016",
                      "06242016",
                      "06252016",
                      "06262016",
                      "06272016",
                      "06282016",
                      "06292016",
                      "07022016",
                      "07032016",
                      "07042016",
                      "07052016",
                      "07062016",
                      "07072016",
                      "07082016",
                      "07132016"
                      ]
    for date in datesToInclude:
        updateFromCloudForDate(oAuthKeyFilename, environment, archivePath, wwMasterFilename, dwMasterFilename, date, includeDailySummaryTables)

"""Recompile all the data from SD16 from the Google drive."""
def recompileSD16FromDrive(oAuthKeyFilename, environment, archivePath, wwMasterFilename, dwMasterFilename, includeDailySummaryTables):
    datesToInclude = ["07252016",
                      "07292016",
                      "07312016",
                      "08032016",
                      "08042016",
                      "08052016",
                      "08082016",
                      "08092016",
                      "08112016",
                      "08122016",
                      "08162016",
                      "08182016",
                      "08252016",
                      "08262016",
                      "08302016"
                      ]                  
    
    for date in datesToInclude:
        updateFromCloudForDate(oAuthKeyFilename, environment, archivePath, wwMasterFilename, dwMasterFilename, date, includeDailySummaryTables)
    
"""Recompile all the data from AZ16 from local summary files."""
def recompileAZ16FromLocalSummaryFiles(phenotypeArchivePath, wwMasterFilename, dwMasterFilename):
    
    localFiles = [#phenotypeArchivePath + os.sep + "MLC_AZ16_06062016_WetWeights_LeafCurl_Summary.csv",
                  #phenotypeArchivePath + os.sep + "MLC_AZ16_06072016_WetWeights_LeafCurl_Summary.csv",
                  #phenotypeArchivePath + os.sep + "MLC_AZ16_06092016_WetWeights_LeafCurl_Summary.csv",
                  phenotypeArchivePath + os.sep + "MLC_AZ16_06112016_WetWeights_LeafCurl_Summary.csv",
                  phenotypeArchivePath + os.sep + "MLC_AZ16_06122016_WetWeights_LeafCurl_Summary.csv",
                  phenotypeArchivePath + os.sep + "MLC_AZ16_06132016_WetWeights_LeafCurl_Summary.csv",
                  phenotypeArchivePath + os.sep + "MLC_AZ16_06142016_WetWeights_LeafCurl_Summary.csv",
                  phenotypeArchivePath + os.sep + "MLC_AZ16_06152016_WetWeights_LeafCurl_Summary.csv",
                  phenotypeArchivePath + os.sep + "MLC_AZ16_06172016_WetWeights_LeafCurl_Summary.csv",
                  phenotypeArchivePath + os.sep + "MLC_AZ16_06182016_WetWeights_LeafCurl_Summary.csv",
                  phenotypeArchivePath + os.sep + "MLC_AZ16_06192016_WetWeights_LeafCurl_Summary.csv",
                  phenotypeArchivePath + os.sep + "MLC_AZ16_06202016_WetWeights_LeafCurl_Summary.csv",
                  phenotypeArchivePath + os.sep + "MLC_AZ16_06212016_WetWeights_LeafCurl_Summary.csv",
                  phenotypeArchivePath + os.sep + "MLC_AZ16_06222016_WetWeights_LeafCurl_Summary.csv",
                  phenotypeArchivePath + os.sep + "MLC_AZ16_06232016_WetWeights_LeafCurl_Summary.csv",
                  phenotypeArchivePath + os.sep + "MLC_AZ16_06242016_WetWeights_LeafCurl_Summary.csv",
                  phenotypeArchivePath + os.sep + "MLC_AZ16_06252016_WetWeights_LeafCurl_Summary.csv",
                  phenotypeArchivePath + os.sep + "MLC_AZ16_06262016_WetWeights_LeafCurl_Summary.csv",
                  phenotypeArchivePath + os.sep + "MLC_AZ16_06272016_WetWeights_LeafCurl_Summary.csv",
                  phenotypeArchivePath + os.sep + "MLC_AZ16_06282016_WetWeights_LeafCurl_Summary.csv",
                  phenotypeArchivePath + os.sep + "MLC_AZ16_06292016_WetWeights_LeafCurl_Summary.csv",
                  phenotypeArchivePath + os.sep + "MLC_AZ16_07032016_WetWeights_LeafCurl_Summary.csv",
                  phenotypeArchivePath + os.sep + "MLC_AZ16_07042016_WetWeights_LeafCurl_Summary.csv",
                  phenotypeArchivePath + os.sep + "MLC_AZ16_07052016_WetWeights_LeafCurl_Summary.csv",
                  phenotypeArchivePath + os.sep + "MLC_AZ16_07072016_WetWeights_LeafCurl_Summary.csv",
                  phenotypeArchivePath + os.sep + "MLC_AZ16_06102016_DryWeights.csv",
                  phenotypeArchivePath + os.sep + "MLC_AZ16_06172016_DryWeights.csv",
                  phenotypeArchivePath + os.sep + "MLC_AZ16_06192016_DryWeights.csv",
                  phenotypeArchivePath + os.sep + "MLC_AZ16_06202016_DryWeights.csv",
                  phenotypeArchivePath + os.sep + "MLC_AZ16_06272016_DryWeights.csv",
                  phenotypeArchivePath + os.sep + "MLC_AZ16_07022016_DryWeights.csv",
                  phenotypeArchivePath + os.sep + "MLC_AZ16_07042016_DryWeights.csv",
                  phenotypeArchivePath + os.sep + "MLC_AZ16_07052016_DryWeights.csv",
                  phenotypeArchivePath + os.sep + "MLC_AZ16_07062016_DryWeights.csv",
                  phenotypeArchivePath + os.sep + "MLC_AZ16_07072016_DryWeights.csv",
                  phenotypeArchivePath + os.sep + "MLC_AZ16_07082016_DryWeights.csv",
                  phenotypeArchivePath + os.sep + "MLC_AZ16_07132016_DryWeights.csv"
                  ]
    
    updateFromLocalSummaryTableFiles(localFiles,wwMasterFilename,dwMasterFilename)

"""Recompile all the data from SD16 from local summary files."""
def recompileSD16FromLocalSummaryFiles(phenotypeArchivePath, wwMasterFilename, dwMasterFilename):
    
    localFiles = [phenotypeArchivePath + os.sep + "MLC_SD16_07252016_WetWeights_LeafCurl_Summary.csv",
                  phenotypeArchivePath + os.sep + "MLC_SD16_07292016_WetWeights_LeafCurl_Summary.csv",
                  phenotypeArchivePath + os.sep + "MLC_SD16_07312016_WetWeights_LeafCurl_Summary.csv",
                  phenotypeArchivePath + os.sep + "MLC_SD16_08032016_WetWeights_LeafCurl_Summary.csv",
                  phenotypeArchivePath + os.sep + "MLC_SD16_08042016_WetWeights_LeafCurl_Summary.csv",
                  phenotypeArchivePath + os.sep + "MLC_SD16_08052016_WetWeights_LeafCurl_Summary.csv",
                  phenotypeArchivePath + os.sep + "MLC_SD16_08082016_WetWeights_LeafCurl_Summary.csv",
                  phenotypeArchivePath + os.sep + "MLC_SD16_08092016_DryWeights.csv",
                  phenotypeArchivePath + os.sep + "MLC_SD16_08092016_WetWeights_LeafCurl_Summary.csv",
                  phenotypeArchivePath + os.sep + "MLC_SD16_08112016_WetWeights_LeafCurl_Summary.csv",
                  phenotypeArchivePath + os.sep + "MLC_SD16_08122016_DryWeights.csv",
                  phenotypeArchivePath + os.sep + "MLC_SD16_08122016_WetWeights_LeafCurl_Summary.csv",
                  phenotypeArchivePath + os.sep + "MLC_SD16_08162016_WetWeights_LeafCurl_Summary.csv",
                  phenotypeArchivePath + os.sep + "MLC_SD16_08182016_WetWeights_LeafCurl_Summary.csv",
                  phenotypeArchivePath + os.sep + "MLC_SD16_08182016_DryWeights.csv",
                  phenotypeArchivePath + os.sep + "MLC_SD16_08252016_WetWeights_LeafCurl_Summary.csv",
                  phenotypeArchivePath + os.sep + "MLC_SD16_08262016_DryWeights.csv",
                  phenotypeArchivePath + os.sep + "MLC_SD16_08302016_DryWeights.csv",
                ]
    
    updateFromLocalSummaryTableFiles(localFiles,wwMasterFilename,dwMasterFilename)

"""Recompile all the data from AZ17 from local summary files."""
def recompileAZ17FromLocalSummaryFiles(phenotypeArchivePath, wwMasterFilename, dwMasterFilename):
    
    localFiles = [phenotypeArchivePath + os.sep + "MLC_AZ17_06012017_WetWeights_LeafCurl_Summary.csv",
                  phenotypeArchivePath + os.sep + "MLC_AZ17_06022017_WetWeights_LeafCurl_Summary.csv"
                ]
    
    updateFromLocalSummaryTableFiles(localFiles,wwMasterFilename,dwMasterFilename)
    
"""Executable"""
if __name__ == "__main__":
    
    """
    ============================================================================
    ==== CONFIGURATION
    ============================================================================
    """
    
    """The environment to use. This variable is used to identify spreadsheets by
    title on the drive, and to specify the output path."""
    #environment = "AZ16"
    #environment = "SD16"
    #environment = "PT17" # pilot
    #environment = "AZ17"
    environment = "SD17"    
    """
    ============================================================================
    ==== SETUP
    ============================================================================
    """

    logger.debug('PhenotypeArchiver launched ' + str(time.asctime()))

    """Filepath to the keyfile authorizing access to the developer account, with
    which all spreadsheets should be shared. This is the same keyfile regardless
    of which account the sheets are on."""
    oAuthKeyFilename = "MaizeCuticleDataManagement-55f89af56893.json"
        
    """Path to the parent directory for the archive for the environment"""
    envArchivePath = Common.archivePath + os.sep + "raw" + os.sep + environment
    
    """Output paths for the environment-specific master tables"""
    wwMasterFilename = envArchivePath + os.sep + "MLC_" + environment + "_Master_WetWeights_LeafCurl.csv"
    dwMasterFilename = envArchivePath + os.sep + "MLC_" + environment + "_Master_DryWeights.csv"
    
    """
    ============================================================================
    ==== ARCHIVE SPECIFIC DATES FROM GOOGLE DRIVE
    ============================================================================
    """
     
    """Set this variable to true if the script should write summary files (in
    addition to updating the local master tables) for each of the given dates.
    """
    includeDailySummaryTables = True
     
    """Specify the list of dates for which to fetch data"""
    datesToInclude = ["06012017","06022017","06032017","06052017","06062017","06072017","06082017","06092017","06102017","06112017","06122017","06142017","06152017","06162017","06172017","06182017","06192017","06202017","06212017","06222017","06232017","06242017"]
    datesToInclude = ["080502017"]
    
    """Archive all the data collected on each of the specified dates. Uncomment
    to use, otherwise leave commented."""
    for date in datesToInclude:
        updateFromCloudForDate(oAuthKeyFilename, environment, envArchivePath, wwMasterFilename, dwMasterFilename, date, includeDailySummaryTables)
    
    """
    ============================================================================
    ==== REARCHIVE ALL DATA FROM GOOGLE DRIVE BY ENVIRONMENT 
    ============================================================================
    """
     
    """These functions will recompile the database and daily summary files for
    for an environment from scratch, using the sheets on the Google drive. Make
    sure the local archive path specified above is empty before running."""
    #recompileAZ16FromDrive(oAuthKeyFilename, environment, envArchivePath, wwMasterFilename, dwMasterFilename, True)
    #recompileSD16FromDrive(oAuthKeyFilename, environment, envArchivePath, wwMasterFilename, dwMasterFilename, True)
    
    """
    ============================================================================
    ==== REBUILD MASTER TABLES FROM LOCAL DAILY SUMMARY TABLES BY ENVIRONMENT
    ============================================================================
    """
     
    """These functions will recompile the master tables for an environment from
    the daily summary files. This is much faster than scraping from the original
    scoring sheets on the Google Drives. Daily summary files can be downloaded
    from the "PhenotypeArchive" folder on the gorefieldbook drive. Place them in
    the local archive path specified above (making sure previous iterations are
    deleted first), then run the following to recompile."""
    #recompileAZ16FromLocalSummaryFiles(envArchivePath, wwMasterFilename, dwMasterFilename)
    #recompileSD16FromLocalSummaryFiles(envArchivePath, wwMasterFilename, dwMasterFilename)
    #recompileAZ17FromLocalSummaryFiles(envArchivePath, wwMasterFilename, dwMasterFilename)
     
    logger.debug("PhenotypeArchiver has finished running!")
    
