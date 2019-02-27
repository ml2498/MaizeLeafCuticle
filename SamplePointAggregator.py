"""Aggregates wet weight readings by sample from the processed archive tables.

This module contains a function which returns a dictionary mapping samples to
values. The dictionary keys are sample IDs, represented as a tuple of the
genotype and leaf number, e.g. ('A15',2). The values are given as a list with
the following form (all atoms are strings):
    
    [Clip_Weight, Scorer, Computer, Scanner, Balance, Sensor, Date,
     (Timepoint,Weight,Timestamp),
     (Timepoint,Weight,Timestamp),
     ...]
    
...where there are as many datapoint tuples as are found for that sample. Note
that this module does not perform any QC, and will not alert the user to
inadvertently duplicated datapoints. 

This module may be called as a standalone script, and its default behavior is to
create a table of all "weight datapoints" (weight & timepoint pairs) for the
given environment, and write this table to file in the /qc folder. There, it is
useful for manual cleanup using R or what have you.

This module is also called by other pipeline modules for the dictionary
aggregating all data (including metadata) by sample.



Should this module also integrate dry weight data? What about surface area
coming later? No, I think keep these separate until QC is done on the weight
data and unscaled CE rate has been calculated.

Author:         James Chamness
Last Modified:  05/15/2017
"""

"""Dependencies"""
from datetime import datetime
import logging
import os
import time
from importlib.machinery import SourceFileLoader
Common = SourceFileLoader("Common","../../pipeline/Common.py").load_module()

""""Set up logger routing to console and to file"""
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
logger.addHandler(logging.StreamHandler())
#logger.addHandler(logging.FileHandler('SamplePointAggregator.log', mode='w'))

"""Build the dictionary described in the header."""
def aggregateSampleDataForEnvironment(environment):

    """Paths to the processed archive table(s) for the environment"""
    wwMasterFilename = Common.archivePath + os.sep + "processed" + os.sep + "MLC_" + environment + "_Master_WetWeights_Processed.csv"
    
    master = Common.readTableFromFile(wwMasterFilename)
    sampleValues = {}
    
    ### build the dictionary for all the values except the dry weight
    for entry in master:
        sampleKey = (entry[0],int(entry[1]))
        ### add the datapoint from the entry to the weight dictionary
        if sampleKey not in list(sampleValues.keys()):
            sampleValues[sampleKey] = entry[5:] + [(entry[4],entry[2],entry[3])]
        else:
            sampleValues[sampleKey] = sampleValues[sampleKey] + [(entry[4],entry[2],entry[3])]
    
    return sampleValues

"""Write simplified points table for the given environment."""
def exportSampleDatapoints(environment):
    exportFilename = Common.aggregatedPath + os.sep + "unfiltered" + os.sep + "MLC_" + environment + "_Unfiltered_Sample_Datapoints.csv"
    exportTable = []
    sampleValues = aggregateSampleDataForEnvironment(environment)
    
    for sampleID in list(sampleValues.keys()):
        vals =  sampleValues[sampleID]
        exEntry = [sampleID[0],sampleID[1]] + vals[:7]
        if sampleID[0] == "A385" and sampleID[1] == 7:
            pass
        for i in range(7, len(vals)):
            exEntry.append(vals[i][1])
            exEntry.append(vals[i][2])
        exportTable.append(exEntry)
    
    maxTimepoints = int((max(list(map(lambda x: len(x), exportTable))) - 9) / 2)
    header = "Genotype\tLeaf_Number\tClipWeight\tScorer\tComputer\tScanner\tBalance\tSensor\tDate"
    for i in range(0,maxTimepoints):
        header += "\tWeight" + str(i+1) + "\tTimestamp" + str(i+1)
    Common.writeTableToFile(exportTable,exportFilename,header,pad=True)
    return

"""Executable"""
if __name__ == "__main__":
     
    """
    ============================================================================
    ==== CONFIGURATION
    ============================================================================
    """
     
    """The environment to use. Used to specify paths."""
    environment = "AZ16"
    environment = "SD16"
    environment = "PT17"
    environment = "AZ17"
    
    """
    ============================================================================
    ==== SCRIPT
    ============================================================================
    """
    
    logger.debug('SamplePointAggregator launched ' + str(time.asctime()))
    
    #sampleValues = aggregateSampleDataForEnvironment(environment)
    
    #for key in sorted(list(sampleValues.keys())):
    #    print(key)
    
    exportSampleDatapoints(environment)
    
    logger.debug('SamplePointAggregator has finished running!')
    