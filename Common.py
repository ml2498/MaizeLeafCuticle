"""Variables and functions shared across all analysis pipelines.

Author:         James Chamness
Last Modified:  06/05/2017
"""

"""Dependencies"""
import os

"""
================================================================================
==== Variables
================================================================================
"""

"""Environment designators"""
environments = ["AZ16", "SD16", "PT17", "AZ17","SD17"]

"""Top-level for entire project"""
#projectTopLevel = "/home/james/GoreLab/MaizeLeafCuticle" # my laptop
#projectTopLevel = "/home/jamesc/MLC_2017" # my server container
#projectTopLevel = "/share/MLC_2017" # server public
projectTopLevel = "/Users/Meng/'Google Drive'/TWAS_2018/Pheno" # my laptop

"""Top-level for archive of phenotypes taken from Google Drive"""
archivePath = projectTopLevel + os.sep + "phenotypes" + os.sep + "archive"

"""Top-level for tables with datapoints aggregated by sample"""
aggregatedPath = projectTopLevel + os.sep + "phenotypes" + os.sep + "aggregatedSamplePoints"

"""Top-level for tables with raw (unscaled) evaporation rates"""
rawRatesPath = projectTopLevel + os.sep + "phenotypes" + os.sep + "rawEvaporationRates"

"""Top-level for surface area data"""
saPath = projectTopLevel + os.sep + "phenotypes" + os.sep + "surfaceAreas"

"""Top-level for scaled evaporation rates"""
scaledRatesPath = projectTopLevel + os.sep + "phenotypes" + os.sep + "scaledEvaporationRates"

"""Microclimate data archive top-level directory"""
microclimateArchivePath = projectTopLevel + os.sep + "microclimateMonitors" + os.sep + "archive"

"""Genotype top-level directory"""
genotypeTopLevel = projectTopLevel + os.sep + "genotypes"

"""Path to the design directory"""
designPath = projectTopLevel + os.sep + "design"

"""
================================================================================
==== Functions
================================================================================
"""

"""Return a dictionary mapping plot barcodes to genotypes for given environment.

    In these dictionaries, the check plots will have the plot number suffixed to
    the genotype name (e.g., Mo17_0009) to that different plots can be
    discriminated. This also allows for genotypes to be used both as checks and
    as experimental lines, in different plots; for the experimental plots, no
    suffix is appended.
    
    The genotype names are converted to project standard.

    Sources:
    AZ16: Maricopa_Field_Design_with_Checks_20160422_nk_FINAL.xlsx, worksheet
    titled "Maricopa_Field_Design_with_Chec"
    SD16: UCSD_Field_Design_with_Checks_BOTH_YEARS_05042016.xlsx,worksheet
    titled "Year 1 Labels".
    PT17: barcode_dict.csv
    
    Arguments:
    environment -- e.g. "AZ16"
"""
def getPlotBarcodeDict(environment):
    dict = {}
    if environment == "AZ16":
        design_filename = designPath + os.sep + "Maricopa_Field_Design_with_Checks_20160422_nk_FINAL.csv"
        with open(design_filename) as f:
            lineNum = 0
            for line in f:
                lineNum += 1
                if lineNum == 1:
                    continue
                vals = line.strip().split("\t")
                plot = "0" * (4 - len(vals[0])) + vals[0]
                plotFull = "16AZ" + plot
                if "Check" in vals[3]:
                    geno = geno + "_" + plot
                elif vals[6] == "ABP15":
                    geno = "ABP15"
                else:
                    geno = translateNameToStandard(vals[6], 1)
                dict[plotFull] = geno
        return dict
    elif environment == "SD16":
        design_filename = designPath + os.sep + "SD16BarcodeMappingSource.csv"
        with open(design_filename) as f:
            lineNum = 0
            for line in f:
                lineNum += 1
                if lineNum == 1:
                    continue
                vals = line.strip().split("\t")
                plot = "0" * (4 - len(vals[0])) + vals[0]
                plotFull = "16AZ" + plot # NOTE that somebody screwed up the first year and we used the same prefix for both locations
                if vals[1] == "Check1":
                    geno = "B73_" + plot
                elif vals[1] == "Check2":
                    geno = "MO17_" + plot
                else:
                    geno = translateNameToStandard(vals[2], 2)
                dict[plotFull] = geno
        return dict
    elif environment == "PT17":
        design_filename = designPath + os.sep + "2017_pilots_barcode_dict.csv"
        with open(design_filename) as f:
            for line in f:
                vals = line.strip().split(",")
                dict[vals[0]] = vals[1]
        return dict
    elif environment == "AZ17":
        design_filename = designPath + os.sep + "Maricopa_Fieldbook_Year2_05312017.csv"
        with open(design_filename) as f:
            lineNum = 0
            for line in f:
                lineNum += 1
                if lineNum == 1:
                    continue
                vals = line.strip().split("\t")
                barcode = vals[5]
                plot = "0" * (3 - len(barcode)) + barcode
                
                if vals[8] == "Check1":
                    geno = "N28HT_" + plot
                elif vals[8] == "Check2":
                    geno = "MO17_" + plot
                else:
                    geno = translateNameToStandard(vals[1], 2)
                dict[plot] = geno
        return dict
    elif environment == "SD17":
        design_filename = designPath + os.sep + "SD17BarcodeMappingSource.csv"
        with open(design_filename) as f:
            lineNum = 0
            for line in f:
                lineNum += 1
                if lineNum == 1:
                    continue
                vals = line.strip().split("\t")
                barcode = vals[5]
                plot = "0" * (3 - len(barcode)) + barcode
                
                if vals[8] == "Check1":
                    geno = "N28HT_" + plot
                elif vals[8] == "Check2":
                    geno = "MO17_" + plot
                else:
                    geno = translateNameToStandard(vals[1], 2)
                dict[plot] = geno
        return dict
                
    return

"""Return a table object read from the given input file.

    Arguments:
    filename -- the filename of the file containing the table to read
    sep -- the delimiter character to expect
    header -- whether or not to skip the first line of the file
"""
def readTableFromFile(filename, sep="\t", header=True):
    try:    
        tableFile = open(filename)
    except IOError:
        print("FILE NOT FOUND: " + filename)
        print("Returning None")
        return None
    table = []
    lineNum = 0
    for line in tableFile:
        lineNum += 1
        if lineNum == 1 and header:
            continue
        vals = line.strip().split(sep)
        table.append(vals)
    return table

"""Write the given table to file.

    No trailing delimiter chars or end of file empty lines.
    Any non-string table entries will be coerced to strings.
    
    Arguments:
    table -- the table to write
    filename -- the filename of the file to write to
    header -- the header for the table to write as the first line
    sep -- delimiter character
    pad -- fill in empty "cells" with NA (based on length of longest entry)
    mode -- file write mode. Default is to overwrite/make new file 
"""
def writeTableToFile(table, filename, header=None, sep="\t", pad=True, mode='w'):
    
    # determine length of longest entry in table:
    maxEntryLength = max(list(map(lambda x: len(x), table)))
    
    out = open(filename, mode)
    if header is not None:
        out.write(header + "\n")
        
    for i in range(0, len(table)-1):
        if len(table[i]) == maxEntryLength:
            for j in range(0,maxEntryLength-1):
                out.write(str(table[i][j]) + sep)
            out.write(str(table[i][maxEntryLength-1]) + "\n")
        else:
            for j in range(0,len(table[i])):
                out.write(str(table[i][j]) + sep)
            for j in range(0,maxEntryLength - len(table[i])-1):
                out.write("NA" + sep)
            out.write("NA" + "\n")
            
    if len(table[-1]) == maxEntryLength:
        for j in range(0,maxEntryLength-1):
            out.write(str(table[-1][j]) + sep)
        out.write(str(table[-1][maxEntryLength-1]) + "\n")
    else:
        for j in range(0,len(table[-1])):
            out.write(str(table[-1][j]) + sep)
        for j in range(0,maxEntryLength - len(table[-1])-1):
            out.write("NA" + sep)
        out.write("NA" + "\n")

    out.close()

"""Deprecated"""

# """Return a dictionary mapping plot barcodes to genotypes for given environment.
# 
#     NOTE: For the check lines, which are repeated, the line name is suffixed
#     with the barcode so that different plots can be discriminated.
#     
#     Arguments:
#     environment -- e.g., "AZ16", "SD16"
# """
# def getGenotypeBarcodeDictionary(environment):
#     genotypeBarcodeKeyFile = designPath + os.sep + "MLC_" + environment + "_Genotypes_to_PlotBarcodes_Key.csv"
#     keyFile = open(genotypeBarcodeKeyFile)
#     coding = {}
#     for line in keyFile:
#         vals = line.strip().split(",")
#         if vals[1] != 0: ### experimental genotypes all have sources listed
#             coding[vals[2]] = vals[0]
#         else: ### but the checks do not...
#             ### so we have to add barcode as suffix to get a unique plant/leaf id
#             coding[vals[2]] = vals[0] + "_" + vals[2]
#     return coding
# 
# """Return a dictionary mapping line names to plot barcodes.
#  
#     NOTE: For the check lines, which are repeated, the line name should already
#     be suffixed with the barcode so that different plots can be discriminated.
#      
#     Arguments:
#     environment -- e.g., "AZ16", "SD16"
# """
# def getBarcodeGenotypeDictionary(environment):
#     genotypeBarcodeKeyFile = designPath + os.sep + "MLC_" + environment + "_Genotypes_to_PlotBarcodes_Key.csv"
#     keyFile = open(genotypeBarcodeKeyFile)
#     coding = {}
#     currLineNum = 0
#     for line in keyFile:
#         currLineNum += 1
#         if currLineNum == 1:
#             continue
#         vals = line.strip().split(",")
# #         if vals[1] != 0: ### experimental genotypes all have sources listed
# #             coding[vals[2]] = vals[0]
# #         else: ### but the checks do not...
# #             ### so we have to add barcode as suffix to get a unique plant/leaf id
# #             coding[vals[2]] = vals[0] + "_" + vals[2]
#         if vals[0] == "B73" or vals[0] == "Mo17":
#             coding[vals[0] + "_" + vals[2][-4:]] = vals[2]
#         else:
#             coding[vals[0]] = vals[2]
#      
# #     for key in coding.keys():
# #         print(key + "\t:\t" + coding[key])        
#     return coding

# """Export a table of sample weights and times from a sampleValues dictionary.
# 
#     Write a table with just the sample IDs, raw weights, and timestamps. The R
#     script "MLC_DataCleanup.R" is used to efficiently plot these data and spot
#     outliers.
# 
#     Arguments:
#     sampleValues -- see SampleDataAggregator module for documentation
#     exportFilename -- the file to write the table to
# """
# def exportSampleDatapoints(sampleValues, exportFilename):
#     exportTable = []
#     for sampleID in sampleValues.keys():
#         entry = []
#         sample = sampleValues[sampleID]
#         entry = entry + [sampleID[0]] + [sampleID[1]]
#         for val in sampleValues[sampleID][1:-1]:
#             entry = entry + [val[0]] + [val[1]]
#             #print(str(val))
#         exportTable.append(entry)
#     exportTable = sorted(exportTable,key=lambda x: (x[0],x[1]))
#     writeTableToFile(exportTable, exportFilename)

"""Provided the source, return the "MLC standard name" for the given taxon.

    This does NOT handle the 2 check varieties, which are not in the name
    conversion table, nor "ABP15", a line added to the SD16 environment but that
    is not one of the 466 "experimental" taxa.
    
    Also note: one accession was repeated under two different genotype names
    across the two environments - "Mo43" and "Mo44". These are left separate by
    this method

    Arguments:
    name -- the taxon name as a string literal
    source -- an integer specifying the originating name scheme, as follows:
                1: AZ16 fieldbook
                2: SD16 fieldbook
                3: Hirsch et al., 2014, Supplementary Table 1
                4: Hirsch et al., 2014, Raw SNP Table
                5: Hirsch et al., 2014, Imputed SNP Table
                6: Hansey et al., 2010, Supplementary Table 1. WARNING: this
                    table has duplicate entries for some taxa in our panel
"""
def translateNameToStandard(name, source):
    
    nameConversionTable = readTableFromFile(designPath + os.sep + "taxa" + os.sep + "taxa_&_accession_name_mappings.csv",sep="\t",header=True)
    
    """One accession was repeated under two different genotype names. I'm
    leaving them separate here, but requies standalone code because I didn't
    include an entry for Mo44 in the taxa name mapping spreadsheet."""
    if name == "Mo43" and (source == 1 or source == 2):
        return "MO43"
    """^^^"""
    
    if source == 1:
        nameList = list(map(lambda x: x[source],nameConversionTable))
    elif source == 2:
        nameList = list(map(lambda x: x[source],nameConversionTable))
    elif source == 3:
        nameList = list(map(lambda x: x[source],nameConversionTable))
    elif source == 4:
        nameList = list(map(lambda x: x[source+2],nameConversionTable))
    elif source == 5:
        nameList = list(map(lambda x: x[source+2],nameConversionTable))
    elif source == 6:
        nameList = list(map(lambda x: x[source+2],nameConversionTable))
    else:
        print("No source with this index!")
        return None
    standardNameList = list(map(lambda x: x[0],nameConversionTable))
    #for x in standardNameList:
    #    print(x)
    standardName = standardNameList[nameList.index(name)]
    
    return standardName

"""Translate a given "MLC standard name" to the"fieldbook-style" name.

    Warning: "MO44" backconverts only to "Mo44", and not "Mo43", which are
    distinct plots in SD16. 

    Arguments:
    name -- a taxon name in the MLC standard
"""
def translateNameToFieldbookScheme(name):
    
    nameConversionTable = readTableFromFile(designPath + os.sep + "taxa" + os.sep + "taxa_&_accession_name_mappings.csv",sep=",",header=True)
    
    standardNameList = list(map(lambda x: x[0], nameConversionTable))
    fieldbookNameList = list(map(lambda x: x[1],nameConversionTable))
    
    fieldbookName = fieldbookNameList[standardNameList.index(name)]

    return fieldbookName

"""Return a lexicographically ordered list of the MLC taxa.

    There are 465.
"""
def getTaxaList():
    nameConversionTable = readTableFromFile(designPath + os.sep + "taxa" + os.sep + "taxa_&_accession_name_mappings.csv",delimChar=",",header=True)
    taxaList = list(map(lambda x: x[0], nameConversionTable))
    taxaList = [x for x in taxaList if not x == "NA"]
    return(sorted(taxaList))

"""Executable"""
if __name__ == "__main__":
    pass
    """here for testing purposes"""
    #getGenotypeBarcodeDictionary("AZ16")
    #getGenotypeBarcodeDictionary("SD16")
    getPlotBarcodeDict("AZ17")
