"""generateSummary

Usage:
  generateSummaryInfo.py -d <driverFile> -I <funSeqOutFile> -O <outSummaryFile>
  generateSummaryInfo.py (-h | --help)

Options:
  -h --help     Show this screen.

"""
from docopt import docopt


import sys
import os
import subprocess


def readDriverInfo(driverInputFile):

    """ extract coding and non-coding PCAWG driver genes based in patient-centric driver paper"""
 
    cdsList = [] #coding driver list
    ncList = []  #non-coding driver list


    fileName = open(driverInputFile,'r')
    if fileName:
       for index in range(1,2):
           fileName.next()

       for line in fileName:
           line = line.strip()

           if line.split()[9] == "snv":

              geneInfo = line.split()[2] #gene corresponding to driver element
              snpInfo = line.split()[4]+":"+"chr"+line.split()[5]+":"+line.split()[6]+":"+line.split()[7]+":"+line.split()[8] #driver mutation 
 
              """ store PCAWG coding driver gene Name """
              if line.split()[1] == "CDS" or line.split()[1] == "splice_sites":
                 cdsList.append(geneInfo)

              #""" store PCAWG non-coding driver gene Name """
              elif line.split()[1] not in ['CDS','splice_sites']:
                 ncList.append(geneInfo)



    return list(set(cdsList)),list(set(ncList))



def extractFunSeqScore(funSeqOut,cdsDriverList,ncDriverList):

    """ extract coding and non-coding mutation attributes including drive/non-driver status and coding/promoter/intronic/other-non-coding status"""

    funSeqCodingScore = 0.0
    funSeqNonCodingScore = 0.0

   
    driverInfo = 0
    promoterStatus = 0
    intronStatus = 0
 
    cdsTup = ()
    noncdsTup = ()

    scoreList = []

    fileName = open(funSeqOut,'r')
    if fileName:
       for index in range(1,2):
           fileName.next()

       for line in fileName:
           line = line.strip()

           cds = line.split()[6].split(';')[1] # coding SNV status
           gene = line.split()[6].split(';')[2] # geneId of the coding SNV
           sampleId = line.split()[5] #sampleId

           chromNum = line.split()[0]  #chromosome number
           snvPos = line.split()[2]    #snv position
           refAllele = line.split()[3] #ref allele
           altAllele = line.split()[4] #alternate allele

           if len(line.split()) >= 7:

              """ extract coding mutations info and thier attributes """

              if cds == 'Yes' and gene != ".":

                 annotInfo = 0 #coding(0) or non-coding(1) status 

                 target_gene = line.split()[6].split(';')[11]

                 targetGene = target_gene.split('(')[0].split('[')[0] #gene name containing the mutation

                 """ If the target gene in the PCAWG coding driver gene list, driver status as 1 or 0 otherwise"""

                 if targetGene in cdsDriverList:
                    driverInfo = 1

                 else:
                    driverInfo = 0

                 """ extract funseq impact score for each coding mutation"""

                 if line.split()[6].split(';')[-2:][0] != '.' and line.split()[6].split(';')[-2:][1] == '.':
                    funSeqCodingScore = float(str(line.split()[6].split(';')[-4:][0])) - 1.0

                 elif line.split()[6].split(';')[-2:][0] == '.' and line.split()[6].split(';')[-2:][1] != '.':
                      funSeqCodingScore = float(str(line.split()[6].split(';')[-4:][0])) + 1.0
                       
                 else:
                    funSeqCodingScore = float(str(line.split()[6].split(';')[-4:][0]))


                 snpAtt = chromNum+':'+snvPos+':'+refAllele+':'+altAllele #mutation information

                 cdsTup = (snpAtt,annotInfo,driverInfo,funSeqCodingScore,targetGene,promoterStatus,intronStatus) #tuple containing mutation attributes including functional impact and driver status                  
                 scoreList.append(cdsTup) #store each mutation and it's attribute into a list



              """ extract non-coding mutations info and thier attributes """

              if line.split()[6].split(';')[1] == "No":
             
                 annotInfo = 1 #coding(0) or non-coding(1) status

                 target_gene = line.split()[6]

                 """ extract funseq impact score for each non-coding mutation"""

                 tmpList =  [float(i) for i in line.split()[6].split(';')[-3:][0].split(":")]
                 minScore = min(tmpList)
                 maxScore = max(tmpList)
                 
             
                 if line.split()[6].split(';')[-2:][0] != '.' and line.split()[6].split(';')[-2:][1] == '.':
                    funSeqNonCodingScore = minScore

                 elif line.split()[6].split(';')[-2:][0] == '.' and line.split()[6].split(';')[-2:][1] != '.':
                    funSeqNonCodingScore = minScore + 1.0

                 else:
                    funSeqNonCodingScore = maxScore


                 if target_gene != '.':

                    targetGene = target_gene.split(';')[11].split('(')[0].split('[')[0] #gene name containing the non-coding mutation


                    """ If the target gene in the PCAWG non-coding driver gene list, driver status as 1 or 0 otherwise"""

                    if targetGene in ncDriverList:
                       driverInfo = 1

                    else:
                       driverInfo = 0


                    #""" If the mutation fall within promoter region then promoter status as 1 or 0 otherwise"""

                    if [s for s in target_gene.split(',') if 'Promoter' in s]:

                       for element in target_gene.split(','):
                           if  element[0:4] != 'ENSG' and element != '.' and 'Promoter' in element and 'Intron&Promoter' not in element:

                               promoterStatus = 1

                           else: 
                               promoterStatus = 0


                    #""" If the non-coding mutation fall within intron then Intron status as 1 or 0 otherwise"""
 
                    elif [s for s in target_gene.split(',') if 'Intron' in s]:

                       for element in target_gene.split(','):
                           if  element[0:4] != 'ENSG' and element != '.' and 'Intron' in element:

                               intronStatus = 1

                           else:
                               intronStatus = 0


                 snpAtt = chromNum+':'+snvPos+':'+refAllele+':'+altAllele #mutation information

                 noncdsTup = (snpAtt,annotInfo,driverInfo,funSeqNonCodingScore,targetGene,promoterStatus,intronStatus) #tuple containing mutation attributes including functional impact and driver status                  
                 scoreList.append(noncdsTup) #store each mutation and it's attribute into a list
    
                 driverInfo = 0 
                 promoterStatus = 0
                 intronStatus = 0

 

    return scoreList                       




def main (arguments):

    codingDriver = []
    ncDriver = []

    funSeqScoreList = []


    funSeqOutFile = arguments["<funSeqOutFile>"]     #directory containing original funseq output 

    driverFile = arguments["<driverFile>"] #PCAGW driver information

    outSummaryFile = arguments["<outSummaryFile>"] #mutation summary output file


    ''' extract PCAWG coding and non-coding driver gene list '''
    codingDriver,ncDriver = readDriverInfo(driverFile)
   
    ''' extract funseq score and annotation status for coding and non-coding mutations '''
    funSeqScoreList = extractFunSeqScore(funSeqOutFile,codingDriver,ncDriver)


    fileOutput = open(outSummaryFile,'w')

    for item in funSeqScoreList:
        fileOutput.write(str(item[0])+'\t'+str(item[1])+'\t'+str(item[2])+'\t'+str(item[3])+'\t'+str(item[5])+'\t'+str(item[6])+'\t'+str(item[4]))
        fileOutput.write("\n")



if __name__ == '__main__':

    arguments = docopt(__doc__, version='generateSummaryInfo 1.0')
    main(arguments)

