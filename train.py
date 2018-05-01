#!/usr/bin/python

__author__ = 'Chuang Ma'
__date__ = "$2015-9-1$"

import sys
import os
import smtplib
from email.mime.text import MIMEText
from email.header import Header

args = sys.argv

def psource(module):
    file = os.path.basename(module)
    dir = os.path.dirname(module)

    toks = file.split('.')
    modname = toks[0]

    # Check if dirrectory is really a directory
    if (os.path.exists(dir)):

        # Check if the file directory already exists in the sys.path array
        paths = sys.path
        pathfound = 0
        for path in paths:
            if (dir == path):
                pathfound = 1

        # If the dirrectory is not part of sys.path add it
        if not pathfound:
            sys.path.append(dir)

    # exec works like MEL's eval but you need to add in globals()
    # at the end to make sure the file is imported into the global
    # namespace else it will only be in the scope of this function
    exec ('import ' + modname) in globals()

    # reload the file to make sure its up to date
    exec ('reload( ' + modname + ' )') in globals()

    # This returns the namespace of the file imported
    return modname


# end


def train(fullname,email):
    ##################################################################################################
    ##############################parameters section##################################################
    ##################################################################################################
    ##the directory of RNAfold                                                                     ###
    RNAFoldDic = "/usr/bin/"  ###
    ###
    ##file directory of source.py                                                                  ###
    sourceDic = os.getcwd()+"/"  ###
    ###
    ##result directory, include trainDataFileName, predDataFileName, predDataAnnotFileName         ###
    resultDic = sourceDic + "results/"  ###
    ###
    ###
    ##a logical parameter indicate whether run cross validation on training dataset
    cross_validation_flag = False
    # if args[1]== 'true' or args[1] == 'True':
    #     cross_validation_flag = True  ###

    ###
    ##miRNAs and pre-miRNAs for training                                                           ###
    trainDataFileName = fullname.split('/')[1]  ###
    ##pre-miRNAs for prediction                                                                    ###
    predDataFileName = "predictionData.txt"  ###
    ##the annotated miRNAs and pre-miRNAs for prediction dataset                                   ###
    ###
    ###
    ##################################optional parameters ############################################
    ##the full path of trained prediction model. If a full path is given,                          ###
    ##miRLocator will locate it directly, otherwise, miRLocator will run the training program.     ###
    predModelFileDir = resultDic + "trained_prediction_model"  ###
    ##If the file name is defined, miRLocator will evaluate the prediction results                 ###
    ##based on annotation infomation in this file                                                  ###
    predDataAnnotFileName = "predictionData_Annotated.txt"  ###
    ##################################################################################################


    #####################################################################################################################################
    #####################################################################################################################################
    # keep this command, temp file used to record temporary results
    tempFileDir = resultDic + "tempResult.txt"
    psource(sourceDic + "source.py")

    trainDataFileDir = sourceDic + email + '/' + trainDataFileName
    predDataFileDir = resultDic + predDataFileName
    predDataAnnotFileDir = resultDic + predDataAnnotFileName

    predResultFileDir = resultDic + "miRLocator_predResults.txt"
    evalResultFileDir = resultDic + "miRLocator_evalResults.txt"

    # create file directories for recording dp_ss files for training and prediction
    dpSSFileDic_train = sourceDic + email + '/' + "dp_ss_train/"
    dpSSFileDic_pred = sourceDic +  email + '/' + "dp_ss_pred/"
    source.createDict(dpSSFileDic_train)
    source.createDict(dpSSFileDic_pred)

    #####################################################################################################################################
    #####################################################################################################################################
    # check file for training and prediction, and fold pre-miRNA sequences
    ##InputFileDir: the full path of input file including pre-miRNA and mature miRNA sequences (and structures)
    ##RNAFoldDic: the full path of RNAfold
    ##sourceDic:the directory of source.py, miRLocator.py
    ##dpSSFileDic:the directory of dp.ps, ss.ps files
    ##dpSSFlag: True: dp_ss files already exists, not need to be generated again; False: no dp_ss files, need to be generated in this run
    ##trainDataFileDir_Ref: refined data for training.

    # trainDataFileDir_Ref = source.checkFileForTraining(trainDataFileDir, RNAFoldDic, sourceDic, dpSSFileDic_train, dpSSFlag = True )


    ##cross_validation_test
    if (cross_validation_flag == True):
        print "start cross validation.\n"
        source.cross_validation_function(source.cv, trainDataFileDir_Ref, source.minCandidateMiRNALen,
                                         source.maxCandidateMiRNALen, source.upOffSet, source.downOffSet,
                                         source.minPredScore,
                                         resultDic, dpSSFileDic_train, RNAFoldDic, tempFileDir)
    # end if


    ##train prediction model on all miRNAs in input file
    ##trainDataFileDir_Ref: refined data for training
    ##dpSSFileDic_train: dp_ss files for training data
    ##RNAFoldDic: the directory of RNAfold
    ##tempFileDir: tempoary file
    ##source.upOffSet: 5 in default
    ##source.downOffSet: 5 in default
    # predModelFileDir: the full path of trained prediction model
    ##resultDic: result directory
    predModel = ""

    trainDataFileDir_Ref = source.checkFileForTraining(trainDataFileDir, RNAFoldDic, sourceDic, dpSSFileDic_train,
                                                           dpSSFlag=True)
    newdir = os.getcwd()+'/'+'result_'+email
    isExits = os.path.exists(newdir)
    if not isExits:
        os.mkdir(newdir)
    newdir = newdir+'/'
    print "start to train prediction models.\n"
    predModelFileDir = newdir + "trained_prediction_model"
    predModel = source.train_prediction_model(trainDataFileDir_Ref, dpSSFileDic_train, RNAFoldDic, tempFileDir,
                                                  source.upOffSet, source.downOffSet, predModelFileDir, resultDic)
    return predModelFileDir
    # end else


