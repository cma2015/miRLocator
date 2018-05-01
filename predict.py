#!/usr/bin/python

__author__ = 'Chuang Ma'
__date__ = "$2015-9-1$"

import sys
import os
import smtplib
from email.mime.text import MIMEText
from email.header import Header

def predict(fullname,modelname,evaname,email):
    RNAFoldDic = "/usr/bin/"  ###
    ###
    ##file directory of source.py                                                                  ###
    sourceDic = os.getcwd() + "/"  ###
    ###
    ##result directory, include trainDataFileName, predDataFileName, predDataAnnotFileName         ###
    resultDic = sourceDic + "results/"  ###
    ###
    ###
    ##a logical parameter indicate whether run cross validation on training dataset
    cross_validation_flag = False  ###
    ###
    ##miRNAs and pre-miRNAs for training                                                           ###
    trainDataFileName = "trainingData.txt"  ###
    ##pre-miRNAs for prediction                                                                    ###
    predDataFileName = fullname  ###
    ##the annotated miRNAs and pre-miRNAs for prediction dataset                                   ###
    ###
    ###
    ##################################optional parameters ############################################
    ##the full path of trained prediction model. If a full path is given,                          ###
    ##miRLocator will locate it directly, otherwise, miRLocator will run the training program.
    #     ###
    if modelname == '':
        predModelFileDir = resultDic + "trained_prediction_model"  ###
    else:
        predModelFileDir = modelname
    ##If the file name is defined, miRLocator will evaluate the prediction results                 ###
    ##based on annotation infomation in this file
    #                                              ###

    predDataAnnotFileName = "predictionData_Annotated.txt"  ###

    ##################################################################################################

    #####################################################################################################################################
    #####################################################################################################################################
    # keep this command, temp file used to record temporary results
    tempFileDir = resultDic + "tempResult.txt"
    psource(sourceDic + "source.py")
    eval = False
    trainDataFileDir = resultDic + trainDataFileName
    predDataFileDir = sourceDic  + predDataFileName
    if evaname=='':
        eval = False
    else:
        predDataAnnotFileDir = sourceDic  + predDataAnnotFileName
        eval = True

    predResultFileDir = sourceDic + email + '/' + "miRLocator_predResults.txt"
    evalResultFileDir = sourceDic + email + '/' + "miRLocator_evalResults.txt"

    # create file directories for recording dp_ss files for training and prediction
    dpSSFileDic_train = sourceDic + email + '/' + "dp_ss_train/"
    dpSSFileDic_pred = sourceDic + email + '/' + "dp_ss_pred/"
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
    # if (len(predModelFileDir) == 0):
    #     trainDataFileDir_Ref = source.checkFileForTraining(trainDataFileDir, RNAFoldDic, sourceDic, dpSSFileDic_train,
    #                                                        dpSSFlag=True)
    #     print "start to train prediction models.\n"
    #     predModelFileDir = resultDic + "trained_prediction_model"
    #     predModel = source.train_prediction_model(trainDataFileDir_Ref, dpSSFileDic_train, RNAFoldDic, tempFileDir,
    #                                               source.upOffSet, source.downOffSet, predModelFileDir, resultDic)
    # else:
    #     print "start to load prediction models.\n"
    #     predModel = source.loadPredModel(predModelFileDir)
    # # end else

    # prediction mature miRNAs within pre-miRNA sequences
    ##predDataFileDir_Ref:refined data for prediction
    ##dpSSFileDic_pred: dp_ss files for prediction
    ##RNAFoldDic:the directory of RNAfold
    ##tempFileDir:tempoary file
    ##source.minCandidateMiRNALen: 16 in default
    ##source.maxCandidateMiRNALen: 30 in default
    ##source.upOffSet: 5 in default
    ##source.downOffSet: 5 in default
    ##source.minPredScore: 0.25 in default, but not used in this version
    ##predModel: trained prediction model
    ##predResultFileDir: the full path of prediction results
    predModel = source.loadPredModel(predModelFileDir)
    predDataFileDir_Ref = ""
    if (len(predDataFileName) != 0):
        print "start to predict miRNAs.\n"
        predDataFileDir_Ref = source.checkFileForPrediction(predDataFileDir, RNAFoldDic, sourceDic, dpSSFileDic_pred,
                                                            dpSSFlag=True)
        source.prediction(predDataFileDir_Ref, dpSSFileDic_pred, RNAFoldDic, tempFileDir,
                          source.minCandidateMiRNALen, source.maxCandidateMiRNALen, source.upOffSet, source.downOffSet,
                          source.minPredScore,
                          predModel, predResultFileDir)

    # evaluate the prediction results
    ##predDataAnnotFileDir: the annotation of pre-miRNAs used for prediction. For each pre-miRNA, the annotation information include: miRNA_ID, pre-miRNA_ID, miRNA_seq, pre-miRNA_seq, pre-miRNA_structure
    ##predResultFileDir: the full path of predicted results output by miRLocator
    ##evalResultFileDir: the full path of evaluation result
    if eval:
        print "start to evaluate the prediction performance at different resolutions\n"
        source.eval_dif_resolutions(predDataAnnotFileDir, predResultFileDir, evalResultFileDir)
    # end if

    print "miRLocator run successfully..."
    return predResultFileDir

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


##################################################################################################
##############################parameters section##################################################
##################################################################################################
##the directory of RNAfold                                                                     ###

