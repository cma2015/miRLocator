#!/usr/bin/python

import sys
import os
import argparse


def psource(module):
    file = os.path.basename( module )
    dir = os.path.dirname( module )

    toks = file.split( '.' )
    modname = toks[0]

    # Check if dirrectory is really a directory
    if( os.path.exists( dir ) ):

        # Check if the file directory already exists in the sys.path array
        paths = sys.path
        pathfound = 0
        for path in paths:
            if(dir == path):
                pathfound = 1

        # If the dirrectory is not part of sys.path add it
        if not pathfound:
            sys.path.append( dir )

    # exec works like MEL's eval but you need to add in globals()
    # at the end to make sure the file is imported into the global
    # namespace else it will only be in the scope of this function
    exec ('import ' + modname) in globals()

    # reload the file to make sure its up to date
    exec( 'reload( ' + modname + ' )' ) in globals()

    # This returns the namespace of the file imported
    return modname
#end


def train():
    #create file directories for recording dp_ss files for training and prediction
    trainDataFileDir_Ref = source.checkFileForTraining(trainDataFileDir, RNAFoldDic, sourceDic, dpSSFileDic_train, dpSSFlag = True )
    ##cross_validation_test
    if cross_validation_flag != 0:
        print "start cross validation\n"
        source.cross_validation_function( cross_validation_flag, trainDataFileDir_Ref, source.minCandidateMiRNALen,
                                          source.maxCandidateMiRNALen, source.upOffSet, source.downOffSet, source.minPredScore,
                                          resultDic, dpSSFileDic_train, RNAFoldDic, tempFileDir )
    #end if
    ##train prediction model on all miRNAs in input file
    ##trainDataFileDir_Ref: refined data for training
    ##dpSSFileDic_train: dp_ss files for training data
    ##RNAFoldDic: the directory of RNAfold
    ##tempFileDir: tempoary file
    ##source.upOffSet: 5 in default
    ##source.downOffSet: 5 in default
    #predModelFileDir: the full path of trained prediction model
    ##resultDic: result directory
    predModel = ""
    print "start to train prediction models.\n"
    predModel = source.train_prediction_model( trainDataFileDir_Ref, dpSSFileDic_train, RNAFoldDic, tempFileDir, source.upOffSet, source.downOffSet, predModelFileDir, resultDic )


def predict():
    #prediction mature miRNAs within pre-miRNA sequences
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
    predDataFileDir_Ref = ""
    print "== start to predict miRNAs ==\n"
    predDataFileDir_Ref = source.checkFileForPrediction(predDataFileDir, RNAFoldDic, sourceDic, dpSSFileDic_pred, dpSSFlag = True )
    source.prediction( predDataFileDir_Ref, dpSSFileDic_pred, RNAFoldDic, tempFileDir,
                       source.minCandidateMiRNALen, source.maxCandidateMiRNALen, source.upOffSet, source.downOffSet, source.minPredScore,
                       predModel, predResultFileDir )
    #evaluate the prediction results
    ##predDataAnnotFileDir: the annotation of pre-miRNAs used for prediction. For each pre-miRNA, the annotation information include: miRNA_ID, pre-miRNA_ID, miRNA_seq, pre-miRNA_seq, pre-miRNA_structure
    ##predResultFileDir: the full path of predicted results output by miRLocator
    ##evalResultFileDir: the full path of evaluation result
    if predDataAnnotFileDir is not None:
        print "start to evaluate the prediction performance at different resolutions\n"
    #source.eval_dif_resolutions( predDataAnnotFileDir, predResultFileDir, evalResultFileDir )
        source.eval_dif_resolutions( predDataAnnotFileDir, predResultFileDir, evalResultFileDir )
    #end if


parse=argparse.ArgumentParser(description="Predict  miRNAs from pre-miRNA sequences")

parse.add_argument('-p', dest="type", help="Select the process to run: training or prediction")
parse.add_argument('-i', dest="inputfile", help="The input file for training or prediction")
parse.add_argument('-o', dest="outDir", help="The folder in which the results will be output")
parse.add_argument('-m', dest="model", help="The model file. For the training process, the model file with the specified location will be generated. While for the prediction process, the model file at the specified location is used")
parse.add_argument('-a', dest="annotation", default= None, help="The annotation of prediction file")
parse.add_argument('-k', type = int, default= 0, dest="cross_validation", help="Default: 0 (no cross validation); otherwise, k-fold cross validation is performed")

options=parse.parse_args()

sourceDic ="./"
psource(sourceDic + "source.py")
RNAFoldDic = "/usr/bin/"

if options.outDir:
    newFile = os.path.join("./", options.outDir)
    if os.path.isdir(newFile):
        os.system("rm -rf "+newFile)
    os.mkdir(newFile)
    resultDic = "./"+options.outDir+"/"
    tempFileDir = resultDic + "tempResult.txt"
    dpSSFileDic_train = resultDic + "dp_ss_train/"
    dpSSFileDic_pred = resultDic + "dp_ss_pred/"
    source.createDict(dpSSFileDic_train)
    source.createDict(dpSSFileDic_pred)

    predResultFileDir = resultDic + "miRLocator_predResults.txt"
    evalResultFileDir = resultDic + "miRLocator_evalResults.txt"


if options.type not in ["training", "prediction"]:
    parse.print_help()
    sys.exit()

print (options.type)
if options.type == "training":
    trainDataFileDir = options.inputfile
    predModelFileDir = options.model
    cross_validation_flag = options.cross_validation
    train()

elif options.type == "prediction":
    predDataFileDir = options.inputfile
    predModelFileDir = options.model
    predModel = source.loadPredModel(predModelFileDir)
    predDataAnnotFileDir = options.annotation
    predict()



