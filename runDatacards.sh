#!/bin/sh
root -l -b -q 'ProduceDataCards.C(-1,2017,0,0,21,"/nfs/dust/cms/user/cardinia/HtoTauTau/HiggsCP/DNN/CMSSW_10_2_16/src/HiggsCP/Outputs/test_1S_3B/NTuples_mt_2017/predictions_2017/")'
root -l -b -q 'ProduceDataCards.C(-1,2017,1,1,21,"/nfs/dust/cms/user/cardinia/HtoTauTau/HiggsCP/DNN/CMSSW_10_2_16/src/HiggsCP/Outputs/test_1S_3B/NTuples_mt_2017/predictions_2017/")'
root -l -b -q AddDataCards.C
mv CombinedDataCard.root shapes/DESY/CP/test_cp/2017/htt_mt.inputs-sm-13TeV-2D.root