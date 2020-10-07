import FWCore.ParameterSet.Config as cms

#
# Sequence to add lepton MVA
#


def leptonMvaSequence(process, options, tnpVars):
    #
    # One difficulty with lepton mva's is that their input variables are dependent on jet variables, so we need JEC etc... to be in sync
    # By default we simply re-run the JEC and needed b-tag algorithms, to be sure they are in sync with the used global tag, assuming the training was also up to date
    #
    if(options['isMC']): jetCorrectorLevels = ['L1FastJet', 'L2Relative', 'L3Absolute']
    else:                jetCorrectorLevels = ['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual']
    from RecoBTag.MXNet.pfDeepBoostedJet_cff import pfDeepBoostedJetTags, pfMassDecorrelatedDeepBoostedJetTags
    from RecoBTag.MXNet.Parameters.V02.pfDeepBoostedJetPreprocessParams_cfi import pfDeepBoostedJetPreprocessParams as pfDeepBoostedJetPreprocessParamsV02
    from RecoBTag.MXNet.Parameters.V02.pfMassDecorrelatedDeepBoostedJetPreprocessParams_cfi import pfMassDecorrelatedDeepBoostedJetPreprocessParams as pfMassDecorrelatedDeepBoostedJetPreprocessParamsV02
    pfDeepBoostedJetTags.preprocessParams = pfDeepBoostedJetPreprocessParamsV02
    pfDeepBoostedJetTags.model_path = 'RecoBTag/Combined/data/DeepBoostedJet/V02/full/resnet-symbol.json'
    pfDeepBoostedJetTags.param_path = 'RecoBTag/Combined/data/DeepBoostedJet/V02/full/resnet-0000.params'
    pfMassDecorrelatedDeepBoostedJetTags.preprocessParams = pfMassDecorrelatedDeepBoostedJetPreprocessParamsV02
    pfMassDecorrelatedDeepBoostedJetTags.model_path = 'RecoBTag/Combined/data/DeepBoostedJet/V02/decorrelated/resnet-symbol.json'
    pfMassDecorrelatedDeepBoostedJetTags.param_path = 'RecoBTag/Combined/data/DeepBoostedJet/V02/decorrelated/resnet-0000.params'
    from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
    
    updateJetCollection(
       process,
       jetSource = cms.InputTag('slimmedJets'),
       labelName = 'Updated',
       jetCorrections = ('AK4PFchs', cms.vstring(jetCorrectorLevels), 'None'),
       btagDiscriminators = [
        'pfDeepFlavourJetTags:probb',
        'pfDeepFlavourJetTags:probbb',
        'pfDeepFlavourJetTags:problepb',
        ],
        printWarning = False
    )
    leptonMva_sequence = cms.Sequence(process.patAlgosToolsTask)

    #
    # For the calculation of isolations and jet-lep variables we rely on the NanoAOD modules
    # Because in some PAGs people like to use prehistoric effective areas (or some lepton mva developer found 0.00001% better discrimination with those),
    # we need to have the PFIso and MiniIso in all its variations (at least they all use them relative to the lepton pt)
    #
    from PhysicsTools.NanoAOD.electrons_cff import isoForEle, ptRatioRelForEle
    process.ptRatioRelForEle        = ptRatioRelForEle
    process.ptRatioRelForEle.srcJet = cms.InputTag('selectedUpdatedPatJetsUpdated')
    leptonMva_sequence += cms.Sequence(process.ptRatioRelForEle)

    def makeIsoForEle(leptonMva_sequence, name, effAreas):
      isoForEleModule = isoForEle.clone(relative = cms.bool(False))
      setattr(isoForEleModule, 'EAFile_MiniIso', cms.FileInPath(effAreas))
      setattr(isoForEleModule, 'EAFile_PFIso',   cms.FileInPath(effAreas))
      setattr(process, name, isoForEleModule)
      leptonMva_sequence += cms.Sequence(getattr(process, name))

    makeIsoForEle(leptonMva_sequence, 'isoForEleFall17',   'RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt')
    makeIsoForEle(leptonMva_sequence, 'isoForEleSummer16', 'RecoEgamma/ElectronIdentification/data/Summer16/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_80X.txt')
    makeIsoForEle(leptonMva_sequence, 'isoForEleSpring15', 'RecoEgamma/ElectronIdentification/data/Spring15/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_25ns.txt')


    #
    # Calculate the lepton mva's
    #   (at some point we can clean up the older TTH and Ghent ones, keeping only the TOP)
    #https://github.com/cms-data/PhysicsTools-NanoAOD
    
    process.leptonMvaTTH = cms.EDProducer('LeptonMvaProducer',
      leptonMvaType        = cms.string("leptonMvaTTH"),
      weightFile           = cms.FileInPath('EgammaAnalysis/TnPTreeProducer/data/el_BDTG_20%s.weights.xml' % ('16' if '2016' in options['era'] else '17')),
      probes               = cms.InputTag('slimmedElectrons'),
      miniIsoChg           = cms.InputTag('isoForEle%s:miniIsoChg' % ('Spring15'  if '2016' in options['era'] else 'Fall17')),
      miniIsoAll           = cms.InputTag('isoForEle%s:miniIsoAll' % ('Spring15'  if '2016' in options['era'] else 'Fall17')),
      PFIsoAll04           = cms.InputTag('isoForEle%s:PFIsoAll04' % ('Summer16'  if '2016' in options['era'] else 'Fall17')),
      ptRatio              = cms.InputTag('ptRatioRelForEle:ptRatio'),
      ptRel                = cms.InputTag('ptRatioRelForEle:ptRel'),
      jetNDauChargedMVASel = cms.InputTag('ptRatioRelForEle:jetNDauChargedMVASel'),
      closestJet           = cms.InputTag('ptRatioRelForEle:jetForLepJetVar'),                                          
      mvas                 = cms.InputTag('electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV2Values'),
      debug                = cms.bool(False), # set to True if you want to sync with your analysis
    )

    process.leptonMvaGhent = cms.EDProducer('LeptonMvaProducer',
      leptonMvaType        = cms.string("leptonMvaGhent"),
      weightFile           = cms.FileInPath('EgammaAnalysis/TnPTreeProducer/data/el_tZqTTV%s_BDTG.weights.xml' % ('16' if '2016' in options['era'] else '17')),
      probes               = cms.InputTag('slimmedElectrons'),
      miniIsoChg           = cms.InputTag('isoForEle%s:miniIsoChg' % ('Spring15' if '2016' in options['era'] else 'Fall17')),
      miniIsoAll           = cms.InputTag('isoForEle%s:miniIsoAll' % ('Spring15' if '2016' in options['era'] else 'Fall17')),
      PFIsoAll             = cms.InputTag('isoForEle%s:PFIsoAll' % ('Summer16' if '2016' in options['era'] else 'Fall17')),
      ptRatio              = cms.InputTag('ptRatioRelForEle:ptRatio'),
      ptRel                = cms.InputTag('ptRatioRelForEle:ptRel'),
      jetNDauChargedMVASel = cms.InputTag('ptRatioRelForEle:jetNDauChargedMVASel'),
      closestJet           = cms.InputTag('ptRatioRelForEle:jetForLepJetVar'),
      mvas                 = cms.InputTag('electronMVAValueMapProducer:ElectronMVAEstimator%sValues' % ('Run2Spring16GeneralPurposeV1' if '2016' in options['era'] else 'Run2Fall17NoIsoV1')),
      debug                = cms.bool(False), # set to True if you want to sync with your analysis
    )

    process.leptonMvaTOP = cms.EDProducer('LeptonMvaProducer',
      leptonMvaType        = cms.string("leptonMvaTOP"),
      weightFile           = cms.FileInPath('EgammaAnalysis/TnPTreeProducer/data/el_TOP%s_BDTG.weights.xml' % (options['era'].replace('20', '').replace('UL', ''))),
      probes               = cms.InputTag('slimmedElectrons'),
      miniIsoChg           = cms.InputTag('isoForEle%s:miniIsoChg' % ('Spring15' if '2016' in options['era'] else 'Fall17')),
      miniIsoAll           = cms.InputTag('isoForEle%s:miniIsoAll' % ('Spring15' if '2016' in options['era'] else 'Fall17')),
      PFIsoAll             = cms.InputTag('isoForEle%s:PFIsoAll' % ('Summer16' if '2016' in options['era'] else 'Fall17')),
      ptRatio              = cms.InputTag('ptRatioRelForEle:ptRatio'),
      ptRel                = cms.InputTag('ptRatioRelForEle:ptRel'),
      jetNDauChargedMVASel = cms.InputTag('ptRatioRelForEle:jetNDauChargedMVASel'),
      closestJet           = cms.InputTag('ptRatioRelForEle:jetForLepJetVar'),
      mvas                 = cms.InputTag('electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV2Values'),
      debug                = cms.bool(False), # set to True if you want to sync with your analysis
    )

    leptonMva_sequence += cms.Sequence(
      process.leptonMvaTTH #+
      #process.leptonMvaGhent +
      #process.leptonMvaTOP
    )

    #
    # Adding the new variables to the trees
    #   (currently only adding most recent version of miniIso)
    #
    newVariables = {
        'el_leptonMva_ttH'         : cms.InputTag('leptonMvaTTH:leptonMvaTTH'),
        'el_ptAM'                  : cms.InputTag('leptonMvaTTH:ptAM'),
        'el_etaAM'     	      	   : cms.InputTag('leptonMvaTTH:etaAM'),
        'el_miniRelIsoChargedAM'   : cms.InputTag('leptonMvaTTH:miniRelIsoChargedAM'),
        'el_miniRelIsoNeutralAM'   : cms.InputTag('leptonMvaTTH:miniRelIsoNeutralAM'),
        'el_jetPtRelv2AM'     	   : cms.InputTag('leptonMvaTTH:jetPtRelv2AM'),
        'el_jetDFAM'     	   : cms.InputTag('leptonMvaTTH:jetDFAM'),
        'el_jetPtRatioAM'     	   : cms.InputTag('leptonMvaTTH:jetPtRatioAM'),
        'el_dxyAM'     	      	   : cms.InputTag('leptonMvaTTH:dxyAM'),
        'el_sip3dAM'     	   : cms.InputTag('leptonMvaTTH:sip3dAM'),
        'el_dzAM'     	      	   : cms.InputTag('leptonMvaTTH:dzAM'),
        'el_mvaFall17V2noIsoAM'    : cms.InputTag('leptonMvaTTH:mvaFall17V2noIsoAM'), 
        'el_jetPtAM'     	   : cms.InputTag('leptonMvaTTH:jetPtAM'),
        'el_jetEtaAM'     	   : cms.InputTag('leptonMvaTTH:jetEtaAM'),
        'el_jetRelIsoAM'           : cms.InputTag('leptonMvaTTH:jetRelIsoAM'),

      #'el_leptonMva_ghent'       : cms.InputTag('leptonMvaGhent:leptonMvaGhent'),
      #'el_leptonMva_TOP'         : cms.InputTag('leptonMvaTOP:leptonMvaTOP'),
     ## 'el_miniIsoAll_fall17'     : cms.InputTag('isoForEleFall17:miniIsoAll'),
     ## 'el_miniIsoChg_fall17'     : cms.InputTag('isoForEleFall17:miniIsoChg'),
     ## 'el_miniIsoAll_Spring15'   : cms.InputTag('isoForEleSpring15:miniIsoAll'),
     ## 'el_miniIsoChg_Spring15'   : cms.InputTag('isoForEleSpring15:miniIsoChg'),
     ## 'el_relIso_fall17'         : cms.InputTag('isoForEleFall17:PFIsoAll'),
     ## 'el_PFelIso04_Summer16'    : cms.InputTag('isoForEleSummer16:PFIsoAll04'),
     ## 'el_ptRatio'               : cms.InputTag('ptRatioRelForEle:ptRatio'),
     ## 'el_ptRel'                 : cms.InputTag('ptRatioRelForEle:ptRel'),
     ## 'el_closestJetDeepFlavour' : cms.InputTag('leptonMvaTTH:closestJetDeepFlavour'), # For those crazy people who want to add even more cuts on top of their leptonMva but can't tell why they need it
     ## 'el_closestJetpt'          : cms.InputTag('leptonMvaTTH:closestJetpt'),
     ## 'el_closestJeteta'         : cms.InputTag('leptonMvaTTH:closestJeteta'),
    ##'el_jetNDauChargedMVASel'  : cms.InputTag('ptRatioRelForEle:jetNDauChargedMVASel'),
     ## 'el_mvas'                  : cms.InputTag('electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV2Values'),

    }
    for i, j in newVariables.iteritems():
      setattr(tnpVars.CommonStuffForGsfElectronProbe.variables, i, j)

    return leptonMva_sequence
