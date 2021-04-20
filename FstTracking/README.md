    ### this is a code package for 2020 FST Cosmic Test Stand @ UIC

#### list
- run list for FST cosmic test stand    
> - FstPed\_*Mod*\_*HV*.list: raw ped scan list
> - FstData\_*Mod*\_*HV*.list: raw data list
> - FstCluster\_*Mod*\_*HV*.list: list for QA & tracking study

---------------

#### data
- directory for raw data file

---------------

#### output
- directory for output file

---------------

#### src

##### IstTracking
- extract pedestal & noise
- find hit cluster
- extract position dependent tracking resolution and efficiency

##### FstUtil
- provide utility class for analysis
- `FstRawHit`: class of raw hit
- `FstCluster`: class of cluster
- `FstTrack`: class of track of raw hit and cluster
- `FstEvent`: class of all the information saved in one event
- `FstCons`: all constants used in the analysis => geometry & alignment & threshold etc.

##### FstClusterMaker
- identify raw hits (4.5 for FST & 6.0 for IST)
- group raw hits into clusters
- find tracking based on hits and clusters in IST1 & IST3
- save all information to FstEvent
- Thershold and number of used time bins can be adjusted as an input argument in `macros/FstClusterMaker/makeFstCluster.C`

##### FstTracking
- read in FstEvent
- generate histograms for resoluiton and efficiency

##### FstQAStudy
- read in FstEvent
- generate histgram for QA study
- generate Hit Maps 
- generate Cluster Maps
- generate Track Hit Maps

---------------

#### macro
- `IstTracking` extract IST position dependent tracking resolution and efficiency.
- `FstUtil` test FstRawHit, FstCluster, FstTrack and FstEvent classes.
- `FstClusterMaker` analyze raw data and produce RawHits, Clusters and Tracks.
- `FstTracking` do the alignment and generate track resolution and efficiency histograms
- `FstQAStudy`: generate Event Display histograms
- `FstPlotMacro`: produce all QA plots for one specific configuration => Threshold & NumOfTimeBins

---------------

#### compile the code
- mkdir build
- cd ./build
- cmake ..
- make install 
- a shared lib can be found at ${PROJECT\_SOURCE\_DIR}/lib

---------------

#### analysis work flow
1. compile the code
2. generate data list for the module and HV
3. `FstClusterMaker`: root -l -b -q makeFstCluster.C
>  - produce TTree with FstEvent
>  - loopAll.sh to loop over all configurations
>  - output file can be sepcified in the macro
4. generate cluster list for Tracking and QA
5. `FstTracking`
>  - use getAlignment\_FSTClusters\_3Layer.C to extract alignment parameters for each sensor
6. update alignment paramters in `src/FstUtil/FstCons.h`
7. re-compile the code and run `FstClusterMaker` with updated alignment parameters 
>  - loopAll.sh to loop over all configurations
8. `FstTracking`: root -l -b -q doFstTracking.C
>  - produce histograms for residual and efficiency
>  - output file canbe sepcified in the macro
>  - loopAll.sh to loop over all configurations
9. `FstQAStudy`: root -l -b -q doFstQAStudy.C
>  - generate TTree for EventDisplay
>  - output file canbe sepcified in the macro
>  - loopAll.sh to loop over all configurations
10. `FstPlotMacro`: . ./plotQA.sh
>  - change **hv** and **config** in `plotQA.sh` for different HV and configuration
>  - all QA plots can be found in **figures/module\_hv\_config**
>  - loopAll.sh to loop over all configurations

