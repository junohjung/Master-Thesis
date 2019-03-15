#ifndef ZFSSTRCTRDFQVARIABLES_H
#define ZFSSTRCTRDFQVARIABLES_H

class ZFSStrctrdFQVariables {
 public:

  ZFSStrctrdFQVariables()
    {
      //this is the start configuration,
      //depending on which fq fields are needed
      //it will be altered later
      AVG_U = 0;
      AVG_V = 1;
      AVG_W = 2;
      AVG_P = 3;
      AVG_RHO = 4;
      NU_T = 5;
      SPONGE_RHO = 6;
      SPONGE_RHO_E = 7;
      PRESSURE = 8;
      LAMBDA2 =9;
      VORTX= 10;
      VORTY= 11;
      VORTZ= 12;
      VORTICITY = new ZFSId[3];
      VORTICITY[0] = VORTX;
      VORTICITY[1] = VORTY;
      VORTICITY[2] = VORTZ;
      CELLID=13;
      BLOCKID=14;
      SPONGE_FACTOR =15;
      SLOPEX = 16;
      SLOPEY = 17;
      SLOPEZ = 18;
      SLOPE = new ZFSId[3];
      SLOPE[0] = SLOPEX;
      SLOPE[1] = SLOPEY;
      SLOPE[2] = SLOPEZ;
      MU_T = 19;
      WALLDISTANCE = 20;

      /* FLUC_U = 21; */
      /* FLUC_V = 22; */
      /* FLUC_W = 23; */
      /* FLUC_UU = 24; */
      /* FLUC_VV = 25; */
      /* FLUC_WW = 26; */
      /* FLUC_UV = 27; */
      /* FLUC_UW = 28; */
      /* FLUC_VW = 29; */
      
      /* RECONST_NUT= 30; */
      /* RECONST_NUTILDE = 31; */
      /* NUTILDE = 32; */




      noFQVariables = 0;
      /* maxNoFQVariables = 33; */
      maxNoFQVariables = 21;
      noFQPlaneOutput = 0;
      noFQBoxOutput = 0;
    }

  ~ZFSStrctrdFQVariables()
    {
      delete[] VORTICITY;
      delete[] SLOPE;
    }

  ZFSId AVG_U;
  ZFSId AVG_V;
  ZFSId AVG_W;
  ZFSId AVG_P;
  ZFSId AVG_RHO;
  ZFSId NU_T;
  ZFSId SPONGE_RHO;
  ZFSId SPONGE_RHO_E;
  ZFSId PRESSURE;
  ZFSId LAMBDA2;
  ZFSId VORTX;
  ZFSId VORTY;
  ZFSId VORTZ;
  ZFSId BLOCKID;
  ZFSId CELLID;
  ZFSId SPONGE_FACTOR;
  ZFSId* VORTICITY;
  ZFSId SLOPEX;
  ZFSId SLOPEY;
  ZFSId SLOPEZ;
  ZFSId* SLOPE;
  ZFSId MU_T;
  ZFSId WALLDISTANCE;

  /* ZFSId FLUC_U ; */
  /* ZFSId FLUC_V ; */
  /* ZFSId FLUC_W ; */
  /* ZFSId FLUC_UU; */
  /* ZFSId FLUC_VV ; */
  /* ZFSId FLUC_WW ; */
  /* ZFSId FLUC_UV ; */
  /* ZFSId FLUC_UW ; */
  /* ZFSId FLUC_VW ; */
      
  /* ZFSId RECONST_NUT; */
  /* ZFSId RECONST_NUTILDE; */
  /* ZFSId NUTILDE ; */


  ZFSId noFQVariables;
  ZFSId noFQPlaneOutput;
  ZFSId noFQBoxOutput;
  ZFSId maxNoFQVariables;

  ZFSId* neededFQVariables;
  ZFSBool* outputFQVariables;  //write output only for certain FQ variables
  ZFSBool* planeOutputFQVariables;
  ZFSBool* boxOutputFQVariables;
  ZFSBool* loadedFromRestartFile;
  std::vector<ZFSString> fqNames;
  std::vector<ZFSBool> fqWriteOutput;
  std::vector<ZFSBool> fqWriteOutputPlanes;
  std::vector<ZFSBool> fqWriteOutputBoxes;

  void activateFQField(ZFSId fqIndex, ZFSId currentPos, ZFSBool writeOutput, ZFSBool writeOutputPlanes, ZFSBool writeOutputBoxes)
  {
    fqWriteOutput.push_back(writeOutput);
    fqWriteOutputPlanes.push_back(writeOutputPlanes);
    fqWriteOutputBoxes.push_back(writeOutputBoxes);

    switch(fqIndex)
      {
      case 0:
        AVG_U = currentPos;
        fqNames.push_back("avgU");;
        zfs_log << "FQ-Field: AVG_U activated - Write Output: " << fqWriteOutput[fqWriteOutput.size()-1]
                << " , write PlaneOutput: " << fqWriteOutputPlanes[fqWriteOutputPlanes.size()-1]
                << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size()-1]  << endl;
        break;
      case 1:
        AVG_V = currentPos;
        fqNames.push_back("avgV");
        zfs_log << "FQ-Field: AVG_V activated - Write Output: " << fqWriteOutput[fqWriteOutput.size()-1]
                << " , write PlaneOutput: " << fqWriteOutputPlanes[fqWriteOutputPlanes.size()-1]
                << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size()-1]  << endl;
        break;
      case 2:
        AVG_W = currentPos;
        fqNames.push_back("avgW");
        zfs_log << "FQ-Field: AVG_P activated - Write Output: " << fqWriteOutput[fqWriteOutput.size()-1]
                << " , write PlaneOutput: " << fqWriteOutputPlanes[fqWriteOutputPlanes.size()-1]
                << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size()-1]  << endl;
        break;
      case 3:
        AVG_P = currentPos;
        fqNames.push_back("avgP");
        zfs_log << "FQ-Field: AVG_P activated - Write Output: " << fqWriteOutput[fqWriteOutput.size()-1]
                << " , write PlaneOutput: " << fqWriteOutputPlanes[fqWriteOutputPlanes.size()-1]
                << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size()-1]  << endl;
        break;
      case 4:
        AVG_RHO = currentPos;
        fqNames.push_back("avgRho");
        zfs_log << "FQ-Field: AVG_RHO activated - Write Output: " << fqWriteOutput[fqWriteOutput.size()-1]
                << " , write PlaneOutput: " << fqWriteOutputPlanes[fqWriteOutputPlanes.size()-1]
                << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size()-1] << endl;
        break;
      case 5:
        NU_T = currentPos;
        fqNames.push_back("nu_t");
        zfs_log << "FQ-Field: NU_T activated - Write Output: " << fqWriteOutput[fqWriteOutput.size()-1]
                << " , write PlaneOutput: " << fqWriteOutputPlanes[fqWriteOutputPlanes.size()-1]
                << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size()-1]  << endl;
        break;
      case 6:
        SPONGE_RHO = currentPos;
        fqNames.push_back("spongeRho");
        zfs_log << "FQ-Field: SPONGE_RHO activated - Write Output: " << fqWriteOutput[fqWriteOutput.size()-1]
                << " , write PlaneOutput: " << fqWriteOutputPlanes[fqWriteOutputPlanes.size()-1]
                << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size()-1]  << endl;
        break;
      case 7:
        SPONGE_RHO_E = currentPos;
        fqNames.push_back("spongeRhoE");
        zfs_log << "FQ-Field: SPONGE_RHO_E activated - Write Output: " << fqWriteOutput[fqWriteOutput.size()-1]
                << " , write PlaneOutput: " << fqWriteOutputPlanes[fqWriteOutputPlanes.size()-1]
                << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size()-1]  << endl;
        break;
      case 8:
        PRESSURE = currentPos;
        fqNames.push_back("pressure");
        zfs_log << "FQ-Field: PRESSURE activated - Write Output: " << fqWriteOutput[fqWriteOutput.size()-1]
                << " , write PlaneOutput: " << fqWriteOutputPlanes[fqWriteOutputPlanes.size()-1]
                << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size()-1]  << endl;
        break;
      case 9:
        LAMBDA2 = currentPos;
        fqNames.push_back("lambda2");
        zfs_log << "FQ-Field: LAMBDA_2 activated - Write Output: " << fqWriteOutput[fqWriteOutput.size()-1]
                << " , write PlaneOutput: " << fqWriteOutputPlanes[fqWriteOutputPlanes.size()-1]
                << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size()-1]  << endl;
        break;
      case 10:
        VORTX = currentPos;
        VORTICITY[0] = currentPos;
        fqNames.push_back("vortX");
        zfs_log << "FQ-Field: x-Vorticity activated - Write Output: " << fqWriteOutput[fqWriteOutput.size()-1]
                << " , write PlaneOutput: " << fqWriteOutputPlanes[fqWriteOutputPlanes.size()-1]
                << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size()-1]  << endl;
        break;
      case 11:
        VORTY = currentPos;
        VORTICITY[1] = currentPos;
        fqNames.push_back("vortY");
        zfs_log << "FQ-Field: y-Vorticity activated - Write Output: " << fqWriteOutput[fqWriteOutput.size()-1]
                << " , write PlaneOutput: " << fqWriteOutputPlanes[fqWriteOutputPlanes.size()-1]
                << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size()-1]  << endl;
        break;
      case 12:
        VORTZ = currentPos;
        VORTICITY[2] = currentPos;
        fqNames.push_back("vortZ");
        zfs_log << "FQ-Field: z-Vorticity activated - Write Output: " << fqWriteOutput[fqWriteOutput.size()-1]
                << " , write PlaneOutput: " << fqWriteOutputPlanes[fqWriteOutputPlanes.size()-1]
                << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size()-1]  << endl;
        break;
      case 13:
        CELLID = currentPos;
        fqNames.push_back("cellId");
        zfs_log << "FQ-Field: cellId activated - Write Output: " << fqWriteOutput[fqWriteOutput.size()-1]
                << " , write PlaneOutput: " << fqWriteOutputPlanes[fqWriteOutputPlanes.size()-1]
                << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size()-1]  << endl;
        break;
      case 14:
        BLOCKID = currentPos;
        fqNames.push_back("blockId");
        zfs_log << "FQ-Field: blockId activated - Write Output: " << fqWriteOutput[fqWriteOutput.size()-1]
                << " , write PlaneOutput: " << fqWriteOutputPlanes[fqWriteOutputPlanes.size()-1]
                << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size()-1]  << endl;
        break;
      case 15:
        SPONGE_FACTOR = currentPos;
        fqNames.push_back("spongeFactor");
        zfs_log << "FQ-Field: SPONGE_FACTOR activated - Write Output: " << fqWriteOutput[fqWriteOutput.size()-1]
                << " , write PlaneOutput: " << fqWriteOutputPlanes[fqWriteOutputPlanes.size()-1]
                << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size()-1]  << endl;
        break;
      case 16:
        SLOPEX = currentPos;
        SLOPE[0] = currentPos;
        fqNames.push_back("slopeX");
        zfs_log << "FQ-Field: x-Slope activated - Write Output: " << fqWriteOutput[fqWriteOutput.size()-1]
                << " , write PlaneOutput: " << fqWriteOutputPlanes[fqWriteOutputPlanes.size()-1]
                << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size()-1]  << endl;
        break;
      case 17:
        SLOPEY = currentPos;
        SLOPE[1] = currentPos;
        fqNames.push_back("slopeY");
        zfs_log << "FQ-Field: y-Slope activated - Write Output: " << fqWriteOutput[fqWriteOutput.size()-1]
                << " , write PlaneOutput: " << fqWriteOutputPlanes[fqWriteOutputPlanes.size()-1]
                << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size()-1]  << endl;
        break;
      case 18:
        SLOPEZ = currentPos;
        SLOPE[2] = currentPos;
        fqNames.push_back("slopeZ");
        zfs_log << "FQ-Field: z-Slope activated - Write Output: " << fqWriteOutput[fqWriteOutput.size()-1]
                << " , write PlaneOutput: " << fqWriteOutputPlanes[fqWriteOutputPlanes.size()-1]
                << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size()-1]  << endl;
        break;
      case 19:
        MU_T = currentPos;
        fqNames.push_back("mu_t");
        zfs_log << "FQ-Field: mu_t activated - Write Output: " << fqWriteOutput[fqWriteOutput.size()-1]
                << " , write PlaneOutput: " << fqWriteOutputPlanes[fqWriteOutputPlanes.size()-1]
                << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size()-1]  << endl;
        break;
      case 20:
        WALLDISTANCE = currentPos;
        fqNames.push_back("wallDistance");
        zfs_log << "FQ-Field: wallDistance activated - Write Output: " << fqWriteOutput[fqWriteOutput.size()-1]
                << " , write PlaneOutput: " << fqWriteOutputPlanes[fqWriteOutputPlanes.size()-1]
                << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size()-1]  << endl;
        break;
	/////////junoh/////////
		
      /* case 21: */
      /*   FLUC_U = currentPos; */
      /*   fqNames.push_back("FLUC_U"); */
      /*   zfs_log << "FQ-Field: FLUC_U activated - Write Output: " << fqWriteOutput[fqWriteOutput.size()-1] */
      /*           << " , write PlaneOutput: " << fqWriteOutputPlanes[fqWriteOutputPlanes.size()-1] */
      /*           << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size()-1]  << endl; */
      /*   break; */
		
      /* case 22: */
      /*   FLUC_V = currentPos; */
      /*   fqNames.push_back("FLUC_V"); */
      /*   zfs_log << "FQ-Field: FLUC_V activated - Write Output: " << fqWriteOutput[fqWriteOutput.size()-1] */
      /*           << " , write PlaneOutput: " << fqWriteOutputPlanes[fqWriteOutputPlanes.size()-1] */
      /*           << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size()-1]  << endl; */
      /*   break; */
		
      /* case 23: */
      /*   FLUC_W = currentPos; */
      /*   fqNames.push_back("FLUC_W"); */
      /*   zfs_log << "FQ-Field: FLUC_W activated - Write Output: " << fqWriteOutput[fqWriteOutput.size()-1] */
      /*           << " , write PlaneOutput: " << fqWriteOutputPlanes[fqWriteOutputPlanes.size()-1] */
      /*           << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size()-1]  << endl; */
      /*   break; */
		
      /* case 24: */
      /*   FLUC_UU = currentPos; */
      /*   fqNames.push_back("FLUC_UU"); */
      /*   zfs_log << "FQ-Field: FLUC_UU activated - Write Output: " << fqWriteOutput[fqWriteOutput.size()-1] */
      /*           << " , write PlaneOutput: " << fqWriteOutputPlanes[fqWriteOutputPlanes.size()-1] */
      /*           << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size()-1]  << endl; */
      /*   break; */
		
      /* case 25: */
      /*   FLUC_VV = currentPos; */
      /*   fqNames.push_back("FLUC_VV"); */
      /*   zfs_log << "FQ-Field: FLUC_VV activated - Write Output: " << fqWriteOutput[fqWriteOutput.size()-1] */
      /*           << " , write PlaneOutput: " << fqWriteOutputPlanes[fqWriteOutputPlanes.size()-1] */
      /*           << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size()-1]  << endl; */
      /*   break; */
     
      /* case 26: */
      /*   FLUC_WW = currentPos; */
      /*   fqNames.push_back("FLUC_WW"); */
      /*   zfs_log << "FQ-Field: FLUC_WW activated - Write Output: " << fqWriteOutput[fqWriteOutput.size()-1] */
      /*           << " , write PlaneOutput: " << fqWriteOutputPlanes[fqWriteOutputPlanes.size()-1] */
      /*           << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size()-1]  << endl; */
      /*   break; */
     
      /* case 27: */
      /*   FLUC_UV = currentPos; */
      /*   fqNames.push_back("FLUC_UV"); */
      /*   zfs_log << "FQ-Field: FLUC_UV activated - Write Output: " << fqWriteOutput[fqWriteOutput.size()-1] */
      /*           << " , write PlaneOutput: " << fqWriteOutputPlanes[fqWriteOutputPlanes.size()-1] */
      /*           << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size()-1]  << endl; */
      /*   break; */
     
      /* case 28: */
      /*   FLUC_UW = currentPos; */
      /*   fqNames.push_back("FLUC_UW"); */
      /*   zfs_log << "FQ-Field: FLUC_UW activated - Write Output: " << fqWriteOutput[fqWriteOutput.size()-1] */
      /*           << " , write PlaneOutput: " << fqWriteOutputPlanes[fqWriteOutputPlanes.size()-1] */
      /*           << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size()-1]  << endl; */
      /*   break; */
	

      /* case 29: */
      /*   FLUC_VW = currentPos; */
      /*   fqNames.push_back("FLUC_VW"); */
      /*   zfs_log << "FQ-Field: FLUC_VW activated - Write Output: " << fqWriteOutput[fqWriteOutput.size()-1] */
      /*           << " , write PlaneOutput: " << fqWriteOutputPlanes[fqWriteOutputPlanes.size()-1] */
      /*           << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size()-1]  << endl; */
      /*   break; */
	
      /* case 30: */
      /*   RECONST_NUT = currentPos; */
      /*   fqNames.push_back("RECONST_NUT"); */
      /*   zfs_log << "FQ-Field: RECONST_NUT activated - Write Output: " << fqWriteOutput[fqWriteOutput.size()-1] */
      /*           << " , write PlaneOutput: " << fqWriteOutputPlanes[fqWriteOutputPlanes.size()-1] */
      /*           << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size()-1]  << endl; */
      /*   break; */
	
	
      /* case 31: */
      /*   RECONST_NUTILDE = currentPos; */
      /*   fqNames.push_back("RECONST_NUTILDE"); */
      /*   zfs_log << "FQ-Field: RECONST_NUTILDE activated - Write Output: " << fqWriteOutput[fqWriteOutput.size()-1] */
      /*           << " , write PlaneOutput: " << fqWriteOutputPlanes[fqWriteOutputPlanes.size()-1] */
      /*           << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size()-1]  << endl; */
      /*   break; */
	
      /* case 32: */
      /*   NUTILDE = currentPos; */
      /*   fqNames.push_back("NUTILDE"); */
      /*   zfs_log << "FQ-Field: NUTILDE activated - Write Output: " << fqWriteOutput[fqWriteOutput.size()-1] */
      /*           << " , write PlaneOutput: " << fqWriteOutputPlanes[fqWriteOutputPlanes.size()-1] */
      /*           << " , write BoxOutput: " << fqWriteOutputBoxes[fqWriteOutputBoxes.size()-1]  << endl; */
      /*   break; */



	//////////////junoh//////

      default:
        zfsTerm(1, __CALLING_FUNCTION__, "Can't find the FQ field!");
        break;
      }
  }
};

#endif
