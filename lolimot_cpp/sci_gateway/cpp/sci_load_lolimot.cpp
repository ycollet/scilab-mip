#include <LL_Dimension.h>
#include <LL_Mesure.h>
#include <LL_Lolimot.h>

// lolimot_struct = lolimot_load(filename);

extern "C" int sci_load_lolimot(char * fname)
{
  unsigned int Load_NbCoeff, Load_NbDimension, Load_NbPartition;
  string       Load_Name;
  
  FileName.str("");
  FileName << Load_Model_Filename;
  
  if (Display)
    {
      sciprint("Loading the Lolimot network : %s\n", FileName.str().c_str());
    }
  
  ModelFile_Load = new ifstream;
  ModelFile_Load->open(FileName.str().c_str());
  
  if (!ModelFile_Load->is_open())
    {
      Scierror(999, "Load Model Error: can't open file %s\n", FileName.str().c_str());
      
      return 0;
    }
  
  (*ModelFile_Load) >> Load_NbPartition;
  (*ModelFile_Load) >> Load_NbDimension;
  (*ModelFile_Load) >> Load_NbCoeff;
  
  Lolimot->cleanPartitions();
  Lolimot->cleanDimensions();
  Lolimot->cleanMesures();
  
  for(i=0; i<Load_NbPartition; i++)
    {
      Lolimot->addPartition();
    }
  
  for(i=0; i<Load_NbPartition; i++)
    {
      Lolimot->getPartition(i)->setNbDimension(Load_NbDimension);
      
      // We load the boundaries of the partitions
      for(j=0; j<Load_NbDimension; j++)
        {
          (*ModelFile_Load) >> StdResidu;
          Lolimot->getPartition(i)->setDimensionMin(j, StdResidu);
          (*ModelFile_Load) >> StdResidu;
          Lolimot->getPartition(i)->setDimensionMax(j, StdResidu);
          (*ModelFile_Load) >> StdResidu;
          Lolimot->getPartition(i)->setInhibe(j, (bool)StdResidu);
        }
      
      // We store the model corresponding to the partition
      (*ModelFile_Load) >> StdResidu;
      Lolimot->getPartition(i)->setCoeff0(StdResidu);
      (*ModelFile_Load) >> StdResidu;
      Lolimot->getPartition(i)->setFreezeCoeff0((bool)StdResidu);
      
      for(k=0; k<Load_NbCoeff; k++)
        {
          (*ModelFile_Load) >> StdResidu;
          Lolimot->getPartition(i)->setCoeff(k, StdResidu);
          (*ModelFile_Load) >> StdResidu;
          Lolimot->getPartition(i)->setFreezeCoeff(k, (bool)StdResidu);
        }
      
      // We store some complementary informations related to the partitions
      (*ModelFile_Load) >> StdResidu;
      Lolimot->getPartition(i)->setResidu(StdResidu);
      (*ModelFile_Load) >> StdResidu;
      Lolimot->getPartition(i)->SetDoNotCut((bool)StdResidu);
      (*ModelFile_Load) >> Load_Name;
      Lolimot->getPartition(i)->setName(Load_Name);
    }
  
  // We store some informations related to the dimensions
  
  for(i=0; i<Load_NbDimension; i++)
    {
      Lolimot->addDimension("", 0.0, 0.0, 0.0, 0);
    }
  
  for(i=0; i<Load_NbDimension; i++)
    {
      (*ModelFile_Load) >> StdResidu;
      Lolimot->getDimension(i)->setMin(StdResidu);
      (*ModelFile_Load) >> StdResidu;
      Lolimot->getDimension(i)->setMax(StdResidu);
      (*ModelFile_Load) >> StdResidu;
      Lolimot->getDimension(i)->setEcartType(StdResidu);
      (*ModelFile_Load) >> StdResidu;
      Lolimot->getDimension(i)->setNbDiscretisations((int)StdResidu);
      (*ModelFile_Load) >> Load_Name;
      Lolimot->getDimension(i)->setName(Load_Name);
      (*ModelFile_Load) >> StdResidu;
      Lolimot->setDontCutVar(i, (bool)StdResidu);
    }
  
  // We store some informations related to the lolimot network
  (*ModelFile_Load) >> StdResidu;
  Lolimot->setNbMaxPartitions((int)StdResidu);
  (*ModelFile_Load) >> StdResidu;
  Lolimot->setResiduGapPercentage(StdResidu);
  (*ModelFile_Load) >> StdResidu;
  Lolimot->setSigma(StdResidu);
  
  (*ModelFile_Load) >> Load_Name;
  if (Load_Name=="none")
    Lolimot->setCustomArguments("");
  else
    Lolimot->setCustomArguments(Load_Name);
  
  (*ModelFile_Load) >> Load_Name;
  if (Load_Name=="none")
    Lolimot->setCustomReturn_Before("");
  else
    Lolimot->setCustomReturn_Before(Load_Name);
  
  (*ModelFile_Load) >> Load_Name;
  if (Load_Name=="none")
    Lolimot->setCustomReturn_After("");
  else
    Lolimot->setCustomReturn_After(Load_Name);
  
  (*ModelFile_Load) >> StdResidu;
  Lolimot->setUseTransformLog((bool)StdResidu);
  (*ModelFile_Load) >> StdResidu;
  Lolimot->setTransformLogAlpha((double)StdResidu);
  (*ModelFile_Load) >> StdResidu;
  Lolimot->setMembershipThreshold((double)StdResidu);
  
  (*ModelFile_Load) >> StdResidu;
  Lolimot->setTransformLogMin(StdResidu);
  (*ModelFile_Load) >> StdResidu;
  Lolimot->setTransformLogMax(StdResidu);
  (*ModelFile_Load) >> StdResidu;
  Lolimot->setTransformLogEps(StdResidu);
  
  ModelFile_Load->close();
  
  delete ModelFile_Load;

  return 0;
}
