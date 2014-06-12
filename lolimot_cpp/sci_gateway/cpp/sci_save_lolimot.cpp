#include <LL_Dimension.h>
#include <LL_Mesure.h>
#include <LL_Lolimot.h>

// lolimot_save(filename, lolimot_struct);

extern "C" int sci_save_lolimot(char * fname)
{
  FileName.str("");
  FileName << Save_Model_Filename;
  
  if (Display)
    {
      sciprint("Saving the LOLIMOT network : %s\n", FileName.str().c_str());
    }
  
  ModelFile_Save = new ofstream;
  ModelFile_Save->open(FileName.str().c_str());
  
  if (!ModelFile_Save->is_open())
    {
      Scierror(999,"Save Model Error: can't open file %s\n", FileName.str().c_str());
      
      return 0;
    }
  
  (*ModelFile_Save) << Lolimot->getPartitionSet().size() << endl;
  (*ModelFile_Save) << Lolimot->getDimensionSet().size() << endl;
  (*ModelFile_Save) << Lolimot->getPartition(0)->getCoeffSet().size() << endl;
  
  for(i=0; i<Lolimot->getPartitionSet().size(); i++)
    {
      // We store the bounds of the partitions
      for(j=0; j<Lolimot->getDimensionSet().size(); j++)
        {
          (*ModelFile_Save) << Lolimot->getPartition(i)->getDimensionMin(j) << endl;
          (*ModelFile_Save) << Lolimot->getPartition(i)->getDimensionMax(j) << endl;
          (*ModelFile_Save) << Lolimot->getPartition(i)->getInhibe(j) << endl;
        }
      
      // We store the model corresponding to the partition
      (*ModelFile_Save) << Lolimot->getPartition(i)->getCoeff0() << endl;
      (*ModelFile_Save) << Lolimot->getPartition(i)->getFreezeCoeff0() << endl;
      
      for(k=0; k<Lolimot->getPartition(i)->getCoeffSet().size(); k++)
        {
          (*ModelFile_Save) << Lolimot->getPartition(i)->getCoeff(k) << endl;
          (*ModelFile_Save) << Lolimot->getPartition(i)->getFreezeCoeff(k) << endl;
        }
      
      // We store some complementary informations related to the partitions
      (*ModelFile_Save) << Lolimot->getPartition(i)->getResidu() << endl;
      (*ModelFile_Save) << Lolimot->getPartition(i)->DoNotCut() << endl;
      (*ModelFile_Save) << Lolimot->getPartition(i)->getName() << endl;
    }
  
  // We store some informations related to the dimensions
  for(i=0; i<Lolimot->getDimensionSet().size(); i++)
    {
      (*ModelFile_Save) << Lolimot->getDimension(i)->getMin() << endl;
      (*ModelFile_Save) << Lolimot->getDimension(i)->getMax() << endl;
      (*ModelFile_Save) << Lolimot->getDimension(i)->getEcartType() << endl;
      (*ModelFile_Save) << Lolimot->getDimension(i)->getNbDiscretisations() << endl;
      (*ModelFile_Save) << Lolimot->getDimension(i)->getName() << endl;
      (*ModelFile_Save) << Lolimot->getDontCutVar(i) << endl;
      
    }
  
  // We store some informations related to the LOLIMOT network
  (*ModelFile_Save) << Lolimot->getNbMaxPartitions() << endl;
  (*ModelFile_Save) << Lolimot->getResiduGapPercentage() << endl;
  (*ModelFile_Save) << Lolimot->getSigma() << endl;
  
  // We store the customisation
  if (Lolimot->getCustomArguments()!="")
    (*ModelFile_Save) << Lolimot->getCustomArguments() << endl;
  else
    (*ModelFile_Save) << "none" << endl;
  
  if (Lolimot->getCustomReturn_Before()!="")
    (*ModelFile_Save) << Lolimot->getCustomReturn_Before() << endl;
  else
    (*ModelFile_Save) << "none" << endl;
  
  if (Lolimot->getCustomReturn_After()!="")
    (*ModelFile_Save) << Lolimot->getCustomReturn_After() << endl;
  else
    (*ModelFile_Save) << "none" << endl;
  
  (*ModelFile_Save) << Lolimot->getUseTransformLog() << endl;
  (*ModelFile_Save) << Lolimot->getTransformLogAlpha() << endl;
  (*ModelFile_Save) << Lolimot->getMembershipThreshold() << endl;
  (*ModelFile_Save) << Lolimot->getTransformLogMin() << endl;
  (*ModelFile_Save) << Lolimot->getTransformLogMax() << endl;
  (*ModelFile_Save) << Lolimot->getTransformLogEps() << endl;
  
  ModelFile_Save->close();
  
  delete ModelFile_Save;

  return 0;
}
