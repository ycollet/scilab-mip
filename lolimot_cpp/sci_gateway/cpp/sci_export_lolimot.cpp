#include <LL_Dimension.h>
#include <LL_Mesure.h>
#include <LL_Lolimot.h>

#include <ExportC.h>
#include <ExportCpp.h>
#include <ExportMatlab.h>

int sci_export_lolimot(char * fname)
{
  bool ExportC = false, ExportCpp = false, ExportMatlab = false, NbColSet = false;
  bool D_ExportC = false, D_ExportCpp = false, D_ExportMatlab = false, Echant_AddMinMax = false;
  bool DD_ExportC = false, DD_ExportCpp = false, DD_ExportMatlab = false;
  bool D_ExportC2 = false, D_ExportCpp2 = false, D_ExportMatlab2 = false, UseLasso = false;
  bool DD_ExportC2 = false, DD_ExportCpp2 = false, DD_ExportMatlab2 = false;

  GET_PARAM_DOUBLE("exportc", tmp_double, 0, tmp_res);
  if ((tmp_res!=-1)&&(tmp_double!=1)) ExportC = true;

  GET_PARAM_DOUBLE("exportcpp", tmp_double, 0, tmp_res);
  if ((tmp_res!=-1)&&(tmp_double!=1)) ExportCpp = true;

  GET_PARAM_DOUBLE("exportmatlab", tmp_double, 0, tmp_res);
  if ((tmp_res!=-1)&&(tmp_double!=1)) ExportMatlab = true;

  GET_PARAM_DOUBLE("dexportc", tmp_double, 0, tmp_res);
  if ((tmp_res!=-1)&&(tmp_double!=1)) D_ExportC = true;

  GET_PARAM_DOUBLE("dexportc2", tmp_double, 0, tmp_res);
  if ((tmp_res!=-1)&&(tmp_double!=1)) D_ExportC2 = true;

  GET_PARAM_DOUBLE("dexportcpp", tmp_double, 0, tmp_res);
  if ((tmp_res!=-1)&&(tmp_double!=1)) D_ExportCpp = true;

  GET_PARAM_DOUBLE("dexportcpp2", tmp_double, 0, tmp_res);
  if ((tmp_res!=-1)&&(tmp_double!=1)) D_ExportCpp2 = true;

  GET_PARAM_DOUBLE("dexportmatlab", tmp_double, 0, tmp_res);
  if ((tmp_res!=-1)&&(tmp_double!=1)) D_ExportMatlab = true;

  GET_PARAM_DOUBLE("dexportmatlab2", tmp_double, 0, tmp_res);
  if ((tmp_res!=-1)&&(tmp_double!=1)) D_ExportMatlab2 = true;

  GET_PARAM_DOUBLE("ddexportc", tmp_double, 0, tmp_res);
  if ((tmp_res!=-1)&&(tmp_double!=1)) DD_ExportC = true;

  GET_PARAM_DOUBLE("ddexportc2", tmp_double, 0, tmp_res);
  if ((tmp_res!=-1)&&(tmp_double!=1)) DD_ExportC2 = true;

  GET_PARAM_DOUBLE("ddexportcpp", tmp_double, 0, tmp_res);
  if ((tmp_res!=-1)&&(tmp_double!=1)) DD_ExportCpp = true;

  GET_PARAM_DOUBLE("ddexportcpp2", tmp_double, 0, tmp_res);
  if ((tmp_res!=-1)&&(tmp_double!=1)) DD_ExportCpp2 = true;

  GET_PARAM_DOUBLE("ddexportmatlab", tmp_double, 0, tmp_res);
  if ((tmp_res!=-1)&&(tmp_double!=1)) DD_ExportMatlab = true;

  GET_PARAM_DOUBLE("ddexportmatlab2", tmp_double, 0, tmp_res);
  if ((tmp_res!=-1)&&(tmp_double!=1)) DD_ExportMatlab2 = true;

  if (ExportMatlab)
    {
      FileName.str("");
      FileName << ModelName << List_Model_Output[Measure].Name;
      if (Display)
	{
	  sciprint("Export LOLIMOT model to matlab : %s\n", FileName.str().c_str());
	}
      exportFunctionInMatlab(Lolimot, FileName.str(), List_Model_Input);
    }
  
  if (ExportCpp)
    {
      FileName.str("");
      FileName << ModelName << List_Model_Output[Measure].Name;
      if (Display)
	{
	  sciprint("Export LOLIMOT model to c++ : %s\n", FileName.str().c_str());
	}
      exportFunctionInCpp(Lolimot, FileName.str(), List_Model_Input);
      exportHeaderInCpp(Lolimot, FileName.str(), List_Model_Input);
    }
  
  if (ExportC)
    {
      FileName.str("");
      FileName << ModelName << List_Model_Output[Measure].Name;
      if (Display)
	{
	  sciprint("Export LOLIMOT model to c : %s\n", FileName.str().c_str());
	}
      exportFunctionInC(Lolimot, FileName.str(), List_Model_Input);
      exportHeaderInC(Lolimot, FileName.str(), List_Model_Input);
    }
  
  if (D_ExportMatlab2)
    {
      FileName.str("");
      FileName << "D2_" << ModelName << List_Model_Output[Measure].Name;
      if (Display)
	{
	  sciprint("Export type 2 partial derivatives LOLIMOT model to matlab : %s\n", FileName.str().c_str());
	}
      exportDerivativeFunctionInMatlab2(Lolimot, FileName.str(), List_Model_Input);
    }
  
  if (D_ExportCpp2)
    {
      FileName.str("");
      FileName << "D2_" << ModelName << List_Model_Output[Measure].Name;
      if (Display)
	{
	  sciprint("Export type 2 partial derivatives LOLIMOT model to c++ : %s\n", FileName.str().c_str());
	}
      exportDerivativeFunctionInCpp2(Lolimot, FileName.str(), List_Model_Input);
      exportHeaderInCpp2(Lolimot, FileName.str(), List_Model_Input);
    }
  
  if (D_ExportC2)
    {
      FileName.str("");
      FileName << "D2_" << ModelName << List_Model_Output[Measure].Name;
      if (Display)
	{
	  sciprint("Export type 2 partial derivatives LOLIMOT model to c : %s\n", FileName.str().c_str());
	}
      exportDerivativeFunctionInC2(Lolimot, FileName.str(), List_Model_Input);
      exportHeaderInC2(Lolimot, FileName.str(), List_Model_Input);
    }
  
  if (DD_ExportMatlab2)
    {
      FileName.str("");
      FileName << "DD2_" << ModelName << List_Model_Output[Measure].Name;
      if (Display)
	{
	  sciprint("Export type 2 second partial derivatives LOLIMOT model to matlab : %s\n", FileName.str().c_str());
	}
      exportSecondDerivativeFunctionInMatlab2(Lolimot, FileName.str(), List_Model_Input);
    }
  
  if (DD_ExportCpp2)
    {
      FileName.str("");
      FileName << "DD2_" << ModelName << List_Model_Output[Measure].Name;
      if (Display)
	{
	  sciprint("Export type 2 second partial derivatives LOLIMOT model to c++ : %s\n", FileName.str().c_str());
	}
      exportSecondDerivativeFunctionInCpp2(Lolimot, FileName.str(), List_Model_Input);
      exportHeaderInCpp2(Lolimot, FileName.str(), List_Model_Input);
    }
  
  if (DD_ExportC2)
    {
      FileName.str("");
      FileName << "DD2_" << ModelName << List_Model_Output[Measure].Name;
      if (Display)
	{
	  sciprint("Export type 2 second partial derivatives LOLIMOT model to c : %s\n", FileName.str().c_str());
	}
      exportSecondDerivativeFunctionInC2(Lolimot, FileName.str(), List_Model_Input);
      exportHeaderInC2(Lolimot, FileName.str(), List_Model_Input);
    }
  
  for(i=0; i<Lolimot->getDimensionSet().size(); i++)
    {
      if (D_ExportMatlab)
	{
	  FileName.str("");
	  FileName << "D_" << i << "_" << ModelName << List_Model_Output[Measure].Name;
	  if (Display)
	    {
	      sciprint("Export partial derivatives LOLIMOT model to matlab : %s\n", FileName.str().c_str());
	    }
	  exportDerivativeFunctionInMatlab(Lolimot, FileName.str(), i+1, List_Model_Input);
	}
      
      if (D_ExportCpp)
	{
	  FileName.str("");
	  FileName << "D_" << i << "_" << ModelName << List_Model_Output[Measure].Name;
	  if (Display)
	    {
	      sciprint("Export partial derivatives LOLIMOT model to c++ : %s\n", FileName.str().c_str());
	    }
	  exportDerivativeFunctionInCpp(Lolimot, FileName.str(), i, List_Model_Input);
	  exportHeaderInCpp(Lolimot, FileName.str(), List_Model_Input);
	}
      
      if (D_ExportC)
	{
	  FileName.str("");
	  FileName << "D_" << i << "_" << ModelName << List_Model_Output[Measure].Name;
	  if (Display)
	    {
	      sciprint("Export partial derivatives LOLIMOT model to c : %s\n", FileName.str().c_str());
	    }
	  exportDerivativeFunctionInC(Lolimot, FileName.str(), i, List_Model_Input);
	  exportHeaderInC(Lolimot, FileName.str(), List_Model_Input);
	}
    }
  
  
  for(i=0; i<Lolimot->getDimensionSet().size(); i++)
    {
      for(j=0; j<Lolimot->getDimensionSet().size(); j++)
	{
	  if (DD_ExportMatlab)
	    {
	      FileName.str("");
	      FileName << "DD_" << i << "_" << j << "_" << ModelName << List_Model_Output[Measure].Name;
	      if (Display)
		{
		  sciprint("Export second partial derivatives LOLIMOT model to matlab : %s\n", FileName.str().c_str());
		}
	      exportSecondDerivativeFunctionInMatlab(Lolimot, FileName.str(), i+1, j+1, List_Model_Input);
	    }
	  
	  if (DD_ExportCpp)
	    {
	      FileName.str("");
	      FileName << "DD_" << i << "_" << j << "_" << ModelName << List_Model_Output[Measure].Name;
	      if (Display)
		{
		  sciprint("Export second partial derivatives LOLIMOT model to c++ : %s\n", FileName.str().c_str());
		}
	      exportSecondDerivativeFunctionInCpp(Lolimot, FileName.str(), i, j, List_Model_Input);
	      exportHeaderInCpp(Lolimot, FileName.str(), List_Model_Input);
	    }
	  
	  if (DD_ExportC)
	    {
	      FileName.str("");
	      FileName << "DD_" << i << "_" << j << "_" << ModelName << List_Model_Output[Measure].Name;
	      if (Display)
		{
		  sciprint("Export second partial derivatives LOLIMOT model to c : %s\n", FileName.str().c_str());
		}
	      exportSecondDerivativeFunctionInC(Lolimot, FileName.str(), i, j, List_Model_Input);
	      exportHeaderInC(Lolimot, FileName.str(), List_Model_Input);
	    }
	}
    }

  return 0;
}
