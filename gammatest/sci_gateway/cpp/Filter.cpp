// Filter.cpp
// Author: Y. Collette
// Date: 19/11/2005
// Renault - DR - 64240

#include <vector>
#include <limits>
#include <iostream>
#include <limits>
#include <utility>
#include <algorithm>
#include <fstream>

#include <stdlib.h>
#include <math.h>

//#define LOG 1

using namespace std;

class Compare_Pair
{
public:
  int operator()(const pair<double, unsigned int> & P1,
		 const pair<double, unsigned int> & P2) const
  {
    return (P1.first > P2.first);
  }
};

void Filter_AddMinMax(vector<vector<double> > & InputData,
		      vector<double> & MeasureData,
		      vector<vector<double> > & Selected_InputData,
		      vector<double> & Selected_MeasureData)
{
  unsigned int i, j, k, IndexMin = 0, IndexMax = 0;
  double        Min, Max;
  bool         Same = true;

  for(i=0; i<InputData[0].size(); i++)
    {
      // Recherche du min et du max pour la dimension i

#ifdef LOG
      cout << "Filter_AddMinMax:: adding points corresponding to min and max value of dimension " << i << endl;
#endif

      Min = numeric_limits<double>::max();
      Max = numeric_limits<double>::min();

      for(j=0; j<InputData.size(); j++)
	{
	  if (Min>InputData[j][i])
	    {
	      Min = InputData[j][i];
	      IndexMin = j;
	    } /* End If */
	  if (Max<InputData[j][i])
	    {
	      Max = InputData[j][i];
	      IndexMax = j;
	    } /* End If */
	} /* End For */

      // Verification que ces deux points ne sont pas déjà dans Selected_InputData
      // Vérification du point min
      for(j=0; j<Selected_InputData.size(); j++)
	{
	  Same = true;

	  for(k=0; k<Selected_InputData[j].size(); k++)
	    {
	      Same = Same && (Selected_InputData[j][k]==InputData[IndexMin][k]);
	    } /* End For */

	  if (Same) break;
	} /* End For */

      if (!Same) 
	{
	  Selected_InputData.push_back(InputData[IndexMin]);
	  Selected_MeasureData.push_back(MeasureData[IndexMin]);
	} /* End If */

      // Vérification du point max
      for(j=0; j<Selected_InputData.size(); j++)
	{
	  Same = true;

	  for(k=0; k<Selected_InputData[j].size(); k++)
	    {
	      Same = Same && (Selected_InputData[j][k]==InputData[IndexMax][k]);
	    } /* End For */

	  if (Same) break;
	} /* End For */

      if (!Same) 
	{
	  Selected_InputData.push_back(InputData[IndexMax]);
	  Selected_MeasureData.push_back(MeasureData[IndexMax]);
	} /* End If */
    } /* End For */


  for(i=0; i<MeasureData.size(); i++)
    {
      // Recherche du min et du max pour la dimension i

#ifdef LOG
      cout << "Filter_AddMinMax:: adding points corresponding to min and max value of measure " << i << endl;
#endif

      Min = numeric_limits<double>::max();
      Max = numeric_limits<double>::min();

      if (Min>MeasureData[i])
	{
	  Min = MeasureData[i];
	  IndexMin = i;
	} /* End If */
      if (Max<MeasureData[i])
	{
	  Max = MeasureData[i];
	  IndexMax = i;
	} /* End If */

      // Verification que ces deux points ne sont pas déjà dans Selected_MeasureData
      // Vérification du point min
      for(j=0; j<Selected_MeasureData.size(); j++)
	{
	  if (Selected_MeasureData[j]==MeasureData[IndexMin]) break;
	} /* End For */

      if (!Same) 
	{
	  Selected_InputData.push_back(InputData[IndexMin]);
	  Selected_MeasureData.push_back(MeasureData[IndexMin]);
	} /* End If */

      // Vérification du point max
      for(j=0; j<Selected_MeasureData.size(); j++)
	{
	  if (Selected_MeasureData[j]==MeasureData[IndexMax]) break;
	} /* End For */

      if (!Same) 
	{
	  Selected_InputData.push_back(InputData[IndexMax]);
	  Selected_MeasureData.push_back(MeasureData[IndexMax]);
	} /* End If */
    } /* End For */
}

void Filter_Select(vector<vector<double> > & InputData,
		   vector<double> & MeasureData,
		   unsigned int Select_Step)
{ 
  vector<vector<double> > Selected_InputData;
  vector<double>          Selected_MeasureData;
  vector<double>          OneInput, OneMeasure;
  unsigned int            i, j;

  OneInput.resize(InputData[0].size());
  OneMeasure.resize(MeasureData.size());
  
  Selected_MeasureData.resize(0);
  Selected_InputData.resize(0);
  
  for(i=0; i<InputData.size(); i+=Select_Step)
    {
      for(j=0; j<InputData[i].size(); j++)
	{
	  OneInput[j] = InputData[i][j];
	} /* End For j */
	  
      Selected_InputData.push_back(OneInput);
      Selected_MeasureData.push_back(MeasureData[i]);
    } /* End For */

  InputData   = Selected_InputData;
  MeasureData = Selected_MeasureData;

  Selected_InputData.clear();
  Selected_MeasureData.clear();
  OneInput.clear();
  OneMeasure.clear();
}

void Filter_EchantStat(vector<vector<double> > & InputData,
		       vector<double> & MeasureData,
		       unsigned int NbBar,
		       unsigned int NbPoints)
{ 
  vector<vector<double> > Selected_InputData, TempList_InputData;
  vector<double>          Selected_MeasureData, TempList_MeasureData;
  vector<unsigned int>   TempList_Index;
  vector<double>          Min, Max;
  vector<int>            Histogram, CumSum;
  unsigned int           NbPointsPerVar, Index, countNbPointsPerVar;
  double                 RandValue;
  vector<bool>           Selected;
  unsigned int           i, j;
  unsigned int           Index_Min, Index_Max;

  Selected_MeasureData.resize(0);
  Selected_InputData.resize(0);

  Selected.resize(InputData.size(), false);

  Histogram.resize(NbBar, 0);
  CumSum.resize(NbBar, 0);

  Max.resize(InputData[0].size(), numeric_limits<double>::min());
  Min.resize(InputData[0].size(), numeric_limits<double>::max());

  NbPointsPerVar = (int)(NbPoints / (double)InputData[0].size());

#ifdef LOG
  cout << "Echantillonnage: préparation des listes min et max" << endl;
#endif

  for(i=0; i<InputData.size(); i++)
    {
      for(j=0; j<InputData[i].size(); j++)
	{
	  if (Min[j]>InputData[i][j]) Min[j] = InputData[i][j];
	  if (Max[j]<InputData[i][j]) Max[j] = InputData[i][j];
	} /* End For j */
    } /* End For i */

  for(i=0; i<InputData[0].size(); i++)
    {
      // Construction de l'histogramme suivant la variable i

#ifdef LOG
      cout << "Echantillonnage: préparation de l'histogramme pour sélection suivant la variable " << i << endl;
#endif

      for(j=0; j<Histogram.size(); j++)
	{
	  Histogram[j] = 0;
	} /* End For j */
      
      for(j=0; j<InputData.size(); j++)
	{
	  Index = (unsigned int)((InputData[j][i] - Min[i])/(Max[i] - Min[i]) * (Histogram.size()-1));
	  Histogram[Index]++;
	} /* End For j */
      
      // Construction de la somme cumulée
#ifdef LOG
      cout << "Echantillonnage: préparation de la somme cumulée" << endl;
#endif

      CumSum[0] = Histogram[0];
      for(j=1; j<Histogram.size(); j++)
	{
	  CumSum[j] = CumSum[j-1] + Histogram[j];
	} /* End For j */
	        
      countNbPointsPerVar = 0;
      
      while(countNbPointsPerVar<NbPointsPerVar)
	{
#ifdef LOG
	  cout << "Echantillonnage: tirage du point " << countNbPointsPerVar << " sur " << NbPointsPerVar << endl;
#endif
	  // Tirage suivant la distribution de proba définie par histogram
	  RandValue = (CumSum[CumSum.size()-1] - CumSum[0])*(rand()/(double)RAND_MAX) + CumSum[0];

	  Index = 0;
	  while ((CumSum[Index]<=RandValue)&&(Index<CumSum.size())) Index++;

	  // On décale l'Index de façon à ne pas être embêté par des tests sur Index.
	  //Echantillonnage moins rigoureux, mais plus facile à implémenter.
	  if (Index>=Histogram.size()-1) Index = Histogram.size()-2;
	  
	  // On récupère les valeurs contenues dans cette barre
	  TempList_InputData.resize(0);
	  TempList_MeasureData.resize(0);
	  TempList_Index.resize(0);

	  Index_Min = Index;
	  Index_Max = Index+1;

	  double Bound_0 = (Index_Min/(double)(Histogram.size()-1)*(Max[i] - Min[i]) + Min[i]);
	  double Bound_1 = (Index_Max/(double)(Histogram.size()-1)*(Max[i] - Min[i]) + Min[i]);

	  // Si après le tirage suivant la distribution Histogram l'ensemble TempList_Index
	  // est vide, on élargit les bornes soit en diminuant Index_Min, soit en augmentant Index_Max.
	  // Si on atteint 0 ou Histogram.size()-1 comme limite, alors on laisse tomber et on passe
	  // à la variable suivante.
	  // On procède de la même manière avec EchantDyn.

	  while((TempList_Index.size()==0)&&(Index_Min>=0)&&(Index_Max<=Histogram.size()-1))
	    {
	      for(j=0; j<InputData.size(); j++)
		{
		  if ((InputData[j][i]<=Bound_1)&&(InputData[j][i]>=Bound_0)&&(!Selected[j]))
		    {
		      TempList_InputData.push_back(InputData[j]);
		      TempList_MeasureData.push_back(MeasureData[j]);
		      TempList_Index.push_back(j);
		    } /* End If */
		} /* End For j */
	
	      if (rand()/(double)RAND_MAX>0.5) Index_Min--;
	      else                             Index_Max++;
	      
	      Bound_0 = (Index_Min/(double)(CumSum.size()-1)*(Max[i] - Min[i]) + Min[i]);
	      Bound_1 = (Index_Max/(double)(CumSum.size()-1)*(Max[i] - Min[i]) + Min[i]);	   
	    } /* End While */

	  // Tirage aléatoire d'une valeur dans cette liste
	  if (TempList_InputData.size()>0)
	    {
	      Index = (int)(rand()/(double)RAND_MAX * (TempList_InputData.size()-1));
	      Selected_InputData.push_back(TempList_InputData[Index]);
	      Selected_MeasureData.push_back(TempList_MeasureData[Index]);
	      Selected[TempList_Index[Index]] = true;
	      countNbPointsPerVar++;
	    } /* End If */
	} /* End While */

#ifdef LOG
      cout << "Echantillonnage: variable " << i << " terminée" << endl;
#endif
    } /* End For i */

  // Ajout des points correspondant aux valeurs min et max pour chaque dimension et chaque mesure
  Filter_AddMinMax(InputData, MeasureData, Selected_InputData, Selected_MeasureData);
  
  InputData   = Selected_InputData;
  MeasureData = Selected_MeasureData;
  
  Selected_InputData.clear();
  Selected_MeasureData.clear();
}

void Filter_EchantRand(vector<vector<double> > & InputData,
		       vector<double> & MeasureData,
		       unsigned int NbPoints)
{ 
  vector<vector<double> > Selected_InputData;
  vector<double>          Selected_MeasureData;
  vector<bool>            Selected;
  unsigned int            RandIndex, countNbPoints = 0;

  Selected_MeasureData.resize(0);
  Selected_InputData.resize(0);

  Selected.resize(InputData.size(), false);

  while(countNbPoints<NbPoints)
    {
      RandIndex = (unsigned int)(rand()/(double)RAND_MAX * InputData.size());

      while ((Selected[RandIndex])&&(RandIndex<=InputData.size()-1))
	RandIndex = (unsigned int)(rand()/(double)RAND_MAX * (InputData.size()-1));

      Selected_InputData.push_back(InputData[RandIndex]);
      Selected_MeasureData.push_back(MeasureData[RandIndex]);
      Selected[RandIndex] = true;
      countNbPoints++;
#ifdef LOG
      cout << "Filter_EchantRand:: point " << RandIndex << " sélectionnée - ";
      cout << countNbPoints << " / " << NbPoints << endl;
#endif
    } /* End While */
#ifdef LOG
  cout << "Echantillonnage: terminée" << endl;
#endif

  // Ajout des points correspondant aux valeurs min et max pour chaque dimension et chaque mesure
  Filter_AddMinMax(InputData, MeasureData, Selected_InputData, Selected_MeasureData);

  InputData   = Selected_InputData;
  MeasureData = Selected_MeasureData;
  
  Selected_InputData.clear();
  Selected_MeasureData.clear();
  Selected.clear();
}

void Filter_EchantGreedy(vector<vector<double> > & InputData,
			 vector<double> & MeasureData,
			 vector<double> & InputMin,
			 vector<double> & InputMax,
			 unsigned int NbPtsToSelect)
{ 
  vector<vector<double> > Distance_Matrix;
  vector<bool>            Selected;
  vector<vector<double> > Selected_InputData;
  vector<double>          Selected_MeasureData;
  double                  Distance, D_Aux;
  unsigned int            NbPtsToRemove, Index = 0, NbPtsRemoved = 0;
  bool                    ResultMin, ResultMax;
  unsigned int            i, j, k, l;

  Selected_MeasureData.resize(0);
  Selected_InputData.resize(0);
  
  Selected.resize(InputData.size(), false);

  Distance_Matrix.resize(InputData.size());
  for(i=0; i<InputData.size(); i++)
    {
      Distance_Matrix[i].resize(InputData.size(), 0.0);
    } /* End For */

  // Calcule de la matrice des distances

#ifdef LOG
  cout << "Filter_Greedy:: Calcule de la matrice des distances" << endl;
#endif

  for(i=0; i<InputData.size()-1; i++)
    {
      for(j=i+1; j<InputData.size(); j++)
	{
	  Distance = 0.0;
	  for(k=0; k<InputData[0].size(); k++)
	    {
	      D_Aux = (InputData[i][k] - InputData[j][k]) / (InputMax[k] - InputMin[k]);
	      Distance += D_Aux * D_Aux;
	    } /* End For k */

	  Distance_Matrix[i][j] = sqrt(Distance);
	} /* End For j */
    } /* End For i */

  // Suppression des N points les plus proches

  NbPtsToRemove = InputData.size() - NbPtsToSelect;

#ifdef LOG
  cout << "Filter_Greedy:: Selection de " << NbPtsToRemove << " points à supprimer" << endl;
#endif

  for(i=0; i<NbPtsToRemove; i++)
    {
      Distance = numeric_limits<double>::max();
      for(j=0; j<InputData.size()-1; j++)
	{
	  if (Selected[j]) continue;
	  for(k=j+1; k<InputData.size(); k++)
	    {
	      if (Selected[k]) continue;

	      if (Distance>Distance_Matrix[j][k])
		{
		  ResultMin = false;
		  ResultMax = false;

		  // On regarde si une coordonnée vaut Min ou Max.
		  // Si oui, on ne sélectionne pas ce point car il correspond à une
		  // valeur extrême

		  for(l=0; l<InputData[k].size(); l++)
		    {
		      ResultMin = ResultMin || (InputData[k][l]==InputMin[l]);
		      ResultMax = ResultMax || (InputData[k][l]==InputMax[l]);
		    } /* End For */

		  if (!ResultMin && !ResultMax)
		    {
		      Distance = Distance_Matrix[j][k];
		      Index    = k;
		    } /* End If */
		} /* End If */
	    } /* End For */
	} /* End For */
      Selected[Index] = true;
      NbPtsRemoved++;
#ifdef LOG
      cout << "Filter_Greedy:: Point " << Index << " supprimé - ";
      cout << NbPtsRemoved << " / " << NbPtsToRemove << endl;
#endif
    } /* End For */

  // Suppression de la matrice des distances

  for(i=0; i<InputData.size(); i++)
    {
      Distance_Matrix[i].clear();
    } /* End For */
  Distance_Matrix.clear();

  // Création de la liste des points sélectionnés

#ifdef LOG
  cout << "Filter_Greedy:: Création de la liste des points sélectionnés" << endl;
#endif

  for(i=0; i<InputData.size(); i++)
    {
      if (!Selected[i])
	{
	  Selected_MeasureData.push_back(MeasureData[i]);
	  Selected_InputData.push_back(InputData[i]);
	} /* End If */
    } /* End For */

  // Ajout des points correspondant aux valeurs min et max pour chaque dimension et chaque mesure
  Filter_AddMinMax(InputData, MeasureData, Selected_InputData, Selected_MeasureData);

  InputData   = Selected_InputData;
  MeasureData = Selected_MeasureData;
  
  Selected_InputData.clear();
  Selected_MeasureData.clear();
  Selected.clear();
}

void Filter_EchantGreedy_Fast(vector<vector<double> > & InputData,
			      vector<double> & MeasureData,
			      vector<double> & InputMin,
			      vector<double> & InputMax,
			      unsigned int NbPtsToSelect,
			      unsigned int NbPtsAtTheEnd)
{ 
  vector<double>          Distance_List;
  vector<bool>            Selected;
  vector<vector<double> > Selected_InputData;
  vector<double>          Selected_MeasureData;
  double                  Distance, D_Aux;
  unsigned int            NbPtsToRemove, Index = 0, NbPtsRemoved = 0;
  unsigned int            i, j, k;

  Selected_MeasureData.resize(0);
  Selected_InputData.resize(0);
  
  Selected.resize(InputData.size(), false);

  Distance_List.resize(InputData.size());

  NbPtsToRemove = (InputData.size() - NbPtsAtTheEnd) / NbPtsToSelect;

  for(i=0; i<NbPtsToSelect; i++)
    {
      Index = (unsigned int)(rand()/(double)RAND_MAX * (InputData.size() - 1));
      while (Selected[Index]) Index = (unsigned int)(rand()/(double)RAND_MAX * (InputData.size() - 1));

#ifdef LOG
      cout << "Filter_Greedy_Fast:: Sélection du point " << Index << endl;
#endif
      // Calcule de la liste des distances

#ifdef LOG
      cout << "Filter_Greedy_Fast:: Calcule de la liste des distances" << endl;
#endif

      for(j=0; j<InputData.size(); j++)
	{
	  if (j==Index)
	    {
	      Distance_List[Index] = numeric_limits<double>::max();
	      continue;
	    } /* End If */

	  Distance = 0.0;
	  for(k=0; k<InputData[0].size(); k++)
	    {
	      D_Aux = (InputData[i][k] - InputData[j][k]) / (InputMax[k] - InputMin[i]);
	      Distance += D_Aux * D_Aux;
	    } /* End For k */
	  Distance_List[Index] = sqrt(Distance);
	} /* End For j */

      // Suppression des N points les plus proches

#ifdef LOG
      cout << "Filter_Greedy_Fast:: Selection de " << NbPtsToRemove << " points à supprimer" << endl;
#endif

      NbPtsRemoved = 0;

      for(j=0; j<NbPtsToRemove; j++)
	{
	  Distance = numeric_limits<double>::max();

	  for(k=0; k<InputData.size(); k++)
	    {
	      if (Selected[k]) continue;

	      if (Distance>Distance_List[k])
		{
		  Distance = Distance_List[k];
		  Index    = k;
		} /* End If */
	    } /* End For */

	  Selected[Index] = true;
	  NbPtsRemoved++;
#ifdef LOG
	  cout << "Filter_Greedy_Fast:: Point " << Index << " supprimé - ";
	  cout << NbPtsRemoved << " / " << NbPtsToRemove << endl;
#endif
	} /* End For */
    } /* End For */

  // Suppression de la liste des distances

  Distance_List.clear();

  // Création de la liste des points sélectionnés

#ifdef LOG
  cout << "Filter_Greedy_Fast:: Création de la liste des points sélectionnés" << endl;
#endif

  for(i=0; i<InputData.size(); i++)
    {
      if (!Selected[i])
	{
	  Selected_MeasureData.push_back(MeasureData[i]);
	  Selected_InputData.push_back(InputData[i]);
	} /* End If */
    } /* End For */

  // Ajout des points correspondant aux valeurs min et max pour chaque dimension et chaque mesure
  Filter_AddMinMax(InputData, MeasureData, Selected_InputData, Selected_MeasureData);

  InputData   = Selected_InputData;
  MeasureData = Selected_MeasureData;
  
  Selected_InputData.clear();
  Selected_MeasureData.clear();
  Selected.clear();
}
