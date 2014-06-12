// Author: Y. Collette
// EMail: ycollet@freesurf.fr
// Date: 02/04/2007

/*
This file is part of Lolimot.

    Foobar is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Foobar is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#ifndef TRAINLOLIMOTSTRUCT_H
#define TRAINLOLIMOTSTRUCT_H

enum Model_Input_Type {Input = 0, Output = 1, RetroInput = 2, Weight = 3};

struct Model_Input
{
  string       Name;
  unsigned int Pos;
  unsigned int Tempo;
  unsigned int FileNumber;
  bool         DontCut;
  unsigned int GroupNb;
  Model_Input_Type Type;
};

struct Filter_Entry
{
  float        Min;
  float        Max;
  unsigned int InputNb;
  bool         Add;
};

struct Delta_Entry
{
  vector<unsigned int> ListOfInputs;
  vector<int>          ListOfTempo;
};

struct PartToCut
{
  string       PartName;
  unsigned int Dimension;
  float        CutPosition;
  string       PartName_LowerPart;
  string       PartName_UpperPart;
};

struct Model_Coeff
{
  string PartName;
  int    Dimension;
  float  Coeff;
};

// Temporize function
// This function is used to generate some data temporally offseted

float Temporize(vector<float> & ListTempo, float Value);

#endif
