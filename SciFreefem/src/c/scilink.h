#ifndef SCILINK_H
#define SCILINK_H

#include "machine.h"

#include "list.cpp"

Analyseur * scilabana = 0;
enum SciSyntax { withText = 1, withTab };
int getMatProfil = 0; // Pass matrix profil to Scilab ?

extern "C" 
{
  void C2F(delete_lb_lobj)();
  void C2F(delete_lexp)();
  void C2F(add_lobj)();
  void C2F(delete_ptr)(void * ptr);
  void C2F(delete_scilabana)(int * ierr);
  void C2F(freefem_code)(char * Message, int * err);
  void C2F(ff_problem)(char * Message, int * flag, int * err);
  void C2F(get_ff_result)(char * s, int * lhs,
			  double ** tNode, int * mNode, int * nNode, 
			  double ** tFunc, int * mFunc, int * nFunc, 
			  int ** tTriangle, int * mTriangle, 
			  int * nTriangle,
			  char * Message,int * err);
  void C2F(put_scilab_mesh)(char * s,double * tNode, int * mNode,
			    int * tTriangle, int * mTriangle,
			    char * Message, int * err);
  void C2F(put_scilab_border_1)(char * nameBorder,double * tBorder,
				int * nbData,int * normal,
				char * Message,int * err);
  void C2F(put_scilab_border_2)(char * nameBorder,char * equation,
				double * start, double * stop,
				int * ng,int * nrm, 
				char * Message, int * err);
  void C2F(build_scilab_mesh)(char * nameMesh,
			      char * Message,int * err);
  void C2F(get_matrix)(double ** tMatrix, int * sztMatrix, 
		       int ** tjlow,int ** tjhigh, int * sztVect, 
		       int ** szblock, char * Message, int * err);
};

class Border
{
 public:
  double * tCoor;  // Array of points on a segment of the boundary
  int nbCoorData;  // Number of points on a segment of the boundary
  int normal;      // normal == 1 => mesh on the right side
                   // normal == -1 => mesh on the left side
  int label;

  Analvar * an;

  Border() {}
 Border(Analvar * ann,double * tCoor_fp,int nbCoorData_fp, 
	int normal_fp,int label_fp) : an(ann),tCoor(tCoor_fp),
    nbCoorData(nbCoorData_fp),
    normal(normal_fp),label(label_fp) {}
  void eval() const;
};

void Border::eval() const
{
  frontiere & front = *an -> front;
  int num = 0, oldnum;

  for (int i = 0; i < nbCoorData; i++)
    {
      oldnum = num;
      if (front.step) 
        {
          num = front.addPoint(tCoor[i], tCoor[i + nbCoorData], label);

          if (i) front.addSegment(oldnum - 1, num - 1, label); 
        }
      else if (front.initboundingbox)
        {  
          front.xmin = Min((float) tCoor[i],front.xmin);
          front.xmax = Max((float) tCoor[i],front.xmax);
          front.ymin = Min((float) tCoor[i + nbCoorData],front.ymin);
          front.ymax = Max((float) tCoor[i + nbCoorData],front.ymax);
        }
      else
        {  
          front.initboundingbox =1;
          front.xmin = tCoor[i];
          front.xmax = tCoor[i];
          front.ymin = tCoor[i + nbCoorData];
          front.ymax = tCoor[i + nbCoorData];
        }
    }
  if (front.step) 
    {
      front.hh[num - 1] = front.hh[oldnum - 1]; 
      front.sd[2 * front.nbsd] = normal * (front.nbs - 1); 
      front.sd[2 * front.nbsd + 1] = front.nbsd + 1; 
      front.nbsd++;
    }
}

List<Border> lb;
List<Expr *> lexp;
List<int> lobj;

class ESciMesh: public MeshExpr
{ // scilab mesh construction
  Analvar * an;
        
 public:
  ESciMesh (Analvar * ann) : an(ann) {}

  Grid * eval ()
  {      
    an -> front -> init();

    int i = 0,j = 0;
    Node<int> * p = lobj.Begin();
    while (p != NULL)
      {
	if (p -> obj == withText)
	  lexp[i++] -> eval();
	else
	  lb[j++].eval();  

	p = p -> Next;
      }

    an -> front -> step = 1; 
    frontiere & f(*an -> front);
    assert(f.initboundingbox);
    f.epsilon = Max(f.xmax-f.xmin,f.ymax-f.ymin) * 1.0e-5;                
    p = lobj.Begin(); i = 0; j = 0;
    while (p != NULL)
      {
	if (p -> obj == withText)
	  lexp[i++] -> eval();
	else
	  lb[j++].eval();  

	p = p -> Next;
      }

    Grid * g = new Grid();   
    g -> buildit(an -> front, *an -> wait -> storage);
    return g;
  }
};

Instr * Analyseur::InitScilabMesh(char * s,double * vt,int nbvt, 
				  int * tr, int nbtr)
{
  Iden * id = new Iden(s);
  id -> newVar();
  id -> type = Iden::maillage;
  id -> fg = NULL;

  Instr * res = new MeshCode(id,gen_scilabmesh(s,vt,nbvt,tr,nbtr,an),
			     &an);
  return res;
}

Instr * Analyseur::BuildScilabMesh(char * nameMesh)
{
  GestChar buf(nameMesh);
  buf = buf + "0";

  Iden * id = table.find(buf.Data());
  table[nameMesh].newVar();
  table[nameMesh].type = Iden::maillage;
  table[nameMesh].fg = NULL;

  Instr * res = new MeshCode(&table[nameMesh],new ESciMesh(&an),&an);
  curMesh = &table[nameMesh];

  return res;
}


void Analyseur::InitBorder_1(char * nameBorder, double * tBorder, 
			     int nbPoint, int normal)
{
  GestChar buf(nameBorder);
  buf = buf + "0";

  Iden * id = table.find(buf.Data());
  table[nameBorder].newVar();
  table[nameBorder].type = Iden::courbe;

  *an.ng -> storage = (float) (an.bdyLabel);
  table[nameBorder].fn = (void *) new IB(0,0,0,0,NULL,NULL,rien,
					 an.bdyLabel);

  an.front = front;
  Border b(&an,tBorder,nbPoint,normal,an.bdyLabel);
  if (lb.Insert(b)) 
    throw ErrorMemory("Not enough memory to allocate the border data");
  
  an.bdyLabel++;
}

void Analyseur::InitBorder_2(int ng)
{
  nextSym();
  Iden * id = curIden;
  id = curIden;
  match(iden);
  id -> newVar();
  id -> type = Iden::courbe;
  match(lpar);
  Iden * t = curIden; match(iden); match(becomes);
  if (t -> type == Iden::inconnu) t -> newVar();
  Expr * start = expression(); match(comma);
  Expr * stop = expression(); 
  match(rpar);

  *an.ng -> storage = (float)(an.bdyLabel); 
  id -> fn = (void *) new IB(an.x -> storage,an.y -> storage,
			     an.ng -> storage, t -> storage,start, 
			     stop, instruction(),an.bdyLabel);

  an.front = front;
  Expr * e = new EC(ng);
  Expr * res = new EB(&an,(IB *) id -> fn,e);
  if (lexp.Insert(res)) 
    throw ErrorMemory("Not enough memory to allocate the border data");
  
  an.bdyLabel++;
}

void C2F(delete_lb_lobj)()
{
  lb.DeleteAll();
  lobj.DeleteAll();
}

void C2F(delete_lexp)()
{
  lexp.DeleteAll();
}

void C2F(delete_ptr)(void * ptr)
{
  delete [] ptr;
}

void C2F(delete_scilabana)(int * ierr = 0)
{
  if (scilabana != NULL)
    {
      if (!ierr) ierr = new int;
      *ierr = 0;
      delete scilabana;
      scilabana = 0;
      C2F(delete_lb_lobj)();
      C2F(delete_lexp)();
    }
  else
    *ierr = 1;
}

void C2F(add_lobj)()
{
  if (lobj.Insert(withText))
    throw ErrorMemory("Not enough memory to allocate the border data (Scilink.h : add_lobj)");
}

void C2F(freefem_code)(char * Message, int * err)
{
  
  istrstream is(Message);
  *err = 0;
  toScilab = true;

  try 
    { 
      if (scilabana == NULL) 
	{
	  scilabana = new Analyseur(&is);
	  if (scilabana == NULL)
	    throw ErrorMemory("(Scilink.h : freefemcode)");
	}
      else 
	scilabana -> setBuffer(&is);
      
      cout << "-- Compile --" << endl;
      Instr * it = scilabana -> compile(); 
      cout << "-- Execute --" << endl;
      if (it) it -> execute();
    }
  catch(Error & e)
    { 
      GestChar Text(e.getErrorMessage());
      sprintf(Message,Text.Data());
      
      e.WriteError();
      *err = e.CodeError();

      C2F(delete_scilabana)();
    }
}

void C2F(ff_problem)(char * Message, int * flag, int * err)
{
  getMatProfil = *flag;
  C2F(freefem_code)(Message,err);
}

void C2F(get_ff_result)(char * s,int * lhs,
			double ** tNode, int * mNode, int * nNode, 
			double ** tFunc, int * mFunc, int * nFunc, 
			int ** tTriangle, int * mTriangle, 
			int * nTriangle, char * Message,int * err)
{
  if (scilabana == NULL) 
    {
      *err = 1;
      sprintf(Message,"None data has been capted !!!\r\n");
      return ;
    }

  *err = 0;

  if (*lhs != 2) //Case where we want the result
    {
      GestChar t("plot(");
      t = t + s + ");";
      istrstream is(t.Data());

      toScilab = true;
      scilabana -> setBuffer(&is);

      try 
	{ 
	  Instr * it;
	  cout << "-- Compile --" << endl;
	  it = scilabana -> compile(); 
	  cout << "-- Execute --" <<endl;
	  it -> execute();
	}
      catch(Error & e)
	{ 
	  GestChar Text(e.getErrorMessage());
	  sprintf(Message,Text.Data());

	  e.WriteError();
	  *err = e.CodeError();

	  C2F(delete_scilabana)();
	  
	  return ;
	}
    }

  Grid * gr;
  if (*lhs != 2) 
    gr = dts -> g; //Specified mesh
  else
    gr = (scilabana -> GetAnalvarData()).activeMesh; //Current mesh;

  if (*lhs != 1) *tNode = new double [gr -> nv * 3];
  if (*lhs != 2) *tFunc = new double [gr -> nv];
  if (*lhs != 1) *tTriangle = new int [gr -> nt * 5];
  int k = 0;

  // Remplissage de tNode
  if (*lhs != 1) // Si lhs == 1, on ressort que le resultat 'tFunc'
    {
      *mNode = gr -> nv; *nNode = 3;
      for (int j = 0; j < 3; j++) 
	for (int i = 0; i < gr -> nv; i++)
	  switch (j) 
	    {
	    case 0 : (*tNode)[k++] = gr -> no(&gr -> v[i]); break;
	    case 1 : (*tNode)[k++] = gr -> v[i].x; break;
	    case 2 : (*tNode)[k++] = gr -> v[i].y; break;
	    }
    }
  

  // Remplissage de tFunc
  if (*lhs != 2) // Si lhs == 2, on ressort les donnees sur le 
    {            // maillage, ie 'tNode' et 'tTriangle'
      *mFunc = gr -> nv; *nFunc = 1;
      for (int i = 0; i < gr -> nv; i++)
	(*tFunc)[i] = dts -> an_eval[i];

    }


  // Remplissage de tTriangle
  if (*lhs != 1)
    {
      k = 0;
      *mTriangle = gr -> nt; *nTriangle = 5;
      for (int j = 0; j < 5; j++)
	for (int i = 0; i < gr -> nt; i++)
	  switch (j)
	    {
	    case 0 : (*tTriangle)[k++] = i; break;
	    case 1 : (*tTriangle)[k++] = gr -> no(gr -> t[i].v[0]) + 1; break;
	    case 2 : (*tTriangle)[k++] = gr -> no(gr -> t[i].v[1]) + 1; break;
	    case 3 : (*tTriangle)[k++] = gr -> no(gr -> t[i].v[2]) + 1; break;
	    case 4 : (*tTriangle)[k++] = 0; break;
	    }
    }

  if (*lhs != 2) delete dts;
}

void C2F(put_scilab_mesh)(char * s,double * tNode, int * mNode,
			  int * tTriangle, int * mTriangle,
			  char * Message, int * err)
{ 
  *err = 0;

  try
    {
      if (scilabana == NULL) 
	{
	  scilabana = new Analyseur();
	  if (scilabana == NULL)
	    throw ErrorMemory("(Scilink.h : put_scilab_mesh)");
	}

      Instr * it = scilabana -> InitScilabMesh(s,tNode,*mNode,
					       tTriangle,*mTriangle);
      it -> execute();
    }
  catch(Error & e)
    {
      GestChar Text(e.getErrorMessage());
      sprintf(Message,Text.Data());

      e.WriteError();
      *err = e.CodeError();

      C2F(delete_scilabana)();
      
      return ;
    }
}

void C2F(put_scilab_border_1)(char * nameBorder,double * tBorder,
			      int * nbPoint,int * normal,
			      char * Message,int * err)
{
  try
    {
      if (scilabana == NULL) 
	{
	  scilabana = new Analyseur();
	  if (scilabana == NULL)
	    throw ErrorMemory("(Scilink.h : put_scilab_mesh)");
	}
      if (lobj.Insert(withTab))
	throw ErrorMemory("Not enough memory to allocate the border data (Scilink.h : put_scilab_border_1)");

      scilabana -> InitBorder_1(nameBorder,tBorder,*nbPoint,*normal);
    }
  catch(Error & e)
    {
      GestChar Text(e.getErrorMessage());
      sprintf(Message,Text.Data());

      e.WriteError();
      *err = e.CodeError();

      C2F(delete_scilabana)();
      
      return ;
    }
}

void C2F(put_scilab_border_2)(char * nameBorder,char * equation,
			      double * start, double * stop,int * ng, 
			      int * nrm, char * Message, int * err)
{
  char t1 [10], t2 [10];
  sprintf(t1,"%lf",*start);
  sprintf(t2,"%lf",*stop);
  GestChar t(nameBorder);
  t = t + " (t=" + t1 + "," + t2 + + ") {" + equation + "}";
  istrstream is(t.Data());

  try
    {
      C2F(add_lobj)();

      if (scilabana == NULL) 
	{
	  scilabana = new Analyseur(&is);
	  if (scilabana == NULL)
	    throw ErrorMemory("(Scilink.h : put_scilab_border_2)");
	}
      else 
	scilabana -> setBuffer(&is);

      scilabana -> InitBorder_2(*ng * (*nrm));
    }
  catch(Error & e)
    {
      GestChar Text(e.getErrorMessage());
      sprintf(Message,Text.Data());

      e.WriteError();
      *err = e.CodeError();

      C2F(delete_scilabana)();
      
      return ;
    }
}

void C2F(build_scilab_mesh)(char * nameMesh,char * Message,int * err)
{
  *err = 0;
  toScilab = true;

  try
    {
      cout << "\n-- Build Scilab mesh --" << endl;
      Instr * it = scilabana -> BuildScilabMesh(nameMesh);
      it -> execute();
    }
  catch(Error & e)
    {
      GestChar Text(e.getErrorMessage());
      sprintf(Message,Text.Data());

      e.WriteError();
      *err = e.CodeError();

      C2F(delete_scilabana)();
	
      return ;
    }
}

void C2F(get_matrix)(double ** tMatrix,int * sztMatrix,int ** tjlow, 
		     int ** tjhigh, int * sztVect, int ** szblock,
		     char * Message, int * err)
{
  try 
    {
      if (__gmp == NULL) 
	throw ErrorExec("Any work has been done to get the matrix !");

      *sztMatrix = __gmp -> SizeMatrix;
      *sztVect = __gmp -> SizeVect;
      *szblock = new int;
      **szblock = __gmp -> SizeBlock;
      int dim = *sztMatrix * (**szblock) * (**szblock);
      
      *tMatrix = new double [dim];
      if (!*tMatrix) throw ErrorMemory("(scilink.h : get_matrix)");
      *tjlow = new int [*sztVect];
      if (!*tjlow) throw ErrorMemory("(scilink.h : get_matrix)");
      *tjhigh = new int [*sztVect];
      if (!*tjhigh) throw ErrorMemory("(scilink.h : get_matrix)");

      for (int i = 0; i < dim ; i++)
	(*tMatrix)[i] = __gmp -> Matrix[i];
      for (int i = 0; i < (*sztVect); i++)
	{
	  (*tjlow)[i] = __gmp -> jlow[i];
	  (*tjhigh)[i] = __gmp -> jhigh[i];
	}

      delete __gmp;
      __gmp = NULL;
    }
  catch(Error & e)
    {
      GestChar Text(e.getErrorMessage());
      sprintf(Message,Text.Data());
      
      e.WriteError();
      *err = e.CodeError();

      if (__gmp) { delete __gmp; __gmp = NULL; }
      C2F(delete_scilabana)();
    }

  getMatProfil = false;
}

#endif
