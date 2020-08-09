#include <iostream>

using namespace std;

#include <assert.h>
#include <stdlib.h>

#include "list.h"

template <typename Object> 
List<Object>::~List()
{
  DeleteAll();
}

template <typename Object>
int List<Object>::Insert(const Object & obj)
{
  Node<Object> * p = new Node<Object>;
  if (p == NULL) return -1;
  
  if (First == NULL)
    {
      p -> obj = obj;
      p -> Next = NULL;
      First = p;
      Cur = First;
    }
  else
    {
      p -> Next = NULL;
      p -> obj = obj;
      Cur -> Next = p;
      Cur = Cur -> Next;
    }

  return 0;
}

template <typename Object>
bool List<Object>::Delete( const Object & obj )
{
  if (Empty()) return false;

  bool Find = false;
  Node<Object> * p = First;

  if (First -> obj == obj)
    {
      if (First -> Next != NULL)
	First = First -> Next;
      else
	First = NULL;
      delete p;
      Find = true;
    }
  else
    while (!Find && ((p -> Next) != NULL))
      {
	if (((p -> Next) -> obj) == obj) 
	  {
	    Find = true;
	    Node<Object> * q = p -> Next;
	    p -> Next = q -> Next;
	    delete q;
	  }
	else
	  p = p -> Next;
      }

  return Find;
}

template <typename Object>
void List<Object>::DeleteAll()
{
  if (Empty()) return;

  Node<Object> * pt = First, * p = pt;
  
  while (p != NULL)
    {
      p = pt -> Next;
      delete pt;
      pt = p;
    }

  Cur = First = NULL;
}

template <typename Object>
void List<Object>::Display()
{
  Node<Object> * p = First;
  while (p != NULL)
    {
      cout << p -> obj;
      p = p -> Next;
    }
}


template <typename Object>
const Object & List<Object>::operator [] ( int i )
{
  assert( i >= 0 );

  Node<Object> * p = First;
  int j = 0;
  
  while (i != j) 
    {
      p = p -> Next;
      if (p == NULL) 
	{
	  cerr << "Debordement d'indice List::operator []" << endl;
	  exit(0);
	}
      j++;
    }

  return p -> obj;
}
