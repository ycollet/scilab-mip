#ifndef LIST_H
#define LIST_H

template <typename Object>
class Node
{
 public:
  Object obj;
  Node<Object> * Next;
  Node() {}
};

template <typename Object>
class List
{
  Node<Object> * Cur;
  Node<Object> * First;
 public:
 List() : Cur(new Node<Object>),First(new Node<Object>) { First = NULL; }
  ~List();
  Node<Object> * Begin() const { return First; }
  Node<Object> * Current() const { return Cur; }
  
  bool Empty() { return (First == NULL); }
  int Insert( const Object & );
  bool Delete( const Object & );
  void DeleteAll();
  void Display();
  
  const Object & operator [] (int i);
};
#endif
