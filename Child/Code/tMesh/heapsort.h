//-*-c++-*- 

/***************************************************************************/
/**
**  @file
**  @brief a simple heap sort
*/
/***************************************************************************/

template<class TYPE>
void shunt(int ,int ,TYPE const& ,TYPE *);
template<class TYPE>
void heapsort(int ,TYPE *);

template<class TYPE>
void shunt(int k,int m,TYPE const& temp,TYPE array[]){
  /*Heap building and re-ordering function for heapsort 
    Mike Bithell 03/10/00
  */
  int i=k;
  int j=k+k;
  while (j <= m){
    if (j < m){
      if(array[j-1] < array[j])j=j+1;
    }
    if (temp < array[j-1]){
      array[i-1]=array[j-1];
      i=j;
      j=j+j;
    }
    else{
      j=m+1;
    }
  }
  array[i-1]=temp;
}
template<class TYPE>
void heapsort(int n,TYPE array[]){
  /* A function that sorts the array, of dimension n, into increasing
     order using a heap sort (cribbed and adapted from numerical recipes)
     calls shunt to reorder the array elements. First set of calls orders
     the array into a heap and second set of calls does the sort.
     Performance should be O(NlogN) even for initially sorted arrays, but
     a little slower than quicksort.
     Now altered to include a template - note that class TYPE must have the 
     < operator defined.
     Mike Bithell 03/10/00
  */

  if (n<2) return;

  TYPE temp;

  {
  for (int k=n/2;k>=1;k--){
    temp=array[k-1];
    shunt(k,n,temp,array);
  }
  }
  {
  for(int k=n;k>1;k--) {
    temp=array[k-1];
    array[k-1]=array[0];
    shunt(1,k-1,temp,array);
  }
  }
  array[0]=temp;
}

#if defined(TEST_)

#include <iostream>

int main(){
  /*Short test program for heapsort*/
  int n=2;
  int a[n];

  for (int i=0;i<n;i++)a[i]=(4-i);
  std::cout << "Original array" << std::endl;
  for (int i=0;i<n;i++)std::cout << a[i] << ' ';
  std::cout << std::endl;
  heapsort(n,a);
  std::cout << "Sorted!" << std::endl;
  for (int i=0;i<n;i++)std::cout << a[i] << ' ';
  std::cout << '\n';
}

#endif
