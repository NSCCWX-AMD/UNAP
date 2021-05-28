#include <stdio.h>
#include <stdlib.h>

#include <algorithm>  //- to use find_if
#include <iostream>
#include <vector>

//- test for vector and objects comparison

using namespace std;

const int N = 11;
int loc1[N] = {0, 0, 1, 2, 2, 3, 4, 5, 5, 6, 7};
int loc2[N] = {0, 1, 1, 2, 2, 3, 3, 4, 5, 6, 6};

class labelPair
{
private:
  int first_;
  int second_;
  int faceI_;

public:
  labelPair(int i, int j) : first_(i), second_(j), faceI_(-1) {}

  int first() const { return first_; }

  int second() const { return second_; }

  int faceI() const { return faceI_; }

  void faceI(const int i) { faceI_ = i; }

  bool operator==(const labelPair &a) const
  {
    return ((this->first() == a.first()) && (this->second() == a.second()));
  }
};

int main()
{
  vector<labelPair> cellsToCoarseFace;
  vector<int> dynFaceCells;
  vector<int> dynFaceRestrictAddressing;

  for (int i = 0; i < N; i++)
  {
    labelPair cellPair(loc1[i], loc2[i]);
    vector<labelPair>::iterator it =
        find(cellsToCoarseFace.begin(), cellsToCoarseFace.end(), cellPair);

    if (it == cellsToCoarseFace.end())
    {
      //- new face
      int coarseI = dynFaceCells.size();
      dynFaceRestrictAddressing.push_back(coarseI);
      dynFaceCells.push_back(loc1[i]);
      cellPair.faceI(coarseI);
      cellsToCoarseFace.push_back(cellPair);
    }
    else
    {
      dynFaceRestrictAddressing.push_back(it->faceI());
    }
  }

  cout << "dynFaceCells: " << endl;
  for (int i = 0; i < dynFaceCells.size(); i++)
  {
    cout << "At i = " << i << ", dynFaceCells = " << dynFaceCells[i] << endl;
  }

  cout << "dynFaceRestrictAddressing: " << endl;
  for (int i = 0; i < dynFaceRestrictAddressing.size(); i++)
  {
    cout << "At i = " << i
         << ", dynFaceRestrictAddressing = " << dynFaceRestrictAddressing[i]
         << endl;
  }

  return 0;
}