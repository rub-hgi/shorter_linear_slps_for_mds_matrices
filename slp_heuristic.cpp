/* November 2017
This code is a modified version of code available at

  http://www.imada.sdu.dk/~joan/xor/Improved2.cc

The only differences are:
  - Minor changes in formatting.
  - Input is now read from stdin instead of from a fixed file.
  - The code now stores and keeps track of a string representation of the SLP
    in the variable called Program.
  - The SLP is printed at the end from this variable.
*/


// Improved.cc
// version August 31, 2008

/* The idea of the heuristic is to
  - keep a "base" of available signals (initially the
    base is just the set of variables x0, ..., xn);
  - for each required matrix output (I called them
    "Targets") keep a "distance" from the base to the
    output. e.g. Dist[3]+1 is the smallest number of base
    elements that I need to sum in order to obtain the
    third row of the matrix;
  - greedily pick a new basis element by adding two existing
    basis elements;

  The current criteria for picking the new basis element is
   - if a target is the sum of two basis elements, pick those
   - otherwise, pick the one that minimizes the sum of new distances
   - resolve ties by maximizing the euclidean norm of the
   vector of new distances;

  The tie resolution criteria is kind of counter-intuitive. I haven't
  thought much about it, but the genesis was that I preferred a distance
  vector like 0,0,3,1 to one like 1,1,1,1. In the latter case, it seemed
  like I would need 4 more gates to finish. In the former, 3 or 2 might
  do it.

*/
/* On Dec. 06, 2007, this was run on the matrix describing the
bottom and top linear transformations in our new AES S-box circuit
(maple program on file) . The matrices are r (bottom) =

0 0 0 1 1 0 1 1 0 1 1 0 0 0 0 1 1 0
1 1 0 0 0 0 1 1 0 1 1 0 0 0 0 1 1 0
1 0 1 0 0 0 1 0 1 0 0 0 1 0 1 1 0 1
1 1 0 1 1 0 0 0 0 1 1 0 0 0 0 1 1 0
0 1 1 0 1 1 0 0 0 1 1 0 0 0 0 1 1 0
1 0 1 1 1 0 0 1 1 0 1 1 1 0 1 1 1 0
0 0 0 0 1 1 0 1 1 0 0 0 1 1 0 1 1 0
1 0 1 1 0 1 0 0 0 0 0 0 1 1 0 1 1 0

The program output a circuit with 30 XORs.
*/

#include <math.h>
#include <ctype.h>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <string>

using namespace std;

const int MaxBaseSize=1000;
const bool PRINTROWS=true;

int NumInputs;
int NumTargets;
int XorCount;
long long int Target[MaxBaseSize];
int Dist[MaxBaseSize]; //distance from current base to Target[i]
int NDist[MaxBaseSize]; //what Dist would be if NewBase was added
long long int Base[MaxBaseSize];
string Program[MaxBaseSize];
int BaseSize;
int TargetsFound;

void InitBase();
void ReadTargetMatrix();
bool is_target(long long int x);
bool is_base(long long int x);
int NewDistance(int u); //calculates the distance from the base to Target[u]
int TotalDistance(); //returns the sum of distances to targets
bool reachable(long long int T, int K, int S);
bool EasyMove(); //if any two bases add up to a target, pick them
void PickNewBaseElement();
void binprint(long long int x); //outputs last NumInputs bits of x

int main(int argc, char *argv[]) {
    int NumMatrices;

    cin >> NumMatrices;
    for (int i = 0; i < NumMatrices; i++) {
        ReadTargetMatrix();
        InitBase();
        XorCount = 0;
        int counter = 0;
        while (TargetsFound < NumTargets) {
            counter++;
            if (!EasyMove()) PickNewBaseElement();
        }
        cout << "SLP Heuristic XorCount: " << XorCount << endl;
        cout << "SLP: " << endl << endl;
        for (int j = 0; j < XorCount; j++) {
            cout << Program[NumInputs + j] << endl;
        }
    }

    return 0;
}

void InitBase() {
    TargetsFound = 0;
    Base[0] = 1;
    Program[0] = "x0";
    for (int i = 1; i < NumInputs; i++) {
        Base[i] = 2*Base[i-1];
        Program[i] = "x" + to_string(i);
    }
    BaseSize = NumInputs; //initial base is just the xi's
    for (int i = 0; i < NumTargets; i++) {
        if (Dist[i] == 0) {
            TargetsFound++;
        }
    }
}

int TotalDistance() { //returns the sum of distances to targets
    int D = 0;
    int t;
    for (int i = 0; i < NumTargets; i++) {
        t = NewDistance(i);
        NDist[i] = t;
        D = D + t;
    }
    return D;
}

long long int NewBase; //global variable containing a candidate new base

bool EasyMove() {
    int t;
    bool foundone = false;

    //see if anything in the distance vector is 1
    for(int i = 0; i < NumTargets; i++) {
        if (Dist[i] == 1) {
            foundone = true;
            t = i;
            break;
        }
    }
    if (!foundone) {
        return false;
    }
    //update Dist array
    NewBase = Target[t];
    for (int u = 0; u < NumTargets; u++) {
        Dist[u] = NewDistance(u);
    }
    //update Base with NewBase
    Base[BaseSize] = NewBase;
    //find which lines in Base caused this
    string a,b;
    for (int i = 0; i < BaseSize; i++) {
        for (int j = i + 1; j < BaseSize; j++) {
            if ((Base[i] ^ Base[j]) == Target[t]) {
                a = Program[i].substr(0, Program[i].find(" "));
                b = Program[j].substr(0, Program[j].find(" "));
                break;
            }
        }
    }
    Program[BaseSize] = "y" + to_string(t) + " = " + a + " + " + b;
    BaseSize++;
    XorCount++;
    TargetsFound++;
    return true;
}

// PickNewBaseElement is only called when there are no 1's in Dist[]
void PickNewBaseElement() {
    int MinDistance;
    long long int TheBest;
    int ThisDist;
    int ThisNorm, OldNorm;
    int besti,bestj, d;
    bool easytarget;
    int BestDist[MaxBaseSize];

    MinDistance = BaseSize*NumTargets; //i.e. something big
    OldNorm = 0; //i.e. something small
    //try all pairs of bases
    for (int i = 0; i < BaseSize - 1; i++) {
        for (int j = i+1; j < BaseSize; j++) {
            NewBase = Base[i] ^ Base[j];
            //sanity check
            if (NewBase == 0) { cout << "a base is 0, should't happen " << endl; exit(0); }
            //if NewBase is not new continue
            if (is_base(NewBase)) continue;
            //if NewBase is target then choose it
            easytarget = false;
            if (is_target(NewBase)) {
                cout << "shouldn't find an easy target here " << endl;
                exit(0);
                easytarget = true;
                besti = i;
                bestj = j;
                TheBest = NewBase;
                break;
            }
            ThisDist = TotalDistance(); //this also calculates NDist[]
            if (ThisDist <= MinDistance) {
                //calculate Norm
                ThisNorm = 0;
                for (int k = 0; k < NumTargets; k++) {
                    d = NDist[k];
                    ThisNorm = ThisNorm + d*d;
                }
                //resolve tie in favor of largest norm
                if ((ThisDist < MinDistance) || (ThisNorm > OldNorm) ) {
                    besti = i;
                    bestj = j;
                    TheBest = NewBase;
                    for (int uu = 0; uu < NumTargets; uu++) {
                        BestDist[uu] = NDist[uu];
                    }
                    MinDistance = ThisDist;
                    OldNorm = ThisNorm;
                }
            }
        }
        if (easytarget) {
            break;
        }
    }
    //update Dist array
    NewBase = TheBest;
    for (int i = 0; i < NumTargets; i++) {
        Dist[i] = BestDist[i];
    }
    //update Base with TheBest
    Base[BaseSize] = TheBest;
    //update linear program
    string a = Program[besti].substr(0, Program[besti].find(" "));
    string b = Program[bestj].substr(0, Program[bestj].find(" "));
    Program[BaseSize] = "t" + to_string(XorCount) + " = " + a + " + " + b;
    BaseSize++;
    XorCount++;
    if (is_target(TheBest)) {
        TargetsFound++; //this shouldn't happen
    }
}

void binprint(long long int x) { //outputs last NumInputs bits of x
    long long int t = x;
    for (int i = 0; i < NumInputs; i++) {
        if (t%2) {
            cout << "1 ";
        }
        else {
            cout << "0 ";
        }
        t = t/2;
    }
}

void ReadTargetMatrix() {
    cin >> NumTargets;
    cin >> NumInputs;
    //check that NumInputs is < wordsize
    if (NumInputs >= 8*sizeof(long long int)) {
        cout << "too many inputs" << endl;
        exit(0);
    }

    int bit;
    for (int i = 0; i < NumTargets; i++) { //read row i
        long long int PowerOfTwo  = 1;
        Target[i] = 0;
        Dist[i] = -1; //initial distance from Target[i] is Hamming weight - 1
        for (int j = 0; j < NumInputs; j++) {
            cin >> bit;
            if (bit) {
                Dist[i]++;
                Target[i] = Target[i] + PowerOfTwo;
            }
            PowerOfTwo = PowerOfTwo * 2;
        }
    }
}

bool is_target(long long int x)
{
    for (int i = 0; i < NumTargets; i++) {
        if (x == Target[i]) {
            return true;
        }
    }
    return false;
}

bool is_base(long long int x) {
    //sanity check, shouldn't ask if 0 is base
    if (x==0) {
        cout << "asking if 0 is in Base " << endl;
        exit(0);
    }

    for (int i = 0; i < BaseSize; i++) {
        if (x == Base[i]) {
            return true;
        }
    }
    return false;
}

// Distance is 1 less than the number of elements
// in the base that I need to add in order to get Target[u].
// The next function calculates the distance from the base,
// augmented by NewBase, to Target[u]. Uses the following observations:
// Adding to the base can only decrease distance.
// Also, since NewBase is the sum of two old base
// elements, the distance from the augmented base
// to Target[u] can decrease at most by 1. If the
// the distance decreases, then NewBase must be one
// of the summands.

int NewDistance(int u) {
    //if Target[u] is in augmented base return 0;
    if (is_base(Target[u]) || (NewBase == Target[u])) {
        return 0;
    }

    // Try all combinations of Dist[u]-1 base elements until one sums
    // to Target[u] + NewBase. If this is true, then Target[u] is the
    // sum of Dist[u] elements in the augmented base, and therefore
    // the distance decreases by 1.

    if (reachable(Target[u] ^ NewBase, Dist[u] - 1, 0)) {
        return (Dist[u]-1);
    } else {
        return Dist[u]; //keep old distance
    }
}


//return true if T is the sum of K elements among Base[S..BaseSize-1]
bool reachable(long long int T, int K, int S) {
    if ((BaseSize-S) < K) {
        return false; //not enough base elements
    }

    if (K==0) {
        return false; //this is probably not reached
    }

    if (K==1) {
        for (int i=S; i < BaseSize; i++) if (T == Base[i]) {
            return true;
        }
        return false;
    }

    //consider those sums containing Base[S]
    if (reachable(T^Base[S], K-1, S+1)) {
        return true;
    }

    //consider those sums not containing Base[S]
    if (reachable(T, K, S+1)) {
        return true;
    }

    //not found
    return false;
}

