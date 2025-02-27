#include "utilities.h"
#include <cassert>


int modFactorial(int n, int modulus) {
    if (n == 0 || n == 1) { return 1; } 
    else { return ( (n * modFactorial(n - 1, modulus)) % modulus ); }
}

int gcdExtended(int a, int b, int* x, int* y)
{
    if (a == 0) { *x = 0, *y = 1; return b; }
    int x1, y1; int gcd = gcdExtended(b % a, a, &x1, &y1);
    *x = y1 - (b / a) * x1; *y = x1;
    return gcd;
}

int modInverse(int A, int M)
{
    int x, y; int g = gcdExtended(A, M, &x, &y);
    assert(g = 1); // Otherwise, no inverse
    return (x % M + M) % M; 
}


std::vector<int64_t> repFillSlots(std::vector<int64_t> vecIn, int maxSlots)
{
    int n = vecIn.size();
    std::vector<int64_t> resVec(maxSlots,0);
    int halfSlots = std::floor(maxSlots/2);
    int repNum = std::floor(halfSlots/n);
    for (int slot = 0; slot < repNum*n; slot++){
        resVec[slot] = vecIn[slot%n];
        resVec[maxSlots-repNum*n+slot] = vecIn[slot%n];
    }
    return resVec;
}

std::vector<std::vector<int64_t>> matrixDiagonals(std::vector<std::vector<int64_t>> matIn) 
{
    int d = matIn.size();
    std::vector<std::vector<int64_t>> diagonals;
    for (int l=0;l<d;l++){
        std::vector<int64_t> diagonal(d,0);
        for (int i=0;i<d;i++){diagonal[i] =  matIn[i][(l+i)%d];}
        diagonals.push_back(diagonal);
    }
    return diagonals;
}

CryptoOpsLogger::CryptoOpsLogger() : innerProdOps_(0), innerProdTime_(0.0), 
                                    //  addManyOps_(0), addManyTime_(0.0),
                                     multOps_(0), multTime_(0.0), 
                                     addOps_(0), addTime_(0.0),
                                     rotOps_(0), rotTime_(0.0) {}

void CryptoOpsLogger::logInnerProd(double ms) { innerProdOps_ += 1 ; innerProdTime_ = innerProdTime_ + ms; };
void CryptoOpsLogger::logAddMany(int n, double ms) { addOps_ += n-1 ; addTime_ = addTime_ + ms; };
void CryptoOpsLogger::logMult(double ms) { multOps_ += 1 ; multTime_ = multTime_ + ms; };
void CryptoOpsLogger::logAdd(double ms) { addOps_ += 1 ; addTime_ = addTime_ + ms; };
void CryptoOpsLogger::logRot(double ms) { rotOps_ += 1 ; rotTime_ = rotTime_ + ms; };

int CryptoOpsLogger::innerProdOps() { return innerProdOps_; };  
double CryptoOpsLogger::innerProdTime() { return innerProdTime_; };
// int CryptoOpsLogger::addManyOps() { return addManyOps_; };      
// double CryptoOpsLogger::addManyTime() { return addManyTime_; };
int CryptoOpsLogger::multOps() { return multOps_; };            
double CryptoOpsLogger::multTime() { return multTime_; };
int CryptoOpsLogger::addOps() { return addOps_; };              
double CryptoOpsLogger::addTime() { return addTime_; };
int CryptoOpsLogger::rotOps() { return rotOps_; };              
double CryptoOpsLogger::rotTime() { return rotTime_; };
int CryptoOpsLogger::totalTime() { return innerProdTime_+multTime_+addTime_+rotTime_; }; 


VectorIter::VectorIter(int modulus, int slots) {
    mod_ = modulus;
    for (int i = 0; i < slots; ++i) { vectorIter_.push_back(0); }
};
std::vector<int> VectorIter::value() {
    return vectorIter_;
};
bool VectorIter::iterate() {
    return iterate_(0);
};
bool VectorIter::iterate_(int slot) {
    // Case 0: Highest idx cannot be incremented further.
    if (slot == int(vectorIter_.size())-1 && (vectorIter_[slot] == mod_-1)) { return false; }
    // Case 1: Increment & carry. 
    else if (vectorIter_[slot] == mod_ - 1) {
        vectorIter_[slot] = 0; if (iterate_(slot+1)) { return true; } else { return false; }}
    // Case 2: Increment.
    else { vectorIter_[slot] = vectorIter_[slot]+1; return true; }
};

