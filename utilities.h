#ifndef UTILITIES_H
#define UTILITIES_H

int modFactorial(int n, int modulus);
int gcdExtended(int a, int b, int* x, int* y);
int modInverse(int A, int M);

class CryptoOpsLogger {
public:
    CryptoOpsLogger();
    void logInnerProd(double ms);
    void logAddMany(double ms);
    void logMult(double ms);
    void logAdd(double ms);
    void logRot(double ms);

    int innerProdOps();     double innerProdTime();
    int addManyOps();       double addManyTime();
    int multOps();          double multTime();
    int addOps();           double addTime();
    int rotOps();           double rotTime();
    int totalTime();
private:
    int innerProdOps_; double innerProdTime_;
    int addManyOps_; double addManyTime_;
    int multOps_; double multTime_;
    int addOps_; double addTime_;
    int rotOps_; double rotTime_; // TODO: Rotation degree not logged.
};



#endif

