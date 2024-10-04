#ifndef CRYPTO_MATRIX_OPERATIONS_H
#define CRYPTO_MATRIX_OPERATIONS_H

#include "openfhe.h"
#include "utilities.h"
#include "crypto_utilities.h"
#include "crypto_enc_transform.h"
#include "crypto_prefix_mult.h"
#include <map>
#include <omp.h>

using namespace lbcrypto;

Ciphertext<DCRTPoly> evalDiagMatrixVecMult(std::vector<Ciphertext<DCRTPoly>> &encMatDiagonals, // Must be filled.
                                           Ciphertext<DCRTPoly> encVec,                        // Must be filled.
                                           CryptoContext<DCRTPoly> &cryptoContext);


// Class initializes rotation keys and encrypted masks.
class InitMatrixMult {
public:
    InitMatrixMult(CryptoContext<DCRTPoly> &cryptoContext, KeyPair<DCRTPoly> keyPair, int d);
    std::map<int, Ciphertext<DCRTPoly>> u_sigma();
    std::map<int, Ciphertext<DCRTPoly>> u_tau();
    std::map<int, Ciphertext<DCRTPoly>> v1();
    std::map<int, Ciphertext<DCRTPoly>> v2();
    Ciphertext<DCRTPoly> matrixMask();
    const int d;
private:
    std::map<int, Ciphertext<DCRTPoly>> _u_sigma;
    std::map<int, Ciphertext<DCRTPoly>> _u_tau;
    std::map<int, Ciphertext<DCRTPoly>> _v1;
    std::map<int, Ciphertext<DCRTPoly>> _v2;
    Ciphertext<DCRTPoly> _matrixMask;
};


Ciphertext<DCRTPoly> evalMatrixMult(CryptoContext<DCRTPoly> &cryptoContext,
                                    Ciphertext<DCRTPoly> encA,
                                    Ciphertext<DCRTPoly> encB,
                                    InitMatrixMult &initMatrixMult);

Ciphertext<DCRTPoly> evalMatrixMultParallel(CryptoContext<DCRTPoly> &cryptoContext,
                                            Ciphertext<DCRTPoly> encA,
                                            Ciphertext<DCRTPoly> encB,
                                            InitMatrixMult &initMatrixMult);

#endif
