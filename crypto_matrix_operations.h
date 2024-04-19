#ifndef CRYPTO_MATRIX_OPERATIONS_H
#define CRYPTO_MATRIX_OPERATIONS_H

#include "openfhe.h"
#include "utilities.h"
#include "crypto_utilities.h"
#include "crypto_enc_transform.h"
#include "crypto_prefix_mult.h"
#include <map>

using namespace lbcrypto;

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
                                    InitMatrixMult &initMatrixMult,
                                    KeyPair<DCRTPoly> &keyPair); // TODO: keypair for debugging.

std::vector<std::vector<Ciphertext<DCRTPoly>>> // Element-wise-encrypted output matrix.
    evalMatrixMul2Pow(std::vector<std::vector<std::vector<Ciphertext<DCRTPoly>>>> &encMatsElems, // Elem-wise encrypted input matrices.
                      int packingMode, // 0: row/col | 1: full matrix
                      CryptoContext<DCRTPoly> &cryptoContext,
                      InitRotsMasks &initRotsMasks,
                      CryptoOpsLogger &cryptoOpsLogger);

std::vector<std::vector<std::vector<Ciphertext<DCRTPoly>>>> // Element-wise-encrypted output matrix.
    evalMatSquarings(std::vector<std::vector<Ciphertext<DCRTPoly>>> &encMatElems, // Row-wise-encrypted input matrix.
                     int sqs,
                     int packingMode, // 0: row/col | 1: full matrix
                     CryptoContext<DCRTPoly> &cryptoContext,
                     InitRotsMasks &initRotsMasks,
                     CryptoOpsLogger &cryptoOpsLogger,
                     KeyPair<DCRTPoly> &keyPair);

// Fast, but incurs log(n) additional mul depth. 
std::vector<std::vector<Ciphertext<DCRTPoly>>> 
        evalMatSqMul(std::vector<std::vector<Ciphertext<DCRTPoly>>> &encMatElems, 
                     int exponent,
                     int packingMode, // 0: row/col | 1: full matrix
                     CryptoContext<DCRTPoly> &cryptoContext,
                     InitRotsMasks &initRotsMasks,
                     CryptoOpsLogger &cryptoOpsLogger,
                     KeyPair<DCRTPoly> &keyPair);

// Very slow but optimal mul depth. For example, M^4(r=row, c=col) = sum_{i,j,k} M_ri x M_ij x M_jk x M_kc
std::vector<Ciphertext<DCRTPoly>> evalMatrixExp(std::vector<Ciphertext<DCRTPoly>> &encRows, 
                                                int exponent,
                                                CryptoContext<DCRTPoly> &cryptoContext,
                                                InitRotsMasks &initRotsMasks,
                                                CryptoOpsLogger &CryptoOpsLogger); 

Ciphertext<DCRTPoly> evalDiagMatrixVecMult(std::vector<Ciphertext<DCRTPoly>> &encMatDiagonals, // Must be filled.
                                           Ciphertext<DCRTPoly> encVec,                        // Must be filled.
                                           CryptoContext<DCRTPoly> &cryptoContext,            
                                           CryptoOpsLogger &cryptoOpsLogger); 

Ciphertext<DCRTPoly> evalMatrixVecMult(std::vector<Ciphertext<DCRTPoly>> &encRows, 
                                       Ciphertext<DCRTPoly> &enc_vec,
                                       CryptoContext<DCRTPoly> &cryptoContext,            
                                       InitRotsMasks &initRotsMasks);


Ciphertext<DCRTPoly> evalVecMatrixMult(Ciphertext<DCRTPoly> &enc_vec,
                                       std::vector<Ciphertext<DCRTPoly>> &encRows, 
                                       CryptoContext<DCRTPoly> &cryptoContext,            
                                       InitRotsMasks &initRotsMasks,
                                       CryptoOpsLogger &cryptoOpsLogger);




#endif
