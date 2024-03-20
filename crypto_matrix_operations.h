#ifndef CRYPTO_MATRIX_OPERATIONS_H
#define CRYPTO_MATRIX_OPERATIONS_H

#include "openfhe.h"
#include "utilities.h"
#include "crypto_utilities.h"
#include "crypto_enc_transform.h"

using namespace lbcrypto;



std::vector<std::vector<Ciphertext<DCRTPoly>>> // Element-wise-encrypted output matrix.
    evalMatrixMul2Pow(std::vector<std::vector<std::vector<Ciphertext<DCRTPoly>>>> &encMatsElems, // Elem-wise encrypted input matrices.
                      CryptoContext<DCRTPoly> &cryptoContext,
                      InitRotsMasks &initRotsMasks,
                      CryptoOpsLogger &cryptoOpsLogger);


std::vector<std::vector<std::vector<Ciphertext<DCRTPoly>>>> // Element-wise-encrypted output matrix.
    evalMatSquarings(std::vector<Ciphertext<DCRTPoly>> &encMatRows, // Row-wise-encrypted input matrix.
                     int sqs,
                     CryptoContext<DCRTPoly> &cryptoContext,
                     InitRotsMasks &initRotsMasks,
                     CryptoOpsLogger &cryptoOpsLogger,
                     KeyPair<DCRTPoly> &keyPair);

// Fast, but incurs log(n) additional mul depth. 
std::vector<std::vector<Ciphertext<DCRTPoly>>> 
        evalMatSqMul(std::vector<Ciphertext<DCRTPoly>> &encRows, 
                     int exponent,
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
