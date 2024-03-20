#ifndef CRYPTO_ENC_TRANSFORM_H
#define CRYPTO_ENC_TRANSFORM_H

#include "openfhe.h"
#include "utilities.h"
#include "crypto_utilities.h"

using namespace lbcrypto;

// Helper method for matrix exponentiation: transforms row encryptions to encryptions of columns.
// Multiplicative depth: 1
std::vector<Ciphertext<DCRTPoly>> rowToColEnc(std::vector<Ciphertext<DCRTPoly>> &encRows, 
                                              CryptoContext<DCRTPoly> &cryptoContext,
                                              InitRotsMasks &InitRotsMasks,
                                              CryptoOpsLogger &cryptoOpsLogger);

std::vector<Ciphertext<DCRTPoly>> encElem2Rows(std::vector<std::vector<Ciphertext<DCRTPoly>>> &encMatElems,
                                               CryptoContext<DCRTPoly> &cryptoContext,
                                               InitRotsMasks &initRotsMasks,
                                               CryptoOpsLogger &cryptoOpsLogger);

std::vector<Ciphertext<DCRTPoly>> encElem2Cols(std::vector<std::vector<Ciphertext<DCRTPoly>>> &encMatElems,
                                               CryptoContext<DCRTPoly> &cryptoContext,
                                               InitRotsMasks &initRotsMasks,
                                               CryptoOpsLogger &cryptoOpsLogger);


#endif

