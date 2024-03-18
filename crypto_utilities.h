#ifndef CRYPTO_UTILITIES_H
#define CRYPTO_UTILITIES_H

#include "openfhe.h"

using namespace lbcrypto;

// Helper decrypt and print functions for debugging.
void printEncMatRows(std::vector<Ciphertext<DCRTPoly>> &encMatRows, CryptoContext<DCRTPoly> &cryptoContext, KeyPair<DCRTPoly> keyPair);
void printEncMatElems(std::vector<std::vector<Ciphertext<DCRTPoly>>> &encMatElems, CryptoContext<DCRTPoly> &cryptoContext, KeyPair<DCRTPoly> keyPair);

// Class initializes rotation keys and encrypted masks.
class InitRotsMasks {
public:
    InitRotsMasks(CryptoContext<DCRTPoly> &cryptoContext, KeyPair<DCRTPoly> keyPair, int slots);
    std::vector<Ciphertext<DCRTPoly>> encMasks();
    const int slots;
private:
    std::vector<Ciphertext<DCRTPoly>> encMasks_;
};

#endif

