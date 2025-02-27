#ifndef CRYPTO_UTILITIES_H
#define CRYPTO_UTILITIES_H

#include "openfhe.h"

using namespace lbcrypto;


// Helper decrypt and print functions for debugging.
void printEnc(Ciphertext<DCRTPoly> &ciphertext, int slots, CryptoContext<DCRTPoly> &cryptoContext, KeyPair<DCRTPoly> keyPair);
void printEncMatRows(std::vector<Ciphertext<DCRTPoly>> &encMatRows, CryptoContext<DCRTPoly> &cryptoContext, KeyPair<DCRTPoly> keyPair);
void printEncMatElems(std::vector<std::vector<Ciphertext<DCRTPoly>>> &encMatElems, CryptoContext<DCRTPoly> &cryptoContext, KeyPair<DCRTPoly> keyPair);

std::vector<std::vector<int64_t>> matrixDiagonals(std::vector<std::vector<int64_t>> matIn); 

// Class initializes rotation keys and encrypted masks.
class InitRotsMasks {
public:
    InitRotsMasks(CryptoContext<DCRTPoly> &cryptoContext, KeyPair<DCRTPoly> keyPair, int slots);
    std::vector<Ciphertext<DCRTPoly>> encMasks();
    std::vector<Ciphertext<DCRTPoly>> encMasksFullyPacked();
    Ciphertext<DCRTPoly> encZeroes();

    const int slots;
private:
    std::vector<Ciphertext<DCRTPoly>> encMasks_;
    std::vector<Ciphertext<DCRTPoly>> encMasksFullyPacked_;
    Ciphertext<DCRTPoly> encZeroes_;
};


// Ciphertext exponentiation, via square and multiply. Multiplicative depth: log(exponent)
Ciphertext<DCRTPoly> evalExponentiate(Ciphertext<DCRTPoly> &ciphertext, int exponent, 
                                      CryptoContext<DCRTPoly> &cryptoContext);


// Decrypt and encrypt to reset ciphertext noise.
void refreshInPlace(Ciphertext<DCRTPoly> &ciphertext, int slots, 
                    KeyPair<DCRTPoly> keyPair, CryptoContext<DCRTPoly> &cryptoContext);

// Decrypt and return encryptions of slot elements individually (in the first slot position).
std::vector<Ciphertext<DCRTPoly>> refreshElems(Ciphertext<DCRTPoly> &ciphertext, int slots, 
                                               KeyPair<DCRTPoly> keyPair, CryptoContext<DCRTPoly> &cryptoContext);


#endif

