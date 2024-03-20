#ifndef CRYPTO_NOTEQZERO_H
#define CRYPTO_NOTEQZERO_H

#include "openfhe.h"
#include "utilities.h"

using namespace lbcrypto;


class InitNotEqualZero {
public:
    InitNotEqualZero(CryptoContext<DCRTPoly> &cryptoContext, KeyPair<DCRTPoly> keyPair, int slots, int range);
    Ciphertext<DCRTPoly> encOne();
    Ciphertext<DCRTPoly> encInvFactorial();
    std::vector<Ciphertext<DCRTPoly>> encNegRange();

    const int slots;
    const int range;

private:
    Ciphertext<DCRTPoly> encOne_;
    Ciphertext<DCRTPoly> encInvFactorial_;
    std::vector<Ciphertext<DCRTPoly>> encNegRange_;
};


Ciphertext<DCRTPoly> evalNotEqualZero(Ciphertext<DCRTPoly> &ciphertext,
                                      CryptoContext<DCRTPoly> &cryptoContext,
                                      InitNotEqualZero &initNotEqualZero);


#endif