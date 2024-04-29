#ifndef CRYPTO_PREFIX_MULT_H
#define CRYPTO_PREFIX_MULT_H

#include "openfhe.h"
#include "utilities.h"

using namespace lbcrypto;


Ciphertext<DCRTPoly> evalPrefixMult(Ciphertext<DCRTPoly> &ciphertext,
                                    int slots, 
                                    CryptoContext<DCRTPoly> &cryptoContext);


Ciphertext<DCRTPoly> evalPrefixAdd(Ciphertext<DCRTPoly> &ciphertext,
                                    int slots, 
                                    CryptoContext<DCRTPoly> &cryptoContext);


class InitPreserveLeadOne {
public:
    InitPreserveLeadOne(CryptoContext<DCRTPoly> &cryptoContext, KeyPair<DCRTPoly> keyPair, int slots);
    Ciphertext<DCRTPoly> encOnes();
    Ciphertext<DCRTPoly> encNegOnes();
    Ciphertext<DCRTPoly> encLeadingOne();

    const int slots;

private:
    Ciphertext<DCRTPoly> encOnes_;
    Ciphertext<DCRTPoly> encNegOnes_;
    Ciphertext<DCRTPoly> encLeadingOne_;
};


Ciphertext<DCRTPoly> evalPreserveLeadOne(Ciphertext<DCRTPoly> &ciphertext,
                                         CryptoContext<DCRTPoly> &cryptoContext,
                                         InitPreserveLeadOne &initPreserveLeadOne);


#endif