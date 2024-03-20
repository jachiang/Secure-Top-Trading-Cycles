#include "crypto_noteqzero.h"


InitNotEqualZero::InitNotEqualZero(CryptoContext<DCRTPoly> &cryptoContext, KeyPair<DCRTPoly> keyPair, int slots, int range) :
    slots(slots), range(range)
{
    // TODO: assert range leq plaintext slots in cryptoContext
    // Precompute constants.
    auto plaintxtModulus = cryptoContext->GetCryptoParameters()->GetPlaintextModulus();
    int factorialRange = modFactorial(range, plaintxtModulus); 
    int invFactorialRange = modInverse(factorialRange, plaintxtModulus); 

    // Packed encryption of constants.
    std::vector<int64_t> invfPacked(slots,invFactorialRange);
    encInvFactorial_= cryptoContext->Encrypt(keyPair.publicKey,
                                                cryptoContext->MakePackedPlaintext(invfPacked));
    std::vector<int64_t> onePacked(slots,1); 
    encOne_ = cryptoContext->Encrypt(keyPair.publicKey,
                                        cryptoContext->MakePackedPlaintext(onePacked));    
    for (size_t i=1 ; i <= range ; ++i){ 
    std::vector<int64_t> negIntPacked(slots,plaintxtModulus-i); 
    encNegRange_.push_back(cryptoContext->Encrypt(keyPair.publicKey,
                            cryptoContext->MakePackedPlaintext(negIntPacked)));
    };
}

Ciphertext<DCRTPoly> InitNotEqualZero::encOne() { return encOne_; }
Ciphertext<DCRTPoly> InitNotEqualZero::encInvFactorial() { return encInvFactorial_; }
std::vector<Ciphertext<DCRTPoly>> InitNotEqualZero::encNegRange() { return encNegRange_; }


Ciphertext<DCRTPoly> evalNotEqualZero(Ciphertext<DCRTPoly> &ciphertext,
                                  CryptoContext<DCRTPoly> &cryptoContext,
                                  InitNotEqualZero &initNotEqualZero) {
    // Compute binary map for range r: 1-(x-1)(x-2)...(x-r)/r!
    std::vector<Ciphertext<DCRTPoly>> encDiffs;
    for (size_t i=0 ; i < initNotEqualZero.range ; ++i){ 
        encDiffs.push_back(cryptoContext->EvalAdd(ciphertext, initNotEqualZero.encNegRange()[i]));
    }
    encDiffs.push_back(initNotEqualZero.encInvFactorial());
    auto encMult = cryptoContext->EvalMultMany(encDiffs);
    if (initNotEqualZero.range % 2) { return cryptoContext->EvalAdd(initNotEqualZero.encOne(),encMult); } 
    else { return cryptoContext->EvalAdd(initNotEqualZero.encNegRange()[0],encMult); }
}

