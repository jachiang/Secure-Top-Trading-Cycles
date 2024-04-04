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
    int slotsPadded = std::pow(2, std::ceil(std::log2(slots)));
    std::vector<int64_t> invfPacked(slots*slotsPadded,invFactorialRange);
    encInvFactorial_= cryptoContext->Encrypt(keyPair.publicKey,
                                                cryptoContext->MakePackedPlaintext(invfPacked));
    std::vector<int64_t> onePacked(slots*slotsPadded,1); 
    encOne_ = cryptoContext->Encrypt(keyPair.publicKey,
                                        cryptoContext->MakePackedPlaintext(onePacked));    
    std::vector<int64_t> negOnePacked(slots*slotsPadded,plaintxtModulus-1); 
    encNegOne_ = cryptoContext->Encrypt(keyPair.publicKey,
                                        cryptoContext->MakePackedPlaintext(negOnePacked));  
    for (size_t i=1 ; i <= range ; ++i){ 
    std::vector<int64_t> negIntPacked(slots*slotsPadded,plaintxtModulus-i); 
    encNegRange_.push_back(cryptoContext->Encrypt(keyPair.publicKey,
                           cryptoContext->MakePackedPlaintext(negIntPacked)));
    };
}

Ciphertext<DCRTPoly> InitNotEqualZero::encOne() { return encOne_; }
Ciphertext<DCRTPoly> InitNotEqualZero::encNegOne() { return encNegOne_; }
Ciphertext<DCRTPoly> InitNotEqualZero::encInvFactorial() { return encInvFactorial_; }
std::vector<Ciphertext<DCRTPoly>> InitNotEqualZero::encNegRange() { return encNegRange_; }


Ciphertext<DCRTPoly> evalNotEqualZero(Ciphertext<DCRTPoly> &ciphertext,
                                  CryptoContext<DCRTPoly> &cryptoContext,
                                  InitNotEqualZero &initNotEqualZero) {
    // If x is in range, outputs 1. 
    // 1-(x-1)(x-2)...(x-r)/r! 
    std::vector<Ciphertext<DCRTPoly>> encDiffs;
    for (size_t i=0 ; i < initNotEqualZero.range ; ++i){ 
        encDiffs.push_back(cryptoContext->EvalAdd(ciphertext, initNotEqualZero.encNegRange()[i]));
    }
    encDiffs.push_back(initNotEqualZero.encInvFactorial());
    if (initNotEqualZero.range % 2 - 1) { encDiffs.push_back(initNotEqualZero.encNegOne()); }
    auto encMult = cryptoContext->EvalMultMany(encDiffs);
    // TODO: handle even/odd range.
    // if (initNotEqualZero.range % 2) { return cryptoContext->EvalAdd(initNotEqualZero.encOne(),encMult); } 
    // else { return cryptoContext->EvalAdd(initNotEqualZero.encNegRange()[0],encMult); }
    return cryptoContext->EvalAdd(initNotEqualZero.encOne(),encMult);
}

