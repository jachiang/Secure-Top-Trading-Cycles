#include "crypto_prefix_mult.h"


// TODO: Take InitRotMask as input.
Ciphertext<DCRTPoly> evalPrefixMult(Ciphertext<DCRTPoly> &ciphertext,
                                   int slots, CryptoContext<DCRTPoly> &cryptoContext) {
    // Assert slots leq rotation keys. 
    // assert(initRotsMasks.slots >= slots);

    // Prefix computation: 
    int levels = std::ceil(std::log2(slots));
        
    // Precompute rotation steps and plaintext masks. 
    // TODO: Move to precomputation (low priority, no cryptographic ops).
    std::vector<int32_t> rotSteps; std::vector<Plaintext> leadingOnesPlaintxts;
    for (size_t i = 0; i < levels; i++) { 
        rotSteps.push_back(std::pow(2, i)); 
        std::vector<int64_t> prefixOnes(slots,0); 
        for (size_t j = 0; j < rotSteps.back(); j++) { prefixOnes[j] = 1; }
        leadingOnesPlaintxts.push_back(cryptoContext->MakePackedPlaintext(prefixOnes));
    }

    // Compute prefix multiplications.
    auto ciphertext1 = ciphertext;
    for (size_t i = 0; i < levels; i++) {
        auto ciphertext2 = cryptoContext->EvalRotate(ciphertext1, -rotSteps[i]);
        // Pad ciphertext2 with leading 1's.        
        ciphertext2 = cryptoContext->EvalAdd(ciphertext2, leadingOnesPlaintxts[i]);
        // multiply ciphertext1 and ciphertext 2
        ciphertext1 = cryptoContext->EvalMult(ciphertext1, ciphertext2);
        cryptoContext->ModReduceInPlace(ciphertext1);
    }
    return ciphertext1;
} 

InitPreserveLeadOne::InitPreserveLeadOne(CryptoContext<DCRTPoly> &cryptoContext, KeyPair<DCRTPoly> keyPair, int slots) : slots(slots)
    // : initPrefixMult(cryptoContext, keyPair, slots), slots(slots) 
{
    // TODO: assert slots leq plaintext slots in cryptoContext
    // Generate rotation keys.
    cryptoContext->EvalRotateKeyGen(keyPair.secretKey, {-1});
    // Initialize prefix multiplication.
    // Generate encryption of masks. 
    std::vector<int64_t> ones(slots,1);
    std::vector<int64_t> leadingOne(slots,0); leadingOne[0] = 1;
    std::vector<int64_t> negOnes(slots,cryptoContext->GetCryptoParameters()->GetPlaintextModulus()-1);
    encOnes_ = cryptoContext->Encrypt(keyPair.publicKey,
                                        cryptoContext->MakePackedPlaintext(ones));
    encNegOnes_ = cryptoContext->Encrypt(keyPair.publicKey,
                                            cryptoContext->MakePackedPlaintext(negOnes));
    encLeadingOne_ = cryptoContext->Encrypt(keyPair.publicKey,
                                            cryptoContext->MakePackedPlaintext(leadingOne));
}

Ciphertext<DCRTPoly> InitPreserveLeadOne::encOnes() { return encOnes_; }
Ciphertext<DCRTPoly> InitPreserveLeadOne::encNegOnes() { return encNegOnes_; }
Ciphertext<DCRTPoly> InitPreserveLeadOne::encLeadingOne() { return encLeadingOne_; }


Ciphertext<DCRTPoly> evalPreserveLeadOne(Ciphertext<DCRTPoly> &ciphertext,
                                         CryptoContext<DCRTPoly> &cryptoContext,
                                         InitPreserveLeadOne &initPreserveLeadOne) {
    // (1-x0),(1-x1),...,(1-xn).
    auto encDiffs = cryptoContext->EvalAdd(cryptoContext->EvalMult(ciphertext,initPreserveLeadOne.encNegOnes()),
                                           initPreserveLeadOne.encOnes());
    // y0, y1,..., yn: yi = ith multiplicative prefix.
    auto encPrefix = evalPrefixMult(encDiffs, initPreserveLeadOne.slots, cryptoContext);
    // x0, x1*y0 ,...,   xn*yn-1
    auto result = cryptoContext->EvalMult(ciphertext,
                                   cryptoContext->EvalAdd(initPreserveLeadOne.encLeadingOne(), 
                                                          cryptoContext->EvalRotate(encPrefix,-1)));
    cryptoContext->ModReduceInPlace(result);
    return result;
}