#include "crypto_prefix_mult.h"


Ciphertext<DCRTPoly> evalPrefixMult(Ciphertext<DCRTPoly> &ciphertext,
                                   int n, CryptoContext<DCRTPoly> &cryptoContext) {

    int depth = std::ceil(std::log2(n));
    int slotsPadded = std::pow(2,depth);
    std::vector<int32_t> rotSteps; std::vector<Plaintext> leadingOnesPlaintxts;
    for (int i = 0; i < depth; i++) { 
        rotSteps.push_back(std::pow(2, i)); 
        std::vector<int64_t> prefixOnes;
        for (int copy=0;copy<slotsPadded*n;copy++){ 
            for (int elem=0; elem<rotSteps.back(); elem++) { prefixOnes.push_back(1); }
            for (int elem=0; elem<slotsPadded-rotSteps.back(); elem++) { prefixOnes.push_back(0); }
        }
        leadingOnesPlaintxts.push_back(cryptoContext->MakePackedPlaintext(prefixOnes));
    }
    auto ciphertext1 = ciphertext;
    for (int lvl = 0; lvl < depth; lvl++) {
        auto ciphertext2 = cryptoContext->EvalRotate(ciphertext1, -rotSteps[lvl]);
        ciphertext2 = cryptoContext->EvalAdd(ciphertext2, leadingOnesPlaintxts[lvl]);
        ciphertext1 = cryptoContext->EvalMult(ciphertext1, ciphertext2);
        cryptoContext->ModReduceInPlace(ciphertext1);
    }
    return ciphertext1;
} 


Ciphertext<DCRTPoly> evalPrefixAdd(Ciphertext<DCRTPoly> &ciphertext,
                                   int slots, CryptoContext<DCRTPoly> &cryptoContext) {

    int levels = std::ceil(std::log2(slots));
    std::vector<int32_t> rotSteps; 
    for (int i = 0; i < levels; i++) { 
        rotSteps.push_back(std::pow(2, i)); 
    }
    auto ciphertext1 = ciphertext;
    for (int i = 0; i < levels; i++) {
        auto ciphertext2 = cryptoContext->EvalRotate(ciphertext1, rotSteps[i]);
        ciphertext1 = cryptoContext->EvalAdd(ciphertext1, ciphertext2);
    }
    return ciphertext1;
} 


InitPreserveLeadOne::InitPreserveLeadOne(CryptoContext<DCRTPoly> &cryptoContext, KeyPair<DCRTPoly> keyPair, int slots) : slots(slots)
{

    cryptoContext->EvalRotateKeyGen(keyPair.secretKey, {-1});
    // Generate encryption of masks. 
    std::vector<int64_t> ones(slots,1);
    std::vector<int64_t> negOnes(slots,cryptoContext->GetCryptoParameters()->GetPlaintextModulus()-1);
    std::vector<int64_t> leadingOne(slots,0); leadingOne[0]=1;
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
    auto encPrefix = evalPrefixMult(encDiffs,initPreserveLeadOne.slots,cryptoContext);
    // x0, x1*y0 ,...,   xn*yn-1
    auto result = cryptoContext->EvalMult(ciphertext,
                                   cryptoContext->EvalAdd(initPreserveLeadOne.encLeadingOne(), 
                                                          cryptoContext->EvalRotate(encPrefix,-1)));
    cryptoContext->ModReduceInPlace(result);
    return result;
}