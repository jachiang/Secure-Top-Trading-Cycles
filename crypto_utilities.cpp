#include "crypto_utilities.h"


void printEncMatRows(std::vector<Ciphertext<DCRTPoly>> &encMatRows, CryptoContext<DCRTPoly> &cryptoContext, KeyPair<DCRTPoly> keyPair){
    for (int row=0; row < encMatRows.size(); ++row){
        Plaintext plaintext;
        cryptoContext->Decrypt(keyPair.secretKey, encMatRows[row], &plaintext); 
        plaintext->SetLength(encMatRows.size()); auto payload = plaintext->GetPackedValue();
        std::cout << payload << std::endl;
    }

}


void printEncMatElems(std::vector<std::vector<Ciphertext<DCRTPoly>>> &encMatElems, CryptoContext<DCRTPoly> &cryptoContext, KeyPair<DCRTPoly> keyPair){
    for (int row=0; row < encMatElems.size(); ++row){
        for (int col=0; col < encMatElems.size(); ++col){
            Plaintext plaintext;
            cryptoContext->Decrypt(keyPair.secretKey, encMatElems[row][col], &plaintext); 
            plaintext->SetLength(1); auto payload = plaintext->GetPackedValue();
            std::cout << payload << " ";
        }
        std::cout << std::endl;
    }
}


InitRotsMasks::InitRotsMasks(CryptoContext<DCRTPoly> &cryptoContext, KeyPair<DCRTPoly> keyPair, int slots) :
    slots(slots) {
    std::vector<int32_t> rotIndices;
    // (1a) Generate rotation keys for +/-[slots] number of steps.
    for (size_t i = 0; i <= slots; i++) { rotIndices.push_back(-i); rotIndices.push_back(i);}
    // (1b) Generate rotation keys for -i*|slots+1| number of steps.
    for (size_t i = 1; i < slots*slots; i++) { rotIndices.push_back(-i*(slots+1)); }
    cryptoContext->EvalRotateKeyGen(keyPair.secretKey, rotIndices);
    // (2)Generate Eval Sum Key for EvalInnerProduct.
    cryptoContext->EvalSumKeyGen(keyPair.secretKey);
    // (3a) Generate ciphertext masks for extraction of individual ciphertext slot values.
    for (size_t elem=0 ; elem < slots ; ++elem){ 
        std::vector<int64_t> mask(slots,0); mask[elem] = 1;
        encMasks_.push_back(cryptoContext->Encrypt(keyPair.publicKey,
                                                    cryptoContext->MakePackedPlaintext(mask)));
    }
    // (3b) Generate ciphertext masks for extraction of slot values in fully packed ciphertexts.
    for (size_t elem=0 ; elem < slots*slots ; ++elem){ 
        std::vector<int64_t> mask(slots*slots*(slots+1),0); mask[elem*(slots+1)] = 1;
        encMasksFullyPacked_.push_back(cryptoContext->Encrypt(keyPair.publicKey,
                                                              cryptoContext->MakePackedPlaintext(mask)));
    }
}

std::vector<Ciphertext<DCRTPoly>> InitRotsMasks::encMasks() { return encMasks_; }
std::vector<Ciphertext<DCRTPoly>> InitRotsMasks::encMasksFullyPacked() { return encMasksFullyPacked_; }


Ciphertext<DCRTPoly> evalExponentiate(Ciphertext<DCRTPoly> &ciphertext, int exponent, 
                                      CryptoContext<DCRTPoly> &cryptoContext) {
    // Get msb position of exponent.
    int numBits = sizeof(int) * 8;
    int msbPosition = -1;
    for (int i = numBits - 1; i >= 0; i--) {
        if ((exponent >> i & 1) && (msbPosition == -1)) {
            msbPosition = i + 1;
            break;
        }
    }
    // Compute squarings of p.
    std::vector<Ciphertext<DCRTPoly>> ciphertexts_squarings;
    ciphertexts_squarings.push_back(ciphertext);
    for (int i = 1; i < msbPosition; i++) {
        ciphertexts_squarings.push_back(cryptoContext->EvalMult(ciphertexts_squarings[i-1], 
                                                                ciphertexts_squarings[i-1]));
    }
    // Select required squarings.
    std::vector<Ciphertext<DCRTPoly>> ciphertexts_squarings_container;
    for (int i = 0; i < msbPosition; i++) {
        if (exponent >> i & 1) {
            ciphertexts_squarings_container.push_back(ciphertexts_squarings[i]) ;
        }
    }
    // Multiply selected squarings.
    return cryptoContext->EvalMultMany(ciphertexts_squarings_container);
}


void refreshInPlace(Ciphertext<DCRTPoly> &ciphertext, int slots, 
                    KeyPair<DCRTPoly> keyPair, CryptoContext<DCRTPoly> &cryptoContext){
    Plaintext plaintextExpRes;
    cryptoContext->Decrypt(keyPair.secretKey, ciphertext, &plaintextExpRes); 
    plaintextExpRes->SetLength(slots); auto payload = plaintextExpRes->GetPackedValue();
    ciphertext = cryptoContext->Encrypt(keyPair.publicKey, cryptoContext->MakePackedPlaintext(payload));
}

std::vector<Ciphertext<DCRTPoly>> refreshElems(Ciphertext<DCRTPoly> &ciphertext, int slots, 
                                               KeyPair<DCRTPoly> keyPair, CryptoContext<DCRTPoly> &cryptoContext){
    std::vector<Ciphertext<DCRTPoly>> ciphertexts; 
    Plaintext plaintext;
    cryptoContext->Decrypt(keyPair.secretKey, ciphertext, &plaintext); 
    plaintext->SetLength(slots); auto payload = plaintext->GetPackedValue();
    for (int i = 0; i < slots; i++) {
        std::vector<int64_t> elementPlaintext(slots,0);
        elementPlaintext[0] = payload[i];
        ciphertexts.push_back(cryptoContext->Encrypt(keyPair.publicKey, cryptoContext->MakePackedPlaintext(elementPlaintext)));
    }
    return ciphertexts;
}
