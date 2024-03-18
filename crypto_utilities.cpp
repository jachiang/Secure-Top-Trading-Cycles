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

// Class to initialize rotation keys and masks.
InitRotsMasks::InitRotsMasks(CryptoContext<DCRTPoly> &cryptoContext, KeyPair<DCRTPoly> keyPair, int slots) :
    slots(slots) {
    // (1) Generate rotation keys for |slots| number of steps.
    std::vector<int32_t> rotIndices;
    for (size_t i = 0; i <= slots; i++) { rotIndices.push_back(-i); rotIndices.push_back(i);}
    cryptoContext->EvalRotateKeyGen(keyPair.secretKey, rotIndices);
    // (2)Generate Eval Sum Key for EvalInnerProduct.
    cryptoContext->EvalSumKeyGen(keyPair.secretKey);
    // (3) Generate ciphertext masks for extraction of ciphertext slot values.
    for (size_t elem=0 ; elem < slots ; ++elem){ 
        std::vector<int64_t> mask(slots,0); mask[elem] = 1;
        encMasks_.push_back(cryptoContext->Encrypt(keyPair.publicKey,
                                                    cryptoContext->MakePackedPlaintext(mask)));
    }
}
std::vector<Ciphertext<DCRTPoly>> InitRotsMasks::encMasks() { return encMasks_; }