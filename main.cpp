//==================================================================================
// BSD 2-Clause License
//
// Copyright (c) 2014-2022, NJIT, Duality Technologies Inc. and other contributors
//
// All rights reserved.
//
// Author TPOC: contact@openfhe.org
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//==================================================================================

/*
  Example of a computation circuit of depth 3
  BGVrns demo for a homomorphic multiplication of depth 6 and three different approaches for depth-3 multiplications
 */

#define PROFILE

#include "utilities.h"
#include "crypto_utilities.h"
#include "crypto_enc_transform.h"
#include "crypto_matrix_operations.h"
#include "crypto_prefix_mult.h"
#include "crypto_noteqzero.h"

#include <cassert>
#include <iostream>
#include <iterator>
#include <cmath>
#include <vector>

#include "openfhe.h"

using namespace lbcrypto;





int main(int argc, char* argv[]) {

    ////////////////////////////////////////////////////////////
    // Set-up of parameters
    ////////////////////////////////////////////////////////////

    // TODO: Generate dedicated parameters for phase (1), (2a), (2b) and (3).

    // benchmarking variables
    TimeVar t;
    double processingTime(0.0);

    // Crypto Parameters
    // # of evalMults = 3 (first 3) is used to support the multiplication of 7
    // ciphertexts, i.e., ceiling{log2{7}} Max depth is set to 3 (second 3) to
    // generate homomorphic evaluation multiplication keys for s^2 and s^3
    CCParams<CryptoContextBGVRNS> params1, params2a, params2b, params3;

    int chosen_depth1 = 10;
    int chosen_depth2a = 10;
    int chosen_depth2b = 10;
    int chosen_depth3 = 10;

    int chosen_ptxtmodulus = 65537;
    params1.SetPlaintextModulus(chosen_ptxtmodulus); params1.SetMultiplicativeDepth(chosen_depth1);  params1.SetMaxRelinSkDeg(3);
    params2a.SetPlaintextModulus(chosen_ptxtmodulus);params2a.SetMultiplicativeDepth(chosen_depth2a);params2a.SetMaxRelinSkDeg(3);
    params2b.SetPlaintextModulus(chosen_ptxtmodulus);params2b.SetMultiplicativeDepth(chosen_depth2b);params2b.SetMaxRelinSkDeg(3);
    params3.SetPlaintextModulus(chosen_ptxtmodulus); params3.SetMultiplicativeDepth(chosen_depth3);  params3.SetMaxRelinSkDeg(3);
    // TODO: Fermats thm works for p = 786433, dep = 20.

    CryptoContext<DCRTPoly> cryptoContext1 = GenCryptoContext(params1);
    CryptoContext<DCRTPoly> cryptoContext2a = GenCryptoContext(params2a);
    CryptoContext<DCRTPoly> cryptoContext2b = GenCryptoContext(params2b);
    CryptoContext<DCRTPoly> cryptoContext3 = GenCryptoContext(params3);

    std::cout << "Ring dimension N: " << cryptoContext1->GetRingDimension() << std::endl;

    cryptoContext1->Enable(PKE); cryptoContext1->Enable(KEYSWITCH); cryptoContext1->Enable(LEVELEDSHE); cryptoContext1->Enable(ADVANCEDSHE);
    cryptoContext2a->Enable(PKE); cryptoContext2a->Enable(KEYSWITCH); cryptoContext2a->Enable(LEVELEDSHE); cryptoContext2a->Enable(ADVANCEDSHE);
    cryptoContext2b->Enable(PKE); cryptoContext2b->Enable(KEYSWITCH); cryptoContext2b->Enable(LEVELEDSHE); cryptoContext2b->Enable(ADVANCEDSHE);
    cryptoContext3->Enable(PKE); cryptoContext3->Enable(KEYSWITCH); cryptoContext3->Enable(LEVELEDSHE); cryptoContext3->Enable(ADVANCEDSHE);

    // std::cout << "Plaintext modulus p = " << cryptoContext1->GetCryptoParameters()->GetPlaintextModulus() << std::endl;
    // std::cout << "Cyclotomic order n = " << cryptoContext1->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2 << std::endl;
    // std::cout << "log2 q = "
    //           << log2(cryptoContext1->GetCryptoParameters()->GetElementParams()->GetModulus().ConvertToDouble())
    //           << std::endl;

    // Initialize Public Key Containers
    KeyPair<DCRTPoly> keyPair;

    // Perform Key Generation Operation
    TIC(t);

    keyPair = cryptoContext1->KeyGen();

    processingTime = TOC(t);
    std::cout << "Key generation time: " << processingTime << "ms" << std::endl;

    if (!keyPair.good()) {
        std::cout << "Key generation failed!" << std::endl;
        exit(1);
    }

    std::cout << "Running key generation for homomorphic multiplication "
                 "evaluation keys..."
              << std::endl;

    TIC(t);

    cryptoContext1->EvalMultKeysGen(keyPair.secretKey);

    processingTime = TOC(t);
    std::cout << "Key generation time for homomorphic multiplication evaluation keys: " << processingTime << "ms"
              << std::endl;


    ////////////////////////////////////////////////////////////
    // Top Trading Cycle Algorithm.
    ////////////////////////////////////////////////////////////

    //==========================================================
    // Offline phase.
    //==========================================================

    // Offline: User preferences.
    std::vector<std::vector<int64_t>> userInputs;

    // userInputs.push_back({4, 1, 2, 3, 0});
    // userInputs.push_back({4, 3, 2, 1, 0});
    // userInputs.push_back({4, 1, 0, 2, 3});
    // userInputs.push_back({1, 3, 4, 0, 2});
    // userInputs.push_back({3, 1, 2, 0, 4});

    // Required depth: 8.
    // userInputs.push_back({4, 1, 2, 3, 0, 5});
    // userInputs.push_back({4, 3, 2, 1, 0, 5});
    // userInputs.push_back({4, 1, 0, 2, 3, 5});
    // userInputs.push_back({1, 3, 4, 0, 2, 5});
    // userInputs.push_back({3, 1, 2, 0, 4, 5});
    // userInputs.push_back({3, 1, 2, 0, 4, 5});
    // userInputs.push_back({3, 1, 2, 0, 4, 5});
    
    // Required depth: 9.
    userInputs.push_back({4, 1, 2, 3, 0, 5, 6, 7, 8, 9});
    userInputs.push_back({4, 3, 2, 1, 0, 5, 6, 7, 8, 9});
    userInputs.push_back({4, 1, 0, 2, 3, 5, 6, 7, 8, 9});
    userInputs.push_back({1, 3, 4, 0, 2, 5, 6, 7, 8, 9});
    userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9});
    userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9});
    userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9});
    userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9});
    userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9});
    userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9});

    // Required depth: 10.
    // userInputs.push_back({4, 1, 2, 3, 0, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19});
    // userInputs.push_back({4, 3, 2, 1, 0, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19});
    // userInputs.push_back({4, 1, 0, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19});
    // userInputs.push_back({1, 3, 4, 0, 2, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19});
    // userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19});
    // userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19});
    // userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19});
    // userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19});
    // userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19});
    // userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19});
    // userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19});
    // userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19});
    // userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19});
    // userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19});
    // userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19});
    // userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19});
    // userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19});
    // userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19});
    // userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19});
    // userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19});

    auto n = userInputs.size();

    // Offline: Encryption of user preferences.
    // -----------------------------------------------------------------------
    // Represent user preferences as permutation matrices and their transpose.
    // Generate row-wise encryptions of each user matrix (packing mode 1),
    // or a single encryption of all user matrices (packing mode 2).
    
    std::vector<std::vector<std::vector<int64_t>>> userPrefList;
    std::vector<std::vector<std::vector<int64_t>>> userPrefTransposedList;
    std::vector<std::vector<Ciphertext<DCRTPoly>>> encUserPrefList;             
    std::vector<std::vector<Ciphertext<DCRTPoly>>> encUserPrefTransposedList;   

    int k_ceil = std::ceil(std::log2(n));
    int slotsPadded = std::pow(2,k_ceil);

    for (int user=0; user < n; ++user) { 
        // Build user preference-permutation matrix.
        std::vector<std::vector<int64_t>> userPrefMatrix;
        for (int col=0; col < n; ++col) {
            std::vector<int64_t> row(n,0); row[userInputs[user][col]] = 1;
            userPrefMatrix.push_back(row);
        } 
        userPrefList.push_back(userPrefMatrix);
        // Transpose user preference-permutation matrix.
        std::vector<std::vector<int64_t>> userPrefMatrixTransposed;
        std::vector<int64_t> zeroRow(n,0);
        for (int row=0; row < n; ++row) { userPrefMatrixTransposed.push_back(zeroRow); }
        for (int row=0; row < n; ++row) { for (int col=0; col < n; ++col) {
            userPrefMatrixTransposed[col][row] = userPrefMatrix[row][col]; }}
        userPrefTransposedList.push_back(userPrefMatrixTransposed);
        // For packing mode 1: row-wise encryption of each permutation matrix.
        std::vector<Ciphertext<DCRTPoly>> encUserPref;
        std::vector<Ciphertext<DCRTPoly>> encUserPrefTransposed;
        for (int row=0; row<n ; ++row){ 
            encUserPref.push_back(cryptoContext1->Encrypt(keyPair.publicKey,
                                  cryptoContext1->MakePackedPlaintext(userPrefMatrix[row])));
            encUserPrefTransposed.push_back(cryptoContext1->Encrypt(keyPair.publicKey,
                                            cryptoContext1->MakePackedPlaintext(userPrefMatrixTransposed[row])));
        }
        encUserPrefList.push_back(encUserPref);
        encUserPrefTransposedList.push_back(encUserPrefTransposed); 
    }
    // For packing mode 2: encryption of all user permutation matrices in individual ciphertexts.
    std::vector<int64_t> userPrefPacked(slotsPadded*slotsPadded*n,0);
    std::vector<int64_t> userPrefTransposedPacked(slotsPadded*slotsPadded*n,0);
    for (int user=0; user<n ; ++user){ 
        for (int row=0; row<n ; ++row){ 
            for (int col=0; col<n ; ++col){ 
                int pos = user*slotsPadded*slotsPadded+row*slotsPadded+col;
                userPrefPacked[pos] = userPrefList[user][row][col];
                userPrefTransposedPacked[pos] = userPrefTransposedList[user][row][col];
            }
        }
    }
    auto encUserPrefFullyPacked = cryptoContext1->Encrypt(keyPair.publicKey,
                                  cryptoContext1->MakePackedPlaintext(userPrefPacked));
    auto encUserPrefTransPacked = cryptoContext1->Encrypt(keyPair.publicKey,
                                  cryptoContext1->MakePackedPlaintext(userPrefTransposedPacked));

    // Offline: Init objects.
    // -----------------------------------------------------------------------
    InitRotsMasks initRotsMasks(cryptoContext1,keyPair,n);
    InitNotEqualZero initNotEqualZero(cryptoContext1,keyPair,n,userInputs.size()); // TODO: Over packed ciphertexts?
    InitPreserveLeadOne initPreserveLeadOne(cryptoContext1,keyPair,n);
    
    // Offline: Initialize encrypted constants.
    // -----------------------------------------------------------------------
    auto slotTotal = cryptoContext1->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder()/2;
    // TODO: move to initialization class.
    std::vector<int64_t> zeros(slotTotal,0);
    std::vector<int64_t> range; for (int i=0; i<n; ++i) { range.push_back(i+1); } // TODO: not suited for packing.
    std::vector<int64_t> ones; std::vector<int64_t> negOnes;
    for (int iter=0; iter<slotsPadded*slotsPadded; iter++){ 
        for (int elem = 0; elem <n; elem++){ ones.push_back(1); negOnes.push_back(-1); }
        for (int elem = 0; elem <slotsPadded-n; elem++){ ones.push_back(0); negOnes.push_back(0); }
     }
    auto encZeros = cryptoContext1->Encrypt(keyPair.publicKey,
                                           cryptoContext1->MakePackedPlaintext(zeros));
    auto encOnes = cryptoContext1->Encrypt(keyPair.publicKey,
                                          cryptoContext1->MakePackedPlaintext(ones));
    auto encNegOnes = cryptoContext1->Encrypt(keyPair.publicKey,
                                             cryptoContext1->MakePackedPlaintext(negOnes));
    auto encRange = cryptoContext1->Encrypt(keyPair.publicKey,
                                           cryptoContext1->MakePackedPlaintext(range));


    //==========================================================
    // Online phase.
    //==========================================================

    // 1: Row/Col-wise matrix packing
    // 2: Full matrix packing
    int packingMode = 2;

    // Init output ciphertext.
    auto enc_output = encZeros; 
    auto encUserAvailability = encOnes; 

    // Log crypto operations.
    CryptoOpsLogger cryptoOpsLogger;

    // Loop over all rounds.
    for (int i = 0; i < n ; ++i)
    {
        std::cout << "-------------" << std::endl;
        std::cout << "Round ... " << i+1 << "/" << n << std::endl;
        std::cout << "-------------" << std::endl;

        //----------------------------------------------------------
        // (1) Update adjacency matix.
        //----------------------------------------------------------

        std::vector<Ciphertext<DCRTPoly>> encRowsAdjMatrix;
        Ciphertext<DCRTPoly> encAdjMatrixPacked;

        // Begin: Timer.
        TIC(t);
        if (packingMode == 1){
            for (int user = 0; user < n; ++user){
                auto encUserAvailablePref = evalMatrixVecMult(encUserPrefList[user], encUserAvailability,
                                                              cryptoContext1, initRotsMasks); // 2 mult-depth
                auto encUserFirstAvailablePref = evalPreserveLeadOne(encUserAvailablePref, cryptoContext1, initPreserveLeadOne); // 2+log(n) mult-depth
                // TODO: extract rows in refresh.
                encRowsAdjMatrix.push_back(evalMatrixVecMult(encUserPrefTransposedList[user], encUserFirstAvailablePref,
                                                             cryptoContext1, initRotsMasks)); // 2 mult-depth
            }
        }
        else if (packingMode == 2) {
            // Multiplication: packed user preference matrices * availability matrix (repeated) 
            // padded(user1row1)        | padded(user1row2)        | padded(user1row3) ...    | padded(user2row1)        | padded(user2row2)        ...
            // padded(userAvailability) | padded(userAvailability) | padded(userAvailability) | padded(userAvailability) | padded(userAvailability) ...
            auto encResMult = cryptoContext1->EvalMult(encUserAvailability,encUserPrefFullyPacked); 
            auto encResInnerProd = evalPrefixAdd(encResMult,slotsPadded,cryptoContext1);      

            // Extract resulting slot values.
            std::vector<std::vector<Ciphertext<DCRTPoly>>> encAvailPrefElems;
            for (int user = 0; user < n; user++){
                std::vector<Ciphertext<DCRTPoly>> encElemsRow;
                for (int slot = 0; slot < n; slot++){
                    auto encMaskedElem = cryptoContext1->EvalMult(encResInnerProd,initRotsMasks.encMasksFullyPacked()[user*n+slot]);
                    cryptoContext1->ModReduceInPlace(encMaskedElem);
                    encElemsRow.push_back(cryptoContext1->EvalRotate(encMaskedElem,((user*slotsPadded+slot)*slotsPadded)));
                }
                encAvailPrefElems.push_back(encElemsRow);
            }
            // Repack into single ciphertext (for later matrix multiplication).
            // padded(user1AvailPref) ... | padded(user2AvailPref)| ... | padded(user3AvailPref) | ... (each repeated 2^ceil(log(n)) times)
            std::vector<Ciphertext<DCRTPoly>> encAvailPrefContainer;
            for (int user = 0; user < n; user++){
                for (int slot=0; slot < n; slot++){
                    encAvailPrefContainer.push_back(cryptoContext1->EvalRotate(
                                                    encAvailPrefElems[user][slot],-(user*slotsPadded*slotsPadded+slot)));
                }
            }
            auto encUniqueAvailRows = cryptoContext1->EvalAddMany(encAvailPrefContainer); 
            auto encRepAvailRows = encUniqueAvailRows;
            for (int k = 0; k < k_ceil; k++) {  // Copy log(n) times.
                auto encRepAvailRows_ = cryptoContext1->EvalRotate(encRepAvailRows,-std::pow(2,k)*slotsPadded);
                encRepAvailRows = cryptoContext1->EvalAdd(encRepAvailRows,encRepAvailRows_); 
            } 
            // Preserve leading One.
            auto encFirstAvailPrefPacked = evalPreserveLeadOne(encRepAvailRows, cryptoContext1, initPreserveLeadOne); // 2+log(n) mult-depth

            // Preferred available  * Transposed Preference Matrices
            auto encResMult2 = cryptoContext1->EvalMult(encFirstAvailPrefPacked,encUserPrefTransPacked); 
            encAdjMatrixPacked = evalPrefixAdd(encResMult2,slotsPadded,cryptoContext1);
        }
        else { return 1; }
        processingTime = TOC(t);    // End: Timer.
        std::cout << "Online part 1 - Adjacency matrix update time: " << processingTime << "ms" << std::endl;

        // Refresh after phase (1).
        //----------------------------------------------------------
        std::cout << "Adjacency Matrix: " << std::endl;
        std::vector<std::vector<Ciphertext<DCRTPoly>>> encAdjMatrixElems;
        if (packingMode == 1){
            for (int row=0; row < n; ++row){
                Plaintext plaintext;
                cryptoContext1->Decrypt(keyPair.secretKey, encRowsAdjMatrix[row], &plaintext); 
                plaintext->SetLength(n); auto payload = plaintext->GetPackedValue();
                // Print adjacence matrix.
                std::cout << payload << std::endl;
                // Refresh row-wise encrypted ciphertexts.
                encAdjMatrixElems.push_back(refreshElems(encRowsAdjMatrix[row],n,keyPair,cryptoContext1));
                // WIP: for step 3, refresh adjacency matrix ciphertext in row form.
                refreshInPlace(encRowsAdjMatrix[row],n,keyPair,cryptoContext1); 
            }
        }
        else if (packingMode == 2) {
            Plaintext plaintext;
            cryptoContext1->Decrypt(keyPair.secretKey,encAdjMatrixPacked,&plaintext); 
            plaintext->SetLength(slotsPadded*slotsPadded*n); auto payload = plaintext->GetPackedValue();
            // For phase 2, refresh adjacency matrix and output element wise encryptions.
            encRowsAdjMatrix.clear();
            for (int row=0; row < n; ++row){
                std::vector<int64_t> adjMatrixRow;
                std::vector<Ciphertext<DCRTPoly>> encAdjMatrixRow;
                for (int col=0; col < n; ++col){
                    int pos = row*slotsPadded*slotsPadded + col*slotsPadded;
                    adjMatrixRow.push_back(payload[pos]);
                    std::vector<int64_t> elementFirst(n,0); elementFirst[0] = payload[pos];
                    encAdjMatrixRow.push_back(cryptoContext1->Encrypt(keyPair.publicKey,
                                              cryptoContext1->MakePackedPlaintext(elementFirst)));
                }
                encAdjMatrixElems.push_back(encAdjMatrixRow);
                // Print adjacency matrix.
                std::cout << adjMatrixRow << std::endl; 
                // WIP: for phase 3, refresh adjacency matrix ciphertext in row form.
                encRowsAdjMatrix.push_back(cryptoContext1->Encrypt(keyPair.publicKey,
                                           cryptoContext1->MakePackedPlaintext(adjMatrixRow)));
            }
        }
        else { return 1; }

        //----------------------------------------------------------
        // (2) Cycle finding.
        //----------------------------------------------------------

        // 2a) Matrix exponentiation.
        //----------------------------------------------------------
        // Cycle finding result [res_1, ..., res_n]. On cycle, res_i = 1. Not on cycle: res_i = 0.
        std::vector<std::vector<Ciphertext<DCRTPoly>>> encMatrixExpElems; 

        TIC(t);
        if (packingMode == 1){
            // No packing (mode 0), Row/Col matrix packing (mode 1), Full matrix packing (mode 2) 
            encMatrixExpElems = evalMatSqMul(encAdjMatrixElems,n,
                                             1, // Packing mode.
                                             cryptoContext1,initRotsMasks,cryptoOpsLogger,keyPair); // TODO: keypair param for debugging
        }
        else if (packingMode == 2) {
            encMatrixExpElems = evalMatSqMul(encAdjMatrixElems,n,
                                             2, // Packing mode.
                                             cryptoContext1,initRotsMasks,cryptoOpsLogger,keyPair); // TODO: keypair param for debugging, remove later.   
        }       
        else { return 1; }
        processingTime = TOC(t);
        std::cout << "Online part 2a - Matrix exponentiation: " << processingTime << " ms" << std::endl;

        // Refresh after phase (2a)
        //----------------------------------------------------------
        std::vector<Ciphertext<DCRTPoly>> encMatrixExp; // Refreshed & repacked output in packing mode 1.
        Ciphertext<DCRTPoly> encMatrixExpPacked;        // Refreshed & repacked output in packing mode 2.

        if (packingMode == 1){
            // TODO: perform this elem to column transform in refresh.
            encMatrixExp = encElem2Cols(encMatrixExpElems,cryptoContext1,initRotsMasks,cryptoOpsLogger);
            for (int col=0; col < n; ++col){ refreshInPlace(encMatrixExp[col],n,keyPair,cryptoContext1); } 
            printEncMatRows(encMatrixExp,cryptoContext1,keyPair);
        }
        else if (packingMode == 2) {
            // Refresh and repack ciphertexts
            // packed(col1) | packed(col2) | packed(col3) ...     
            std::vector<int64_t> packedMatrix(slotsPadded*n,0);
            for (int row = 0; row < n; row++){
                std::vector<int64_t> matrixElemsRow;
                for (int col = 0; col < n; col++){
                    Plaintext plaintext;
                    cryptoContext1->Decrypt(keyPair.secretKey,encMatrixExpElems[row][col],&plaintext); 
                    plaintext->SetLength(n); auto payload = plaintext->GetPackedValue();
                    int pos = col*slotsPadded + row;
                    packedMatrix[pos] = payload[0];
                }
            }
            encMatrixExpPacked = cryptoContext1->Encrypt(keyPair.publicKey,
                                 cryptoContext1->MakePackedPlaintext(packedMatrix));   
        }
        else { return 1; }

        // 2b) Cycle computation.
        //----------------------------------------------------------
        Ciphertext<DCRTPoly> enc_u;             // Output in packing mode 1.
        Ciphertext<DCRTPoly> enc_u_unmasked;    // Output in packing mode 2.

        TIC(t); // Begin: Timer.
        if (packingMode == 1){
            enc_u = evalVecMatrixMult(encOnes,encMatrixExp,cryptoContext1,initRotsMasks,cryptoOpsLogger); 
            enc_u = evalNotEqualZero(enc_u,cryptoContext1,initNotEqualZero); 
        }
        else if (packingMode == 2) {
            // packed(col1) | packed(col2) | packed(col3) ...     
            // packed(ones) | packed(ones) | packed(ones) ...
            auto encResMult = cryptoContext1->EvalMult(encMatrixExpPacked,encOnes); 
            auto encResInnerProd = evalPrefixAdd(encResMult,slotsPadded,cryptoContext1);      
            // Not equal zero (SIMD operation).
            enc_u_unmasked = evalNotEqualZero(encResInnerProd,cryptoContext1,initNotEqualZero); 
        }
        else { return 1; }
        processingTime = TOC(t); // End: Timer.
        std::cout << "Online part 2b - Cycle computation: " << processingTime << "ms" << std::endl;

        // Refresh after phase (2b)
        //----------------------------------------------------------
        if (packingMode == 1){
            refreshInPlace(enc_u,n,keyPair, cryptoContext1); 
        }
        // Refresh and extract resulting slot values.
        else if (packingMode == 2) {
            // Generate fresh encryption of [ u_1 | u_2 | ... | u_n ... ] (replicated in packed fashion)
            Plaintext plaintext;
            cryptoContext1->Decrypt(keyPair.secretKey,enc_u_unmasked,&plaintext); 
            plaintext->SetLength(n*slotsPadded); auto payload = plaintext->GetPackedValue();
            std::vector<int64_t> uElems;
            for (int user = 0; user < n; user++){
                uElems.push_back(payload[user*slotsPadded]); // TODO: Pack this for phase 3.
            }
            enc_u = cryptoContext1->Encrypt(keyPair.publicKey,
                    cryptoContext1->MakePackedPlaintext(uElems));    
        }
        else { return 1; }

        //----------------------------------------------------------
        // (3) Update user availability and outputs.
        //----------------------------------------------------------

        TIC(t); // Begin: Timer.
        // Compute current preference index (t) for all users in packed ciphertext.
        std::vector<Ciphertext<DCRTPoly>> enc_elements;
        for (int user=0; user < n; ++user){ 
            // Note: encRowsAdjMatrix must be refreshed after (1)
            auto enc_t_user = cryptoContext1->EvalInnerProduct(encRowsAdjMatrix[user], encRange, 
                                                              encRowsAdjMatrix.size());
            enc_t_user = cryptoContext1->EvalMult(enc_t_user, initRotsMasks.encMasks()[0]); 
            cryptoContext1->ModReduceInPlace(enc_t_user);
            enc_elements.push_back(cryptoContext1->EvalRotate(enc_t_user,-user));
        }
        auto enc_t = cryptoContext1->EvalAddMany(enc_elements);
        // o: Update output for all users in packed ciphertext: o <- t x u + o x (1-u)
        auto enc_t_mult_u = cryptoContext1->EvalMult(enc_t, enc_u); cryptoContext1->ModReduceInPlace(enc_t);
        auto enc_one_min_u = cryptoContext1->EvalAdd(encOnes,
                                                    cryptoContext1->EvalMult(enc_u, encNegOnes));
        // output <- t x u + o x (1-u)
        enc_output = cryptoContext1->EvalAdd(enc_t_mult_u, 
                                            cryptoContext1->EvalMult(enc_output,enc_one_min_u));
        // Update availability: 1-NotEqualZero(output)
        auto enc_output_reduced = evalNotEqualZero(enc_output,cryptoContext1,initNotEqualZero); 
        encUserAvailability = cryptoContext1->EvalAdd(encOnes,
                                                     cryptoContext1->EvalMult(enc_output_reduced, encNegOnes));
        processingTime = TOC(t); // End: Timer.
        std::cout << "Online part 3 - User availability & output update: " << processingTime << "ms" << std::endl;
        

        // Refresh after phase (3).
        //----------------------------------------------------------
        // refreshInPlace(enc_output,n,keyPair, cryptoContext1);
        if (packingMode == 1 || packingMode == 2) {
            // Refresh encrypted output vector.
            Plaintext plaintext; cryptoContext1->Decrypt(keyPair.secretKey, enc_output, &plaintext); 
            plaintext->SetLength(n); auto output = plaintext->GetPackedValue();
            std::cout << "Output vector: " << output << std::endl;
            enc_output = cryptoContext1->Encrypt(keyPair.publicKey,
                                                cryptoContext1->MakePackedPlaintext(output));       
            // For packing mode 2: Refresh & pack copies of user availability vector into single ciphertext.
            cryptoContext1->Decrypt(keyPair.secretKey, encUserAvailability, &plaintext); 
            plaintext->SetLength(n); auto userAvailability = plaintext->GetPackedValue();
            std::cout << "Availability vector: " << userAvailability << std::endl;
            // padded(availabilityVector) | padded(availabilityVector) | ... (2^ceil(log(n)) x 2^ceil(log(n)) times.          
            std::vector<int64_t> fullyPackedPlaintext;
            for (int copy=0; copy < slotsPadded*slotsPadded; ++copy) {
                for (int elem=0; elem < n; ++elem) { fullyPackedPlaintext.push_back(userAvailability[elem]); }
                for (int elem=0; elem<slotsPadded-n; ++elem){ fullyPackedPlaintext.push_back(0); }
            }
            encUserAvailability = cryptoContext1->Encrypt(keyPair.publicKey,
                                                         cryptoContext1->MakePackedPlaintext(fullyPackedPlaintext));            
        }
        else { return 1; }

    // End loop.
    }

    return 0;
}
