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

    // benchmarking variables
    TimeVar t;
    double processingTime(0.0);

    // Crypto Parameters
    // # of evalMults = 3 (first 3) is used to support the multiplication of 7
    // ciphertexts, i.e., ceiling{log2{7}} Max depth is set to 3 (second 3) to
    // generate homomorphic evaluation multiplication keys for s^2 and s^3
    CCParams<CryptoContextBGVRNS> parameters;
    CCParams<CryptoContextBGVRNS> params2;
    CCParams<CryptoContextBGVRNS> params3;

    int chosen_ptxtmodulus = 65537;
    parameters.SetPlaintextModulus(chosen_ptxtmodulus);
    // p = 65537, depth = 13 -> "Please provide a q and a m satisfying: (q-1)/m is an integer. The values of primeModulus = 65537 and m = 131072 do not."
    // Fermats thm works for p = 786433, dep = 20.
    int chosen_depth = 9;
    parameters.SetMultiplicativeDepth(chosen_depth);
    parameters.SetMaxRelinSkDeg(3);

    CryptoContext<DCRTPoly> cryptoContext = GenCryptoContext(parameters);
    std::cout << "Ring dimension N: " << cryptoContext->GetRingDimension() << std::endl;

    cryptoContext->Enable(PKE);
    cryptoContext->Enable(KEYSWITCH);
    cryptoContext->Enable(LEVELEDSHE);
    cryptoContext->Enable(ADVANCEDSHE);

    std::cout << "Plaintext modulus p = " << cryptoContext->GetCryptoParameters()->GetPlaintextModulus() << std::endl;
    std::cout << "Cyclotomic order n = " << cryptoContext->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2
              << std::endl;
    // std::cout << "log2 q = "
    //           << log2(cryptoContext->GetCryptoParameters()->GetElementParams()->GetModulus().ConvertToDouble())
    //           << std::endl;

    // Initialize Public Key Containers
    KeyPair<DCRTPoly> keyPair;

    // Perform Key Generation Operation
    TIC(t);

    keyPair = cryptoContext->KeyGen();

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

    cryptoContext->EvalMultKeysGen(keyPair.secretKey);

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

    userInputs.push_back({4, 1, 2, 3, 0, 5});
    userInputs.push_back({4, 3, 2, 1, 0, 5});
    userInputs.push_back({4, 1, 0, 2, 3, 5});
    userInputs.push_back({1, 3, 4, 0, 2, 5});
    userInputs.push_back({3, 1, 2, 0, 4, 5});
    userInputs.push_back({3, 1, 2, 0, 4, 5});

    auto n = userInputs.size();

    // Offline: User preference as a permutation matrix.
    // Translate user preference into permutation matrices
    // and generate row-wise encryptions of each user matrix (packing mode 1),
    // or a single encryption of all rows of all user matrices (packing mode 2).

    std::vector<std::vector<Ciphertext<DCRTPoly>>> encUserPrefList;             // Packing mode: 1 (row/col-wise packing of each matrices)
    std::vector<std::vector<Ciphertext<DCRTPoly>>> encUserPrefTransposedList;   // Packing mode: 1 (row/col-wise packing of each matrices)
    std::vector<int64_t> userPrefFullyPacked;                                   // Packing mode: 2 (full packing of user all matrices)
    std::vector<int64_t> userPrefTransposedFullyPacked;                         // Packing mode: 2 (full packing of user all matrices)
    
    int k_ceil = std::ceil(std::log2(n));
    int slotsPadded = std::pow(2,k_ceil);

    for (int user=0; user < n; ++user) { 
        // Derive user preference-permutation matrix.
        std::vector<std::vector<int64_t>> userPrefMatrix;
        for (int col=0; col < n; ++col) {
            // For packing mode 1: row-vector of user permutation matrix.
            std::vector<int64_t> row(n,0); row[userInputs[user][col]] = 1;
            userPrefMatrix.push_back(row);
             // For packing mode 2: vector with all transposed user permutation matrices.
            for (auto elem: row) { userPrefFullyPacked.push_back(elem);}
            for (int i=0; i<slotsPadded-n; ++i){ userPrefFullyPacked.push_back(0); } 
        }   
        // Transpose user preference-permutation matrix.
        std::vector<std::vector<int64_t>> userPrefMatrixTransposed;
        std::vector<int64_t> zeroRow(n,0);
        // For packing mode 1: row-vector of transposed user permutation matrix.
        for (int row=0; row < n; ++row) { userPrefMatrixTransposed.push_back(zeroRow); }
        for (int row=0; row < n; ++row) { for (int col=0; col < n; ++col) {
            // if (userPrefMatrix[row][col] == 1) { userPrefMatrixTransposed[col][row] = 1; break; }
            userPrefMatrixTransposed[col][row] = userPrefMatrix[row][col];
        }}
        // For packing mode 2: vector with all transposed user permutation matrices.
        for (int row=0; row < n; ++row) { 
            for (int col=0; col < n; ++col) {
                userPrefTransposedFullyPacked.push_back(userPrefMatrixTransposed[row][col]);
            }
            for (int i=0; i<slotsPadded-n; ++i){ userPrefTransposedFullyPacked.push_back(0); }
        }    
        // For packing mode 1: row-wise encryption of each permutation matrix.
        std::vector<Ciphertext<DCRTPoly>> encUserPref;
        std::vector<Ciphertext<DCRTPoly>> encUserPrefTransposed;
        for (int j=0; j<n ; ++j){ 
            encUserPref.push_back(cryptoContext->Encrypt(keyPair.publicKey,
                                  cryptoContext->MakePackedPlaintext(userPrefMatrix[j])));
            encUserPrefTransposed.push_back(cryptoContext->Encrypt(keyPair.publicKey,
                                            cryptoContext->MakePackedPlaintext(userPrefMatrixTransposed[j])));
        }
        encUserPrefList.push_back(encUserPref);
        encUserPrefTransposedList.push_back(encUserPrefTransposed);        
    }
    // For packing mode 2: encryption of all permutation matrices in one ciphertext.
    auto encUserPrefFullyPacked = cryptoContext->Encrypt(keyPair.publicKey,
                                  cryptoContext->MakePackedPlaintext(userPrefFullyPacked));
    auto encUserPrefTransposedFullyPacked = cryptoContext->Encrypt(keyPair.publicKey,
                                            cryptoContext->MakePackedPlaintext(userPrefTransposedFullyPacked));


    // Server Offline: Precompute encryptions of constants, initialize rotation keys.
    InitRotsMasks initRotsMasks(cryptoContext,keyPair,n);
    InitNotEqualZero initNotEqualZero(cryptoContext,keyPair,cryptoContext->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder()/2,
                                                            userInputs.size()); // (..., slots, range)
    InitPreserveLeadOne initPreserveLeadOne(cryptoContext,keyPair,n);
    
    // Server Offline: Initialize enc(user-availability), enc(ones), enc([0:n])
    // TODO: move to initialization class.
    std::vector<int64_t> ones(n,1); std::vector<int64_t> zeros(n,0);
    std::vector<int64_t> negOnes(n,cryptoContext->GetCryptoParameters()->GetPlaintextModulus()-1);
    std::vector<int64_t> range; for (int i=0; i<n; ++i) { range.push_back(i+1); } // 1, 2, ..., n
    auto encOnes = cryptoContext->Encrypt(keyPair.publicKey,
                                          cryptoContext->MakePackedPlaintext(ones));
    auto encZeros = cryptoContext->Encrypt(keyPair.publicKey,
                                           cryptoContext->MakePackedPlaintext(zeros));
    auto encNegOnes = cryptoContext->Encrypt(keyPair.publicKey,
                                             cryptoContext->MakePackedPlaintext(negOnes));
    auto encRange = cryptoContext->Encrypt(keyPair.publicKey,
                                           cryptoContext->MakePackedPlaintext(range));


    //==========================================================
    // Tests for inner product with full matrix packing.
    //==========================================================
               
    // std::vector<std::vector<int64_t>> testMatRows1;
    // std::vector<std::vector<int64_t>> testMatRows2;
    // testMatRows1.push_back({4, 1, 2, 3, 1});
    // testMatRows1.push_back({4, 3, 2, 1, 1});
    // testMatRows1.push_back({4, 1, 1, 2, 3});
    // testMatRows1.push_back({1, 3, 4, 1, 2});
    // testMatRows1.push_back({3, 1, 2, 1, 4});

    // testMatRows2.push_back({4, 1, 2, 3, 1});
    // testMatRows2.push_back({4, 3, 2, 1, 1});
    // testMatRows2.push_back({4, 1, 1, 2, 3});
    // testMatRows2.push_back({1, 3, 4, 1, 2});
    // testMatRows2.push_back({3, 1, 2, 1, 4});

    // testMatRows1.push_back({4, 1, 2, 3, 0, 1});
    // testMatRows1.push_back({4, 3, 2, 1, 0, 1});
    // testMatRows1.push_back({4, 1, 0, 2, 3, 1});
    // testMatRows1.push_back({1, 3, 4, 0, 2, 1});
    // testMatRows1.push_back({3, 1, 2, 0, 4, 1});
    // testMatRows1.push_back({3, 1, 2, 0, 4, 1});

    // testMatRows2.push_back({4, 1, 2, 3, 0, 1});
    // testMatRows2.push_back({4, 3, 2, 1, 0, 1});
    // testMatRows2.push_back({4, 1, 0, 2, 3, 1});
    // testMatRows2.push_back({1, 3, 4, 0, 2, 1});
    // testMatRows2.push_back({3, 1, 2, 0, 4, 1});
    // testMatRows2.push_back({3, 1, 2, 0, 4, 1});

    // auto dimMat = testMatRows1.size();
    // int k_ceil = std::ceil(std::log2(dimMat));
    // int slotsPadded =std::pow(2, k_ceil);

    // // Fully pack both matrices into individual ciphertexts.
    // std::vector<int64_t> packedRows;
    // std::vector<int64_t> packedCols;
    // for (int row = 0; row < dimMat; row++){
    //     for (int col = 0; col < dimMat; col++){ 
    //         for (int i = 0; i < dimMat; i++) {
    //             packedRows.push_back(testMatRows1[row][i]);
    //             packedCols.push_back(testMatRows2[i][col]);
    //         }
    //         for (int i = 0; i < slotsPadded-dimMat; i++) {
    //             packedRows.push_back(0);
    //             packedCols.push_back(0);
    //         }        
    //     }
    // }
    // auto testEncMat1 = cryptoContext->Encrypt(keyPair.publicKey,cryptoContext->MakePackedPlaintext(packedRows));
    // auto testEncMat2 = cryptoContext->Encrypt(keyPair.publicKey,cryptoContext->MakePackedPlaintext(packedCols));

    // // Matrix multiplication; single inner product.
    // auto encResMult= cryptoContext->EvalMult(testEncMat1, testEncMat2);
    // auto encResInnerProd = evalPrefixAdd(encResMult,slotsPadded,cryptoContext);

    // // Extract individual ciphertexts and rotate to first slot.
    // std::vector<std::vector<Ciphertext<DCRTPoly>>> encElemsMat;
    // for (int row = 0; row < dimMat; row++){
    //     std::vector<Ciphertext<DCRTPoly>> encElemsRow;
    //     for (int col = 0; col < dimMat; col++){
    //         auto encMaskedElem = cryptoContext->EvalMult(encResInnerProd, initRotsMasks.encMasksFullyPacked()[row*dimMat+col]);
    //         cryptoContext->ModReduceInPlace(encMaskedElem);
    //         encElemsRow.push_back(cryptoContext->EvalRotate(encMaskedElem,(row*dimMat*slotsPadded+col*slotsPadded)));
    //     }
    //     encElemsMat.push_back(encElemsRow);
    // }
    // printEncMatElems(encElemsMat,cryptoContext,keyPair);

    // return 0;


    //==========================================================
    // Tests for matrix exponentiation.
    //==========================================================
               
    // std::vector<std::vector<int64_t>> testMatRows;
    // testMatRows.push_back({1, 0, 1, 0, 0});
    // testMatRows.push_back({0, 1, 0, 0, 1});
    // testMatRows.push_back({1, 0, 1, 0, 0});
    // testMatRows.push_back({0, 1, 0, 0, 1});
    // testMatRows.push_back({1, 0, 1, 0, 0});

    // std::vector<Ciphertext<DCRTPoly>> testEncMatRows;
    // for (int i = 0; i < testMatRows.size(); i++){
    // testEncMatRows.push_back(cryptoContext->Encrypt(keyPair.publicKey,
    //                          cryptoContext->MakePackedPlaintext(testMatRows[i])));
    // }
    // Test 1:
    // auto encMatsResElems = evalMatSquarings(testEncMatRows,3,cryptoContext,initRotsMasks, keyPair);
    // auto encMatResElems = evalMatrixMul2Pow(encMatsResElems,cryptoContext,initRotsMasks);
    // auto encMatResRows = encElem2Rows(encMatResElems,cryptoContext,initRotsMasks);

    // Test 2:
    // auto encMatResElems = evalMatSqMul(testEncMatRows, 15, cryptoContext,initRotsMasks, keyPair);
    // auto encMatResRows = encElem2Rows(encMatResElems,cryptoContext,initRotsMasks);

    // std::cout << "Test Matrix: " << std::endl;
    // printEncMatRows(encMatResRows,cryptoContext,keyPair);

    // return 0;

    //==========================================================
    // Runtime test for single operations.
    //==========================================================

    // auto res = encOnes;
    // for (int i = 0; i<chosen_depth; ++i){
    //     TIC(t); res = cryptoContext->EvalRotate(res,std::pow(2,3));
    //     processingTime = TOC(t);
    //     std::cout << "Rotation: " << processingTime << "ms" << std::endl;
    //     TIC(t); res = cryptoContext->EvalMult(res,res);
    //     cryptoContext->ModReduceInPlace(res);
    //     processingTime = TOC(t);
    //     std::cout << "Multiplication + mod reduction: " << processingTime << "ms" << std::endl;

    // }
    // res = encOnes;
    // for (int i = 0; i<chosen_depth; ++i){
    //     TIC(t);
    //     res = cryptoContext->EvalAdd(res,res);
    //     processingTime = TOC(t);
    //     std::cout << "Addition: " << processingTime << "ms" << std::endl;
    // }
    // for (int i = 0; i<n; ++i){
    //     TIC(t);
    //     res = cryptoContext->EvalRotate(res,i);
    //     processingTime = TOC(t);
    //     std::cout << "Rotation: " << processingTime << "ms" << std::endl;
    // }
    // return 0;

    //==========================================================
    // Ciphertext packing test.
    //==========================================================

    // std::vector<int64_t> testInputs1;
    // std::vector<int64_t> testInputs2;

    // auto slots = cryptoContext->GetRingDimension();

    // for (int i = 0; i < slots; ++i) {
    //     testInputs1.push_back(2);
    //     testInputs2.push_back(3);
    // }
    // auto encTestInputs1 = cryptoContext->Encrypt(keyPair.publicKey,
    //                                              cryptoContext->MakePackedPlaintext(testInputs1));
    // auto encTestInputs2 = cryptoContext->Encrypt(keyPair.publicKey,
    //                                              cryptoContext->MakePackedPlaintext(testInputs2));
    
    // auto resTest = cryptoContext->EvalMult(encTestInputs1, encTestInputs2); 
    // cryptoContext->ModReduceInPlace(resTest);
    
    // Plaintext testPlaintext;
    // cryptoContext->Decrypt(keyPair.secretKey, resTest, &testPlaintext); 
    // testPlaintext->SetLength(slots); auto testPayload = testPlaintext->GetPackedValue();
    // std::cout << testPayload[slots-1] << std::endl;

    // return 0;


    //==========================================================
    // Online phase.
    //==========================================================

    // Init output ciphertext.
    // auto enc_output = encNegOnes; // -1 output => not on cycle.
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

        // 0: No packing 
        // 1: Row/Col matrix packing
        // 2: Full matrix packing
        int packingModeStep1 = 1;

        // Begin: Timer.
        TIC(t);

        if (packingModeStep1 == 0){
             // TODO: 
        }
        else if (packingModeStep1 == 1){
            // Generate adjacency matrix for all users.
            // TODO: This step runs over row/col-wise packed ciphertexts.
            for (int user = 0; user < n; ++user){
                auto encUserAvailablePref = evalMatrixVecMult(encUserPrefList[user], encUserAvailability,
                                                              cryptoContext, initRotsMasks); // 2 mult-depth
                auto encUserFirstAvailablePref = evalPreserveLeadOne(encUserAvailablePref, cryptoContext, initPreserveLeadOne); // 2+log(n) mult-depth
                encRowsAdjMatrix.push_back(evalMatrixVecMult(encUserPrefTransposedList[user], encUserFirstAvailablePref,
                                                             cryptoContext, initRotsMasks)); // 2 mult-depth
            }
        }
        else {
            assert(packingModeStep1 == 2);
            // TODO: Availability matrix * Preference Matrices (row_11 | ... | row_1n) | (row_21 | ... | row_2n)
            // padded(available) | padded(available)| ... (nxn repetitions)
            // padded(row1)      | padded(row2)     | padded(row3)

            // TODO: EvalPreserveLeadOne over padded(preferredAvailable1) | padded(preferredAvailable2) | ...
            // TODO: PreferredAvailable matrix * Preference Matrices
        }

        std::cout << "Adjacency Matrix: " << std::endl;
        for (int row=0; row < n; ++row){
            Plaintext plaintext;
            cryptoContext->Decrypt(keyPair.secretKey, encRowsAdjMatrix[row], &plaintext); 
            plaintext->SetLength(n); auto payload = plaintext->GetPackedValue();
            std::cout << payload << std::endl;
        }

        // End: Timer.
        processingTime = TOC(t);
        std::cout << "Online part 1 - Adjacency matrix update time: " << processingTime << "ms" << std::endl;

        // Refresh ciphertexts.
        std::vector<std::vector<Ciphertext<DCRTPoly>>> encAdjMatrixElems;
        for (int row=0; row < n; ++row){
            encAdjMatrixElems.push_back(refreshElems(encRowsAdjMatrix[row],n,keyPair,cryptoContext));
            refreshInPlace(encRowsAdjMatrix[row],n,keyPair,cryptoContext); // This ciphertext is used in part (3).
        }

        //----------------------------------------------------------
        // (2) Matrix exponentiation for cycle finding.
        //----------------------------------------------------------

       // 0: No packing 
        // 1: Row/Col matrix packing
        // 2: Full matrix packing
        int packingModeStep2 = 1;

        // Matrix exponentiation.
        // Begin: Timer.
        TIC(t);

        // Cycle finding result [res_1, ..., res_n]; 
        // On cycle: res_i = 1. Not on cycle: res_i = 0.
        Ciphertext<DCRTPoly> enc_u; 

        if (packingModeStep2 == 0){
            // TODO: ...
        }
        else if (packingModeStep2 == 1){
            // No packing (mode 0), Row/Col matrix packing (mode 1), Full matrix packing (mode 2) 
            auto encMatrixExpElems = evalMatSqMul(encAdjMatrixElems,n,
                                                1, // Packing mode.
                                                cryptoContext,initRotsMasks,cryptoOpsLogger,keyPair); // TODO: keypair param for debugging
            // End: Timer.
            processingTime = TOC(t);
            std::cout << "Online part 2a - Matrix exponentiation: " << processingTime << " ms" << std::endl;

            // TODO: perform this transform in refresh.
            auto encMatrixExp = encElem2Cols(encMatrixExpElems,cryptoContext,initRotsMasks,cryptoOpsLogger);

            // Refresh ciphertexts.
            for (int col=0; col < n; ++col){ refreshInPlace(encMatrixExp[col],n,keyPair,cryptoContext); } 
            printEncMatRows(encMatrixExp,cryptoContext,keyPair);

            // Extract computed cycles.
            // Begin: Timer.
            TIC(t);
            enc_u = evalVecMatrixMult(encOnes,encMatrixExp,cryptoContext,initRotsMasks,cryptoOpsLogger); 
            enc_u = evalNotEqualZero(enc_u,cryptoContext,initNotEqualZero); 
            // End: Timer.
            processingTime = TOC(t);  
            std::cout << "Online part 2b - Cycle computation: " << processingTime << "ms" << std::endl;

            // Refresh ciphertexts.
            refreshInPlace(enc_u,n,keyPair, cryptoContext);
        }
        else {
            assert(packingModeStep2 == 2);
            auto encMatrixExpElems = evalMatSqMul(encAdjMatrixElems,n,
                                                2, // Packing mode.
                                                cryptoContext,initRotsMasks,cryptoOpsLogger,keyPair); // TODO: keypair param for debugging, remove later.
            // TODO: refresh ciphertexts
            std::vector<int64_t> packedMatrix;
            for (int row = 0; row < n; row++){
                for (int col = 0; col < n; col++){
                    // payload<-decrypt mat[row][col]
                    // packedMatrix.push_back(payload);
                }
            }
                // TODO: Pack exp-matrix
                // packed(col1) | packed(col2) | packed(col3) ...     
                // packed(ones) | packed(ones) | packed(ones) ...

            // Extract computed cycles.
                // TODO: ones-vector x exp-matrix 
                // TODO: EvalNotEqualZero (SIMD)
                // TODO: Extract individual elements.
                
            // Refresh and encrypt packed
                // TODO: fresh encryption of [ u_1 | u_2 | ... | u_n ]
        }           

        //----------------------------------------------------------
        // (3) Update user availability and outputs.
        //----------------------------------------------------------
        
        // Begin: Timer.
        TIC(t);

        // Compute current preference index (t) for all users in packed ciphertext.
        std::vector<Ciphertext<DCRTPoly>> enc_elements;
        for (int user=0; user < n; ++user){ 
            // Note: encRowsAdjMatrix must be refreshed after (1)
            auto enc_t_user = cryptoContext->EvalInnerProduct(encRowsAdjMatrix[user], encRange, 
                                                              encRowsAdjMatrix.size());
            enc_t_user = cryptoContext->EvalMult(enc_t_user, initRotsMasks.encMasks()[0]); 
            cryptoContext->ModReduceInPlace(enc_t_user);
            enc_elements.push_back(cryptoContext->EvalRotate(enc_t_user,-user));
        }
        auto enc_t = cryptoContext->EvalAddMany(enc_elements);
        // o: Update output for all users in packed ciphertext: o <- t x u + o x (1-u)
        auto enc_t_mult_u = cryptoContext->EvalMult(enc_t, enc_u); cryptoContext->ModReduceInPlace(enc_t);
        auto enc_one_min_u = cryptoContext->EvalAdd(encOnes,
                                                    cryptoContext->EvalMult(enc_u, encNegOnes));
        // output <- t x u + o x (1-u)
        enc_output = cryptoContext->EvalAdd(enc_t_mult_u, 
                                            cryptoContext->EvalMult(enc_output,enc_one_min_u));
        // Update availability: 1-NotEqualZero(output)
        auto enc_output_reduced = evalNotEqualZero(enc_output,cryptoContext,initNotEqualZero); 
        encUserAvailability = cryptoContext->EvalAdd(encOnes,
                                                     cryptoContext->EvalMult(enc_output_reduced, encNegOnes));

        // End: Timer.
        processingTime = TOC(t);  
        std::cout << "Online part 3 - User availability & output update: " << processingTime << "ms" << std::endl;
        
        // Print availability vector.
        std::cout << "Availability vector: "; printEnc(encUserAvailability,n,cryptoContext,keyPair); 

        // Print output vector.
        std::cout << "Output vector: "; printEnc(enc_output,n,cryptoContext,keyPair); 

        // Refresh ciphertexts.
        // TODO: During refresh, produce (0) elem-wise, (1) row/col-wise or (2) matrix-wise packed ciphertext.
        if (packingModeStep1 == 0){
            // TODO: refresh individually encrypted ciphertexts.
        }
        else if (packingModeStep1 == 1){
            refreshInPlace(encUserAvailability,n,keyPair, cryptoContext);
            refreshInPlace(enc_output,n,keyPair, cryptoContext);
        }
        else {
            assert(packingModeStep1 == 2);
            // TODO: assume refreshed availability and output vectors.
            // padded(availabilityVector) | padded(availabilityVector) | ... (n x n times, each padding up to next 2^k)
        }
    // End loop.
    }

    return 0;
}
