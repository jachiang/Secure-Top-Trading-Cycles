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
    int chosen_depth = 10;
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

    // depth 8.
    // userInputs.push_back({4, 1, 2, 3, 0, 5});
    // userInputs.push_back({4, 3, 2, 1, 0, 5});
    // userInputs.push_back({4, 1, 0, 2, 3, 5});
    // userInputs.push_back({1, 3, 4, 0, 2, 5});
    // userInputs.push_back({3, 1, 2, 0, 4, 5});
    // userInputs.push_back({3, 1, 2, 0, 4, 5});
    // userInputs.push_back({3, 1, 2, 0, 4, 5});
    
    // depth 9.
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

    // depth 10.
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
            encUserPref.push_back(cryptoContext->Encrypt(keyPair.publicKey,
                                  cryptoContext->MakePackedPlaintext(userPrefMatrix[row])));
            encUserPrefTransposed.push_back(cryptoContext->Encrypt(keyPair.publicKey,
                                            cryptoContext->MakePackedPlaintext(userPrefMatrixTransposed[row])));
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
    auto encUserPrefFullyPacked = cryptoContext->Encrypt(keyPair.publicKey,
                                  cryptoContext->MakePackedPlaintext(userPrefPacked));
    auto encUserPrefTransPacked = cryptoContext->Encrypt(keyPair.publicKey,
                                  cryptoContext->MakePackedPlaintext(userPrefTransposedPacked));

    // Offline: Init objects.
    // -----------------------------------------------------------------------
    InitRotsMasks initRotsMasks(cryptoContext,keyPair,n);
    InitNotEqualZero initNotEqualZero(cryptoContext,keyPair,n,userInputs.size()); // TODO: Over packed ciphertexts?
    InitPreserveLeadOne initPreserveLeadOne(cryptoContext,keyPair,n);
    
    // Offline: Initialize encrypted constants.
    // -----------------------------------------------------------------------
    auto slotTotal = cryptoContext->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder()/2;
    // TODO: move to initialization class.
    std::vector<int64_t> zeros(slotTotal,0);
    std::vector<int64_t> range; for (int i=0; i<n; ++i) { range.push_back(i+1); } // TODO: not suited for packing.
    std::vector<int64_t> ones; std::vector<int64_t> negOnes;
    for (int iter=0; iter<slotsPadded*slotsPadded; iter++){ 
        for (int elem = 0; elem <n; elem++){ ones.push_back(1); negOnes.push_back(-1); }
        for (int elem = 0; elem <slotsPadded-n; elem++){ ones.push_back(0); negOnes.push_back(0); }
     }
    auto encZeros = cryptoContext->Encrypt(keyPair.publicKey,
                                           cryptoContext->MakePackedPlaintext(zeros));
    auto encOnes = cryptoContext->Encrypt(keyPair.publicKey,
                                          cryptoContext->MakePackedPlaintext(ones));
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
    // For mode 2: Repeated availability vector in single ciphertext.
    auto encUserAvailability = encOnes; 
    // TODO: need gaps, and nxn repetition of availability vector.
    // auto encUserAvailabilityPacked = encOnes; 

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

        // 0: No packing 
        // 1: Row/Col matrix packing
        // 2: Full matrix packing
        int packingModeStep1 = 2;
        int packingModeStep2 = 2;
        // int packingModeStep3 = 1;

        // Begin: Timer.
        TIC(t);
        if (packingModeStep1 == 0){
             // TODO: 
        }
        else if (packingModeStep1 == 1){
            for (int user = 0; user < n; ++user){
                auto encUserAvailablePref = evalMatrixVecMult(encUserPrefList[user], encUserAvailability,
                                                              cryptoContext, initRotsMasks); // 2 mult-depth
                auto encUserFirstAvailablePref = evalPreserveLeadOne(encUserAvailablePref, cryptoContext, initPreserveLeadOne); // 2+log(n) mult-depth
                // TODO: extract rows in refresh.
                encRowsAdjMatrix.push_back(evalMatrixVecMult(encUserPrefTransposedList[user], encUserFirstAvailablePref,
                                                             cryptoContext, initRotsMasks)); // 2 mult-depth
            }
        }
        else {
            assert(packingModeStep1 == 2);
            // Multiplication: packed user preference matrices * availability matrix (repeated) 
            // padded(user1row1)        | padded(user1row2)        | padded(user1row3) ...    | padded(user2row1)        | padded(user2row2)        ...
            // padded(userAvailability) | padded(userAvailability) | padded(userAvailability) | padded(userAvailability) | padded(userAvailability) ...
            auto encResMult = cryptoContext->EvalMult(encUserAvailability,encUserPrefFullyPacked); 
            auto encResInnerProd = evalPrefixAdd(encResMult,slotsPadded,cryptoContext);      

            // Extract resulting slot values.
            std::vector<std::vector<Ciphertext<DCRTPoly>>> encAvailPrefElems;
            for (int user = 0; user < n; user++){
                std::vector<Ciphertext<DCRTPoly>> encElemsRow;
                for (int slot = 0; slot < n; slot++){
                    auto encMaskedElem = cryptoContext->EvalMult(encResInnerProd,initRotsMasks.encMasksFullyPacked()[user*n+slot]);
                    cryptoContext->ModReduceInPlace(encMaskedElem);
                    encElemsRow.push_back(cryptoContext->EvalRotate(encMaskedElem,((user*slotsPadded+slot)*slotsPadded)));
                }
                encAvailPrefElems.push_back(encElemsRow);
            }
            // Repack into single ciphertext (for later matrix multiplication).
            // padded(user1AvailPref) ... | padded(user2AvailPref)| ... | padded(user3AvailPref) | ... (each repeated 2^ceil(log(n)) times)
            std::vector<Ciphertext<DCRTPoly>> encAvailPrefContainer;
            for (int user = 0; user < n; user++){
                for (int slot=0; slot < n; slot++){
                    encAvailPrefContainer.push_back(cryptoContext->EvalRotate(
                                                    encAvailPrefElems[user][slot],-(user*slotsPadded*slotsPadded+slot)));
                }
            }
            auto encUniqueAvailRows = cryptoContext->EvalAddMany(encAvailPrefContainer); 
            auto encRepAvailRows = encUniqueAvailRows;
            for (int k = 0; k < k_ceil; k++) {  // Copy log(n) times.
                auto encRepAvailRows_ = cryptoContext->EvalRotate(encRepAvailRows,-std::pow(2,k)*slotsPadded);
                encRepAvailRows = cryptoContext->EvalAdd(encRepAvailRows,encRepAvailRows_); 
            } 
            // Preserve leading One.
            auto encFirstAvailPrefPacked = evalPreserveLeadOne(encRepAvailRows, cryptoContext, initPreserveLeadOne); // 2+log(n) mult-depth

            // Preferred available  * Transposed Preference Matrices
            auto encResMult2 = cryptoContext->EvalMult(encFirstAvailPrefPacked,encUserPrefTransPacked); 
            encAdjMatrixPacked = evalPrefixAdd(encResMult2,slotsPadded,cryptoContext);
        }

        // End: Timer.
        processingTime = TOC(t);
        std::cout << "Online part 1 - Adjacency matrix update time: " << processingTime << "ms" << std::endl;

        // Refresh adjacency matrix.
        std::cout << "Adjacency Matrix: " << std::endl;
        std::vector<std::vector<Ciphertext<DCRTPoly>>> encAdjMatrixElems;
        if (packingModeStep1 == 0) {
            // TODO
        }
        else if (packingModeStep1 == 1){
            for (int row=0; row < n; ++row){
                // Print adjacence matrix.
                Plaintext plaintext;
                cryptoContext->Decrypt(keyPair.secretKey, encRowsAdjMatrix[row], &plaintext); 
                plaintext->SetLength(n); auto payload = plaintext->GetPackedValue();
                std::cout << payload << std::endl;
                // Refresh row-wise encrypted ciphertexts.
                encAdjMatrixElems.push_back(refreshElems(encRowsAdjMatrix[row],n,keyPair,cryptoContext));
                // This ciphertext is used in part (3).
                refreshInPlace(encRowsAdjMatrix[row],n,keyPair,cryptoContext); 
            }
        }
        else {
            assert(packingModeStep1 == 2);
            // Refresh and repack ciphertexts.
            Plaintext plaintext;
            cryptoContext->Decrypt(keyPair.secretKey,encAdjMatrixPacked,&plaintext); 
            plaintext->SetLength(slotsPadded*slotsPadded*n); auto payload = plaintext->GetPackedValue();
            // WIP: For step 3, refresh adjacency matrix in row form.
            encRowsAdjMatrix.clear();
            for (int row=0; row < n; ++row){
                std::vector<int64_t> adjMatrixRow;
                std::vector<Ciphertext<DCRTPoly>> encAdjMatrixRow;
                for (int col=0; col < n; ++col){
                    int pos = row*slotsPadded*slotsPadded + col*slotsPadded;
                    adjMatrixRow.push_back(payload[pos]);
                    std::vector<int64_t> elementFirst(n,0); elementFirst[0] = payload[pos];
                    encAdjMatrixRow.push_back(cryptoContext->Encrypt(keyPair.publicKey,
                                              cryptoContext->MakePackedPlaintext(elementFirst)));
                }
                encAdjMatrixElems.push_back(encAdjMatrixRow);
                // Print adjacency matrix.
                std::cout << adjMatrixRow << std::endl; 
                // WIP: for step 3, refresh adjacency matrix ciphertext in row form.
                encRowsAdjMatrix.push_back(cryptoContext->Encrypt(keyPair.publicKey,
                                           cryptoContext->MakePackedPlaintext(adjMatrixRow)));
            }
        }


        //----------------------------------------------------------
        // (2) Matrix exponentiation for cycle finding.
        //----------------------------------------------------------

        // Matrix exponentiation.
  
        // Cycle finding result [res_1, ..., res_n]; 
        // On cycle: res_i = 1. Not on cycle: res_i = 0.
        Ciphertext<DCRTPoly> enc_u; 

        if (packingModeStep2 == 0){
            // TODO: ...
        }
        else if (packingModeStep2 == 1){
            TIC(t);
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
            TIC(t);
            auto encMatrixExpElems = evalMatSqMul(encAdjMatrixElems,n,
                                                  2, // Packing mode.
                                                  cryptoContext,initRotsMasks,cryptoOpsLogger,keyPair); 
                                                  // TODO: keypair param for debugging, remove later.
            processingTime = TOC(t);
            std::cout << "Online part 2a - Matrix exponentiation: " << processingTime << " ms" << std::endl;
                                      
            // Refresh and repack ciphertexts
            // packed(col1) | packed(col2) | packed(col3) ...     
            std::vector<int64_t> packedMatrix(slotsPadded*n,0);
            for (int row = 0; row < n; row++){
                std::vector<int64_t> matrixElemsRow;
                for (int col = 0; col < n; col++){
                    Plaintext plaintext;
                    cryptoContext->Decrypt(keyPair.secretKey,encMatrixExpElems[row][col],&plaintext); 
                    plaintext->SetLength(n); auto payload = plaintext->GetPackedValue();
                    int pos = col*slotsPadded + row;
                    packedMatrix[pos] = payload[0];
                }
            }
            auto encMatrixExpPacked = cryptoContext->Encrypt(keyPair.publicKey,
                                      cryptoContext->MakePackedPlaintext(packedMatrix));    

            // Extract computed cycles.
            TIC(t);
            // packed(col1) | packed(col2) | packed(col3) ...     
            // packed(ones) | packed(ones) | packed(ones) ...
            auto encResMult = cryptoContext->EvalMult(encMatrixExpPacked,encOnes); 
            auto encResInnerProd = evalPrefixAdd(encResMult,slotsPadded,cryptoContext);      
            // Not equal zero (SIMD operation).
            auto enc_u_unmasked = evalNotEqualZero(encResInnerProd,cryptoContext,initNotEqualZero); 
            processingTime = TOC(t);  
            std::cout << "Online part 2b - Cycle computation: " << processingTime << "ms" << std::endl;

            // Refresh and extract resulting slot values.
            // Generate fresh encryption of [ u_1 | u_2 | ... | u_n ... ] (replicated in packed fashion)
            Plaintext plaintext;
            cryptoContext->Decrypt(keyPair.secretKey,enc_u_unmasked,&plaintext); 
            plaintext->SetLength(n*slotsPadded); auto payload = plaintext->GetPackedValue();
            std::vector<int64_t> uElems;
            for (int user = 0; user < n; user++){
                uElems.push_back(payload[user*slotsPadded]); // TODO: Pack this for phase 3.
            }
            enc_u = cryptoContext->Encrypt(keyPair.publicKey,
                    cryptoContext->MakePackedPlaintext(uElems));    
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
        refreshInPlace(enc_output,n,keyPair, cryptoContext);
        refreshInPlace(encUserAvailability,n,keyPair, cryptoContext);
        if (packingModeStep1 == 0){
            // TODO: refresh individually encrypted ciphertexts.
        }
        else {
            assert(packingModeStep1 == 1 || packingModeStep1 == 2);
            // For packing mode 2: Pack copies of user availability vector into single ciphertext.
            // padded(availabilityVector) | padded(availabilityVector) | ... (2^ceil(log(n)) x 2^ceil(log(n)) times.          
            // For packing mode 1: Only first copy will be affect computation.
            Plaintext plaintext; cryptoContext->Decrypt(keyPair.secretKey, encUserAvailability, &plaintext); 
            plaintext->SetLength(n); auto userAvailability = plaintext->GetPackedValue();
            std::vector<int64_t> fullyPackedPlaintext;
            for (int copy=0; copy < slotsPadded*slotsPadded; ++copy) {
                for (int elem=0; elem < n; ++elem) { fullyPackedPlaintext.push_back(userAvailability[elem]); }
                for (int elem=0; elem<slotsPadded-n; ++elem){ fullyPackedPlaintext.push_back(0); }
            }
            encUserAvailability = cryptoContext->Encrypt(keyPair.publicKey,
                                                         cryptoContext->MakePackedPlaintext(fullyPackedPlaintext));            
        }
    // End loop.
    }

    return 0;
}
