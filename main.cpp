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
    // Test parameters.
    ////////////////////////////////////////////////////////////

    // 1: Row/Col-wise matrix packing (naive)
    // 2: Flat matrix packing (low mult depth).
    // 3: Diagonal matrix packing (low runtime).
    // int packingMode = 3;
    // TODO: Set user number n here.

    std::vector<std::vector<int64_t>> userInputs;

    // n = 5
    // userInputs.push_back({4, 1, 2, 3, 0});
    // userInputs.push_back({4, 3, 2, 1, 0});
    // userInputs.push_back({4, 1, 0, 2, 3});
    // userInputs.push_back({1, 3, 4, 0, 2});
    // userInputs.push_back({3, 1, 2, 0, 4});
    // int chosen_depth1 = 8;

    // n = 10
    // userInputs.push_back({4, 1, 2, 3, 0, 5, 6, 7, 8, 9});
    // userInputs.push_back({4, 3, 2, 1, 0, 5, 6, 7, 8, 9});
    // userInputs.push_back({4, 1, 0, 2, 3, 5, 6, 7, 8, 9});
    // userInputs.push_back({1, 3, 4, 0, 2, 5, 6, 7, 8, 9});
    // userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9});
    // userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9});
    // userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9});
    // userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9});
    // userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9});
    // userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9});
    // int chosen_depth1 = 9;

    // n = 15
    // userInputs.push_back({4, 1, 2, 3, 0, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14});
    // userInputs.push_back({4, 3, 2, 1, 0, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14});
    // userInputs.push_back({4, 1, 0, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14});
    // userInputs.push_back({1, 3, 4, 0, 2, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14});
    // userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14});
    // userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14});
    // userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14});
    // userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14});
    // userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14});
    // userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14});
    // userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14});
    // userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14});
    // userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14});
    // userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14});
    // userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14});
    // int chosen_depth1 = 9;

    // n = 20
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
    // int chosen_depth1 = 10;

    // n = 25
    userInputs.push_back({4, 1, 2, 3, 0, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
    userInputs.push_back({4, 3, 2, 1, 0, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
    userInputs.push_back({4, 1, 0, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
    userInputs.push_back({1, 3, 4, 0, 2, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
    userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
    userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
    userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
    userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
    userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
    userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
    userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
    userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
    userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
    userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
    userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
    userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
    userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
    userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
    userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
    userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
    userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
    userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
    userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
    userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
    userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
    int chosen_depth1 = 12;

    int n = userInputs.size();

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

    int chosen_ptxtmodulus = 65537;
    params1.SetPlaintextModulus(chosen_ptxtmodulus); 
    params1.SetMultiplicativeDepth(chosen_depth1); 
    params1.SetMaxRelinSkDeg(3);
    params1.SetSecurityLevel(lbcrypto::HEStd_128_classic);
    // TODO: Fermats thm works for p = 786433, dep = 20.

    // Global context.
    CryptoContext<DCRTPoly> cryptoContext1 = GenCryptoContext(params1);
    cryptoContext1->Enable(PKE); cryptoContext1->Enable(KEYSWITCH); cryptoContext1->Enable(LEVELEDSHE); cryptoContext1->Enable(ADVANCEDSHE);
    int slotTotal = cryptoContext1->GetRingDimension();
    std::cout << "Slots: " << slotTotal << std::endl;
    std::cout << "Ring dimension N: " << cryptoContext1->GetRingDimension() << std::endl;
    std::cout << "Plaintext modulus p = " << cryptoContext1->GetCryptoParameters()->GetPlaintextModulus() << std::endl;

    std::cout << "Cyclotomic order n = " << cryptoContext1->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() << std::endl;
    std::cout << "log2 q = "
              << log2(cryptoContext1->GetCryptoParameters()->GetElementParams()->GetModulus().ConvertToDouble())
              << std::endl;
    // std::cout << "q = "
    //           << cryptoContext1->GetCryptoParameters()->GetElementParams()->GetModulus().ConvertToDouble()
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



    //==========================================================
    // Runtime test for single operations.
    //==========================================================

    // std::vector<int32_t> rotIndicesTest;
    // rotIndicesTest.push_back(std::pow(2,5));
    // cryptoContext1->EvalRotateKeyGen(keyPair.secretKey, rotIndicesTest);

    // std::vector<int64_t> onesTest(slotTotal,1);
    // auto encOnesTest = cryptoContext1->Encrypt(keyPair.publicKey,
    //                                            cryptoContext1->MakePackedPlaintext(onesTest)); 

    // // auto encOnesTestElems = encOnesTest->GetElements();
    // // for (size_t i = 0; i < encOnesTestElems.size(); i++)
    // // 	std::cout << "Polynomial " << i << " " << encOnesTestElems[i] << std::endl;
    // for  (int r = 0; r<1; ++r){
    //     auto res = encOnesTest;
    //     for (int i = 0; i<chosen_depth1; ++i){
    //         // TIC(t); 
    //         // res = cryptoContext1->EvalRotate(res,std::pow(2,5));
    //         // processingTime = TOC(t);
    //         // std::cout << "Rotation: " << processingTime << "ms" << std::endl;
    //         // TIC(t); 
    //         res = cryptoContext1->EvalMult(res,res);
    //         auto resElems = res->GetElements();
    //         auto length = resElems.size();
    //         auto size = length * sizeof(resElems[0]);
    //         std::cout << "Elementsize: " << sizeof(resElems[0]) << std::endl;
    //         std::cout << "Bytes: " << size << std::endl;
    //         // processingTime = TOC(t);
    //         // std::cout << "Multiplication + mod reduction: " << processingTime << "ms" << std::endl;
    //     }
    // }

    // return 1;

    ////////////////////////////////////////////////////////////
    // Top Trading Cycle Algorithm.
    ////////////////////////////////////////////////////////////

    
    // Offline: Encryption of user preferences.
    // -----------------------------------------------------------------------
    // Represent user preferences as permutation matrices and their transpose.
    // Encrypt diagonals of permutation matrices.
    
    std::vector<std::vector<std::vector<int64_t>>> userPrefList;
    std::vector<std::vector<std::vector<int64_t>>> userPrefTransposedList;
 
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
        for (int row=0; row < n; ++row) { 
            for (int col=0; col < n; ++col) {
                userPrefMatrixTransposed[col][row] = userPrefMatrix[row][col]; 
            }
        }
        userPrefTransposedList.push_back(userPrefMatrixTransposed);
    }
    // Build diagonals of pref permutation matrix diagonals.
    std::vector<std::vector<std::vector<int64_t>>> usersPrefMatrixDiagonals;
    std::vector<std::vector<std::vector<int64_t>>> usersPrefMatrixTransposedDiagonals;
    for (int user=0; user<n ; ++user){ 
        usersPrefMatrixDiagonals.push_back(matrixDiagonals(userPrefList[user]));
        usersPrefMatrixTransposedDiagonals.push_back(matrixDiagonals(userPrefTransposedList[user]));
    }
    // Encrypt diagonals of pref permutation matrix diagonals.
    std::vector<std::vector<Ciphertext<DCRTPoly>>> encUsersPrefMatrixDiagonals; 
    std::vector<std::vector<Ciphertext<DCRTPoly>>> encUsersPrefMatrixTransposedDiagonals;
    std::vector<int64_t> zeroVec(n,0);
    for (int user=0; user<n ; ++user){ 
        std::vector<Ciphertext<DCRTPoly>> usersPrefDiagonal;
        std::vector<Ciphertext<DCRTPoly>> usersPrefTransposedDiagonal;
        for (int l=0; l<n ; ++l){ 
            usersPrefDiagonal.push_back(cryptoContext1->Encrypt(keyPair.publicKey,
                                        cryptoContext1->MakePackedPlaintext(repFillSlots(usersPrefMatrixDiagonals[user][l],slotTotal))));
            usersPrefTransposedDiagonal.push_back(cryptoContext1->Encrypt(keyPair.publicKey,
                                                  cryptoContext1->MakePackedPlaintext(repFillSlots(usersPrefMatrixTransposedDiagonals[user][l],slotTotal))));
        }
        encUsersPrefMatrixDiagonals.push_back(usersPrefDiagonal);
        encUsersPrefMatrixTransposedDiagonals.push_back(usersPrefTransposedDiagonal);
    }

    // Offline: Init objects and encrypted constants.
    // -----------------------------------------------------------------------

    TIC(t);
    // Generate rotation keys.
    std::vector<int32_t> rotIndices;
    for (int i = 0; i <= n; i++) { rotIndices.push_back(-i); rotIndices.push_back(i); rotIndices.push_back(n*i); } // TODO: n*i required for matrix multiplication.
    cryptoContext1->EvalRotateKeyGen(keyPair.secretKey, rotIndices);
    cryptoContext1->EvalSumKeyGen(keyPair.secretKey);
    processingTime = TOC(t);
    std::cout << "Rotation key generation: " 
              << processingTime << " ms" << std::endl;
    
    TIC(t);
    InitNotEqualZero initNotEqualZero(cryptoContext1,keyPair,n,userInputs.size()); 
    InitPreserveLeadOne initPreserveLeadOne(cryptoContext1,keyPair,n);
    InitMatrixMult initMatrixMult(cryptoContext1,keyPair,n); // n as dimension of nxn adjacency matrix.

    std::vector<int64_t> zeros(slotTotal,0);
    std::vector<int64_t> ones(slotTotal,1); std::vector<int64_t> negOnes(slotTotal,-1);
    std::vector<int64_t> leadingOne(n,0); leadingOne[0] = 1;
    std::vector<int64_t> onesRow(n,1);
    std::vector<int64_t> range; for (int i=0; i<n; ++i) { range.push_back(i+1); }
    auto encZeros = cryptoContext1->Encrypt(keyPair.publicKey,
                                            cryptoContext1->MakePackedPlaintext(zeros));
    auto encOnes = cryptoContext1->Encrypt(keyPair.publicKey,
                                           cryptoContext1->MakePackedPlaintext(ones));
    auto encNegOnes = cryptoContext1->Encrypt(keyPair.publicKey,
                                              cryptoContext1->MakePackedPlaintext(negOnes));
    auto encLeadingOne = cryptoContext1->Encrypt(keyPair.publicKey,
                                                 cryptoContext1->MakePackedPlaintext(leadingOne));
    auto encRange = cryptoContext1->Encrypt(keyPair.publicKey,
                                            cryptoContext1->MakePackedPlaintext(range));
    auto encOnesRow = cryptoContext1->Encrypt(keyPair.publicKey,
                                              cryptoContext1->MakePackedPlaintext(onesRow));
    processingTime = TOC(t);
    std::cout << "Encryption of constants: " 
              << processingTime << " ms" << std::endl;


    //==========================================================
    // Online phase.
    //==========================================================

    // Global: Log crypto operations.
    CryptoOpsLogger cryptoOpsLogger;

    // Global: init output ciphertext.
    auto enc_output = encZeros; 
    Ciphertext<DCRTPoly> encUserAvailability;
    encUserAvailability = encOnes;

    // Main loop for cycle finding algorithm.
    for (int i = 0; i < n ; ++i)
    {
        std::cout << "-------------" << std::endl;
        std::cout << "Round ... " << i+1 << "/" << n << std::endl;
        std::cout << "-------------" << std::endl;

        //----------------------------------------------------------
        // (1) Update adjacency matix.
        //----------------------------------------------------------

        // Global: Row-wise packed adjacency matrix (test mode 1,2).
        std::vector<Ciphertext<DCRTPoly>> encRowsAdjMatrix;
        // Global: Flat packed adjacency matrix (test mode 3).
        Ciphertext<DCRTPoly> encAdjMatrixPacked;

        // Begin: Timer.
        TIC(t);
        for (int user = 0; user < n; ++user){
            auto encUserAvailablePref = evalDiagMatrixVecMult(encUsersPrefMatrixDiagonals[user], encUserAvailability,
                                                                cryptoContext1, cryptoOpsLogger); // 1 mult-depth                                                  
            auto encUserFirstAvailablePref = evalPreserveLeadOne(encUserAvailablePref, cryptoContext1, initPreserveLeadOne); // 2+log(n) mult-depth
            // Mask and replicate availability row left and right.
            encUserFirstAvailablePref = cryptoContext1->EvalMult(encUserFirstAvailablePref,encOnesRow); // 1 mult-depth
            std::vector<Ciphertext<DCRTPoly>> addContainer;
            addContainer.push_back(encUserFirstAvailablePref);
            addContainer.push_back(cryptoContext1->EvalRotate(encUserFirstAvailablePref,-n));
            addContainer.push_back(cryptoContext1->EvalRotate(encUserFirstAvailablePref,n));
            encUserFirstAvailablePref = cryptoContext1->EvalAddMany(addContainer);
            encRowsAdjMatrix.push_back(evalDiagMatrixVecMult(encUsersPrefMatrixTransposedDiagonals[user], encUserFirstAvailablePref,
                                                                cryptoContext1, cryptoOpsLogger)); // 1 mult-depth
        }
        processingTime = TOC(t); // End: Timer.
        std::cout << "Online part 1 - Adjacency matrix update time: " << processingTime << "ms" << std::endl;

        // Refresh after phase (1).
        //----------------------------------------------------------
        std::cout << "Adjacency Matrix: " << std::endl;
        // Global: Element-wise packed adjacency matrix input to matrix exponentiation in modes 1,2.
        std::vector<std::vector<Ciphertext<DCRTPoly>>> encAdjMatrixElems; 
        // Global: Flat packed adjacency matrix input to matrix exponentiation in modes 3.
        Ciphertext<DCRTPoly> encAdjMatrixFlat; 
        
        // (1) Refresh "encRowsAdjMatrix" as encrypted flat packed matrix.
        std::vector<std::vector<int64_t>> rowsAdjMatrix;
        for (int row=0; row < n; ++row){
            Plaintext plaintext;
            cryptoContext1->Decrypt(keyPair.secretKey, encRowsAdjMatrix[row], &plaintext); 
            plaintext->SetLength(n); auto payload = plaintext->GetPackedValue();
            // Print adjacence matrix.
            std::cout << payload << std::endl;
            rowsAdjMatrix.push_back(payload);
            // (2) Also refresh "encRowsAdjMatrix" in row form for phase (3).
            refreshInPlace(encRowsAdjMatrix[row],n,keyPair,cryptoContext1); 
        }
        std::vector<int64_t> flatMatrix(n*n,0);
        for (int row=0; row < n; ++row){
            for (int col=0; col < n; ++col){
                flatMatrix[row*n+col] = rowsAdjMatrix[row][col];
            }
        }
        encAdjMatrixFlat = cryptoContext1->Encrypt(keyPair.publicKey,
                                                    cryptoContext1->MakePackedPlaintext(repFillSlots(flatMatrix,slotTotal)));            

        //----------------------------------------------------------
        // (2) Cycle finding.
        //----------------------------------------------------------

        // 2a) Matrix exponentiation.
        //----------------------------------------------------------
        // Cycle finding result [r_1, ..., r_n]. On cycle, r_i = 1. Not on cycle: r_i = 0.
        std::vector<std::vector<Ciphertext<DCRTPoly>>> encMatrixExpElems; // Output in mode 1,2.
        Ciphertext<DCRTPoly> encMatrixExpFlat; // Output in mode 3.
        double processingTime2a(0.0);

        bool contFlag = true; int sqs = 1;
        while (contFlag) {
            int exp = 2 * sqs;
            if (exp >= n) { contFlag = false; }
            else { sqs = sqs + 1; }
        }

        // Perform squarings, refresh every 3 squarings (each of d-mult. depth)
        TIC(t);
        encMatrixExpFlat = encAdjMatrixFlat;
        int refreshInterval = std::floor(chosen_depth1/3);
        for (int i=1; i <= sqs; i++){
            encMatrixExpFlat = evalMatrixMult(cryptoContext1,encMatrixExpFlat,encMatrixExpFlat,initMatrixMult,keyPair); // TODO: keyPair for debugging only.
            if (i % refreshInterval == 0) {
                processingTime2a += TOC(t);
                refreshInPlace(encMatrixExpFlat,cryptoContext1->GetRingDimension(),keyPair,cryptoContext1);
                TIC(t);
                std::cout << "Refresh during squaring" << std::endl;
            }
        }
        processingTime2a += TOC(t);
        std::cout << "Online part 2a - Matrix exponentiation: " << processingTime2a << " ms" << std::endl;

        // Refresh after phase (2a)
        //----------------------------------------------------------
        std::vector<Ciphertext<DCRTPoly>> encMatrixExp; // Refreshed & repacked output in packing mode 1.
        Ciphertext<DCRTPoly> encMatrixExpPacked;        // Refreshed & repacked output in packing mode 2/3.
  
        std::vector<int64_t> packedMatrix(slotsPadded*n,0);
        Plaintext plaintext;
        cryptoContext1->Decrypt(keyPair.secretKey,encMatrixExpFlat,&plaintext); 
        plaintext->SetLength(n*n); auto payload = plaintext->GetPackedValue();
        for (int row = 0; row < n; row++){
            std::vector<int64_t> matrixElemsRow;
            for (int col = 0; col < n; col++){
                int pos = col*slotsPadded + row;
                packedMatrix[pos] = payload[row*n+col];
            }
        }
        encMatrixExpPacked = cryptoContext1->Encrypt(keyPair.publicKey,
                                cryptoContext1->MakePackedPlaintext(packedMatrix));   

        // 2b) Cycle computation.
        //----------------------------------------------------------
        Ciphertext<DCRTPoly> enc_u;             // Output in packing mode 1.
        Ciphertext<DCRTPoly> enc_u_unmasked;    // Output in packing mode 2.

        TIC(t); // Begin: Timer.

        auto encResMult = cryptoContext1->EvalMult(encMatrixExpPacked,encOnes); 
        auto encResInnerProd = evalPrefixAdd(encResMult,slotsPadded,cryptoContext1);      
        enc_u_unmasked = evalNotEqualZero(encResInnerProd,cryptoContext1,initNotEqualZero); 

        processingTime = TOC(t); // End: Timer.
        std::cout << "Online part 2b - Cycle computation: " << processingTime << "ms" << std::endl;

        // Refresh after phase (2b)
        //----------------------------------------------------------
        Plaintext plaintext2b;
        cryptoContext1->Decrypt(keyPair.secretKey,enc_u_unmasked,&plaintext2b); 
        plaintext2b->SetLength(n*slotsPadded); auto payload2b = plaintext2b->GetPackedValue();
        std::vector<int64_t> uElems;
        for (int user = 0; user < n; user++){
            uElems.push_back(payload2b[user*slotsPadded]);
        }
        enc_u = cryptoContext1->Encrypt(keyPair.publicKey,
                cryptoContext1->MakePackedPlaintext(uElems));    

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
            // enc_t_user = cryptoContext1->EvalMult(enc_t_user, initRotsMasks.encMasks()[0]); 
            enc_t_user = cryptoContext1->EvalMult(enc_t_user, encLeadingOne); 
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
 
        // (1) Refresh encrypted output vector.
        Plaintext plaintext3; cryptoContext1->Decrypt(keyPair.secretKey, enc_output, &plaintext3); 
        plaintext3->SetLength(n); auto output = plaintext3->GetPackedValue();
        std::cout << "Output vector: " << output << std::endl;
        enc_output = cryptoContext1->Encrypt(keyPair.publicKey,
                                            cryptoContext1->MakePackedPlaintext(output));       
        // (2) For test mode 3: Refresh & pack copies of user availability vector into single ciphertext.
        cryptoContext1->Decrypt(keyPair.secretKey, encUserAvailability, &plaintext3); 
        plaintext3->SetLength(n); auto userAvailability = plaintext3->GetPackedValue();
        std::cout << "Availability vector: " << userAvailability << std::endl;
        encUserAvailability = cryptoContext1->Encrypt(keyPair.publicKey,
                                                        cryptoContext1->MakePackedPlaintext(repFillSlots(userAvailability,slotTotal)));            


    // End loop.
    }

    return 0;
}
