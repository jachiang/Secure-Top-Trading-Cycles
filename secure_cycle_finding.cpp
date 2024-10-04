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

    // omp_set_max_active_levels(2);
    // omp_set_dynamic(omp_get_max_threads());
    std::cout << "Thread count: " << omp_get_max_threads() << std::endl;

    ////////////////////////////////////////////////////////////
    // Client inputs.
    ////////////////////////////////////////////////////////////

    // Uncomment chosen test vector.
    // int numParties = 5;
    int numParties = 10;
    // int numParties = 15;
    // int numParties = 20;
    // int numParties = 25;

    std::vector<std::vector<int64_t>> userInputs;
    int chosen_depth(0);

    if (numParties == 5) {
        userInputs.push_back({4, 1, 2, 3, 0});
        userInputs.push_back({4, 3, 2, 1, 0});
        userInputs.push_back({4, 1, 0, 2, 3});
        userInputs.push_back({1, 3, 4, 0, 2});
        userInputs.push_back({3, 1, 2, 0, 4});
        chosen_depth = 8;
    }
    else if (numParties == 10) {
        userInputs.push_back({4, 1, 2, 3, 0, 5, 6, 7, 8, 9});
        userInputs.push_back({4, 3, 2, 1, 0, 5, 6, 7, 8, 9});
        userInputs.push_back({4, 1, 0, 2, 3, 5, 6, 7, 8, 9});
        userInputs.push_back({1, 3, 4, 0, 2, 5, 6, 7, 8, 9});
        userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9});
        userInputs.push_back({0, 1, 2, 3, 4, 9, 6, 7, 8, 5});
        userInputs.push_back({0, 1, 2, 3, 4, 9, 8, 7, 6, 5});
        userInputs.push_back({0, 1, 2, 3, 4, 9, 6, 5, 7, 8});
        userInputs.push_back({0, 1, 2, 3, 4, 6, 8, 9, 5, 7});
        userInputs.push_back({0, 1, 2, 3, 4, 8, 6, 7, 5, 9});
        chosen_depth = 9;
    }
    else if (numParties == 15) {
        userInputs.push_back({4, 1, 2, 3, 0, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14});
        userInputs.push_back({4, 3, 2, 1, 0, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14});
        userInputs.push_back({4, 1, 0, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14});
        userInputs.push_back({1, 3, 4, 0, 2, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14});
        userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14});
        userInputs.push_back({0, 1, 2, 3, 4, 9, 6, 7, 8, 5, 10, 11, 12, 13, 14});
        userInputs.push_back({0, 1, 2, 3, 4, 9, 8, 7, 6, 5, 10, 11, 12, 13, 14});
        userInputs.push_back({0, 1, 2, 3, 4, 9, 6, 5, 7, 8, 10, 11, 12, 13, 14});
        userInputs.push_back({0, 1, 2, 3, 4, 6, 8, 9, 5, 7, 10, 11, 12, 13, 14});
        userInputs.push_back({0, 1, 2, 3, 4, 8, 6, 7, 5, 9, 10, 11, 12, 13, 14});
        userInputs.push_back({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 14, 11, 12, 13, 10});
        userInputs.push_back({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 14, 13, 12, 11, 10});
        userInputs.push_back({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 14, 11, 10, 12, 13});
        userInputs.push_back({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 14, 10, 12});
        userInputs.push_back({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 11, 12, 10, 14});
        chosen_depth = 9;
    }
    else if (numParties == 20) {
        userInputs.push_back({4, 1, 2, 3, 0, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19});
        userInputs.push_back({4, 3, 2, 1, 0, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19});
        userInputs.push_back({4, 1, 0, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19});
        userInputs.push_back({1, 3, 4, 0, 2, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19});
        userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19});
        userInputs.push_back({0, 1, 2, 3, 4, 9, 6, 7, 8, 5, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19});
        userInputs.push_back({0, 1, 2, 3, 4, 9, 8, 7, 6, 5, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19});
        userInputs.push_back({0, 1, 2, 3, 4, 9, 6, 5, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19});
        userInputs.push_back({0, 1, 2, 3, 4, 6, 8, 9, 5, 7, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19});
        userInputs.push_back({0, 1, 2, 3, 4, 8, 6, 7, 5, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19});
        userInputs.push_back({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 14, 11, 12, 13, 10, 15, 16, 17, 18, 19});
        userInputs.push_back({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 14, 13, 12, 11, 10, 15, 16, 17, 18, 19});
        userInputs.push_back({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 14, 11, 10, 12, 13, 15, 16, 17, 18, 19});
        userInputs.push_back({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 14, 10, 12, 15, 16, 17, 18, 19});
        userInputs.push_back({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 11, 12, 10, 14, 15, 16, 17, 18, 19});
        userInputs.push_back({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 19, 16, 17, 18, 15});
        userInputs.push_back({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 19, 18, 17, 16, 15});
        userInputs.push_back({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 19, 16, 15, 17, 18});
        userInputs.push_back({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 19, 15, 17});
        userInputs.push_back({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 18, 16, 17, 15, 19});
        chosen_depth = 10;
    }
    else if (numParties == 25) {
        userInputs.push_back({4, 1, 2, 3, 0, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
        userInputs.push_back({4, 3, 2, 1, 0, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
        userInputs.push_back({4, 1, 0, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
        userInputs.push_back({1, 3, 4, 0, 2, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
        userInputs.push_back({3, 1, 2, 0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
        userInputs.push_back({0, 1, 2, 3, 4, 9, 6, 7, 8, 5, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
        userInputs.push_back({0, 1, 2, 3, 4, 9, 8, 7, 6, 5, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
        userInputs.push_back({0, 1, 2, 3, 4, 9, 6, 5, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
        userInputs.push_back({0, 1, 2, 3, 4, 6, 8, 9, 5, 7, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
        userInputs.push_back({0, 1, 2, 3, 4, 8, 6, 7, 5, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
        userInputs.push_back({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 14, 11, 12, 13, 10, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
        userInputs.push_back({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 14, 13, 12, 11, 10, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
        userInputs.push_back({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 14, 11, 10, 12, 13, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
        userInputs.push_back({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 14, 10, 12, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
        userInputs.push_back({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 11, 12, 10, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
        userInputs.push_back({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 19, 16, 17, 18, 15, 20, 21, 22, 23, 24});
        userInputs.push_back({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 19, 18, 17, 16, 15, 20, 21, 22, 23, 24});
        userInputs.push_back({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 19, 16, 15, 17, 18, 20, 21, 22, 23, 24});
        userInputs.push_back({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 19, 15, 17, 20, 21, 22, 23, 24});
        userInputs.push_back({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 18, 16, 17, 15, 19, 20, 21, 22, 23, 24});
        userInputs.push_back({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 24, 21, 22, 23, 20});
        userInputs.push_back({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 24, 23, 22, 21, 20});
        userInputs.push_back({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 24, 21, 20, 22, 23});
        userInputs.push_back({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21, 23, 24, 20, 22});
        userInputs.push_back({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 23, 21, 22, 20, 24});
        chosen_depth = 10;
    }
    else { return 1; }

    int n = numParties;

    TimeVar t;
    double runtimePhase(0.0);

    ////////////////////////////////////////////////////////////
    // Set-up of BGV parameters
    ////////////////////////////////////////////////////////////

    // Crypto Parameters
    // # of evalMults = 3 (first 3) is used to support the multiplication of 7
    // ciphertexts, i.e., ceiling{log2{7}} Max depth is set to 3 (second 3) to
    // generate homomorphic evaluation multiplication keys for s^2 and s^3
    CCParams<CryptoContextBGVRNS> params1, params2a, params2b, params3;

    int chosen_ptxtmodulus = 65537;
    params1.SetPlaintextModulus(chosen_ptxtmodulus);
    params1.SetMultiplicativeDepth(chosen_depth);
    params1.SetMaxRelinSkDeg(3);
    params1.SetSecurityLevel(lbcrypto::HEStd_128_classic);

    CryptoContext<DCRTPoly> cc = GenCryptoContext(params1);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);

    int slotTotal = cc->GetRingDimension();
    std::cout << "Ciphertext slots: " << slotTotal << std::endl;
    std::cout << "Ring dimension N: " << cc->GetRingDimension() << std::endl;
    std::cout << "Plaintext modulus p = " << cc->GetCryptoParameters()->GetPlaintextModulus() << std::endl;

    std::cout << "Cyclotomic order n = " << cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder()/2 << std::endl;
    std::cout << "log2 q = "
              << log2(cc->GetCryptoParameters()->GetElementParams()->GetModulus().ConvertToDouble())
              << std::endl;

    // Single encryption/decryption key is generated for simplicity.
    // BGV keys are additive; for fixed parameters (p, q, N), the number of decryption keys
    // does not affect runtimes of ciphertext arithmetic or ciphertext size.
    KeyPair<DCRTPoly> keyPair;

    TIC(t);
    keyPair = cc->KeyGen();
    runtimePhase = TOC(t);

    std::cout << "Key generation time: " << runtimePhase << "ms" << std::endl;

    if (!keyPair.good()) {
        std::cout << "Key generation failed!" << std::endl;
        exit(1);
    }

    std::cout << "Running key generation for homomorphic multiplication "
                 "evaluation keys..."
              << std::endl;

    TIC(t);
    cc->EvalMultKeysGen(keyPair.secretKey);
    runtimePhase = TOC(t);

    std::cout << "Key generation time for homomorphic multiplication evaluation keys: " << runtimePhase << "ms"
              << std::endl;


    ////////////////////////////////////////////////////////////
    // Top Trading Cycle Algorithm.
    ////////////////////////////////////////////////////////////

    // Offline: Init objects and encrypted constants.
    // -----------------------------------------------------------------------

    TIC(t);
    // Generate rotation keys.
    std::vector<int32_t> rotIndices;
    for (int i = 0; i <= n; i++) { rotIndices.push_back(-i); rotIndices.push_back(i); rotIndices.push_back(n*i); }
    cc->EvalRotateKeyGen(keyPair.secretKey, rotIndices);
    cc->EvalSumKeyGen(keyPair.secretKey);
    runtimePhase = TOC(t);
    std::cout << "Rotation key generation: "
              << runtimePhase << " ms" << std::endl;

    TIC(t);
    InitNotEqualZero initNotEqualZero(cc,keyPair,n,userInputs.size());
    InitPreserveLeadOne initPreserveLeadOne(cc,keyPair,n);
    InitMatrixMult initMatrixMult(cc,keyPair,n); // n in of nxn matrix.

    std::vector<int64_t> zeros(slotTotal,0);
    std::vector<int64_t> ones(slotTotal,1); std::vector<int64_t> negOnes(slotTotal,-1);
    std::vector<int64_t> leadingOne(n,0); leadingOne[0] = 1;
    std::vector<int64_t> onesRow(n,1);
    std::vector<int64_t> range; for (int i=0; i<n; ++i) { range.push_back(i+1); }
    auto encZeros = cc->Encrypt(keyPair.publicKey,
                                            cc->MakePackedPlaintext(zeros));
    auto encOnes = cc->Encrypt(keyPair.publicKey,
                                           cc->MakePackedPlaintext(ones));
    auto encNegOnes = cc->Encrypt(keyPair.publicKey,
                                              cc->MakePackedPlaintext(negOnes));
    auto encLeadingOne = cc->Encrypt(keyPair.publicKey,
                                                 cc->MakePackedPlaintext(leadingOne));
    auto encRange = cc->Encrypt(keyPair.publicKey,
                                            cc->MakePackedPlaintext(range));
    auto encOnesRow = cc->Encrypt(keyPair.publicKey,
                                              cc->MakePackedPlaintext(onesRow));
    runtimePhase = TOC(t);
    std::cout << "Encryption of constants: "
              << runtimePhase << " ms" << std::endl;


    // Online: Encryption of user preferences.
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
            usersPrefDiagonal.push_back(cc->Encrypt(keyPair.publicKey,
                                        cc->MakePackedPlaintext(repFillSlots(usersPrefMatrixDiagonals[user][l],slotTotal))));
            usersPrefTransposedDiagonal.push_back(cc->Encrypt(keyPair.publicKey,
                                                  cc->MakePackedPlaintext(repFillSlots(usersPrefMatrixTransposedDiagonals[user][l],slotTotal))));
        }
        encUsersPrefMatrixDiagonals.push_back(usersPrefDiagonal);
        encUsersPrefMatrixTransposedDiagonals.push_back(usersPrefTransposedDiagonal);
    }


    // Online: Top Trading Cycle
    // -----------------------------------------------------------------------

    // Log runtime.
    double runtimePhase1Total(0.0);
    double runtimePhase2aTotal(0.0);
    double runtimePhase2bTotal(0.0);
    double runtimePhase3Total(0.0);

    // Initialize availability and output variables.
    Ciphertext<DCRTPoly> encUserAvailability;
    encUserAvailability = encOnes;
    auto enc_output = encZeros;

    // Main loop for cycle finding algorithm.
    for (int i = 0; i < n ; ++i)
    {
        std::cout << "--------------" << std::endl;
        std::cout << "Round ... " << i+1 << "/" << n << std::endl;
        std::cout << "--------------" << std::endl;

        //----------------------------------------------------------
        // (1) Update adjacency matix.
        //----------------------------------------------------------

        std::vector<Ciphertext<DCRTPoly>> encRowsAdjMatrix;
        encRowsAdjMatrix.resize(n);
        Ciphertext<DCRTPoly> encAdjMatrixPacked;
        double runtimePhase1(0.0);

        TIC(t);
        #pragma omp parallel for
        for (int user = 0; user < n; ++user){
            auto encUserAvailablePref = evalDiagMatrixVecMult(encUsersPrefMatrixDiagonals[user],
                                                              encUserAvailability, cc);
            auto encUserFirstAvailablePref = evalPreserveLeadOne(encUserAvailablePref,
                                                                 cc, initPreserveLeadOne);
            // Mask and replicate availability row left and right.
            encUserFirstAvailablePref = cc->EvalMult(encUserFirstAvailablePref,encOnesRow);
            std::vector<Ciphertext<DCRTPoly>> addContainer;
            addContainer.push_back(encUserFirstAvailablePref);
            addContainer.push_back(cc->EvalRotate(encUserFirstAvailablePref,-n));
            addContainer.push_back(cc->EvalRotate(encUserFirstAvailablePref,n));
            encUserFirstAvailablePref = cc->EvalAddMany(addContainer);
            encRowsAdjMatrix[user] = evalDiagMatrixVecMult(encUsersPrefMatrixTransposedDiagonals[user],
                                                             encUserFirstAvailablePref,cc);
        }
        runtimePhase1 = TOC(t);
        runtimePhase1Total += runtimePhase1;
        std::cout << "Online part 1 - Adjacency matrix update time: " << runtimePhase1 << "ms" << std::endl;

        // Refresh after (1) update adjacency matrix.
        //----------------------------------------------------------
        std::cout << "Adjacency Matrix: " << std::endl;

        // Flat encoded adjacency matrix for matrix exponentiation.
        Ciphertext<DCRTPoly> encAdjMatrixFlat;

        // Refresh "encRowsAdjMatrix" as encrypted flat packed matrix.
        std::vector<std::vector<int64_t>> rowsAdjMatrix;
        for (int row=0; row < n; ++row){
            Plaintext plaintext;
            cc->Decrypt(keyPair.secretKey, encRowsAdjMatrix[row], &plaintext);
            plaintext->SetLength(n); auto payload = plaintext->GetPackedValue();
            // Print adjacence matrix.
            std::cout << payload << std::endl;
            rowsAdjMatrix.push_back(payload);
            // Also refresh "encRowsAdjMatrix" in row form for phase (3).
            refreshInPlace(encRowsAdjMatrix[row],n,keyPair,cc);
        }
        std::vector<int64_t> flatMatrix(n*n,0);
        for (int row=0; row < n; ++row){
            for (int col=0; col < n; ++col){
                flatMatrix[row*n+col] = rowsAdjMatrix[row][col];
            }
        }
        encAdjMatrixFlat = cc->Encrypt(keyPair.publicKey,
                                       cc->MakePackedPlaintext(repFillSlots(flatMatrix,slotTotal)));

        //----------------------------------------------------------
        // (2) Cycle finding.
        //----------------------------------------------------------

        // 2a) Matrix exponentiation.
        //----------------------------------------------------------
        // Cycle finding result [r_1, ..., r_n]. On cycle, r_i = 1. Not on cycle: r_i = 0.
        Ciphertext<DCRTPoly> encMatrixExpFlat;
        double runtimePhase2a(0.0);

        bool contFlag = true; int sqs = 1;
        while (contFlag) {
            int exp = 2 * sqs;
            if (exp >= n) { contFlag = false; }
            else { sqs = sqs + 1; }
        }

        TIC(t);
        encMatrixExpFlat = encAdjMatrixFlat;
        int refreshInterval = std::floor(chosen_depth/3);
        for (int i=1; i <= sqs; i++){
            encMatrixExpFlat = evalMatrixMult(cc,encMatrixExpFlat,encMatrixExpFlat,initMatrixMult);
            if (i % refreshInterval == 0) {
                runtimePhase2a += TOC(t);
                refreshInPlace(encMatrixExpFlat,cc->GetRingDimension(),keyPair,cc);
                TIC(t);
            }
        }
        runtimePhase2a += TOC(t);
        runtimePhase2aTotal += runtimePhase2a;
        std::cout << "Online part 2a - Matrix exponentiation: " << runtimePhase2a << " ms" << std::endl;

        // Refresh after (2a) matrix squaring.
        //----------------------------------------------------------
        Ciphertext<DCRTPoly> encMatrixExpPacked;

        std::vector<int64_t> packedMatrix(slotsPadded*n,0);
        Plaintext plaintext;
        cc->Decrypt(keyPair.secretKey,encMatrixExpFlat,&plaintext);
        plaintext->SetLength(n*n); auto payload = plaintext->GetPackedValue();
        for (int row = 0; row < n; row++){
            // std::vector<int64_t> matrixElemsRow;
            for (int col = 0; col < n; col++){
                int pos = col*slotsPadded + row;
                packedMatrix[pos] = payload[row*n+col];
            }
        }
        encMatrixExpPacked = cc->Encrypt(keyPair.publicKey,
                                         cc->MakePackedPlaintext(packedMatrix));

        // 2b) Cycle computation.
        //----------------------------------------------------------
        Ciphertext<DCRTPoly> enc_u_unmasked;
        double runtimePhase2b(0.0);

        TIC(t);

        auto encResMult = cc->EvalMult(encMatrixExpPacked,encOnes);
        auto encResInnerProd = evalPrefixAdd(encResMult,slotsPadded,cc); // TODO: cryptoOpsLogger
        enc_u_unmasked = evalNotEqualZero(encResInnerProd,cc,initNotEqualZero); // TODO: cryptoOpsLogger

        runtimePhase2b = TOC(t);
        runtimePhase2bTotal += runtimePhase2b;
        std::cout << "Online part 2b - Cycle computation: " << runtimePhase2b << "ms" << std::endl;

        // Refresh after (2b) cycle computation.
        //----------------------------------------------------------
        Ciphertext<DCRTPoly> enc_u;

        Plaintext plaintext2b;
        cc->Decrypt(keyPair.secretKey,enc_u_unmasked,&plaintext2b);
        plaintext2b->SetLength(n*slotsPadded); auto payload2b = plaintext2b->GetPackedValue();
        std::vector<int64_t> uElems;
        for (int user = 0; user < n; user++){
            uElems.push_back(payload2b[user*slotsPadded]);
        }
        enc_u = cc->Encrypt(keyPair.publicKey,
                cc->MakePackedPlaintext(uElems));

        //----------------------------------------------------------
        // (3) Update user availability and outputs.
        //----------------------------------------------------------
        double runtimePhase3(0.0);

        TIC(t);

        // Compute current preference index (t) for all users in packed ciphertext.
        std::vector<Ciphertext<DCRTPoly>> enc_elements;
        enc_elements.resize(n);
        #pragma omp parallel for
        for (int user=0; user < n; ++user){
            // Note: encRowsAdjMatrix must be refreshed after (1)
            auto enc_t_user = cc->EvalInnerProduct(encRowsAdjMatrix[user], encRange,
                                                   encRowsAdjMatrix.size());
            // enc_t_user = cc->EvalMult(enc_t_user, initRotsMasks.encMasks()[0]);
            enc_t_user = cc->EvalMult(enc_t_user, encLeadingOne);
            cc->ModReduceInPlace(enc_t_user);
            enc_elements[user] = cc->EvalRotate(enc_t_user,-user);
        }
        auto enc_t = cc->EvalAddMany(enc_elements);
        // o: Update output for all users in packed ciphertext: o <- t x u + o x (1-u)
        auto enc_t_mult_u = cc->EvalMult(enc_t, enc_u); cc->ModReduceInPlace(enc_t);
        auto enc_one_min_u = cc->EvalAdd(encOnes, cc->EvalMult(enc_u, encNegOnes));
        // output <- t x u + o x (1-u)
        enc_output = cc->EvalAdd(enc_t_mult_u, cc->EvalMult(enc_output,enc_one_min_u));
        // Update availability: 1-NotEqualZero(output)
        auto enc_output_reduced = evalNotEqualZero(enc_output,cc,initNotEqualZero);
        encUserAvailability = cc->EvalAdd(encOnes, cc->EvalMult(enc_output_reduced, encNegOnes));

        runtimePhase3 = TOC(t);
        runtimePhase3Total += runtimePhase3;
        std::cout << "Online part 3 - User availability & output update: " << runtimePhase3 << "ms" << std::endl;

        // Refresh after (3) update availability.
        //----------------------------------------------------------

        // Refresh encrypted output vector.
        Plaintext plaintext3; cc->Decrypt(keyPair.secretKey, enc_output, &plaintext3);
        plaintext3->SetLength(n); auto output = plaintext3->GetPackedValue();
        std::cout << "Output vector: " << output << std::endl;
        enc_output = cc->Encrypt(keyPair.publicKey,cc->MakePackedPlaintext(output));
        // Refresh & pack copies of user availability vector into single ciphertext.
        cc->Decrypt(keyPair.secretKey, encUserAvailability, &plaintext3);
        plaintext3->SetLength(n); auto userAvailability = plaintext3->GetPackedValue();
        std::cout << "Availability vector: " << userAvailability << std::endl;
        encUserAvailability = cc->Encrypt(keyPair.publicKey,
                                          cc->MakePackedPlaintext(repFillSlots(userAvailability,slotTotal)));

    // End loop.
    }
    std::cout << "-----------------------------------------" << std::endl;
    std::cout << "Online part 1 - Total runtime: " << runtimePhase1Total << "ms" << std::endl;
    std::cout << "Online part 2a - Total runtime: " << runtimePhase2aTotal << "ms" << std::endl;
    std::cout << "Online part 2b - Total runtime: " << runtimePhase2bTotal << "ms" << std::endl;
    std::cout << "Online part 3 - Total runtime: " << runtimePhase3Total << "ms" << std::endl;
    std::cout << "Online all - Total runtime: " << runtimePhase1Total+runtimePhase2aTotal+runtimePhase2bTotal+runtimePhase3Total << "ms" << std::endl;

    return 0;
}
