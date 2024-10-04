#include "crypto_matrix_operations.h"

#include <cassert>

Ciphertext<DCRTPoly> evalDiagMatrixVecMult(std::vector<Ciphertext<DCRTPoly>> &encMatDiagonals, // Output of repFillSlots()
                                           Ciphertext<DCRTPoly> encVec,                        // Output of repFillSlots()
                                           CryptoContext<DCRTPoly> &cryptoContext) {
    int d = encMatDiagonals.size();
    std::vector<Ciphertext<DCRTPoly>> addContainer;
    addContainer.resize(d);

    #pragma omp parallel for
    for (int l = 0; l < d; l++) {
        auto encVecRot = cryptoContext->EvalRotate(encVec,l);
        auto encVecRotMult = cryptoContext->EvalMult(encMatDiagonals[l],encVecRot);
        #pragma omp critical
        addContainer[l] = encVecRotMult;
    }

    auto res = cryptoContext->EvalAddMany(addContainer);
    return res;
}

InitMatrixMult::InitMatrixMult(CryptoContext<DCRTPoly> &cryptoContext, KeyPair<DCRTPoly> keyPair, int d) :
    d(d) {
        auto maxSlots = cryptoContext->GetRingDimension();
        auto n = d*d;
        // STEP 1-1
        std::vector<int> iterRange;
        for (int k = -d; k <= d; k++){ iterRange.push_back(k); }
         // Pre-process encryption of u_sigma.
        for (int k : iterRange) {
            std::vector<int64_t> u_sigma_k(n,0);
            if (k < 0) {
                for (int l = 0; l < n; l++){
                    if (-k <= (l-(d+k)*d) && (l-(d+k)*d) < d){ u_sigma_k[l] = 1; }
                }
            }
            if (k >= 0) {
                for (int l = 0; l < n; l++){
                    if (0<=(l-d*k) && (l-d*k) < (d-k)){ u_sigma_k[l] = 1; }
                }
            }
            _u_sigma[k] = cryptoContext->Encrypt(keyPair.publicKey,
                                                 cryptoContext->MakePackedPlaintext(repFillSlots(u_sigma_k,maxSlots)));
        }
        // STEP 1-2
         // Pre-process encryption of u_tau.
        for (int k = 0; k < d; k++) {
            std::vector<int64_t> u_tau_k(n,0);
            for (int i = 0; i < d; i++){
                u_tau_k[k+d*i]=1;
            }
            _u_tau[d*k] = cryptoContext->Encrypt(keyPair.publicKey,
                                                 cryptoContext->MakePackedPlaintext(repFillSlots(u_tau_k,maxSlots)));
        }
        // STEP 2
        for (int k = 1; k < d; k++) {
            // Pre-process encryption of v1, v2.
            std::vector<int64_t> v1_k(n,0);
            std::vector<int64_t> v2_k_d(n,0);
            for (int l = 0; l < n; l++){
                if (0 <= l % d && l % d < d-k) { v1_k[l] = 1; }
                if (d-k <= l % d && l % d < d) { v2_k_d[l] = 1; }
            }
            _v1[k] = cryptoContext->Encrypt(keyPair.publicKey,
                                           cryptoContext->MakePackedPlaintext(repFillSlots(v1_k,maxSlots)));
            _v2[k-d] = cryptoContext->Encrypt(keyPair.publicKey,
                                             cryptoContext->MakePackedPlaintext(repFillSlots(v2_k_d,maxSlots)));
        }
        std::vector<int64_t> matrixMask(n,1);
        _matrixMask = cryptoContext->Encrypt(keyPair.publicKey,
                                             cryptoContext->MakePackedPlaintext(matrixMask));
    }

    std::map<int, Ciphertext<DCRTPoly>> InitMatrixMult::u_sigma() { return _u_sigma; }
    std::map<int, Ciphertext<DCRTPoly>> InitMatrixMult::u_tau() { return _u_tau; }
    std::map<int, Ciphertext<DCRTPoly>> InitMatrixMult::v1() { return _v1; }
    std::map<int, Ciphertext<DCRTPoly>> InitMatrixMult::v2() { return _v2; }
    Ciphertext<DCRTPoly> InitMatrixMult::matrixMask() { return _matrixMask; }


Ciphertext<DCRTPoly> evalMatrixMult(CryptoContext<DCRTPoly> &cryptoContext,
                                    Ciphertext<DCRTPoly> encA,
                                    Ciphertext<DCRTPoly> encB,
                                    InitMatrixMult &initMatrixMult) {
        // Note: Encrypted matrix must be consistent with initMatrixMult dimension (d).
        auto d = initMatrixMult.d;
        // STEP 1-1
        std::vector<int> iterRange;
        for (int k = -d; k <= d; k++){ iterRange.push_back(k); }
        std::vector<Ciphertext<DCRTPoly>> A_0_container;
        A_0_container.resize(iterRange.size());
        int container_idx = 0;
        // #pragma omp parallel for
        for (int k : iterRange) {
            Ciphertext<DCRTPoly> A_rot;
            Ciphertext<DCRTPoly> A_rot_mult;
            if (k != 0) {
                A_rot = cryptoContext->EvalRotate(encA,k);
            }
            else {
                A_rot = encA;
            }
            A_rot_mult = cryptoContext->EvalMult(A_rot, initMatrixMult.u_sigma()[k]);
            A_0_container[container_idx] = A_rot_mult;
            container_idx++;
        }
        auto A_0 = cryptoContext->EvalAddMany(A_0_container);
        // STEP 1-2
        std::vector<Ciphertext<DCRTPoly>> B_0_container;
        B_0_container.resize(d);
        // #pragma omp parallel for
        for (int k = 0; k < d; k++) {
            auto B_rot = cryptoContext->EvalRotate(encB,d*k);
            auto B_rot_mult = cryptoContext->EvalMult(B_rot,initMatrixMult.u_tau()[d*k]);
            B_0_container[k] = B_rot_mult;
        }
        auto B_0 = cryptoContext->EvalAddMany(B_0_container);
        // STEP 2
        // std::map<int, Ciphertext<DCRTPoly>> A;
        // std::map<int, Ciphertext<DCRTPoly>> B;
        std::vector<Ciphertext<DCRTPoly>> A;
        std::vector<Ciphertext<DCRTPoly>> B;
        A.resize(d);
        B.resize(d);
        // #pragma omp parallel for
        for (int k = 1; k < d; k++) {
            auto A_k = cryptoContext->EvalMult(initMatrixMult.v1()[k],
                                               cryptoContext->EvalRotate(A_0,k));
            auto A_k_d = cryptoContext->EvalMult(initMatrixMult.v2()[k-d],
                                                 cryptoContext->EvalRotate(A_0,k-d));
            A[k] = cryptoContext->EvalAdd(A_k,A_k_d);
            B[k] = cryptoContext->EvalRotate(B_0,d*k);
        }
        // STEP 3
        std::vector<Ciphertext<DCRTPoly>> AB_container;
        AB_container.resize(d);
        AB_container[0] =cryptoContext->EvalMult(A_0,B_0);
        // #pragma omp parallel for
        for (int k = 1; k < d; k++) {
            AB_container [k] = cryptoContext->EvalMult(A[k],B[k]);
        }
        auto AB =  cryptoContext->EvalAddMany(AB_container);
        return AB;
    }


Ciphertext<DCRTPoly> evalMatrixMultParallel(CryptoContext<DCRTPoly> &cryptoContext,
                                            Ciphertext<DCRTPoly> encA,
                                            Ciphertext<DCRTPoly> encB,
                                            InitMatrixMult &initMatrixMult) {
        // Note: Encrypted matrix must be consistent with initMatrixMult dimension (d).
        auto d = initMatrixMult.d;
        // STEP 1-1
        std::vector<int> iterRange;
        for (int k = -d; k <= d; k++){ iterRange.push_back(k); }

        std::vector<Ciphertext<DCRTPoly>> A_0_container;
        A_0_container.resize(2*d+1);

        int container_idx = 0;
        // std::cout << "Before loop 1" << std::endl;
        #pragma omp parallel for
        for (int k : iterRange) {
            Ciphertext<DCRTPoly> A_rot;
            Ciphertext<DCRTPoly> A_rot_mult;
            if (k != 0) {
                A_rot = cryptoContext->EvalRotate(encA,k);
            }
            else {
                A_rot = encA;
            }
            A_rot_mult = cryptoContext->EvalMult(A_rot, initMatrixMult.u_sigma()[k]);
            #pragma omp critical
            A_0_container[container_idx] = A_rot_mult;
            container_idx++;
        }
        // std::cout << "After loop 1" << std::endl;
        auto A_0 = cryptoContext->EvalAddMany(A_0_container);

        // STEP 1-2
        std::vector<Ciphertext<DCRTPoly>> B_0_container;
        B_0_container.resize(d);
        // std::cout << "Before loop 2" << std::endl;
        #pragma omp parallel for
        for (int k = 0; k < d; k++) {
            auto B_rot = cryptoContext->EvalRotate(encB,d*k);
            auto B_rot_mult = cryptoContext->EvalMult(B_rot,initMatrixMult.u_tau()[d*k]);
            #pragma omp critical
            B_0_container[k] = B_rot_mult;
        }
        // std::cout << "After loop 2" << std::endl;
        auto B_0 = cryptoContext->EvalAddMany(B_0_container);

        // STEP 2
        // std::map<int, Ciphertext<DCRTPoly>> A;
        // std::map<int, Ciphertext<DCRTPoly>> B;
        std::vector<Ciphertext<DCRTPoly>> A;
        std::vector<Ciphertext<DCRTPoly>> B;
        A.resize(d-1);
        B.resize(d-1);
        // std::cout << "Before loop 3" << std::endl;
        #pragma omp parallel for
        for (int k = 1; k < d; k++) {
            auto A_k = cryptoContext->EvalMult(initMatrixMult.v1()[k],
                                               cryptoContext->EvalRotate(A_0,k));
            auto A_k_d = cryptoContext->EvalMult(initMatrixMult.v2()[k-d],
                                                 cryptoContext->EvalRotate(A_0,k-d));
            #pragma omp critical
            A[k-1] = cryptoContext->EvalAdd(A_k,A_k_d);
            #pragma omp critical
            B[k-1] = cryptoContext->EvalRotate(B_0,d*k);
        }
        // std::cout << "After loop 3" << std::endl;
        // STEP 3
        std::vector<Ciphertext<DCRTPoly>> AB_container;
        AB_container.resize(d);
        AB_container[0] = cryptoContext->EvalMult(A_0,B_0);
        // std::cout << "Before loop 4" << std::endl;
        #pragma omp parallel for
        for (int k = 1; k < d; k++) {
            auto res = cryptoContext->EvalMult(A[k-1],B[k-1]);
            #pragma omp critical
            AB_container[k] = res;
        }
        // std::cout << "After loop 4" << std::endl;
        return cryptoContext->EvalAddMany(AB_container);
    }
