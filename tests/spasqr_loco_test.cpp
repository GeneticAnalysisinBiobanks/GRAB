#include "spasqr/null_model.hpp"
#include "tinytest.hpp"

TEST(transform_selection_is_per_phenotype) {
    LocoData loco;
    loco.regeniePhenos.insert("regenie");
    using spasqr_null::resolvePhenoTransform;
    CHECK(resolvePhenoTransform("", nullptr, "plain") == "int");
    CHECK(resolvePhenoTransform("", &loco, "ldak") == "raw");
    CHECK(resolvePhenoTransform("", &loco, "regenie") == "standardize");
    for (const std::string mode : {"raw", "int", "standardize"}) {
        CHECK(resolvePhenoTransform(mode, nullptr, "plain") == mode);
        CHECK(resolvePhenoTransform(mode, &loco, "ldak") == mode);
        CHECK(resolvePhenoTransform(mode, &loco, "regenie") == "standardize");
    }
}

TEST(regenie_scale_uses_covariate_residual_degrees_of_freedom) {
    spasqr_null::PhenoWork pk;
    pk.X.resize(6, 1);
    pk.X << -2, -1, 0, 0, 1, 2;
    Eigen::VectorXd r(6);
    r << 1, -2, 1, 1, -2, 1; // Orthogonal to the intercept and X.
    pk.Y = 10.0 + 3.0 * pk.X.col(0).array() + 2.0 * r.array();
    const Eigen::VectorXd original = pk.Y;
    auto rescaled = pk;
    rescaled.Y *= 7.0;

    spasqr_null::scaleRegenieResponse(pk);
    spasqr_null::scaleRegenieResponse(rescaled);
    // Residual sum of squares = 48; residual degrees of freedom = 6 - 2.
    for (Eigen::Index i = 0; i < pk.Y.size(); ++i) {
        CHECK_CLOSE(pk.Y[i], original[i] / std::sqrt(12.0), 1e-13, 1e-13);
        CHECK_CLOSE(pk.Y[i], rescaled.Y[i], 1e-13, 1e-13);
    }
}

TEST(regenie_scale_without_covariates_uses_sample_sd) {
    spasqr_null::PhenoWork pk;
    pk.X.resize(4, 0);
    pk.Y.resize(4);
    pk.Y << 1, 3, 5, 7;
    const Eigen::VectorXd original = pk.Y;
    spasqr_null::scaleRegenieResponse(pk);
    for (Eigen::Index i = 0; i < pk.Y.size(); ++i)
        CHECK_CLOSE(pk.Y[i], original[i] / std::sqrt(20.0 / 3.0), 1e-13, 1e-13);
}

TEST(regenie_scale_rejects_zero_variance_and_no_residual_degrees_of_freedom) {
    for (Eigen::Index n : {1, 4}) {
        spasqr_null::PhenoWork pk;
        pk.X.resize(n, 0);
        pk.Y = Eigen::VectorXd::Zero(n);
        bool rejected = false;
        try {
            spasqr_null::scaleRegenieResponse(pk);
        } catch (const std::runtime_error &) {
            rejected = true;
        }
        CHECK(rejected);
    }
}

TINYTEST_MAIN
