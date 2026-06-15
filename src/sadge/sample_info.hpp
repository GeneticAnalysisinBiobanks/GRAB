// sample_info.hpp — SADGE pedigree parsing + family-structure model
//
// Port of generateSampleInfo() / extract.familyID() from the R reference
// (example_sadge_zouyu/impute_func.R).  A PLINK .fam supplies, per row, an
// index individual ("sibling") with FID, IID, FatherID, MotherID.  Families
// are classified by the number of GENOTYPED parents (n_par ∈ {0,1,2}) and the
// number of GENOTYPED siblings (n_sib ∈ {1,2}).  This "stage-1" (genotyped)
// structure drives the parental-genotype imputation kernels.  The "stage-2"
// (phenotyped) structure used by the per-phenotype variance/score is rebuilt
// later in precond.{hpp,cpp} from the analysis subjects of each phenotype.
#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace sadge {

// Family-structure index.  The integer order MUST match the R name order in
// var_geno_par_func (score_func_vec.R) so per-structure vectors line up by
// position: [0par1sib, 1par1sib, 2par1sib, 0par2sib, 1par2sib, 2par2sib].
enum class FamStruct : uint8_t {
    s0par1sib = 0,
    s1par1sib = 1,
    s2par1sib = 2,
    s0par2sib = 3,
    s1par2sib = 4,
    s2par2sib = 5
};
constexpr int kNStruct = 6;

// Compose a FamStruct from genotyped parent count (0/1/2) and sib count (1/2).
inline FamStruct famStruct(int nPar, int nSib) {
    if (nSib == 1) {
        return static_cast<FamStruct>(nPar); // 0/1/2 → s{0,1,2}par1sib
    }
    return static_cast<FamStruct>(3 + nPar); // 0/1/2 → s{0,1,2}par2sib
}

inline bool isSingleton(FamStruct fs) {
    return fs == FamStruct::s0par1sib;
}

// One genotyped sibling (an index .fam row in a valid n_sib∈{1,2} family).
struct SibInfo {
    std::string iid;
    long fid;          // numeric FID (grouping key)
    FamStruct fs1;     // genotyped (stage-1) family structure
    int genoNPar;      // 0/1/2 genotyped parents
    int genoNSib;      // 1/2 genotyped siblings
    int32_t pedSelf;   // pedigree-decode row of this sibling
    int32_t pedCoSib;  // pedigree-decode row of the genotyped co-sibling, or -1
    int32_t pedPar1;   // genotyped parent #1 ped row, or -1 (the lone parent for n_par==1)
    int32_t pedPar2;   // genotyped parent #2 ped row, or -1
};

struct SampleInfo {
    // Subject-ID list for SADGE's own genotype cursor: every genotyped
    // sibling plus every genotyped referenced parent, in genotype-file order.
    std::vector<std::string> pedigreeIIDs;
    // All genotyped siblings (valid families), in genotype-file order.
    std::vector<SibInfo> siblings;
    // IID → index into siblings[].
    std::unordered_map<std::string, int> sibIndexByIID;
    // IID → row in the pedigree decode (rank within pedigreeIIDs).
    std::unordered_map<std::string, int32_t> pedRowByIID;
};

// Parse a 4-column PLINK .fam (FID IID FatherID MotherID), retain only rows
// whose IID is in genoIIDs and whose non-"0" parents are in genoIIDs
// (matching the R filters in pipeline1.R), classify families with
// n_sib ∈ {1,2}, and build the pedigree subject set + decode-row maps.
SampleInfo buildSampleInfo(
    const std::string &famFile,
    const std::vector<std::string> &genoIIDs);

} // namespace sadge
