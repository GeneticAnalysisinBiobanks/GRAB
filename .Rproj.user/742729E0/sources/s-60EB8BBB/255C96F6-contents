// // [[Rcpp::depends(RcppArmadillo)]]
// #include <RcppArmadillo.h>
// 
// #include "../thirdParty/cget/include/savvy/reader.hpp"
// #include "../thirdParty/cget/include/savvy/varint.hpp"
// #include "../thirdParty/cget/include/savvy/sav_reader.hpp"
// #include "../thirdParty/cget/include/savvy/variant_group_iterator.hpp"
// 
// #include <iostream>
// #include <fstream>
// #include <vector>
// #include <string>
// #include <Rcpp.h>
// #include <stdlib.h>
// #include <cstring>
// 
// #include "VCF.hpp"
// 
// 
// namespace VCF {
// 
//   VcfClass::VcfClass(std::string t_vcfFileName,
//            std::string t_vcfFileIndex,
//            std::string t_vcfField,
//            std::string t_chr,
//            int32_t t_start,
//            int32_t t_end,
//            bool t_isDropMissingDosages,
//            double t_minMAF,
//            double t_minInfo,
//            double t_maxMAF,
//            std::vector<std::string> t_SampleInModel)
//   {
//     setVcfObj(t_vcfFileName,
//              t_vcfFileIndex,
//              t_vcfField,
//              t_chr,
//              t_start,
//              t_end,
//              t_isDropMissingDosages,
//              t_minMAF,
//              t_minInfo,
//              t_maxMAF);
// 
//     setPosSampleInVcf(t_SampleInModel);
// 
//   }
// 
// 
// 
// 
// 
//   void VcfClass::setVcfObj(const std::string t_vcfFileName,
//                  const std::string t_vcfFileIndex,
//                  const std::string t_vcfField,
//                  const std::string t_chr,
//                  int32_t t_start,
//                  int32_t t_end,
// 		 bool t_isDropMissingDosages,
// 		 double t_minMAF,
//                  double t_minInfo,
// 		 double t_maxMAF)
//   {
//     if(t_vcfField == "DS"){
//       m_reader = savvy::indexed_reader(t_vcfFileName, {t_chr, std::uint32_t(t_start), std::uint32_t(t_end)}, savvy::fmt::ds);
// 
//     }else if(t_vcfField == "GT"){
//       m_reader = savvy::indexed_reader(t_vcfFileName, {t_chr, std::uint32_t(t_start), std::uint32_t(t_end)}, savvy::fmt::ac);	    
//     }
// 
//     bool isVcfOpen = m_reader.good();
// 
//     if(isVcfOpen){
//       std::cout << "Open VCF done" << std::endl;
//       m_vcfField = t_vcfField;
//       std::cout << "To read the field " << m_vcfField << std::endl;
//       std::cout << "Number of meta lines in the vcf file (lines starting with ##): " << m_reader.headers().size() << std::endl;
//       std::cout << "Number of samples in in the vcf file: " << m_reader.samples().size() << std::endl;
//       m_SampleInVcf = m_reader.samples();
//       m_N0 = m_SampleInVcf.size();
//     }else{
//       std::cout << "WARNING: Open VCF failed" << std::endl;
//     }
// 
//     m_isDropMissingDosages = t_isDropMissingDosages;
//     m_minMAF = t_minMAF;
//     m_minInfo = t_minInfo;
//     m_maxMAF = t_maxMAF;
//   }
// 
//   void VcfClass::setPosSampleInVcf(std::vector<std::string> & t_SampleInModel)
//   {	  
//      m_N = t_SampleInModel.size();
//    
//      for(uint32_t i = 0; i < m_N; i++){
//        std::string sample = t_SampleInModel.at(i);
//        auto pos = std::find(m_SampleInVcf.begin(), m_SampleInVcf.end(), sample);
//        if(pos == m_SampleInVcf.end()){
//          Rcpp::stop("At least one subject requested is not in Vcf file.");
//        }
//      }
// 
//      m_posSampleInModel.clear();
//      for(uint32_t i = 0; i < m_N0; i++){
//        std::string sample = m_SampleInVcf.at(i);
//        auto pos = std::find(t_SampleInModel.begin(), t_SampleInModel.end(), sample);
//        if(pos != t_SampleInModel.end()){
//          m_posSampleInModel.push_back(pos - t_SampleInModel.begin());
//        }else{
// 	 m_posSampleInModel.push_back(-1);      
//        }
//      }
// 
//   }
// 
// 
//     // get dosages/genotypes of one marker
//   bool VcfClass::getOneMarker(uint32_t t_posMarker,
// 		         double& t_freq,
//                          double& t_missingRate,
//                          std::vector<uint32_t>& t_posMissingGeno,
//                          std::string& t_a1,
//                          std::string& t_a2,
//                          std::string& t_marker,
//                          uint32_t& t_pd,
// 			 uint8_t& t_chr,
//                          std::string& t_chrMarker,
// 			 arma::vec & t_dosage)
//   {	  
//     m_isReadVariant = m_reader >> m_record; 
//     if(m_isReadVariant){
//       t_marker = m_record.prop("ID");
//       t_chrMarker = m_record.chromosome();
//       t_pd = m_record.position();
//       int numAlt = 1;
//       t_a1 = m_record.ref();
//       t_a2 = m_record.alt();
//       double dosage;
//       double AC = 0;
//       t_dosage.set_size(m_N);
//       std::size_t missing_cnt = 0;
//       std::vector< size_t > indexforMissing;
//       int i;
//       for (std::vector<float>::iterator it = m_record.data().begin(); it != m_record.data().end(); ++it) {
//         i = std::distance(m_record.data().begin(), it);
//         if(m_posSampleInModel.at(i) != -1) {
//           if (std::isnan(*it)) {
//             t_dosage(m_posSampleInModel.at(i)) = -1;
//             ++missing_cnt;
//             indexforMissing.push_back(m_posSampleInModel.at(i));
//           }else {
//             t_dosage(m_posSampleInModel.at(i)) = *it;
//             AC = AC + *it;
//           }
//         }
//       }
// 
//       if(missing_cnt > 0){
//         if(missing_cnt == m_N){
//           t_freq = 0;
//         }else{
//           t_freq = AC / 2 / (double)(m_N - missing_cnt);
//         }
//         if(!m_isDropMissingDosages){	
// 	  double imputeDosage = 2*t_freq;
//           for (unsigned int i = 0; i < indexforMissing.size(); i++){
//             t_dosage(indexforMissing.at(i)) = imputeDosage;
// 	    AC = AC + imputeDosage;
//           }         	
// 	}  
//       }
// 
//     }else{
//       std::cout << "Reach the end of the vcf file" << std::endl;
//     }
//     return(m_isReadVariant);    	  
//   }
// 
//   bool VcfClass::getMultiMarker(std::string& t_markerline,
//                   std::vector<uint32_t>& t_posMissingGeno,
//                   std::vector<double> & t_freqMulti,
//                   std::vector<double> & t_acMulti,
//                   std::vector<std::string> & t_a1Multi,
//                   std::vector<std::string> & t_a2Multi,
//                   std::vector<std::string> & t_markerMulti,
//                   std::vector<uint32_t> & t_pdMulti,
//                   std::vector<uint8_t> & t_chrMarkerMulti,
//                   std::vector<double> & t_dosageMulti,
//                   std::vector<uint32_t> & t_iIndexMulti,
//                   std::vector<uint32_t> & t_jIndexMulti
//                   ){
//     savvy::variant_group_iterator<savvy::compressed_vector<double>> it(m_reader, t_markerline);
//     savvy::variant_group_iterator<savvy::compressed_vector<double>> end{};
//     t_dosageMulti.clear();
//     std::vector<double> MACs;
//     std::vector<double> dosagesforOneMarker;
//     std::vector<uint32_t> iIndexforOneMarker;
//     std::vector<uint32_t> jIndexforOneMarker;
//     std::vector<uint32_t> indexforMissing;
//     double freq, maf;
// 
//     if (it != end){
//       m_isReadVariant = true;	    
//       int missing_cnt = 0;
//       std::vector< int > indexforMissing;
//       int cnt = 0;
//       double AC = 0;
// 
//       MACs.clear();
//       for ( ; it != end; ++it)
//       {
//         AC = 0;
//         missing_cnt = 0;
//         indexforMissing.clear();
//         dosagesforOneMarker.clear();
//         iIndexforOneMarker.clear();
//         jIndexforOneMarker.clear();
//         it.group_id();
//         it.sites();
//         std::string marker_id = it->chromosome() + ":" + std::to_string(it->position()) + "_" + it->ref() + "/" + it->alt();
//         std::string markerInfo_str = it->prop("R2");
//         double markerInfo = strtod((markerInfo_str).c_str(),0);
//         for (auto dose_it = it->data().begin(); dose_it != it->data().end(); ++dose_it){
//           int lengthi = std::distance(it->data().begin(), it->data().end());
//           int i = dose_it.offset();
//  
//           if(m_posSampleInModel.at(i) != -1) {
//             if (std::isnan(*dose_it)) {
//               //dosagesforOneMarker(m_posSampleInModel.at(i)) = -1;
//               ++missing_cnt;
// 	      indexforMissing.push_back(m_posSampleInModel.at(i));
//               t_posMissingGeno.push_back(m_posSampleInModel.at(i));
//             }else {
//               if(*dose_it > 0){   		    
//                 dosagesforOneMarker.push_back(*dose_it);
// 		jIndexforOneMarker.push_back(cnt+1);
// 		iIndexforOneMarker.push_back(m_posSampleInModel.at(i)+1); //1-based for R
//                 AC = AC + *dose_it;
//               }
// 	    }
// 	  }    
//         }
// 
//         if(missing_cnt > 0){
//           if(missing_cnt == m_N){
//             freq = 0;
//           }else{
//             freq = AC / 2 / (double)(m_N - missing_cnt);
//           }
//         }
//         maf = freq;
//         if(freq > 0.5){maf = 1-maf;}	
// 
// 	if(maf >= m_minMAF && maf <= m_maxMAF && markerInfo >= m_minInfo){
//           if(missing_cnt > 0 && !m_isDropMissingDosages){
//             //std::cout << "missing_cnt > 0!" << std::endl;
//             double imputeDosage = 2*freq;
//             for (unsigned int i = 0; i < missing_cnt; i++){
//               dosagesforOneMarker.push_back(imputeDosage);
//               jIndexforOneMarker.push_back(cnt+1);
//               iIndexforOneMarker.push_back(indexforMissing.at(i)+1); //1-based for R
// 	      AC = AC + imputeDosage;
//             }
// 	  }
// 
//           t_dosageMulti.insert(std::end(t_dosageMulti), std::begin(dosagesforOneMarker), std::end(dosagesforOneMarker));
//           t_iIndexMulti.insert(std::end(t_iIndexMulti), std::begin(iIndexforOneMarker), std::end(iIndexforOneMarker));
// 	  t_jIndexMulti.insert(std::end(t_jIndexMulti), std::begin(jIndexforOneMarker), std::end(jIndexforOneMarker));
//           cnt = cnt + 1;
//           t_markerMulti.push_back(marker_id);
//           t_freqMulti.push_back(freq);
// 	  t_acMulti.push_back(AC);
//           t_pdMulti.push_back(it->position());
//        }      
//      }
// 
//     }else{
//       m_isReadVariant = false;
//     }	  
//     return(m_isReadVariant);
//   }
// 
// }
