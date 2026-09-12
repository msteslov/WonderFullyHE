#include "m2424/sparse_bootstrap.hpp"
#include <stdexcept>
#include <iostream>
using namespace m2424;
using S=BootstrapCertificationStatus;
void check(bool b,const char* why) { if(!b) throw std::runtime_error(why); }
int main() {
 try {
    std::vector<PublicRlweFamily> families(2);
    families[0].id="ordinary"; families[0].statement="exact ordinary parameters and relations";
    families[1].id="sparse"; families[1].statement="exact sparse parameters and reciprocal key relation";
    families[1].searchSpaceCeilingBits=180;
    auto make=[](const PublicRlweFamily& f,double bits) { return RlweSecurityEvidence{f.id,f.statement,"test-only mock estimator","fixture-v1","test fixture artifact (not a production estimate)","mock assumptions, test only",bits}; };
    std::vector<RlweSecurityEvidence> evidence{make(families[0],160),make(families[1],140)};
    auto good=certifyPublicRlweFamilies(families,evidence,128);
    check(good.result.status==S::Certified&&good.minimumSecurityBits==140,"min across all families");
    auto missing=evidence; missing.pop_back();
    auto unknown=certifyPublicRlweFamilies(families,missing,128);
    check(!unknown.minimumSecurityBits&&unknown.result.status==S::SecurityBudgetExceeded,"one unknown fails global security");
    auto stale=evidence; stale[1].statement+="modified weight";
    check(certifyPublicRlweFamilies(families,stale,128).result.status==S::SecurityBudgetExceeded,"binding to exact statement");
    auto duplicate=evidence; duplicate.push_back(evidence[0]);
    check(certifyPublicRlweFamilies(families,duplicate,128).result.status==S::SecurityBudgetExceeded,"ambiguous evidence rejected");
    auto lower=evidence; lower[1].lowerSecurityBits=100;
    check(certifyPublicRlweFamilies(families,lower,128).result.status==S::SecurityBudgetExceeded,"below target");
    families[1].searchSpaceCeilingBits=20;
    check(certifyPublicRlweFamilies(families,evidence,128).result.status==S::SecurityBudgetExceeded,"low weight search-space ceiling is not security lower bound");
    evidence[1].lowerSecurityBits.reset();
    check(!certifyPublicRlweFamilies(families,evidence,128).minimumSecurityBits,"Unknown != zero");
    check(certifyPublicRlweFamilies({}, {},128).result.status==S::SecurityBudgetExceeded,"empty inventory rejected");
    std::cout<<"PASS sparse security contract (mock evidence only)\n";
 } catch(const std::exception& e) {std::cerr<<e.what()<<"\n";return 1;}
}
