#include "../core/sparse_bootstrap_internal.hpp"
#include "certified_arithmetic_internal.hpp"
#include <mpfr.h>
#include <sstream>
namespace m2424 {
namespace {
using namespace experimental::arithmetic;
using S=BootstrapCertificationStatus;
std::string statement(const PublicRlweFamily& f,std::uint64_t generation) {
    std::ostringstream o; o<<"PR6/family-v1|"<<generation<<'|'<<f.id<<'|'<<f.degree<<'|';
    for(auto p:f.modulus) o<<p<<',';
    o<<'|'<<f.secretDistribution<<'|'<<f.errorDistribution<<'|'<<f.samples<<'|'<<f.components<<'|'
     <<f.relations<<'|'<<f.purpose<<'|'<<f.provenance; return o.str();
}
BootstrapSecurityReport inventory(const SealAdapter& a,const SparseBootstrapInput& input,const BootstrapInputContext& source) {
    auto key=a.sparseKeyMetadata(); auto QP=a.dataModulusValues(); QP.push_back(a.specialKeyModulusValue());
    mpz_class space; mpz_bin_uiui(space.get_mpz_t(),key.degree,key.weight); space<<=key.weight;
    mpfr_t log; mpfr_init2(log,256); mpfr_set_z(log,space.get_mpz_t(),MPFR_RNDU); mpfr_log2(log,log,MPFR_RNDU);
    const double ceiling=mpfr_get_d(log,MPFR_RNDU); mpfr_clear(log);
    const std::string dense=key.ordinaryDistributionKnown?"iid uniform ternary {-1,0,1}":"Unknown imported secret/key distribution";
    const std::string sparse=key.distribution+"; h="+std::to_string(key.weight);
    const std::string error="CBD difference HW21-HW21; coefficient support [-21,21]";
    std::vector<PublicRlweFamily> families;
    auto add=[&](std::string id,std::string secret,std::size_t count,std::string relation,std::string purpose) {
        PublicRlweFamily f; f.id=id; f.degree=key.degree; f.modulus=QP; f.secretDistribution=secret;
        f.errorDistribution=error; f.samples=count; f.components=2*count; f.relations=relation;
        f.purpose=purpose; f.provenance="Adapter-generated material; SEAL util/rlwe.cpp and keygenerator.cpp; all families share s and s_b as stated; no independent-family security assumption";
        // Every family belongs to the coupled key graph. Recovering s_b allows
        // attacking its encapsulation messages (s); the ceiling is not a lower bound.
        f.searchSpaceCeilingBits=ceiling; f.statement=statement(f,key.generation); families.push_back(f);
    };
    add("ordinary.public",dense,key.publicKeySamples,"(-a*s-e,a); message zero; same s as evaluation and restoration", "Ordinary public key");
    for(std::size_t i=0;i<key.encryptionModuli.size();++i) {
        const auto& entry=key.encryptionModuli[i];
        add("ordinary.encryption/"+std::to_string(i),"fresh iid ternary u per encryption",entry.second*2,
            "(pk0*u+e0,pk1*u+e1) at the exact predecessor modulus; shared public key; one divide-round then message addition at requested level; ciphertext components share u", "Public encryption RLWE-like systems before deterministic ModDown");
        families.back().modulus=entry.first; families.back().components=entry.second*2;
        families.back().statement=statement(families.back(),key.generation);
    }
    add("ordinary.relinearization",dense,key.relinSamples,"Gadget messages P*s^2 under s; decomposition limbs share source secret; circular/KDM relation", "Relinearization keys");
    std::string automorphisms; for(auto g:key.galoisElements) automorphisms+=std::to_string(g)+",";
    add("ordinary.galois",dense,key.galoisSamples,"Gadget messages P*tau_g(s) under s; g="+automorphisms+"; shared/circular relation", "Rotation and conjugation keys");
    add("sparse.encapsulation",sparse,a.dataModulusValues().size(),"Gadget messages P*s under s_b; s_b independent of s at generation", "s -> s_b before ModRaise");
    add("sparse.restoration",dense,a.dataModulusValues().size(),"Gadget messages P*s_b under s; reciprocal relation to encapsulation (two-key cycle)", "s_b -> s after ModRaise");
    for(const auto& stage:std::vector<std::string>{"encapsulated","raised","restored"}) {
        add("derived."+stage,stage=="restored"?dense:sparse,1,
            "Deterministic public transform of the input and the above public switch keys; no independent RLWE sample or error draw; all known-message and reciprocal-key relations retained", "One bootstrap invocation: "+stage);
        auto& f=families.back(); f.modulus=stage=="encapsulated"?source.sourcePrimes:source.raisedPrimes;
        f.errorDistribution="Derived convolution of upstream noise, key CBD errors and componentwise divide-round; not an independent CBD sample";
        f.statement=statement(f,key.generation);
    }
    // Every estimator result is bound to the entire published key graph, not
    // just its own marginal family. Another publication invalidates old evidence.
    std::string graph="|joint-inventory:";
    for(const auto& f:families) graph+=std::to_string(f.statement.size())+":"+f.statement;
    for(auto& f:families) f.statement+=graph;
    return certifyPublicRlweFamilies(std::move(families),input.securityEvidence,input.targetSecurityBits);
}
}
SparseBootstrapPlan prepareSparseBootstrap(const SealAdapter& a,const Cipher& input,const SparseBootstrapInput& request) {
    using namespace experimental::arithmetic;
    SparseBootstrapPlan plan; auto p=std::make_shared<SparseBootstrapPlan::Data>(); plan.data_=p;
    auto& c=p->certificate; c.keys=a.sparseKeyMetadata();
    auto reject=[&](S s,const std::string& gate,const std::string& why) { c.result={s,gate,why}; return plan; };
    auto resolved=resolveBootstrapInput(a,input); if(!resolved.context) { c.result=resolved.result; return plan; }
    c.input=*resolved.context;
    if(!a.hasSparseEncapsulationKey()||!a.hasSparseRestorationKey()) return reject(S::MissingEvaluationKeys,"sparse.keys","Both actual directional switching keys required");
    if(c.keys.weight<2||c.keys.weight>c.keys.degree||!c.keys.ordinaryDistributionKnown||a.info(input).ciphertextSize!=2)
        return reject(S::SparseSecretCertificateUnavailable,"sparse.secret","Generated fixed-weight h>=2 and ordinary ternary keys required");
    c.security=inventory(a,request,c.input);
    const auto noise=request.evaluationKeyNoiseSupport.value_or(experimental::finiteSupportBackendKeyNoise());
    const auto backend=experimental::finiteSupportBackendKeyNoise();
    if(!known(noise)||!known(backend)||noise.upperBound<backend.upperBound)
        return reject(S::KeySwitchBoundUnavailable,"sparse.arithmetic","Unknown or insufficient actual backend evaluation-key error support");
    if(!known(request.messageMagnitude)||!known(request.sourceNoiseMagnitude))
        return reject(S::RequiredBoundUnavailable,"sparse.upstream","Message and original noise coefficient bounds required; Unknown is not zero");
    const auto N=c.keys.degree; const mpq_class scale(a.scale(input));
    mpz_class qsrc=1,Q=1; for(auto prime:c.input.sourcePrimes) qsrc*=integer(prime);
    for(auto prime:c.input.raisedPrimes) Q*=integer(prime);
    auto ks=[&](const std::vector<std::uint64_t>& primes) {
        return mpq_class(experimental::finiteSupportKeyNoise(N,q(noise.upperBound),primes,c.input.specialPrime,scale)
             +experimental::finiteSupportDivideRound(N,1,2,scale));
    };
    const auto enc=ks(c.input.sourcePrimes),restore=ks(c.input.raisedPrimes);
    c.encapsulationError=bound(enc,"s -> s_b: shared K(active source primes)+R(N,H=1,2); canonical error / Delta0; H=1 conservative also for fixed weight");
    c.restorationError=bound(restore,"s_b -> s: shared K(all raised data primes)+R(N,H=1,2); canonical error / Delta0; distinct from nu_b");
    const mpq_class nu=q(request.sourceNoiseMagnitude.upperBound)+enc*scale;
    c.sourceNoiseMagnitude=bound(nu,"nu_b coefficient bound: upstream coefficient noise + Delta0 * canonical encapsulation error (inverse embedding norm <=1)");
    if(q(request.messageMagnitude.upperBound)+nu>=mpq_class(qsrc)/2)
        return reject(S::HeadroomViolation,"sparse.source","m+nu_b not proven strictly inside (-q_src/2,q_src/2)");
    // Each centered component coefficient <=(q-1)/2. Convolution with exactly
    // h signed units has l1=h. Together with |m+nu_b|<q/2 this implies
    // |I|<(h+2)/2<=h for h>=2. We retain the v9 baseline K=h.
    c.K=static_cast<std::uint32_t>(c.keys.weight);
    c.liftProvenance="v9 section 3.1; exact fixed weight h; centered |c_i|<q_src/2; |m+nu_b|<q_src/2; |q_src I|<(h+2)q_src/2; integer I => |I|<=h. No decryption or observed maximum.";
    const mpq_class raised=mpq_class(integer(N)*(integer(c.keys.weight)+1)*(qsrc-1))/2/scale;
    c.raisedMagnitude=bound(raised,"N*(h+1)*(q_src-1)/(2 Delta0): centered component convolution, canonical embedding <=N*coefficient norm");
    c.restoredMagnitude=bound(raised+restore,"Sparse ModRaise ideal magnitude plus separately certified restoration canonical error");
    c.rhoBeforeCoeffToSlot=bound((q(request.messageMagnitude.upperBound)+nu+restore*scale)/qsrc,
        "(M_m+M_nu_b+Delta0*E_restore)/q_src; excludes CtS operator/local error, so not final rho_cert");
    if(2*scale*(raised+restore)>=mpq_class(Q))
        return reject(S::HeadroomViolation,"sparse.raised","Raised/restoration semantic envelope lacks centered headroom in actual Q");
    c.result={S::Certified,"sparse.arithmetic","Encapsulation, centered lift and restoration finite-support certificates; security is a separate gate"};
    return plan;
}
}
