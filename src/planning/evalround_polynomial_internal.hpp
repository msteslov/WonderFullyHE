#pragma once
#include "certified_arithmetic_internal.hpp"
#include "m2424/experimental/evalmod_analysis/exact_decimal.hpp"
#include <map>
namespace m2424::experimental::arithmetic {
// Generic power-basis polynomial compiler. Integer numerators and one exact
// common denominator avoid encoding every rational coefficient at scale one.
class PolynomialCompiler {
    Builder& b; std::map<std::size_t,std::size_t> powers;
    std::size_t power(std::size_t n) {
        if(powers.count(n)) return powers[n];
        auto a=power(n/2),c=power(n-n/2);
        if(b.states[a].components>2) a=b.reduce(a,"polynomial power");
        if(b.states[c].components>2) c=b.reduce(c,"polynomial power");
        if(b.states[a].level<b.states[c].level)a=b.alignLevel(a,c,"polynomial power");
        if(b.states[c].level<b.states[a].level)c=b.alignLevel(c,a,"polynomial power");
        return powers[n]=b.mul(a,c,"polynomial power");
    }
public:
    PolynomialCompiler(Builder& builder,std::size_t input):b(builder),powers{{1,input}} {}
    std::size_t compile(const EvalModPolynomial& polynomial,const std::string& stage) {
        if(polynomial.basis!=PolynomialBasis::Monomial||polynomial.decimalCoefficients.empty()||polynomial.decimalCoefficients.size()>257)
            throw Failure{Status::ExtractionNotCertified,"Known monomial coefficients (degree <=256) required"};
        std::vector<mpq_class> c;for(const auto& v:polynomial.decimalCoefficients)c.push_back(parseExactDecimal(v));
        while(c.size()>1&&c.back()==0)c.pop_back();
        if(c.size()<2) throw Failure{Status::ExtractionNotCertified,"Constant-only extraction has no nontransparent input-dependent ciphertext path"};
        mpz_class denominator=1;for(std::size_t i=1;i<c.size();++i)mpz_lcm(denominator.get_mpz_t(),denominator.get_mpz_t(),c[i].get_den_mpz_t());
        std::optional<std::size_t> sum;
        for(std::size_t i=c.size();i-->1;) if(c[i]!=0) {
            auto term=power(i);const mpq_class numerator=c[i]*denominator;
            if(numerator!=1)term=b.scalar(term,numerator,1,stage);
            if(!sum) {sum=term;continue;}
            auto a=*sum;
            if(b.states[a].level<b.states[term].level)a=b.alignLevel(a,term,stage);
            if(b.states[term].level<b.states[a].level)term=b.alignLevel(term,a,stage);
            const double sa=b.states[a].scale,st=b.states[term].scale;
            if(bits(sa)!=bits(st)) {
                const double multiplier=st*st==sa?st:sa/st;
                // This is an actual multiplyPlain(1), with coefficient encoding
                // and exact-to-runtime scale-ratio errors charged by Builder.
                if(std::isfinite(multiplier)&&multiplier>0&&st*multiplier==sa) term=b.scalar(term,1,multiplier,stage);
                else {a=b.scalar(a,1,st,stage);term=b.scalar(term,1,sa,stage);}
            }
            sum=b.add(Op::Add,{a,term},stage);
        }
        auto out=*sum;
        if(denominator!=1)out=b.scalar(out,mpq_class(1,denominator),up(mpq_class(denominator)),stage);
        if(c.size()>2)out=b.reduce(out,stage);
        if(c[0]!=0)out=b.plus(out,c[0],stage);
        return out;
    }
};
}
