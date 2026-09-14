#pragma once
#include "certified_arithmetic_internal.hpp"
#include "m2424/experimental/evalmod_analysis/approximation_lab.hpp"
#include "m2424/experimental/evalmod_analysis/exact_decimal.hpp"
#include <map>
namespace m2424::experimental::arithmetic {
inline std::vector<mpq_class> exactCoefficients(const EvalModPolynomial& polynomial) {
    std::vector<mpq_class> result;
    result.reserve(polynomial.decimalCoefficients.size());
    for(const auto& coefficient:polynomial.decimalCoefficients)
        result.push_back(parseExactDecimal(coefficient));
    while(result.size()>1&&result.back()==0)result.pop_back();
    return result;
}
inline bool exactPolynomialEqual(const EvalModPolynomial& left,const EvalModPolynomial& right) {
    return left.basis==right.basis&&exactCoefficients(left)==exactCoefficients(right);
}
// Generic power-basis polynomial compiler. Prefer the exact common-denominator
// path; use directly encoded rational coefficients only when its runtime scale
// is not finite/valid.
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
    void accumulate(std::optional<std::size_t>& sum,std::size_t term,const std::string& stage) {
        if(!sum) {sum=term;return;}
        auto [a,alignedTerm]=b.alignForAdd(*sum,term,stage);
        term=alignedTerm;
        sum=b.add(Op::Add,{a,term},stage);
    }
    std::size_t compileCommonDenominator(const std::vector<mpq_class>& c,
                                         const mpz_class& denominator,double denominatorScale,
                                         const std::string& stage) {
        std::optional<std::size_t> sum;
        for(std::size_t i=c.size();i-->1;) if(c[i]!=0) {
            auto term=power(i);const mpq_class numerator=c[i]*denominator;
            if(numerator!=1)term=b.scalar(term,numerator,1,stage);
            accumulate(sum,term,stage);
        }
        auto out=*sum;
        if(denominator!=1)out=b.scalar(out,mpq_class(1,denominator),denominatorScale,stage);
        if(c.size()>2)out=b.reduce(out,stage);
        if(c[0]!=0)out=b.plus(out,c[0],stage);
        return out;
    }
    std::size_t compileFiniteScales(const std::vector<mpq_class>& c,const std::string& stage) {
        std::optional<std::size_t> sum;
        for(std::size_t i=c.size();i-->1;) if(c[i]!=0) {
            auto term=power(i);
            const std::string coefficientStage=stage+" finite coefficient "+std::to_string(i);
            const double scale=b.scheduleScalarScale(term,c[i],coefficientStage);
            term=b.scalar(term,c[i],scale,coefficientStage);
            accumulate(sum,term,stage+" finite coefficient accumulation");
        }
        auto out=*sum;
        if(c.size()>2)out=b.reduce(out,stage);
        if(c[0]!=0)out=b.plus(out,c[0],stage+" finite coefficient 0");
        return out;
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
        const auto denominatorScale=projectUp(mpq_class(denominator));
        if(!denominatorScale) return compileFiniteScales(c,stage);
        const auto savedBuilder=b;
        const auto savedPowers=powers;
        try {
            return compileCommonDenominator(c,denominator,*denominatorScale,stage);
        } catch(const Failure& failure) {
            if(failure.status!=Status::ScaleScheduleInfeasible
               &&failure.status!=Status::HeadroomViolation) throw;
            b=savedBuilder;powers=savedPowers;
            return compileFiniteScales(c,stage);
        }
    }
};

// Exact scaled-Chebyshev execution for an already-certified canonical
// polynomial. T_k nodes use memoized fast doubling rather than Clenshaw.
class ScaledChebyshevCompiler {
    Builder& b;
    std::size_t t;
    mpq_class radius;
    std::map<std::size_t,std::size_t> values;
    std::vector<mpq_class> magnitudeBounds;
    double coefficientOutputScale;
    std::size_t multiply(std::size_t a,std::size_t c,const std::string& stage) {
        if(b.states[a].level<b.states[c].level)a=b.alignLevel(a,c,stage);
        if(b.states[c].level<b.states[a].level)c=b.alignLevel(c,a,stage);
        return b.stabilizeScale(b.reduce(b.mul(a,c,stage),stage),stage);
    }
    mpq_class magnitude(std::size_t k) {
        while(magnitudeBounds.size()<=k) {
            const auto n=magnitudeBounds.size();
            mpq_class next=2*radius*magnitudeBounds[n-1]-magnitudeBounds[n-2];
            if(next<0) throw Failure{Status::ExtractionNotCertified,
                "Negative exact Chebyshev magnitude recurrence outside R>=1 contract"};
            magnitudeBounds.push_back(std::move(next));
        }
        return magnitudeBounds[k];
    }
    std::size_t chebyshev(std::size_t k) {
        if(const auto found=values.find(k);found!=values.end())return found->second;
        const std::string stage="Chebyshev T_"+std::to_string(k);
        std::size_t result;
        if(k%2==0) {
            const auto half=chebyshev(k/2);
            result=multiply(half,half,stage);
            result=b.scalar(result,2,1,stage);
            result=b.plus(result,-1,stage);
        } else {
            auto left=chebyshev(k/2),right=chebyshev(k/2+1);
            result=multiply(left,right,stage);
            result=b.scalar(result,2,1,stage);
            auto [product,linear]=b.alignForAddByProduct(result,t,stage);
            result=b.add(Op::Subtract,{product,linear},stage);
            result=b.stabilizeScale(result,stage);
        }
        b.tightenMagnitude(result,magnitude(k),
            stage+": exact |T_k(t)|<=T_k(R), R=(K+rho)/64 recurrence");
        values[k]=result;
        return result;
    }
    void accumulate(std::optional<std::size_t>& sum,std::size_t term,const std::string& stage) {
        if(!sum){sum=term;return;}
        auto [a,c]=b.alignForAdd(*sum,term,stage);
        sum=b.add(Op::Add,{a,c},stage);
    }
public:
    ScaledChebyshevCompiler(Builder& builder,std::size_t normalized,const mpq_class& R):
        b(builder),t(normalized),radius(R),values{{1,normalized}},magnitudeBounds{1,R},
        coefficientOutputScale(builder.baselineDyadicScale(builder.states[normalized].level)
                               *builder.baselineDyadicScale(builder.states[normalized].level)) {
        if(radius<1)throw Failure{Status::ExtractionNotCertified,
            "Chebyshev magnitude recurrence requires exact R>=1"};
        b.tightenMagnitude(t,radius,
            "Normalized Chebyshev input: exact R=(K+rho)/64");
    }
    mpq_class exactMagnitudeBound(std::size_t k) { return magnitude(k); }
    // Exposed to arithmetic-layer tests and other certified polynomial
    // compilers; this is still a real memoized ciphertext node.
    std::size_t compileBasis(std::size_t k) {
        if(k==0)throw Failure{Status::ExtractionNotCertified,
            "T_0 is a plaintext constant, not a ciphertext basis node"};
        return chebyshev(k);
    }
    std::size_t compile(const EvalModPolynomial& polynomial,const std::string& stage) {
        if(polynomial.basis!=PolynomialBasis::Chebyshev
           ||polynomial.decimalCoefficients.empty()||polynomial.decimalCoefficients.size()>257)
            throw Failure{Status::ExtractionNotCertified,
                "Known scaled-Chebyshev coefficients (degree <=256) required"};
        auto coefficients=exactCoefficients(polynomial);
        if(coefficients.size()<2)throw Failure{Status::ExtractionNotCertified,
            "Constant-only Chebyshev extraction has no input-dependent ciphertext path"};
        std::optional<std::size_t> sum;
        for(std::size_t k=coefficients.size();k-->1;)if(coefficients[k]!=0) {
            auto term=chebyshev(k);
            const auto coefficientStage=stage+" Chebyshev coefficient "+std::to_string(k);
            const double scale=b.scheduleScalarScaleForOutput(
                term,coefficients[k],coefficientOutputScale,coefficientStage);
            term=b.scalar(term,coefficients[k],scale,coefficientStage);
            accumulate(sum,term,stage+" Chebyshev accumulation");
        }
        auto out=*sum;
        if(coefficients[0]!=0)out=b.plus(out,coefficients[0],stage+" Chebyshev coefficient 0");
        return out;
    }
};
}
