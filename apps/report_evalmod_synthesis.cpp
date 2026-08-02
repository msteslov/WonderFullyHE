#include "m2424/experimental/evalmod_analysis/synthesis.hpp"

#include <cstdio>
#include <string>

int main(int argc, char** argv) {
    namespace em = m2424::experimental;
    if (argc > 2 || (argc == 2 && std::string(argv[1]) != "--csv")) {
        std::fprintf(stderr, "usage: report_evalmod_synthesis [--csv]\n");
        return 2;
    }
    const em::EvalModProblem problem{
        65537, em::ExactScale::rational(1 << 20, 1), em::ExactScale::rational(1 << 20, 1),
        {em::TailModel::Subgaussian, 4096, 0.0, 0.01, 0.08, 1e-5, 1e-5, 1e-7, 1e-6,
         -128.0, 384, {"component-wise subgaussian composition",
                       "encryption,key-switching,CtS", "centered independent noise sources"}},
        12, 40, 1.0, {1.0, 0.2, 0.05, 0.2, 0.1, 4096}, 1.0, 0.0, 0.0
    };
    const auto result = em::synthesizeEvalMod(problem);
    const std::string report = argc == 2 ? em::evalModSynthesisCsv(result)
                                         : em::evalModSynthesisJson(result) + "\n";
    std::fwrite(report.data(), 1, report.size(), stdout);
    return result.selectedCandidate ? 0 : 1;
}
