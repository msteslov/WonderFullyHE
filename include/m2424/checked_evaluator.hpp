#pragma once

#include "m2424/accuracy.hpp"
#include "m2424/linear_transform.hpp"
#include "m2424/seal_adapter.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace m2424 {

struct CheckedResult {
    std::string operation;
    Cipher cipher;
    CipherInfo info;
    AccuracyReport accuracy;
    bool ok{};
};

class CheckedEvaluator {
public:
    CheckedEvaluator(SealAdapter& adapter, std::size_t payload_size, double tolerance);

    CheckedResult add(const Cipher& lhs, const Cipher& rhs, const std::vector<double>& expected);
    CheckedResult sub(const Cipher& lhs, const Cipher& rhs, const std::vector<double>& expected);
    CheckedResult mul(const Cipher& lhs, const Cipher& rhs, const std::vector<double>& expected);
    CheckedResult rotate(const Cipher& input, int steps, const std::vector<double>& expected);
    CheckedResult sum_slots(const Cipher& input, std::size_t slot_count, const std::vector<double>& expected);
    CheckedResult linear_transform(const Cipher& input, const LinearTransform& transform,
                                   const std::vector<double>& expected);

private:
    CheckedResult finalize(std::string operation, Cipher cipher, const std::vector<double>& expected);

    SealAdapter& adapter_;
    std::size_t payload_size_{};
    double tolerance_{};
};

} // namespace m2424
