#include "m2424/accuracy.hpp"
#include "m2424/profiles.hpp"
#include "m2424/seal_adapter.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

template <class Fn>
double elapsed_ms(Fn&& fn) {
    const auto start = Clock::now();
    fn();
    const auto finish = Clock::now();
    return std::chrono::duration<double, std::milli>(finish - start).count();
}

std::vector<double> head(const std::vector<double>& values, std::size_t n) {
    return std::vector<double>(values.begin(), values.begin() + std::min(values.size(), n));
}

struct Task {
    std::vector<double> input;
    m2424::Cipher encrypted;
    m2424::Cipher result;
};

struct WorkerState {
    m2424::SealAdapter adapter;
    std::vector<Task> tasks;
};

std::vector<double> make_input(std::size_t payload_size, std::size_t task_index) {
    std::vector<double> input;
    input.reserve(payload_size);
    const double shift = static_cast<double>(task_index) / 17.0;
    for (std::size_t i = 0; i < payload_size; ++i) {
        const double x = 0.2 + 0.1 * std::sin(static_cast<double>(i) / 11.0 + shift);
        input.push_back(x);
    }
    return input;
}

std::vector<double> square_ref(const std::vector<double>& input) {
    std::vector<double> ref;
    ref.reserve(input.size());
    for (double x : input) {
        ref.push_back(x * x);
    }
    return ref;
}

std::size_t positiveArgument(const char* text, const char* name) {
    const std::string input(text);
    std::size_t consumed = 0;
    const auto value = std::stoull(input, &consumed);
    if (value == 0 || consumed != input.size() || input.front() == '-') {
        throw std::invalid_argument(std::string(name) + " must be positive");
    }
    return static_cast<std::size_t>(value);
}

std::vector<std::size_t> threadCounts(std::size_t maximum) {
    std::vector<std::size_t> counts;
    for (std::size_t count = 1; count <= maximum; count *= 2) {
        counts.push_back(count);
        if (count > maximum / 2) {
            break;
        }
    }
    if (counts.back() != maximum) {
        counts.push_back(maximum);
    }
    return counts;
}

} // namespace

int main(int argc, char** argv) try {
    if (argc > 4) {
        std::fprintf(stderr, "usage: %s [profile] [tasks] [max_threads]\n", argv[0]);
        return 2;
    }
    const std::string profile_name = argc > 1 ? argv[1] : "basic_ckks";
    const std::size_t task_count = argc > 2 ? positiveArgument(argv[2], "tasks") : 32;
    const std::size_t maximum_threads = argc > 3
        ? positiveArgument(argv[3], "max_threads") : 8;
    const std::size_t payload_size = 64;
    const auto profile = m2424::profiles::by_name(profile_name);

    const auto thread_counts = threadCounts(maximum_threads);

    std::printf("profile,tolerance,threads,tasks,payload_size,setup_ms,runtime_ms,"
                "throughput_tasks_per_sec,max_abs_error,mean_abs_error,status\n");

    bool all_ok = true;
    for (std::size_t thread_count : thread_counts) {
        const std::size_t worker_count = std::min(thread_count, task_count);
        std::vector<std::unique_ptr<WorkerState>> workers;
        workers.reserve(worker_count);

        const double setup_ms = elapsed_ms([&] {
            for (std::size_t worker = 0; worker < worker_count; ++worker) {
                auto state = std::make_unique<WorkerState>();
                state->adapter = m2424::SealAdapter::create(profile);
                state->adapter.generateKeys(true, true);
                workers.push_back(std::move(state));
            }

            for (std::size_t task = 0; task < task_count; ++task) {
                auto& state = *workers[task % worker_count];
                Task item;
                item.input = make_input(payload_size, task);
                item.encrypted = state.adapter.encrypt(state.adapter.encode(item.input));
                state.tasks.push_back(std::move(item));
            }
        });

        const double runtime_ms = elapsed_ms([&] {
            std::vector<std::thread> threads;
            threads.reserve(worker_count);
            for (std::size_t worker = 0; worker < worker_count; ++worker) {
                threads.emplace_back([&workers, worker] {
                    auto& state = *workers[worker];
                    for (auto& task : state.tasks) {
                        task.result = state.adapter.rescaleToNext(
                            state.adapter.relinearize(state.adapter.multiply(task.encrypted, task.encrypted)));
                    }
                });
            }
            for (auto& thread : threads) {
                thread.join();
            }
        });

        double max_abs_error = 0.0;
        double mean_abs_error_sum = 0.0;
        bool ok = true;
        std::size_t checked = 0;

        for (auto& state_ptr : workers) {
            auto& state = *state_ptr;
            for (const auto& task : state.tasks) {
                const auto actual = head(state.adapter.decode(state.adapter.decrypt(task.result)), payload_size);
                const auto expected = square_ref(task.input);
                const auto accuracy = m2424::compare(expected, actual, m2424::kTargetAbsoluteError);
                max_abs_error = std::max(max_abs_error, accuracy.max_abs_error);
                mean_abs_error_sum += accuracy.mean_abs_error;
                ok = ok && accuracy.ok;
                ++checked;
            }
        }

        const double throughput = runtime_ms > 0.0
            ? (1000.0 * static_cast<double>(task_count)) / runtime_ms
            : 0.0;
        const double mean_abs_error = checked > 0
            ? mean_abs_error_sum / static_cast<double>(checked)
            : 0.0;

        std::printf("%s,%.6e,%zu,%zu,%zu,%.6f,%.6f,%.6f,%.6e,%.6e,%s\n",
                    profile_name.c_str(),
                    m2424::kTargetAbsoluteError,
                    worker_count,
                    task_count,
                    payload_size,
                    setup_ms,
                    runtime_ms,
                    throughput,
                    max_abs_error,
                    mean_abs_error,
                    ok ? "PASS" : "FAIL");
        all_ok = all_ok && ok;
    }

    return all_ok ? 0 : 1;
} catch (const std::exception& error) {
    std::fprintf(stderr, "%s\n", error.what());
    return 2;
}
