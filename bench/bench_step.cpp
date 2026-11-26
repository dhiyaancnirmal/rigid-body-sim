// Step-time benchmark: the reference pile at several sizes, with each
// broadphase. Prints JSON. Each configuration is rebuilt and replayed
// `repeats` times over the same window of steps; the fastest window is
// reported because noise on a shared machine only ever adds time.
#include <chrono>
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "rbsim/scenes.hpp"

using namespace rbsim;

namespace {

struct Result {
    double step_ms = 0.0;
    StepTimings phases;
    std::uint64_t hash = 0;
    double pairs = 0.0;
    double contacts = 0.0;
};

Result measure(int n, Broadphase bp, int warmup, int window, int repeats) {
    Result best;
    best.step_ms = 1e300;
    for (int r = 0; r < repeats; ++r) {
        World w = scenes::pile(n, bp);
        for (int i = 0; i < warmup; ++i) w.step(1.0 / 60.0);
        w.timings = StepTimings{};
        double pairs = 0.0, contacts = 0.0;
        const auto t0 = std::chrono::steady_clock::now();
        for (int i = 0; i < window; ++i) {
            w.step(1.0 / 60.0);
            pairs += static_cast<double>(w.pair_count());
            contacts += static_cast<double>(w.manifolds().size());
        }
        const auto t1 = std::chrono::steady_clock::now();
        const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count() / window;
        if (ms < best.step_ms) {
            best.step_ms = ms;
            best.phases = w.timings;
            best.hash = w.state_hash();
            best.pairs = pairs / window;
            best.contacts = contacts / window;
        }
    }
    return best;
}

}  // namespace

int main(int argc, char** argv) {
    std::vector<int> sizes = {250, 500, 1000, 2000, 4000};
    int warmup = 120, window = 120, repeats = 5;
    if (argc > 1 && std::string(argv[1]) == "--quick") {
        sizes = {250, 1000};
        repeats = 1;
        window = 30;
    }
    std::printf("{\n  \"warmup_steps\": %d,\n  \"window_steps\": %d,\n  \"repeats\": %d,\n", warmup,
                window, repeats);
    std::printf("  \"runs\": [\n");
    bool first = true;
    for (int n : sizes) {
        for (Broadphase bp : {Broadphase::BruteForce, Broadphase::SweepAndPrune}) {
            const Result r = measure(n, bp, warmup, window, repeats);
            const double steps = static_cast<double>(r.phases.steps);
            std::printf("%s    {\"bodies\": %d, \"broadphase\": \"%s\", \"step_ms\": %.4f, "
                        "\"broadphase_ms\": %.4f, \"narrowphase_ms\": %.4f, \"solve_ms\": %.4f, "
                        "\"integrate_ms\": %.4f, \"pairs\": %.1f, \"manifolds\": %.1f, "
                        "\"final_hash\": \"0x%016" PRIx64 "\"}",
                        first ? "" : ",\n", n,
                        bp == Broadphase::BruteForce ? "brute_force" : "sweep_and_prune",
                        r.step_ms, 1e3 * r.phases.broadphase_s / steps,
                        1e3 * r.phases.narrowphase_s / steps, 1e3 * r.phases.solve_s / steps,
                        1e3 * r.phases.integrate_s / steps, r.pairs, r.contacts, r.hash);
            std::fflush(stdout);
            first = false;
        }
    }
    std::printf("\n  ]\n}\n");
    return 0;
}
