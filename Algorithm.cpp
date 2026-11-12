#include <vector>
#include <unordered_map>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <cstdint>
#include <limits>
#include <string>
#include <functional>
#include <iomanip>
using namespace std;
using u64 = uint64_t;
using ld = long double;

constexpr ld GAMMA = 0.57721566490153286060651209L;
constexpr ld E_GAMMA = 1.7810724179901979852365041L;

vector<int> sieve(int limit)
{
    if (limit < 2)
        return {};
    vector<bool> isPrime(limit + 1, true);
    isPrime[0] = isPrime[1] = false;
    for (int i = 2; 1LL * i * i <= limit; ++i)
        if (isPrime[i])
            for (int j = i * i; j <= limit; j += i)
                isPrime[j] = false;
    vector<int> primes;
    for (int i = 2; i <= limit; ++i)
        if (isPrime[i])
            primes.push_back(i);
    return primes;
}

bool isPrime64(u64 x)
{
    if (x < 2)
        return false;
    for (u64 d = 2; d * d <= x; ++d)
        if (x % d == 0)
            return false;
    return true;
}

u64 nextPrime(u64 x)
{
    ++x;
    while (!isPrime64(x))
        ++x;
    return x;
}

vector<u64> distinctFactors(u64 x, const vector<int> &primes)
{
    vector<u64> result;
    for (int p : primes)
    {
        if (1ULL * p * p > x)
            break;
        if (x % p == 0)
        {
            result.push_back(p);
            do
                x /= p;
            while (x % p == 0);
        }
    }
    if (x > 1)
        result.push_back(x);
    return result;
}

ld radicalRatio(u64 x, const vector<int> &primes)
{
    ld r = 1.0L;
    for (u64 p : distinctFactors(x, primes))
        r *= (ld)p / (ld)(p - 1);
    return r;
}

ld objective(u64 n, const vector<int> &primes)
{
    ld t = logl(logl((ld)n));
    return (radicalRatio(n, primes) - E_GAMMA * t) * t;
}

struct Best
{
    u64 n = 0;
    ld value = -numeric_limits<ld>::infinity();
    void consider(u64 x, const vector<int> &primes)
    {
        ld v = objective(x, primes);
        if (n == 0 || v > value)
            n = x, value = v;
    }
};

int main(int argc, char **argv)
{
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    u64 R = (1ULL << 32);

    if (R < 3)
    {
        cerr << "R must be >= 3\n";
        return 1;
    }

    vector<int> smallPrimes = sieve(100);
    u64 primorial = 1;
    vector<int> smoothPrimes;
    for (int p : smallPrimes)
    {
        if (primorial > R / p)
            break;
        primorial *= p;
        smoothPrimes.push_back(p);
    }

    const u64 pmax = smoothPrimes.back();
    int sieveLimit = (int)sqrtl((ld)R) + 2;
    sieveLimit = max<int>(sieveLimit, (int)pmax + 2);
    vector<int> primes = sieve(sieveLimit);
    ;

    ld rmax = 1.0L;
    for (u64 p : smoothPrimes)
        rmax *= (ld)p / (ld)(p - 1);
    ld t0 = rmax / (2.0L * E_GAMMA);
    ld n0Real = expl(expl(t0));
    u64 n0 = (n0Real >= (ld)R) ? R : (u64)ceill(n0Real);

    Best best;

    if (n0 >= R)
    {
        for (u64 n = 3; n <= R; ++n)
            best.consider(n, primes);
        cout << setprecision(18) << "pmax=" << pmax << " primorial=" << primorial
             << " n0=" << n0 << "\nmax_n=" << best.n
             << "\nmax_f=" << best.value << '\n';
        return 0;
    }

    u64 K = n0 / 2;
    ld rExc = 0.0L;
    for (u64 k = 1; k <= K; ++k)
    {
        ld rk = radicalRatio(k, primes);
        u64 q = nextPrime(pmax);
        while (k % q == 0)
            q = nextPrime(q);
        rExc = max(rExc, rk * (ld)q / (ld)(q - 1));
    }

    ld tZero = rExc / E_GAMMA;
    ld logB = expl(tZero);
    u64 B;
    if (logB >= logl((ld)R))
        B = R;
    else
        B = min<u64>(R, (u64)ceill(expl(logB)));

    const u64 low = n0 + 1;

    for (u64 n = low; n <= B; ++n)
        best.consider(n, primes);

    unordered_map<u64, u64> minByRadical;
    minByRadical.reserve(300000);

    function<void(size_t, u64, u64)> generate = [&](size_t idx, u64 value, u64 rad)
    {
        if (idx == smoothPrimes.size())
        {
            if (value >= low)
            {
                auto [it, inserted] = minByRadical.emplace(rad, value);
                if (!inserted && value < it->second)
                    it->second = value;
            }
            return;
        }

        generate(idx + 1, value, rad);

        u64 p = (u64)smoothPrimes[idx];
        u64 cur = value;
        while (cur <= R / p && rad <= R / p)
        {
            cur *= p;
            generate(idx + 1, cur, rad * p);
        }
    };
    generate(0, 1, 1);

    size_t afterLemma2 = 0;
    for (const auto &[rad, n] : minByRadical)
    {
        vector<u64> factors = distinctFactors(rad, primes);
        bool reject = false;

        for (u64 q : factors)
        {
            for (int p : primes)
            {
                if ((u64)p >= q)
                    break;
                if (binary_search(factors.begin(), factors.end(), (u64)p))
                    continue;

                u64 replaced = (n / q) * (u64)p;
                if (replaced >= low)
                {
                    reject = true;
                    break;
                }
            }
            if (reject)
                break;
        }

        if (!reject)
        {
            ++afterLemma2;
            best.consider(n, primes);
        }
    }

    cout << setprecision(18);
    cout << "pmax=" << pmax << "\nprimorial=" << primorial
         << "\nrmax=" << rmax << "\nt0=" << t0
         << "\nn0=" << n0 << "\nr_exc=" << rExc
         << "\nB=" << B
         << "\nradicals_after_lemma1=" << minByRadical.size() << "\ncandidates_after_lemma2=" << afterLemma2
         << "\nmax_n=" << best.n
         << "\nmax_f=" << best.value << '\n';
    ;
}
